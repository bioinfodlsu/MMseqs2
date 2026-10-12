/* Vendored subset of src/algo/blast/core/blast_stat.c
 * (github.com/ncbi/ncbi-cxx-toolkit-public HEAD 024628c, 2026-09-20).
 *
 * Extracted Karlin subset only (~330 of 5270 lines); function bodies below
 * are verbatim upstream, including shifted-pointer indexing
 * (sprob -= score_min; sprob[low..high]) and Int4 types:
 *   BlastScoreChk, Blast_ScoreFreqNew/Free,
 *   BlastKarlinLHtoK, NlmKarlinLambdaNR,
 *   Blast_KarlinLambdaNR, BlastKarlinLtoH,
 *   Blast_KarlinBlkUngappedCalc.
 * The bulk of blast_stat.c (score tables, gapped statistics, PSI, matrix
 * readers) is intentionally not vendored. Local includes stand in for
 * <algo/blast/core/...> + boost_erf.h + blast_psi_priv.h, which only the
 * excluded bulk needs.
 *
 * One marked MMSEQS PATCH: lattice compression in BlastKarlinLHtoK for
 * divisor > 1 (rare: needs all observed scores even). Upstream keeps the
 * original array and indexes it with scaled offsets, misreading for d>1
 * (synthetic even-only -4..4: K=3.29 vs Monte-Carlo ~0.08; fixed K=0.086).
 * divisor==1 paths are byte-identical logic. See lib/ncbi-blast/README.md.
 *
 * PUBLIC DOMAIN (US Government work) — see upstream header.
 */

#include "blast_stat.h"
#include "ncbi_math.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>

#define BLAST_SCORE_RANGE_MAX   (BLAST_SCORE_MAX - BLAST_SCORE_MIN) /**< maximum allowed range of BLAST scores. */

/****************************************************************************
For more accuracy in the calculation of K, set K_SUMLIMIT to 0.00001.
For high speed in the calculation of K, use a K_SUMLIMIT of 0.001
Note:  statistical significance is often not greatly affected by the value
of K, so high accuracy is generally unwarranted.
*****************************************************************************/
#define BLAST_KARLIN_K_SUMLIMIT_DEFAULT 0.0001 /**< K_SUMLIMIT_DEFAULT == sumlimit used in BlastKarlinLHtoK() */

#define BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT    (1.e-5) /**< LAMBDA_ACCURACY_DEFAULT == accuracy to which Lambda should be calc'd */

#define BLAST_KARLIN_LAMBDA_ITER_DEFAULT        17 /**< LAMBDA_ITER_DEFAULT == no. of iterations in LambdaBis = ln(accuracy)/ln(2)*/

#define BLAST_KARLIN_LAMBDA0_DEFAULT    0.5 /**< Initial guess for the value of Lambda in BlastKarlinLambdaNR */

#define BLAST_KARLIN_K_ITER_MAX 100 /**< upper limit on iterations for BlastKarlinLHtoK */

/** Check that the lo and hi score are within the allowed ranges
 * @param lo the lowest permitted value [in]
 * @param hi the highest permitted value [in]
 * @return zero on success, 1 otherwise
 */

static Int2
BlastScoreChk(Int4 lo, Int4 hi)
{
   if (lo >= 0 || hi <= 0 ||
         lo < BLAST_SCORE_MIN || hi > BLAST_SCORE_MAX)
      return 1;

   if (hi - lo > BLAST_SCORE_RANGE_MAX)
      return 1;

   return 0;
}

Blast_ScoreFreq*
Blast_ScoreFreqNew(Int4 score_min, Int4 score_max)
{
   Blast_ScoreFreq*  sfp;
   Int4  range;

   if (BlastScoreChk(score_min, score_max) != 0)
      return NULL;

   sfp = (Blast_ScoreFreq*) calloc(1, sizeof(Blast_ScoreFreq));
   if (sfp == NULL)
      return NULL;

   range = score_max - score_min + 1;
   sfp->sprob = (double*) calloc(range, sizeof(double));
   if (sfp->sprob == NULL)
   {
      Blast_ScoreFreqFree(sfp);
      return NULL;
   }

   sfp->sprob0 = sfp->sprob;
   sfp->sprob -= score_min;        /* center around 0 */
   sfp->score_min = score_min;
   sfp->score_max = score_max;
   sfp->obs_min = sfp->obs_max = 0;
   sfp->score_avg = 0.0;
   return sfp;
}

Blast_ScoreFreq*
Blast_ScoreFreqFree(Blast_ScoreFreq* sfp)
{
   if (sfp == NULL)
      return NULL;

   if (sfp->sprob0 != NULL)
      sfree(sfp->sprob0);
   sfree(sfp);
   return sfp;
}

/** The following procedure computes K. The input includes Lambda, H,
 *  and an array of probabilities for each score.
 *  There are distinct closed form for three cases:
 *  1. high score is 1 low score is -1
 *  2. high score is 1 low score is not -1
 *  3. low score is -1, high score is not 1
 *
 * Otherwise, in most cases the value is computed as:
 * -exp(-2.0*outerSum) / ((H/lambda)*(exp(-lambda) - 1)
 * The last term (exp(-lambda) - 1) can be computed in two different
 * ways depending on whether lambda is small or not.
 * outerSum is a sum of the terms
 * innerSum/j, where j is denoted by iterCounter in the code.
 * The sum is truncated when the new term innersum/j i sufficiently small.
 * innerSum is a weighted sum of the probabilities of
 * of achieving a total score i in a gapless alignment,
 * which we denote by P(i,j).
 * of exactly j characters. innerSum(j) has two parts
 * Sum over i < 0  P(i,j)exp(-i * lambda) +
 * Sum over i >=0  P(i,j)
 * The terms P(i,j) are computed by dynamic programming.
 * An earlier version was flawed in that ignored the special case 1
 * and tried to replace the tail of the computation of outerSum
 * by a geometric series, but the base of the geometric series
 * was not accurately estimated in some cases.
 *
 * @param sfp object holding scoring frequency information [in]
 * @param lambda a Karlin-Altschul parameter [in]
 * @param H a Karlin-Altschul parameter [in]
 * @return K, another Karlin-Altschul parameter
 */

static double
BlastKarlinLHtoK(Blast_ScoreFreq* sfp, double lambda, double H)
{
    /*The next array stores the probabilities of getting each possible
      score in an alignment of fixed length; the array is shifted
      during part of the computation, so that
      entry 0 is for score 0.  */
    double         *alignmentScoreProbabilities = NULL;
    Int4            low;    /* Lowest score (must be negative) */
    Int4            high;   /* Highest score (must be positive) */
    Int4            range;  /* range of scores, computed as high - low*/
    double          K;      /* local copy of K  to return*/
    int             i;   /*loop index*/
    int             iterCounter; /*counter on iterations*/
    Int4            divisor; /*candidate divisor of all scores with
                               non-zero probabilities*/
    /*highest and lowest possible alignment scores for current length*/
    Int4            lowAlignmentScore, highAlignmentScore;
    Int4            first, last; /*loop indices for dynamic program*/
    register double innerSum;
    double          oldsum, oldsum2;  /* values of innerSum on previous
                                         iterations*/
    double          outerSum;        /* holds sum over j of (innerSum
                                        for iteration j/j)*/

    double          score_avg; /*average score*/
    /*first term to use in the closed form for the case where
      high == 1 or low == -1, but not both*/
    double          firstTermClosedForm;  /*usually store H/lambda*/
    int             iterlimit; /*upper limit on iterations*/
    double          sumlimit; /*lower limit on contributions
                                to sum over scores*/

    /*array of score probabilities reindexed so that low is at index 0*/
    double         *probArrayStartLow;

    /*pointers used in dynamic program*/
    double         *ptrP, *ptr1, *ptr2, *ptr1e;
    double          expMinusLambda; /*e^^(-Lambda) */

    if (lambda <= 0. || H <= 0.) {
        /* Theory dictates that H and lambda must be positive, so
         * return -1 to indicate an error */
        return -1.;
    }

    /*Karlin-Altschul theory works only if the expected score
      is negative*/
    if (sfp->score_avg >= 0.0) {
        return -1.;
    }

    low   = sfp->obs_min;
    high  = sfp->obs_max;
    range = high - low;

    probArrayStartLow = &sfp->sprob[low];
    /* Look for the greatest common divisor ("delta" in Appendix of PNAS 87 of
       Karlin&Altschul (1990) */
    for (i = 1, divisor = -low; i <= range && divisor > 1; ++i) {
        if (probArrayStartLow[i] != 0.0)
            divisor = BLAST_Gcd(divisor, i);
    }

    high   /= divisor;
    low    /= divisor;
    lambda *= divisor;

    range = high - low;

    /* MMSEQS PATCH (lattice divisor > 1, rare: needs all observed scores
       even): upstream keeps the original array and indexes it with scaled
       offsets (probArrayStartLow[first] with scaled first; special case
       excepted since it uses absolute sprob[low*divisor]). That misreads
       for divisor > 1, so compress to the unit lattice here:
       P_scaled[k] = P_orig[k*divisor]. Off-lattice entries are zero by
       construction of divisor, so compression is lossless. divisor==1
       paths are unaffected (pointer rebased to the same array). */
    double *scaledBase0 = NULL;
    double *scaledBase = NULL;
    {
        Int4 origLow = sfp->obs_min;
        double *origProb = &sfp->sprob[origLow];
        scaledBase0 = (double*) calloc((size_t)range + 1, sizeof(double));
        if (scaledBase0 == NULL)
            return -1.;
        scaledBase = scaledBase0;
        for (i = 0; i <= range; ++i)
            scaledBase[i] = origProb[(size_t)i * (size_t)divisor];
        probArrayStartLow = scaledBase;
    }

    firstTermClosedForm = H/lambda;
    expMinusLambda      = exp((double) -lambda);

    if (low == -1 && high == 1) {
        K = (probArrayStartLow[0] - probArrayStartLow[range]) *
            (probArrayStartLow[0] - probArrayStartLow[range]) / probArrayStartLow[0];
        sfree(scaledBase0);
        return(K);
    }

    if (low == -1 || high == 1) {
        if (high != 1) {
            score_avg = sfp->score_avg / divisor;
            firstTermClosedForm
                = (score_avg * score_avg) / firstTermClosedForm;
        }
        sfree(scaledBase0);
        return firstTermClosedForm * (1.0 - expMinusLambda);
    }

    sumlimit  = BLAST_KARLIN_K_SUMLIMIT_DEFAULT;
    iterlimit = BLAST_KARLIN_K_ITER_MAX;

    alignmentScoreProbabilities =
        (double *)calloc((iterlimit*range + 1), sizeof(*alignmentScoreProbabilities));
    if (alignmentScoreProbabilities == NULL) {
        sfree(scaledBase0);
        return -1.;
    }

    outerSum = 0.;
    lowAlignmentScore = highAlignmentScore = 0;
    alignmentScoreProbabilities[0] = innerSum = oldsum = oldsum2 = 1.;

    for (iterCounter = 0;
         ((iterCounter < iterlimit) && (innerSum > sumlimit));
         outerSum += innerSum /= ++iterCounter) {
        first = last = range;
        lowAlignmentScore  += low;
        highAlignmentScore += high;
        /*dynamic program to compute P(i,j)*/
        for (ptrP = alignmentScoreProbabilities +
                 (highAlignmentScore-lowAlignmentScore);
             ptrP >= alignmentScoreProbabilities;
             *ptrP-- =innerSum) {
            ptr1  = ptrP - first;
            ptr1e = ptrP - last;
            ptr2  = probArrayStartLow + first;
            for (innerSum = 0.; ptr1 >= ptr1e; ) {
                innerSum += *ptr1  *  *ptr2;
		ptr1--;
		ptr2++;
            }
            if (first)
                --first;
            if (ptrP - alignmentScoreProbabilities <= range)
                --last;
        }
        /* Horner's rule */
        innerSum = *++ptrP;
        for( i = lowAlignmentScore + 1; i < 0; i++ ) {
            innerSum = *++ptrP + innerSum * expMinusLambda;
        }
        innerSum *= expMinusLambda;

        for (; i <= highAlignmentScore; ++i)
            innerSum += *++ptrP;
        oldsum2 = oldsum;
        oldsum  = innerSum;
    }

#ifdef ADD_GEOMETRIC_TERMS_TO_K
    /*old code assumed that the later terms in sum were
      asymptotically comparable to those of a geometric
      progression, and tried to speed up convergence by
      guessing the estimated ratio between sucessive terms
      and using the explicit terms of a geometric progression
      to speed up convergence. However, the assumption does not
      always hold, and convergenece of the above code is fast
      enough in practice*/
    /* Terms of geometric progression added for correction */
    {
        double     ratio;  /* fraction used to generate the
                                   geometric progression */

        ratio = oldsum / oldsum2;
        if (ratio >= (1.0 - sumlimit*0.001)) {
            K = -1.;
            if (alignmentScoreProbabilities != NULL)
                sfree(alignmentScoreProbabilities);
            return K;
        }
        sumlimit *= 0.01;
        while (innerSum > sumlimit) {
            oldsum   *= ratio;
            outerSum += innerSum = oldsum / ++iterCounter;
        }
    }
#endif

    K = -exp((double)-2.0*outerSum) /
             (firstTermClosedForm*BLAST_Expm1(-(double)lambda));

    if (alignmentScoreProbabilities != NULL)
        sfree(alignmentScoreProbabilities);
    sfree(scaledBase0);

    return K;
}


/**
 * Find positive solution to
 *
 *     sum_{i=low}^{high} exp(i lambda) * probs[i] = 1.
 *
 * Note that this solution does not exist unless the average score is
 * negative and the largest score that occurs with nonzero probability
 * is positive.
 *
 * @param probs         probabilities of a score occurring
 * @param d             the gcd of the possible scores. This equals 1 if
 *                      the scores are not a lattice
 * @param low           the lowest possible score that occurs with
 *                      nonzero probability
 * @param high          the highest possible score that occurs with
 *                      nonzero probability.
 * @param lambda0       an initial guess for lambda
 * @param tolx          the tolerance to which lambda must be computed
 * @param itmax         the maximum number of times the function may be
 *                      evaluated
 * @param maxNewton     the maximum permissible number of Newton
 *                      iterations; after that the computation will proceed
 *                      by bisection.
 * @param *itn          the number of iterations needed to compute Lambda,
 *                      or itmax if Lambda could not be computed.
 *
 * Let phi(lambda) =  sum_{i=low}^{high} exp(i lambda) - 1. Then
 * phi(lambda) may be written
 *
 *     phi(lambda) = exp(u lambda) f( exp(-lambda) )
 *
 * where f(x) is a polynomial that has exactly two zeros, one at x = 1
 * and one at x = exp(-lambda).  It is simpler to solve this problem
 * in x = exp(-lambda) than it is to solve it in lambda, because we
 * know that for x, a solution lies in [0,1], and because Newton's
 * method is generally more stable and efficient for polynomials than
 * it is for exponentials.
 *
 * For the most part, this function is a standard safeguarded Newton
 * iteration: define an interval of uncertainty [a,b] with f(a) > 0
 * and f(b) < 0 (except for the initial value b = 1, where f(b) = 0);
 * evaluate the function and use the sign of that value to shrink the
 * interval of uncertainty; compute a Newton step; and if the Newton
 * step suggests a point outside the interval of uncertainty or fails
 * to decrease the function sufficiently, then bisect.  There are
 * three further details needed to understand the algorithm:
 *
 * 1)  If y the unique solution in [0,1], then f is positive to the left of
 *     y, and negative to the right.  Therefore, we may determine whether
 *     the Newton step -f(x)/f'(x) is moving toward, or away from, y by
 *     examining the sign of f'(x).  If f'(x) >= 0, we bisect instead
 *     of taking the Newton step.
 * 2)  There is a neighborhood around x = 1 for which f'(x) >= 0, so
 *     (1) prevents convergence to x = 1 (and for a similar reason
 *     prevents convergence to x = 0, if the function is incorrectly
 *     called with probs[high] == 0).
 * 3)  Conditions like  fabs(p) < lambda_tolerance * x * (1-x) are used in
 *     convergence criteria because these values translate to a bound
 *     on the relative error in lambda.  This is proved in the
 *     "Blast Scoring Parameters" document that accompanies the BLAST
 *     code.
 *
 * The iteration on f(x) is robust and doesn't overflow; defining a
 * robust safeguarded Newton iteration on phi(lambda) that cannot
 * converge to lambda = 0 and that is protected against overflow is
 * more difficult.  So (despite the length of this comment) the Newton
 * iteration on f(x) is the simpler solution.
 */
static double
NlmKarlinLambdaNR(double* probs, Int4 d, Int4 low, Int4 high, double lambda0,
                  double tolx, Int4 itmax, Int4 maxNewton, Int4 * itn )
{
  Int4 k;
  double x0, x, a = 0, b = 1;
  double f = 4;  /* Larger than any possible value of the poly in [0,1] */
  Int4 isNewton = 0; /* we haven't yet taken a Newton step. */

  assert( d > 0 );

   x0 = exp( -lambda0 );
  x = ( 0 < x0 && x0 < 1 ) ? x0 : .5;

  for( k = 0; k < itmax; k++ ) { /* all iteration indices k */
    Int4 i;
    double g, fold = f;
    Int4 wasNewton = isNewton; /* If true, then the previous step was a */
                              /* Newton step */
    isNewton  = 0;            /* Assume that this step is not */

    /* Horner's rule for evaluating a polynomial and its derivative */
    g = 0;
    f = probs[low];
    for( i = low + d; i < 0; i += d ) {
      g = x * g + f;
      f = f * x + probs[i];
    }
    g = x * g + f;
    f = f * x + probs[0] - 1;
    for( i = d; i <= high; i += d ) {
      g = x * g + f;
      f = f * x + probs[i];
    }
    /* End Horner's rule */

    if( f > 0 ) {
      a = x; /* move the left endpoint */
    } else if( f < 0 ) {
      b = x; /* move the right endpoint */
    } else { /* f == 0 */
      break; /* x is an exact solution */
    }
    if( b - a < 2 * a * ( 1 - b ) * tolx ) {
      /* The midpoint of the interval converged */
      x = (a + b) / 2; break;
    }

    if( k >= maxNewton ||
        /* If convergence of Newton's method appears to be failing; or */
            ( wasNewton && fabs( f ) > .9 * fabs(fold) ) ||
        /* if the previous iteration was a Newton step but didn't decrease
         * f sufficiently; or */
        g >= 0
        /* if a Newton step will move us away from the desired solution */
        ) { /* then */
      /* bisect */
      x = (a + b)/2;
    } else {
      /* try a Newton step */
      double p = - f/g;
      double y = x + p;
      if( y <= a || y >= b ) { /* The proposed iterate is not in (a,b) */
        x = (a + b)/2;
      } else { /* The proposed iterate is in (a,b). Accept it. */
        isNewton = 1;
        x = y;
        if( fabs( p ) < tolx * x * (1-x) ) break; /* Converged */
      } /* else the proposed iterate is in (a,b) */
    } /* else try a Newton step. */
  } /* end for all iteration indices k */
   *itn = k;
  return -log(x)/d;
}


double
Blast_KarlinLambdaNR(Blast_ScoreFreq* sfp, double initialLambdaGuess)
{
   Int4  low;        /* Lowest score (must be negative)  */
   Int4  high;       /* Highest score (must be positive) */
   Int4     itn;
   Int4  i, d;
   double*  sprob;
   double   returnValue;

   low = sfp->obs_min;
   high = sfp->obs_max;
   if (sfp->score_avg >= 0.) {   /* Expected score must be negative */
      return -1.0;
   }
   if (BlastScoreChk(low, high) != 0) return -1.;

   sprob = sfp->sprob;
   /* Find greatest common divisor of all scores */
   for (i = 1, d = -low; i <= high-low && d > 1; ++i) {
      if (sprob[i+low] != 0.0) {
         d = BLAST_Gcd(d, i);
      }
   }
   returnValue =
      NlmKarlinLambdaNR( sprob, d, low, high,
                           initialLambdaGuess,
                           BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT,
                     20, 20 + BLAST_KARLIN_LAMBDA_ITER_DEFAULT, &itn );


   return returnValue;
}

/** Calculate H, the relative entropy of the p's and q's
 *
 * @param sfp object containing scoring frequency information [in]
 * @param lambda a Karlin-Altschul parameter [in]
 * @return H, a Karlin-Altschul parameter
 */
static double
BlastKarlinLtoH(Blast_ScoreFreq* sfp, double lambda)
{
   Int4  score;
   double   H, etonlam, sum, scale;

   double *probs = sfp->sprob;
   Int4 low   = sfp->obs_min,  high  = sfp->obs_max;

   if (lambda < 0.) {
      return -1.;
   }
   if (BlastScoreChk(low, high) != 0) return -1.;

   etonlam = exp( - lambda );
  sum = low * probs[low];
  for( score = low + 1; score <= high; score++ ) {
    sum = score * probs[score] + etonlam * sum;
  }

  scale = BLAST_Powi( etonlam, high );
  if( scale > 0.0 ) {
    H = lambda * sum/scale;
  } else { /* Underflow of exp( -lambda * high ) */
    H = lambda * exp( lambda * high + log(sum) );
  }
   return H;
}

/**************** Statistical Significance Parameter Subroutine ****************

    Version 1.0     February 2, 1990
    Version 1.2     July 6,     1990

    Program by:     Stephen Altschul

    Address:        National Center for Biotechnology Information
                    National Library of Medicine
                    National Institutes of Health
                    Bethesda, MD  20894

    Internet:       altschul@ncbi.nlm.nih.gov

See:  Karlin, S. & Altschul, S.F. "Methods for Assessing the Statistical
    Significance of Molecular Sequence Features by Using General Scoring
    Schemes,"  Proc. Natl. Acad. Sci. USA 87 (1990), 2264-2268.

    Computes the parameters lambda and K for use in calculating the
    statistical significance of high-scoring segments or subalignments.

    The scoring scheme must be integer valued.  A positive score must be
    possible, but the expected (mean) score must be negative.

    A program that calls this routine must provide the value of the lowest
    possible score, the value of the greatest possible score, and a pointer
    to an array of probabilities for the occurrence of all scores between
    these two extreme scores.  For example, if score -2 occurs with
    probability 0.7, score 0 occurs with probability 0.1, and score 3
    occurs with probability 0.2, then the subroutine must be called with
    low = -2, high = 3, and pr pointing to the array of values
    { 0.7, 0.0, 0.1, 0.0, 0.0, 0.2 }.  The calling program must also provide
    pointers to lambda and K; the subroutine will then calculate the values
    of these two parameters.  In this example, lambda=0.330 and K=0.154.

    The parameters lambda and K can be used as follows.  Suppose we are
    given a length N random sequence of independent letters.  Associated
    with each letter is a score, and the probabilities of the letters
    determine the probability for each score.  Let S be the aggregate score
    of the highest scoring contiguous segment of this sequence.  Then if N
    is sufficiently large (greater than 100), the following bound on the
    probability that S is greater than or equal to x applies:

            P( S >= x )   <=   1 - exp [ - KN exp ( - lambda * x ) ].

    In other words, the p-value for this segment can be written as
    1-exp[-KN*exp(-lambda*S)].

    This formula can be applied to pairwise sequence comparison by assigning
    scores to pairs of letters (e.g. amino acids), and by replacing N in the
    formula with N*M, where N and M are the lengths of the two sequences
    being compared.

    In addition, letting y = KN*exp(-lambda*S), the p-value for finding m
    distinct segments all with score >= S is given by:

                           2             m-1           -y
            1 - [ 1 + y + y /2! + ... + y   /(m-1)! ] e

    Notice that for m=1 this formula reduces to 1-exp(-y), which is the same
    as the previous formula.

*******************************************************************************/
Int2
Blast_KarlinBlkUngappedCalc(Blast_KarlinBlk* kbp, Blast_ScoreFreq* sfp)
{


   if (kbp == NULL || sfp == NULL)
      return 1;

   /* Calculate the parameter Lambda */

   kbp->Lambda = Blast_KarlinLambdaNR(sfp, BLAST_KARLIN_LAMBDA0_DEFAULT);
   if (kbp->Lambda < 0.)
      goto ErrExit;


   /* Calculate H */

   kbp->H = BlastKarlinLtoH(sfp, kbp->Lambda);
   if (kbp->H < 0.)
      goto ErrExit;


   /* Calculate K and log(K) */

   kbp->K = BlastKarlinLHtoK(sfp, kbp->Lambda, kbp->H);
   if (kbp->K < 0.)
      goto ErrExit;
   kbp->logK = log(kbp->K);

   /* Normal return */
   return 0;

ErrExit:
   kbp->Lambda = kbp->H = kbp->K = -1.;
   kbp->logK = HUGE_VAL;
   return 1;
}
