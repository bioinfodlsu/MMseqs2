/* Vendored subset of include/algo/blast/core/blast_stat.h
 * (github.com/ncbi/ncbi-cxx-toolkit-public HEAD 024628c, 2026-09-20).
 * Blast_KarlinBlk / Blast_ScoreFreq layouts, score bounds and the
 * Blast_KarlinBlkUngappedCalc prototype below are verbatim; the
 * gapped/matrix/table declarations of the full header are intentionally not
 * vendored (not needed by the ungapped path). The upstream
 * NCBI_XBLAST_EXPORT attribute is omitted (static object library).
 */
#ifndef MMSEQS_VENDOR_BLAST_STAT_H
#define MMSEQS_VENDOR_BLAST_STAT_H

#include "ncbi_std.h"
#include "ncbi_math.h"

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

#define BLAST_SCORE_MIN INT2_MIN   /**< minimum allowed score (for one letter comparison). */
#define BLAST_SCORE_MAX INT2_MAX   /**< maximum allowed score (for one letter comparison). */

/**
  Structure to hold the Karlin-Altschul parameters.
*/
typedef struct Blast_KarlinBlk {
      double   Lambda; /**< Lambda value used in statistics */
      double   K; /**< K value used in statistics */
      double   logK; /**< natural log of K value used in statistics */
      double   H; /**< H value used in statistics */
      double   paramC;  /**< for use in seed. */
   } Blast_KarlinBlk;

/** Holds score frequencies used in calculation
of Karlin-Altschul parameters for an ungapped search.
*/
typedef struct Blast_ScoreFreq {
    Int4         score_min; /**< lowest allowed scores */
    Int4         score_max; /**< highest allowed scores */
    Int4         obs_min;   /**< lowest observed (actual) scores */
    Int4         obs_max;   /**< highest observed (actual) scores */
    double       score_avg; /**< average score, must be negative for local alignment. */
    double*      sprob0;    /**< arrays for frequency of given score */
    double*      sprob;     /**< arrays for frequency of given score, shifted down by score_min. */
} Blast_ScoreFreq;

/* Plain free() stands in for upstream sfree() (which also nulls). */
#ifndef sfree
#define sfree(p) (free(p))
#endif

Blast_ScoreFreq* Blast_ScoreFreqNew(Int4 score_min, Int4 score_max);
Blast_ScoreFreq* Blast_ScoreFreqFree(Blast_ScoreFreq* sfp);
double Blast_KarlinLambdaNR(Blast_ScoreFreq* sfp, double initialLambdaGuess);

/** Computes the parameters lambda, H K for use in calculating the
 * statistical significance of high-scoring segments or subalignments (see
 * comment on blast_stat.c for more details).
 * @param kbp object containing Lambda, H, and K as well as scoring information [in|out]
 * @param sfp array of probabilities for all scores [in]
 * @return zero on success, 1 on error.
 */
Int2
Blast_KarlinBlkUngappedCalc(Blast_KarlinBlk* kbp, Blast_ScoreFreq* sfp);

#ifdef __cplusplus
}
#endif

#endif /*MMSEQS_VENDOR_BLAST_STAT_H*/
