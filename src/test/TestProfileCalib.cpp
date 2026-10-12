// Correctness tests for ProfileCalibrator (analytic NCBI path).
//
// Strategy: a synthetic profile whose columns cycle through residues
// proportionally to the background reproduces the BLOSUM62 joint column
// distribution, so the vendored NCBI Karlin code must recover the hardcoded
// BLOSUM62 gapless statistics (EvalueComputation.h: lambda=0.3207378152604,
// K=0.13904657125294345). Deterministic: no RNG in the calibration itself.
// Exit 0 on pass, EXIT_FAILURE with diagnostics otherwise.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>

#include <unistd.h>

#include "Parameters.h"
#include "Sequence.h"
#include "StripedSmithWaterman.h"
#include "SubstitutionMatrix.h"
#include "EvalueComputation.h"
#include "ProfileCalibrator.h"
#include "NcbiKarlin.h"

const char* binary_name = "test_profilecalib";
DEFAULT_PARAMETER_SINGLETON_INIT

// Hardcoded BLOSUM62 gapless reference (EvalueComputation.h defaultParameter).
static const double REF_LAMBDA = 0.3207378152604042354;
static const double REF_K = 0.13904657125294345166;

static int failures = 0;
#define CHECK(cond, ...) do { \
    if (!(cond)) { \
        printf("FAIL %s:%d: ", __FILE__, __LINE__); \
        printf(__VA_ARGS__); \
        printf("\n"); \
        ++failures; \
    } else { \
        printf("ok %s:%d: ", __FILE__, __LINE__); \
        printf(__VA_ARGS__); \
        printf("\n"); \
    } \
} while (0)

// Synthetic profile: column i scores like BLOSUM row colResidue[i].
static void buildSyntheticProfile(Sequence &q, SubstitutionMatrix &subMat,
                                  const std::vector<int> &colResidue) {
    const int L = (int)colResidue.size();
    q.L = L;
    for (int i = 0; i < L; ++i) {
        q.numSequence[i] = (unsigned char)colResidue[i];
    }
    for (int a = 0; a < subMat.alphabetSize - 1; ++a) {
        int8_t *row = q.profile_for_alignment + (size_t)a * (size_t)L;
        for (int i = 0; i < L; ++i) {
            row[i] = (int8_t)subMat.subMatrix[a][colResidue[i]];
        }
    }
    memset(q.profile_for_alignment
           + (size_t)(subMat.alphabetSize - 1) * (size_t)L, 0, (size_t)L);
}

// Columns drawn i.i.d. ~ background: emulates the BLOSUM joint column
// distribution. Order does not affect the aggregate distribution (hence
// lambda/K), only the validation simulation.
static std::vector<int> backgroundColumns(SubstitutionMatrix &subMat, int L) {
    std::mt19937_64 rng(20260918ULL);
    double total = 0.0;
    for (int a = 0; a < subMat.alphabetSize - 1; ++a) {
        total += subMat.pBack[a];
    }
    std::vector<int> cols;
    cols.reserve((size_t)L);
    for (int i = 0; i < L; ++i) {
        double u = (double)(rng() - rng.min()) / (double)(rng.max() - rng.min())
                * total;
        double acc = 0.0;
        int pick = subMat.alphabetSize - 2;
        for (int a = 0; a < subMat.alphabetSize - 1; ++a) {
            acc += subMat.pBack[a];
            if (u < acc) { pick = a; break; }
        }
        cols.push_back(pick);
    }
    return cols;
}

// Null target draw for uniformity validation only (test helper, not used by
// the calibrator itself which is fully analytic).
static void drawNullTarget(std::mt19937_64 &rng, unsigned char *target, int n,
                           const double *pBack, int alphabetSize) {
    const int states = alphabetSize - 1;
    double total = 0.0;
    for (int a = 0; a < states; ++a) total += pBack[a];
    for (int i = 0; i < n; ++i) {
        double u = (double)(rng() - rng.min()) / (double)(rng.max() - rng.min())
                * total;
        double acc = 0.0;
        int pick = states - 1;
        for (int a = 0; a < states; ++a) {
            acc += pBack[a];
            if (u < acc) { pick = a; break; }
        }
        target[i] = (unsigned char)pick;
    }
}

int main(int, const char**) {
    Parameters &par = Parameters::getInstance();
    par.initMatrices();
    SubstitutionMatrix subMat(par.scoringMatrixFile.values.aminoacid().c_str(),
                              2.0, 0.0);
    const int alph = subMat.alphabetSize;  // 21 incl. X
    printf("alphabetSize=%d\n", alph);

    const int L = 200;
    const int nLen = 200;

    Sequence qSeq(2048, Parameters::DBTYPE_HMM_PROFILE, &subMat, 0, false,
                  false);
    buildSyntheticProfile(qSeq, subMat, backgroundColumns(subMat, L));
    SmithWaterman sw(2048, alph, false, 0.0f, NULL);
    sw.ssw_init(&qSeq, qSeq.profile_for_alignment, &subMat);

    // 1. Recovery: analytic (lambda, K) close to hardcoded BLOSUM62 gapless.
    ProfileCalibrator::Calib cal = ProfileCalibrator::calibrate(
            qSeq.profile_for_alignment, L, alph, subMat.pBack);
    CHECK(cal.ok, "calibration ok (H=%.4f)", cal.H);
    if (cal.ok) {
        const double lambdaErr = fabs(cal.lambda - REF_LAMBDA) / REF_LAMBDA;
        const double kErr = fabs(cal.K - REF_K) / REF_K;
        printf("analytic lambda=%.6f (ref %.6f) K=%.6f (ref %.6f) H=%.4f\n",
               cal.lambda, REF_LAMBDA, cal.K, REF_K, cal.H);
        CHECK(lambdaErr < 0.02, "lambda within 2%% (err %.4f)", lambdaErr);
        CHECK(kErr < 0.20, "K within 20%% (err %.4f)", kErr);
        // Deterministic: recalibration is bit-identical.
        ProfileCalibrator::Calib cal2 = ProfileCalibrator::calibrate(
                qSeq.profile_for_alignment, L, alph, subMat.pBack);
        CHECK(cal2.ok && cal2.lambda == cal.lambda && cal2.K == cal.K,
              "recalibration bit-identical");
        // 2. Operational agreement vs SLS reference at mid-tail scores.
        EvalueComputation refEval(1000000, &subMat);  // cached gapless params
        for (int s = 0; s < 3; ++s) {
            const double score = (s == 0) ? 25.0 : (s == 1) ? 40.0 : 60.0;
            const double pRef = refEval.computePvalue(score, L, nLen);
            const double pFit = EvalueComputation::computeProfilePvalue(score, L, nLen, cal.lambda, cal.K);
            const double ratio = (pRef > 0 && pFit > 0) ? pFit / pRef : -1.0;
            printf("score=%.0f pRef=%.4g pFit=%.4g ratio=%.3f\n",
                   score, pRef, pFit, ratio);
            if (s > 0) {
                CHECK(ratio > 0.50 && ratio < 2.0,
                      "p agreement within 2x at score %.0f (ratio %.3f)",
                      score, ratio);
            }
        }
        // 3. Uniformity under the null at small p (operating range), scored
        // with the production kernel.
        std::mt19937_64 rng(12345);
        std::vector<unsigned char> tgt((size_t)nLen);
        int below01 = 0, below05 = 0;
        const int NU = 3000;
        std::vector<int> nullScores;
        nullScores.reserve(NU);
        for (int i = 0; i < NU; ++i) {
            drawNullTarget(rng, tgt.data(), nLen, subMat.pBack, alph);
            const int sc = sw.ungapped_alignment(tgt.data(), nLen);
            nullScores.push_back(sc);
            const double p = EvalueComputation::computeProfilePvalue(sc, L, nLen, cal.lambda, cal.K);
            if (p < 0.01) ++below01;
            if (p < 0.05) ++below05;
        }
        std::sort(nullScores.begin(), nullScores.end());
        printf("null score quantiles: min=%d p50=%d p95=%d p99=%d max=%d\n",
               nullScores.front(), nullScores[NU/2], nullScores[(NU*95)/100],
               nullScores[(NU*99)/100], nullScores.back());
        // Textbook control: random SEQUENCE query vs random targets must
        // follow SLS BLOSUM theory. Isolates harness bugs from fixed-query
        // null effects.
        {
            int8_t tiny[21 * 21];
            for (int a = 0; a < 21; ++a)
                for (int b = 0; b < 21; ++b)
                    tiny[a * 21 + b] = (int8_t)subMat.subMatrix[a][b];
            Sequence qs(2048, Parameters::DBTYPE_AMINO_ACIDS, &subMat, 0,
                        false, false);
            std::mt19937_64 rq(999);
            std::vector<unsigned char> qnum((size_t)L);
            drawNullTarget(rq, qnum.data(), L, subMat.pBack, alph);
            for (int i = 0; i < L; ++i) qs.numSequence[i] = qnum[i];
            qs.L = L;
            SmithWaterman sws(2048, alph, false, 0.0f, NULL);
            sws.ssw_init(&qs, tiny, &subMat);
            EvalueComputation refEval2(1000000, &subMat);
            int cBelow05 = 0;
            std::vector<int> cs;
            for (int i = 0; i < 1500; ++i) {
                drawNullTarget(rng, tgt.data(), nLen, subMat.pBack, alph);
                const int sc = sws.ungapped_alignment(tgt.data(), nLen);
                cs.push_back(sc);
                if (refEval2.computePvalue(sc, L, nLen) < 0.05) ++cBelow05;
            }
            std::sort(cs.begin(), cs.end());
            printf("seq-control: f(SLS p<0.05)=%.4f quant min=%d p50=%d p95=%d max=%d\n",
                   (double)cBelow05 / 1500, cs.front(), cs[750],
                   cs[(1500*95)/100], cs.back());
            CHECK(fabs((double)cBelow05 / 1500 - 0.05) < 0.03,
                  "seq control follows SLS theory");
        }
        const double f01 = (double)below01 / NU, f05 = (double)below05 / NU;
        printf("uniformity: f(p<0.01)=%.4f f(p<0.05)=%.4f\n", f01, f05);
        CHECK(f01 > 0.003 && f01 < 0.03, "null f(p<0.01) in [0.003,0.03]");
        CHECK(f05 > 0.02 && f05 < 0.08, "null f(p<0.05) in [0.02,0.08]");
        // 4. Representation equivalence: identical query content scored as
        // profile must equal sequence-query scoring (guards packing bugs).
        {
            std::vector<int> cols = backgroundColumns(subMat, 120);
            Sequence qp(2048, Parameters::DBTYPE_HMM_PROFILE, &subMat, 0,
                        false, false);
            buildSyntheticProfile(qp, subMat, cols);
            SmithWaterman swp(2048, alph, false, 0.0f, NULL);
            swp.ssw_init(&qp, qp.profile_for_alignment, &subMat);
            int8_t tiny[21 * 21];
            for (int a = 0; a < 21; ++a)
                for (int b = 0; b < 21; ++b)
                    tiny[a * 21 + b] = (int8_t)subMat.subMatrix[a][b];
            Sequence qs(2048, Parameters::DBTYPE_AMINO_ACIDS, &subMat, 0,
                        false, false);
            for (int i = 0; i < 120; ++i) qs.numSequence[i] = (unsigned char)cols[i];
            qs.L = 120;
            SmithWaterman sws(2048, alph, false, 0.0f, NULL);
            sws.ssw_init(&qs, tiny, &subMat);
            std::mt19937_64 rq(31337);
            int mism = 0, compared = 0;
            long sumAbsDiff = 0;
            for (int i = 0; i < 300; ++i) {
                drawNullTarget(rq, tgt.data(), nLen, subMat.pBack, alph);
                const int sp = swp.ungapped_alignment(tgt.data(), nLen);
                const int ss = sws.ungapped_alignment(tgt.data(), nLen);
                ++compared;
                sumAbsDiff += labs((long)sp - (long)ss);
                if (sp != ss) {
                    if (mism < 3) {
                        printf("repr mismatch #%d: profile=%d seq=%d\n",
                               mism, sp, ss);
                    }
                    ++mism;
                }
            }
            printf("representation: %d/%d mismatched meanAbsDiff=%.2f\n",
                   mism, compared, (double)sumAbsDiff / compared);
            CHECK(mism == 0, "profile == sequence representation");
        }
        // 5. Cache round-trip exact.
        const std::string tmp = std::string("/tmp/test_profilecalib_")
                + std::to_string((long long)getpid()) + ".tsv";
        std::map<uint64_t, ProfileCalibrator::Calib> out, back;
        const uint64_t h = ProfileCalibrator::hashProfile(
                qSeq.profile_for_alignment, L, alph);
        out[h] = cal;
        CHECK(ProfileCalibrator::saveCache(tmp, out), "saveCache");
        CHECK(ProfileCalibrator::loadCache(tmp, back), "loadCache");
        CHECK(back.count(h) == 1 && back[h].lambda == cal.lambda
              && back[h].K == cal.K && back[h].H == cal.H,
              "round-trip values identical");
        remove(tmp.c_str());
        // 6. Composition dependence: poly-W calibrates fine but differently
        // (rare residue, self-score 11) -- why the gate is per-profile.
        {
            const int rW = subMat.aa2num[(int)'W'];
            Sequence qW(2048, Parameters::DBTYPE_HMM_PROFILE, &subMat, 0,
                        false, false);
            buildSyntheticProfile(qW, subMat, std::vector<int>(120, rW));
            ProfileCalibrator::Calib calW = ProfileCalibrator::calibrate(
                    qW.profile_for_alignment, 120, alph, subMat.pBack);
            CHECK(calW.ok, "poly-W calibration ok");
            if (calW.ok) {
                printf("poly-W lambda=%.6f K=%.6f vs mix lambda=%.6f K=%.6f\n",
                       calW.lambda, calW.K, cal.lambda, cal.K);
                CHECK(fabs(calW.K - cal.K) / cal.K > 0.2,
                      "poly-W K differs from background mix");
            }
        }
    }

    // 7. Degenerate inputs never crash and are classified fail-open
    // (degenerate=true), so the prefilter keeps their pairs instead of
    // applying inapplicable statistics. Zero lengths never reach SLS calc().
    {
        ProfileCalibrator::Calib bad1 = ProfileCalibrator::calibrate(
                NULL, 100, alph, subMat.pBack);
        ProfileCalibrator::Calib bad2 = ProfileCalibrator::calibrate(
                qSeq.profile_for_alignment, 0, alph, subMat.pBack);
        std::vector<int8_t> zeros((size_t)20 * 100, 0);
        ProfileCalibrator::Calib bad3 = ProfileCalibrator::calibrate(
                zeros.data(), 100, alph, subMat.pBack);
        CHECK(!bad1.ok && !bad2.ok && !bad3.ok, "degenerate -> !ok");
        CHECK(bad1.degenerate && bad2.degenerate && bad3.degenerate,
              "degenerate inputs -> degenerate (fail-open), not internal");
        CHECK(cal.ok && !cal.degenerate, "normal profile -> ok, not degenerate");
        CHECK(EvalueComputation::computeProfilePvalue(50, 100, 200, bad1.lambda, bad1.K) == 1.0,
              "failed calib -> p=1");
        CHECK(EvalueComputation::computeProfilePvalue(-5, 100, 200, cal.lambda, cal.K) == 1.0,
              "non-positive score -> p=1");
        EvalueComputation seqEval(1000000, &subMat);
        CHECK(seqEval.computePvalue(50.0, 0.0, 200.0) == 1.0,
              "zero query length -> p=1, no throw");
        CHECK(seqEval.computePvalue(50.0, 100.0, 0.0) == 1.0,
              "zero target length -> p=1, no throw");
    }

    // 8. Lattice (divisor > 1) correctness: lib/ncbi-blast carries the
    // single upstream patch (compressed P_scaled[k]=P_orig[k*d], rare: needs
    // all observed scores even).
    {
        NcbiScoreDist dA;
        dA.lo = -2; dA.hi = 2;
        dA.prob.assign(5, 0.0);
        dA.prob[0] = 0.7; dA.prob[4] = 0.3;
        dA.mean = -2*0.7 + 2*0.3;
        NcbiKarlinParams pA = NcbiKarlinUngapped(dA);
        CHECK(pA.ok, "divisor special-case ok");
        if (pA.ok) {
            const double KexpA = (0.7-0.3)*(0.7-0.3)/0.7;
            printf("divisor-A lambda=%.6f K=%.6f (exp %.6f)\n",
                   pA.lambda, pA.K, KexpA);
            CHECK(fabs(pA.K - KexpA)/KexpA < 0.01,
                  "divisor-A K correct (got %.6f exp %.6f)", pA.K, KexpA);
            CHECK(fabs(pA.lambda - 0.423649) < 0.005,
                  "divisor-A lambda (got %.6f)", pA.lambda);
        }
        NcbiScoreDist dB;
        dB.lo = -4; dB.hi = 4;
        dB.prob.assign(9, 0.0);
        dB.prob[0] = 0.4; dB.prob[2] = 0.2;
        dB.prob[6] = 0.2; dB.prob[8] = 0.2;
        dB.mean = -4*0.4 + -2*0.2 + 2*0.2 + 4*0.2;
        NcbiKarlinParams pB = NcbiKarlinUngapped(dB);
        CHECK(pB.ok, "divisor general-case ok");
        if (pB.ok) {
            printf("divisor-B lambda=%.6f K=%.6f H=%.6f\n",
                   pB.lambda, pB.K, pB.H);
            CHECK(pB.K > 0.03 && pB.K < 0.20,
                  "divisor-B K in [0.03,0.20] (got %.6f)", pB.K);
            CHECK(fabs(pB.lambda - 0.147424) < 0.005,
                  "divisor-B lambda (got %.6f)", pB.lambda);
        }
        // End-to-end: even-only int8 profile calibrates (divisor=2).
        std::vector<int8_t> evenProf((size_t)20 * 50, 0);
        for (int i = 0; i < 50; ++i) {
            for (int a = 0; a < 20; ++a) {
                evenProf[(size_t)a * 50 + (size_t)i] =
                        (int8_t)(((a * 7 + i * 3) % 2 == 0) ? 2 : -4);
            }
        }
        ProfileCalibrator::Calib calEven = ProfileCalibrator::calibrate(
                evenProf.data(), 50, alph, subMat.pBack);
        CHECK(calEven.ok, "even-profile calibration ok");
        if (calEven.ok) {
            printf("even-profile lambda=%.6f K=%.6f H=%.6f\n",
                   calEven.lambda, calEven.K, calEven.H);
            CHECK(calEven.lambda > 0 && calEven.K > 0
                  && calEven.K < 1.0, "even-profile sane");
        }
    }

    if (failures == 0) {
        printf("ALL TESTS PASSED\n");
        return 0;
    }
    printf("%d FAILURES\n", failures);
    return EXIT_FAILURE;
}
