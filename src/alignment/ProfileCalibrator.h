#ifndef MMSEQS_PROFILE_CALIBRATOR_H
#define MMSEQS_PROFILE_CALIBRATOR_H

#include <cstddef>
#include <cstdint>
#include <map>
#include <string>

// Per-query ungapped calibration for profile (PSSM) queries.
//
// Single-matrix Karlin-Altschul statistics do not apply to position-specific
// scores. Following PSI-BLAST's method (NCBI C++ toolkit,
// src/algo/blast/composition_adjustment/composition_adjustment.c,
// s_GetPssmScoreProbs: average the PSSM columns' score distributions against
// the subject background into one aggregate distribution), the aggregate
// distribution is fed to the vendored NCBI ungapped Karlin routines
// (lib/ncbi-blast, upstream Blast_KarlinBlkUngappedCalc).
// Deterministic, millisecond-scale, no simulation or RNG involved.
//
// CRITICAL: the aggregate distribution is built from the PRODUCTION integer
// scores (int8 profile_for_alignment), so the resulting (lambda, K) describe
// the actual scorer, including its rounding.
class ProfileCalibrator {
public:
    struct Calib {
        double lambda;      // fitted Gumbel lambda (> 0 when ok)
        double K;           // fitted Gumbel K (> 0 when ok)
        double H;           // relative entropy (> 0 when ok, diagnostic)
        int queryLen;       // query columns m aggregated
        bool ok;            // statistics exist and the fit succeeded
        bool degenerate;    // !ok and the input has no Karlin-Altschul
                            // statistics (empty/flat/one-sided distribution or
                            // non-negative expectation). The prefilter keeps
                            // the query's pairs (fail-open). !ok with
                            // degenerate == false is an internal failure and
                            // must not be masked.
        Calib() : lambda(0), K(0), H(0), queryLen(0), ok(false), degenerate(false) {}
    };

    // Calibrate a profile query. profile holds alphabetSize rows of L int8
    // production scores (row-major: profile[a * L + i]); the X row
    // (index alphabetSize - 1) is ignored (production forces it to 0).
    // pBack is the subject background over 0..alphabetSize-2 (X excluded).
    // Never throws. Sets ok on success, degenerate for inputs without
    // Karlin-Altschul statistics (invalid profile/length, no negative or no
    // positive score, non-negative mean). Any other failure is internal
    // (ok == false, degenerate == false).
    static Calib calibrate(const int8_t *profile, int L, int alphabetSize,
                           const double *pBack);

    // FNV-1a hash over profile_for_alignment bytes (alphabetSize x L).
    // X row (index alphabetSize-1) is skipped: ssw_init forces it to 0.
    static uint64_t hashProfile(const int8_t *profile, int L, int alphabetSize);

    // Cache version. Bump when the file layout or fit semantics change.
    static const int CACHE_VERSION = 3;

    // TSV sidecar: header "#v3", rows:
    //   key L lambda K H
    // key is the profile content hash (reorder-proof across DB rebuilds).
    // Only successful (ok) calibrations are stored. Returns false (never
    // exits) when the file is missing/unreadable; valid rows are still
    // loaded. Entries with a version mismatch are ignored. Saves are written
    // to a unique temporary file and atomically renamed, so concurrent
    // writers cannot interleave within one cache file.
    static bool loadCache(const std::string &path,
                          std::map<uint64_t, Calib> &out);
    static bool saveCache(const std::string &path,
                          const std::map<uint64_t, Calib> &entries);
};

#endif //MMSEQS_PROFILE_CALIBRATOR_H
