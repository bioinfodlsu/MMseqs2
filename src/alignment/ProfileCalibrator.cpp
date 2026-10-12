#include "ProfileCalibrator.h"
#include "NcbiKarlin.h"
#include "Debug.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include <vector>

ProfileCalibrator::Calib ProfileCalibrator::calibrate(const int8_t *profile,
        int L, int alphabetSize, const double *pBack) {
    Calib result;
    if (profile == NULL || pBack == NULL || L <= 0 || alphabetSize < 2) {
        // Empty or malformed profile entry: no statistics can exist.
        Debug(Debug::WARNING) << "ProfileCalibrator: empty/invalid profile (L=" << L
                              << "); ungapped p-value gate will keep its pairs\n";
        result.degenerate = true;
        return result;
    }
    result.queryLen = L;
    const int states = alphabetSize - 1;  // exclude X
    // Background over the true residues, renormalised (X excluded, mirroring
    // production targets which carry no X mass in the null).
    double total = 0.0;
    for (int a = 0; a < states; ++a) {
        total += (pBack[a] > 0.0) ? pBack[a] : 0.0;
    }
    if (total <= 0.0) {
        // Broken background is a configuration problem, not a per-profile
        // degeneracy: leave degenerate=false so the prefilter aborts.
        return result;
    }
    // Aggregate score distribution over positions and background residues,
    // following NCBI s_GetPssmScoreProbs: each of the L positions contributes
    // weight 1/L, each residue its background probability. Scores are the
    // production int8 values, so rounding is accounted for exactly.
    int obsMin = 0, obsMax = 0;
    bool first = true;
    for (int i = 0; i < L; ++i) {
        for (int a = 0; a < states; ++a) {
            const int s = (int)profile[(size_t)a * (size_t)L + (size_t)i];
            if (first || s < obsMin) obsMin = s;
            if (first || s > obsMax) obsMax = s;
            first = false;
        }
    }
    std::vector<double> prob((size_t)(obsMax - obsMin + 1), 0.0);
    for (int i = 0; i < L; ++i) {
        for (int a = 0; a < states; ++a) {
            const double w = ((pBack[a] > 0.0) ? pBack[a] : 0.0) / total
                    / (double)L;
            const int s = (int)profile[(size_t)a * (size_t)L + (size_t)i];
            prob[(size_t)(s - obsMin)] += w;
        }
    }
    double mean = 0.0;
    for (int s = obsMin; s <= obsMax; ++s) {
        mean += (double)s * prob[(size_t)(s - obsMin)];
    }
    if (obsMin >= 0 || obsMax <= 0 || !(mean < 0.0)) {
        // No negative or no positive score, or a non-negative expectation:
        // Karlin-Altschul statistics do not exist for this distribution.
        Debug(Debug::WARNING) << "ProfileCalibrator: profile has no Karlin-Altschul statistics"
                              << " (range [" << obsMin << "," << obsMax << "] mean=" << mean
                              << "); ungapped p-value gate will keep its pairs\n";
        result.degenerate = true;
        return result;
    }
    NcbiScoreDist dist;
    dist.lo = obsMin;
    dist.hi = obsMax;
    dist.prob = prob;
    dist.mean = mean;
    const NcbiKarlinParams params = NcbiKarlinUngapped(dist);
    if (params.ok == false) {
        // Degeneracy was ruled out above, so this is an internal failure
        // (allocation, non-convergence, non-finite lambda/K/H). The prefilter
        // aborts rather than masking it.
        return result;
    }
    result.lambda = params.lambda;
    result.K = params.K;
    result.H = params.H;
    result.ok = true;
    return result;
}

uint64_t ProfileCalibrator::hashProfile(const int8_t *profile, int L,
                                        int alphabetSize) {
    // FNV-1a 64 over the 20 residue rows; the X row is ssw_init-forced to 0.
    uint64_t h = 1469598103934665603ULL;
    if (profile == NULL || L <= 0 || alphabetSize < 2) {
        return h;
    }
    const int rows = alphabetSize - 1;
    for (int a = 0; a < rows; ++a) {
        const int8_t *row = profile + (size_t)a * (size_t)L;
        for (int i = 0; i < L; ++i) {
            h ^= (uint64_t)(uint8_t)row[i];
            h *= 1099511628211ULL;
        }
    }
    h ^= (uint64_t)L;
    h *= 1099511628211ULL;
    return h;
}

bool ProfileCalibrator::loadCache(const std::string &path,
                                  std::map<uint64_t, Calib> &out) {
    if (path.empty()) {
        return false;
    }
    FILE *f = fopen(path.c_str(), "r");
    if (f == NULL) {
        return false;
    }
    char line[1024];
    char expectedHeader[16];
    snprintf(expectedHeader, sizeof(expectedHeader), "#v%d", CACHE_VERSION);
    const size_t headerLen = strlen(expectedHeader);
    bool headerSeen = false;
    while (fgets(line, sizeof(line), f) != NULL) {
        if (line[0] == '#') {
            if (strncmp(line, expectedHeader, headerLen) == 0) {
                headerSeen = true;
            }
            continue;
        }
        if (headerSeen == false) {
            continue;  // version mismatch: ignore body
        }
        unsigned long long key = 0;
        int L = 0;
        double lambda = 0, K = 0, H = 0;
        const int cols = sscanf(line, "%llu %d %lf %lf %lf",
                &key, &L, &lambda, &K, &H);
        if (cols != 5 || L <= 0 || lambda <= 0.0 || K <= 0.0) {
            continue;
        }
        Calib c;
        c.lambda = lambda;
        c.K = K;
        c.H = H;
        c.queryLen = L;
        c.ok = true;
        out[(uint64_t)key] = c;
    }
    fclose(f);
    return headerSeen;
}

bool ProfileCalibrator::saveCache(const std::string &path,
                                  const std::map<uint64_t, Calib> &entries) {
    if (path.empty() || entries.empty()) {
        return false;
    }
    // Unique temporary file per writer + atomic rename: concurrent writers to
    // the same cache path cannot interleave within one file.
    std::string tmpl = path + ".tmp.XXXXXX";
    std::vector<char> tmpBuf(tmpl.begin(), tmpl.end());
    tmpBuf.push_back('\0');
    const int fd = mkstemp(tmpBuf.data());
    if (fd < 0) {
        return false;
    }
    FILE *f = fdopen(fd, "w");
    if (f == NULL) {
        close(fd);
        unlink(tmpBuf.data());
        return false;
    }
    fprintf(f, "#v%d key L lambda K H\n", CACHE_VERSION);
    for (std::map<uint64_t, Calib>::const_iterator it = entries.begin();
            it != entries.end(); ++it) {
        const Calib &c = it->second;
        if (c.ok == false) {
            continue;
        }
        fprintf(f, "%llu %d %.17g %.17g %.17g\n",
                (unsigned long long)it->first, c.queryLen,
                c.lambda, c.K, c.H);
    }
    if (fclose(f) != 0) {
        unlink(tmpBuf.data());
        return false;
    }
    if (rename(tmpBuf.data(), path.c_str()) != 0) {
        unlink(tmpBuf.data());
        return false;
    }
    return true;
}
