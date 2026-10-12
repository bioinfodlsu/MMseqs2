/* Thin MMseqs adapter: NcbiScoreDist (0-based vector) -> vendor
 * Blast_ScoreFreq (shifted sprob) -> lib/ncbi-blast
 * Blast_KarlinBlkUngappedCalc -> params. The vendored function is used
 * verbatim; the stricter positivity/finiteness validation lives here. All
 * logging lives here; vendor files stay logging-free.
 */

#include "NcbiKarlin.h"
#include "Debug.h"

#include <cmath>

extern "C" {
#include "blast_stat.h"
}

NcbiKarlinParams NcbiKarlinUngapped(const NcbiScoreDist &dist) {
    NcbiKarlinParams out;
    if (dist.lo >= 0 || dist.hi <= 0 || dist.prob.size() !=
            (size_t)(dist.hi - dist.lo + 1)) {
        Debug(Debug::WARNING) << "NcbiKarlin: no ungapped statistics for degenerate score range ["
                            << dist.lo << "," << dist.hi << "]\n";
        return out;
    }
    Blast_ScoreFreq *sfp = Blast_ScoreFreqNew(dist.lo, dist.hi);
    if (sfp == NULL) {
        Debug(Debug::WARNING) << "NcbiKarlin: score range allocation failed for ["
                            << dist.lo << "," << dist.hi << "]\n";
        return out;
    }
    for (int s = dist.lo; s <= dist.hi; ++s) {
        sfp->sprob[s] = dist.prob[(size_t)(s - dist.lo)];
    }
    sfp->obs_min = dist.lo;
    sfp->obs_max = dist.hi;
    sfp->score_avg = dist.mean;
    Blast_KarlinBlk kbp;
    const int rc = Blast_KarlinBlkUngappedCalc(&kbp, sfp);
    Blast_ScoreFreqFree(sfp);
    if (rc != 0 || !(kbp.Lambda > 0.) || !std::isfinite(kbp.Lambda)
            || !(kbp.H > 0.) || !std::isfinite(kbp.H) || !(kbp.K > 0.)
            || !std::isfinite(kbp.K)) {
        // Vendor uses < 0. checks and -1 returns without logging; validate
        // strictly here so callers never store junk (this is an internal
        // failure: the degeneracy checks happened before the call).
        Debug(Debug::WARNING) << "NcbiKarlin: ungapped fit failed for range ["
                            << dist.lo << "," << dist.hi << "] mean=" << dist.mean << "\n";
        return out;
    }
    out.lambda = kbp.Lambda;
    out.K = kbp.K;
    out.H = kbp.H;
    out.ok = true;
    return out;
}
