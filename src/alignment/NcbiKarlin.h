/*   ==========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ==========================================================================*/

/** @file NcbiKarlin.h
 * Thin MMseqs2 adapter over the vendored upstream statistics in
 * lib/ncbi-blast/ (original filenames blast_stat.c / ncbi_math.c,
 * github.com/ncbi/ncbi-cxx-toolkit-public HEAD 024628c, 2026-09-20).
 *
 * The adapter maps NcbiScoreDist/NcbiKarlinParams onto the vendor types and
 * calls upstream Blast_KarlinBlkUngappedCalc; all mathematics lives in the
 * vendor files with shifted-pointer indexing. The single upstream patch
 * (lattice divisor > 1 compression, rare: needs all observed scores even) is
 * documented in lib/ncbi-blast/README.md. Failures log here (never silent);
 * the vendor itself stays logging-free. Callers distinguish "no statistics
 * exist for this distribution" (degenerate input, handled fail-open by the
 * prefilter) from an unexpected internal fit failure (aborted by the
 * prefilter); this header only reports ok / !ok.
 */

#ifndef MMSEQS_NCBI_KARLIN_H
#define MMSEQS_NCBI_KARLIN_H

#include <vector>

// One-column score distribution: prob[(s) - lo] for s in [lo, hi].
struct NcbiScoreDist {
    int lo;                    // lowest score with nonzero probability (< 0)
    int hi;                    // highest score with nonzero probability (> 0)
    std::vector<double> prob;  // size hi - lo + 1, sums to 1
    double mean;               // expected score (must be < 0)
};

struct NcbiKarlinParams {
    double lambda;  // > 0 on success
    double K;       // > 0 on success
    double H;       // relative entropy, > 0 on success
    bool ok;        // false when statistics do not exist for this dist
    NcbiKarlinParams() : lambda(0), K(0), H(0), ok(false) {}
};

// lambda/K/H for the ungapped local-alignment score maximum under dist.
// Returns ok=false (never throws/aborts) when the expectation is
// non-negative, the range is degenerate, or the iteration fails.
NcbiKarlinParams NcbiKarlinUngapped(const NcbiScoreDist &dist);

#endif //MMSEQS_NCBI_KARLIN_H
