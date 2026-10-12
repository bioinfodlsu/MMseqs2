# ncbi-blast (vendored from the NCBI C++ toolkit)

> **AI-generated.** This document was written by AI tools; the vendored code and
> its integration have been manually reviewed and tested.

MMseqs2's profile prefilter needs ungapped Karlin-Altschul parameters
(`lambda`, `K`, `H`) for a per-query score distribution. Rather than
reimplementing the statistics, this directory vendors the exact upstream
functions needed.

## Provenance

* Upstream: `https://github.com/ncbi/ncbi-cxx-toolkit-public`
* Pin: `024628c277cbd1a996e2627d651fd8ac9312300b` (2026-09-20)
* Sources: `src/algo/blast/core/blast_stat.c`,
  `src/algo/blast/core/ncbi_math.c`,
  `include/algo/blast/core/blast_stat.h`
* Extraction: manual subset of the functions and constants below, taken at the
  pin above.

## What is vendored

| file | content | upstream lines at the pin |
| --- | --- | --- |
| `blast_stat.c` | `BlastScoreChk`, `Blast_ScoreFreqNew/Free`, `BlastKarlinLHtoK`, `NlmKarlinLambdaNR`, `Blast_KarlinLambdaNR`, `BlastKarlinLtoH`, `Blast_KarlinBlkUngappedCalc`, `BLAST_SCORE_RANGE_MAX`, `BLAST_KARLIN_*_DEFAULT` | 56-72, 940-950, 2100-2141, 2247-2420, 2491-2599, 2607-2634, 2699-2735 |
| `ncbi_math.c` | `BLAST_Expm1`, `BLAST_Gcd`, `BLAST_Powi` | 33-56, 405-419, 444-471 |
| `blast_stat.h` | `BLAST_SCORE_MIN/MAX`, `Blast_KarlinBlk`, `Blast_ScoreFreq`, `Blast_KarlinBlkUngappedCalc` prototype | 63-72, 121-136, 843-852 |
| `ncbi_std.h` | shim: `Int4`/`Int2`, `INT2_MIN/MAX`, `ABS` (stands in for upstream toolkit headers such as `ncbitype.h`) | — |
| `CMakeLists.txt` | object library `ncbi_blast` | — |

Deliberately excluded (not needed by the ungapped path): the score tables,
gapped statistics, PSI-BLAST and matrix readers, and the `boost_erf.h` /
`blast_psi_priv.h` includes of the upstream `.c` file.

## Design rules

1. **Function bodies are verbatim upstream**, including shifted-pointer indexing
   (`sprob -= score_min`, `sprob[low..high]`) and the `Int4`/`Int2` types. Keep
   it that way when re-syncing: a diff against a fresh extraction must show only
   the one marked patch below.
2. **No MMseqs2 logic in these files.** The adapter is
   `src/alignment/NcbiKarlin.{h,cpp}` (builds `Blast_ScoreFreq`, calls
   `Blast_KarlinBlkUngappedCalc`, validates `lambda`/`K`/`H`); the distribution
   builder and cache is `src/alignment/ProfileCalibrator.{h,cpp}`.
3. The upstream `NCBI_XBLAST_EXPORT` attribute is omitted (static object
   library), and `ncbi_std.h` replaces upstream toolkit includes instead of
   copying them.

## The single modification: `MMSEQS PATCH` in `BlastKarlinLHtoK`

`BlastKarlinLHtoK` computes `K` with a dynamic program over the score
probabilities. Before the DP it finds the greatest common divisor `divisor` of
the observed score differences, divides the scores by it and scales `lambda`, so
the DP runs on the unit lattice. Upstream then keeps the original probability
array but indexes it with the *scaled* offsets (`probArrayStartLow[first]`, where
`first` is a lattice index). That is correct for `divisor == 1`, but for
`divisor > 1` it reads `P_orig[first]` instead of `P_orig[first * divisor]`, so
the DP convolves with the wrong distribution and `K` comes out wrong.

Worked example: an even-only distribution on `-4..4` (`divisor = 2`) gives
upstream `K = 3.29`; the correct value is `K ~= 0.086` (Monte-Carlo `~0.08`).

The vendored copy compresses the array onto the unit lattice first:

    P_scaled[k] = P_orig[k * divisor],   k = 0 .. high - low

This is lossless: by construction of the divisor, every off-lattice entry is
zero. The largest index read is `obs_max`, and the scratch array is freed on
every return path. The `low == -1 && high == 1` special case is
index-equivalent to upstream; `divisor == 1` is unchanged; the
`ADD_GEOMETRIC_TERMS_TO_K` block is verbatim and never defined. The patch is
exercised by `TestProfileCalib` (lattice `divisor-A`, `divisor-B` and the
even-only profile). It looks like a genuine upstream bug for `divisor > 1`;
reporting it upstream is worthwhile.

## Build integration

`CMakeLists.txt` builds the `ncbi_blast` object library as C with warnings
suppressed; `src/CMakeLists.txt` adds its objects to `mmseqs-framework` and
`lib/ncbi-blast` is on the include path (top-level `CMakeLists.txt`).

## Re-sync procedure

1. Clone `https://github.com/ncbi/ncbi-cxx-toolkit-public` and check out the new
   pin.
2. Re-extract the functions and constants listed under "What is vendored" from
   `blast_stat.c`, `ncbi_math.c` and `blast_stat.h` (keep upstream names and
   bodies).
3. `diff -u` the fresh extraction against this directory. The only expected
   difference is the `MMSEQS PATCH` hunk in `BlastKarlinLHtoK`, plus the
   provenance headers and this README's pin/line table. Anything else must be
   reverted or documented here.
4. Re-run `test_profilecalib` and `TestTinyProfileCpuFilter.sh`.
5. Update the pin, date and line numbers here.

## License

PUBLIC DOMAIN (US Government work) — see the file headers (NCBI).
