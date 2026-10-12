/* Minimal NCBI stdint/compat shim for the vendored Karlin subset.
 * Full toolkit: include/corelib/ncbitype.h + src/algo/blast/core/ncbi_std.h
 * (github.com/ncbi/ncbi-cxx-toolkit-public HEAD 024628c, 2026-09-20).
 * Only what blast_stat.c / ncbi_math.c need below is provided.
 */
#ifndef MMSEQS_VENDOR_NCBI_STD_H
#define MMSEQS_VENDOR_NCBI_STD_H

#include <stdint.h>

typedef int32_t Int4;
typedef int16_t Int2;

#define INT2_MIN (-32768)
#define INT2_MAX 32767

#define ABS(a) ((a) >= 0 ? (a) : -(a))

#endif /*MMSEQS_VENDOR_NCBI_STD_H*/
