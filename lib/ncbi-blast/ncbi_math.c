/* Vendored subset of src/algo/blast/core/ncbi_math.c
 * (github.com/ncbi/ncbi-cxx-toolkit-public HEAD 024628c, 2026-09-20).
 * Function bodies of BLAST_Expm1 / BLAST_Gcd / BLAST_Powi are verbatim;
 * the rest of ncbi_math.c (Log1p, Gamma, Romberg, Erf, ...) is not needed
 * by the ungapped Karlin path and is intentionally not vendored.
 *
 * PUBLIC DOMAIN (US Government work) — see upstream header.
 */

#include "ncbi_math.h"

#include <math.h>

double BLAST_Expm1(double	x)
{
  double	absx = ABS(x);

  if (absx > .33)
    return exp(x) - 1.;

  if (absx < 1.e-16)
    return x;

  return x * (1. + x *
             (1./2. + x *
             (1./6. + x *
             (1./24. + x *
             (1./120. + x *
             (1./720. + x *
             (1./5040. + x *
             (1./40320. + x *
             (1./362880. + x *
             (1./3628800. + x *
             (1./39916800. + x *
             (1./479001600. +
              x/6227020800.))))))))))));
}

Int4 BLAST_Gcd(Int4 a, Int4 b)
{
   Int4   c;

   b = ABS(b);
   if (b > a)
      c=a, a=b, b=c;

   while (b != 0) {
      c = a%b;
      a = b;
      b = c;
   }
   return a;
}

double BLAST_Powi(double x, Int4 n)
{
   double   y;

   if (n == 0)
      return 1.;

   if (x == 0.) {
      if (n < 0) {
         return HUGE_VAL;
      }
      return 0.;
   }

   if (n < 0) {
      x = 1./x;
      n = -n;
   }

   y = 1.;
   while (n > 0) {
      if (n & 1)
         y *= x;
      n /= 2;
      x *= x;
   }
   return y;
}
