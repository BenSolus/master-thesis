/* ------------------------------------------------------------
 * This is the file "kernels.c" of this master thesis.
 * All rights reserved, Bennet Carstensen 2017
 * ------------------------------------------------------------ */

/**
 * @file      kernels.c
 * @author    Bennet Carstensen
 * @date      2017
 * @copyright All rights reserved, Bennet Carstensen 2017
 */

#include "kernels.h"

/* ------------------------------------------------------------
 * Kernel functions
 * ------------------------------------------------------------ */

real
laplace_kernel(const field* x, const field* y, const uint dim)
{
 assert((dim == 2) || (dim == 3));

 real norm = dim == 2 ? REAL_NORM2(x[0] - y[0], x[1] - y[1])
                      : REAL_NORM3(x[0] - y[0], x[1] - y[1], x[2] - y[2]);

 if(norm < 1e-15)
   return r_zero;

 return dim == 2 ? r_minus_two_pi * REAL_LOG(norm) : r_four_pi / norm;
}

/* ------------------------------------------------------------
 * Partial derivatives of kernel functions
 * ------------------------------------------------------------ */

real
pdx_laplace_kernel(const field *x, const field *y, const uint dim, const uint i)
{
  assert(i < dim);

  real rnorm = dim == 2 ? REAL_NORMSQR2(x[0] - y[0], x[1] - y[1])
                        : REAL_NORMSQR3(x[0] - y[0], x[1] - y[1], x[2] - y[2]);

  if(rnorm < 1e-16)
    return r_zero;

  rnorm = dim == 2 ? r_one / rnorm : REAL_RSQRT(rnorm * rnorm * rnorm);

  return dim == 2 ? (x[i] - y[i]) * r_minus_two_pi * rnorm
                  : - (x[i] - y[i]) * r_four_pi * rnorm;
}

real
pdy_laplace_kernel(const field *x, const field *y, const uint dim, const uint i)
{
  return -pdx_laplace_kernel(x, y, dim, i);
}
