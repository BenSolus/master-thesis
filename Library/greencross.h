/* ------------------------------------------------------------
 * This is the file "greencross.h" of this master thesis.
 * All rights reserved, Bennet Carstensen 2017
 * ------------------------------------------------------------ */

/** @file   greencross.h
 *  @author Bennet Carstensen
 */

#ifndef GREENCROSS_H
#define GREENCROSS_H

/** @defgroup greencross greencross
 *  @brief Implementation of the green cross approximation.
 *
 *  @{ */

typedef struct _greencross_ocl greencross_ocl;

typedef greencross_ocl *pgreencross_ocl;

typedef const greencross_ocl *pcgreencross_ocl;

#include "clusterbasis.h"
#include "h2matrix.h"

#include "CL/cl.h"

struct _greencross_ocl
{
  cl_uint num_kernels;

  cl_kernel *kernels;

  ph2matrix H2;

  cl_mem *mem_coefs;
};

HEADER_PREFIX void
init_greencross_ocl();

HEADER_PREFIX void
uninit_greencross_ocl();

HEADER_PREFIX void
forward_greencross(pclusterbasis cb, pcavector x, pavector xt);

HEADER_PREFIX void
eval_leafs_forward(pclusterbasis cb,
                   pcavector     x,
                   pavector      xt);

/** @} */

#endif //GREENCROSS_H
