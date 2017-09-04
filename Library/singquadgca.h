/* ------------------------------------------------------------
 * This is the file "singquadgca.h" of this master thesis.
 * All rights reserved, Bennet Carstensen 2017
 * ------------------------------------------------------------ */

/**
 * @file      singquadgca.h
 * @author    Bennet Carstensen
 * @date      2017
 * @copyright All rights reserved, Bennet Carstensen 2017
 */

#include "singquad2d.h"

#ifdef __APPLE__
#include "OpenCL/opencl.h"
#else
#include "CL/cl.h"
#endif

/** @brief @ref singquadgca is just an abbreviation for the struct @ref
                _singquadgca. */
typedef struct _singquadgca singquadgca;

/** @brief Pointer to a @ref _singquadgca "singquadgca" object. */
typedef singquadgca *psingquadgca;

/** @brief Pointer to a constant @ref _singquadgca "singquadgca" object. */
typedef const singquadgca *pcsingquadgca;

struct _singquadgca
{
  uint nq;

  real *xqs;
  real *yqs;
  real *wqs;
  real *bases;

  /** @brief OpenCL buffer objects for the first quadratur points. */
  cl_mem *buf_xqs;

  /** @brief OpenCL buffer objects for the second quadratur points. */
  cl_mem *buf_yqs;

  /** @brief OpenCL buffer objects for the quadratur weights. */
  cl_mem *buf_wqs;

  cl_mem *buf_bases;
};

HEADER_PREFIX void
init_singquadgca(psingquadgca sq_gca);

HEADER_PREFIX void
uninit_singquadgca(psingquadgca sq_gca);

HEADER_PREFIX void
del_singquadgca(psingquadgca sq_gca);

HEADER_PREFIX psingquadgca
build_from_singquad2d(psingquad2d sq);

HEADER_PREFIX void
select_quadrature_singquadgca(psingquadgca sq_gca,
                              const uint   *tv,
                              const uint   *sv,
                              uint         *tp,
                              uint         *sp,
                              real         **x,
                              real         **y,
                              real         **w,
                              real         *base);