/* ------------------------------------------------------------
 * This is the file "fastaddevalgca.h" of this master thesis.
 * All rights reserved, Bennet Carstensen 2017
 * ------------------------------------------------------------ */

/**
 * @file      fastaddevalgca.h
 * @author    Bennet Carstensen
 * @date      2017
 * @copyright All rights reserved, Bennet Carstensen 2017
 */

#ifndef FASTADDEVALGCA_H
#define FASTADDEVALGCA_H

#include "clusterbasis.h"

#ifdef __APPLE__
#include "OpenCL/opencl.h"
#else
#include "CL/cl.h"
#endif

/** @brief @ref fastaddevalgca is just an abbreviation for the struct @ref
                _fastaddevalgca. */
typedef struct _fastaddevalgca fastaddevalgca;

/** @brief Pointer to a @ref _fastaddevalgca "fastaddevalgca" object. */
typedef fastaddevalgca *pfastaddevalgca;

/** @brief Pointer to a constant @ref _fastaddevalgca "fastaddevalgca"
 * object. */
typedef const fastaddevalgca *pcfastaddevalgca;

struct _fastaddevalgca
{
  uint     num_wrk_pkgs;

  cl_mem   *buf_xt;

  cl_event *events_xt;

//#ifndef USE_OPENMP
  cl_mem   *buf_yt;
//#else
//  cl_mem   *buf_yt_ff;
//
//  cl_mem   *buf_yt_nf_common;
//
//  cl_mem   *buf_yt_nf_min_vert;
//
//  cl_mem   *buf_yt_nf_min_edge;
//#endif

  cl_event *events_yt;
};

HEADER_PREFIX void
init_fastaddevalgca(pfastaddevalgca feval);

HEADER_PREFIX void
uninit_fastaddevalgca(pfastaddevalgca feval);

HEADER_PREFIX pfastaddevalgca
new_fastaddevalgca(pcclusterbasis rb,
                   pcclusterbasis cb,
                   const uint     num_wrk_pkgs);

HEADER_PREFIX void
del_fastaddevalgca(pfastaddevalgca feval);

#endif // FASTADDEVALGCA_H
