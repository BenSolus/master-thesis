/* ------------------------------------------------------------
 * This is the file "oclintegralinfos.h" of this master thesis.
 * All rights reserved, Bennet Carstensen 2017
 * ------------------------------------------------------------ */

/**
 * @file      oclintegralinfos.h
 * @author    Bennet Carstensen
 * @date      2017
 * @copyright All rights reserved, Bennet Carstensen 2017
 */

#ifndef OCLINTEGRALINFOS_H
#define OCLINTEGRALINFOS_H

#include "gcopencl.h"
#include "h2matrix.h"

#ifdef __APPLE__
#include "OpenCL/opencl.h"
#else
#include "CL/cl.h"
#endif

/** @brief @ref integralinfos is just an abbreviation for the struct @ref
                _integralinfos. */
typedef struct _integralinfos integralinfos;

/** @brief Pointer to a @ref integralinfos object. */
typedef integralinfos *pintegralinfos;

/** @brief Pointer to a constant @ref gcopencl object. */
typedef const integralinfos *pcintegralinfos;

struct _integralinfos
{
  uint      num_integral_grps;

  uint      *num_integrals;

  cl_mem    *buf_num_integrals;

  uint      *idx_off;

  cl_mem    *buf_idx_off;

  uint      **rows;

  uint      *host_rows;

  cl_mem    *buf_rows;

  uint      **cols;

  uint      *host_cols;

  cl_mem    *buf_cols;

  uint      **cidx;

  uint      *host_cidx;

  cl_mem    *buf_cidx;
};

HEADER_PREFIX void
init_integralinfos(pintegralinfos iinfos);

HEADER_PREFIX void
uninit_integralinfos(pintegralinfos iinfos);

HEADER_PREFIX pintegralinfos
new_integralinfos();

HEADER_PREFIX void
del_integralinfos(pintegralinfos iinfos);

HEADER_PREFIX pintegralinfos
build_from_idxinfos_integralinfos(pcgcopencl idxinfos,
                                  pch2matrix H2,
                                  const uint dim,
                                  void       *poly_idxs,
                                  const uint min_num_id_verts);

HEADER_PREFIX pintegralinfos
build_from_idxinfos(pcgcopencl idxinfos,
                    pch2matrix H2,
                    const uint dim,
                    void       *poly_idxs);

#endif // OCLINTEGRALINFOS_H
