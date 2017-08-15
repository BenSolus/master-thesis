/* ------------------------------------------------------------
 * This is the file "gcopencl.h" of this master thesis.
 * All rights reserved, Bennet Carstensen 2017
 * ------------------------------------------------------------ */

/**
 * @file      gcopencl.h
 * @author    Bennet Carstensen
 * @date      2017
 * @copyright All rights reserved, Bennet Carstensen 2017
 */

#ifndef GCOPENCL_H
#define GCOPENCL_H

#include "h2matrix.h"

#ifdef __APPLE__
#include "OpenCL/opencl.h"
#else
#include "CL/cl.h"
#endif

static const size_t num_kernels = 1;

/** @brief @ref gcopencl is just an abbreviation for the struct @ref
                _gcopencl. */
typedef struct _gcopencl gcopencl;

/** @brief Pointer to a @ref gcopencl object. */
typedef gcopencl *pgcopencl;

/** @brief Pointer to a constant @ref gcopencl object. */
typedef const gcopencl *pcgcopencl;

struct _gcopencl
{
  uint      num_row_leafs;

  uint      *names_row_leafs;

  uint      *num_h2_leafs_per_row;

  cl_mem    *buf_num_h2_leafs_per_row;

  uint      *workload_per_row;

  uint      *idx_off;

  cl_mem    *buf_idx_off;

  ph2matrix **h2_leafs_per_row;

  uint      **col_names_per_row;

  uint      size_ridx;

  uint      *ridx_sizes;

  cl_mem    *buf_ridx_sizes;

  uint      *ridx_off;

  cl_mem    *buf_ridx_off;

  uint      *host_ridx;

  cl_mem    *buf_ridx;

  uint      size_cidx;

  uint      *cidx_sizes;

  cl_mem    *buf_cidx_sizes;

  uint      **cidx_off;

  uint      *host_cidx_off;

  cl_mem    *buf_cidx_off;

  uint      *host_cidx;

  cl_mem    *buf_cidx;

  uint      **xtoffs;

  uint      *host_xtoffs;

  cl_mem    *buf_xtoffs;

  cl_mem    *buf_xt;

  uint      *ytoffs;

  cl_mem    *buf_ytoffs;

  cl_mem    *buf_yt;

  cl_kernel *kernels;

  pcluster  *row_clusters;

  uint      *cur_h2matrix_leafs;

  uint roff;
  uint coff;
};

HEADER_PREFIX void
init_gcopencl(pgcopencl gcocl);

HEADER_PREFIX void
uninit_gcopencl(pgcopencl gcocl);

#endif // GCOPENCL_H