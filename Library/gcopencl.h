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

/** @brief @ref gcopencl is just an abbreviation for the struct @ref
                _gcopencl. */
typedef struct _gcopencl gcopencl;

/** @brief Pointer to a @ref gcopencl object. */
typedef gcopencl *pgcopencl;

/** @brief Pointer to a constant @ref gcopencl object. */
typedef const gcopencl *pcgcopencl;

/** @brief An object containing all informations an OpenCL device needs to
  *        calculate the MVM part of the listed clusters. */
struct _gcopencl
{
  /** @brief Number of row cluster this object manages. */
  uint      num_row_leafs;

  /** @brief Name of the row clusters usable in iterating functions like
    *        @p iterate_h2matrix. */
  uint      *names_row_leafs;

  /** @brief Number of @f$\mathcal{H}^2@f$-matrices a row cluster is part of. */
  uint      *num_h2_leafs_per_row;

  /** @brief The maximum number of @f$\mathcal{H}^2@f$-matrices a row cluster
   *         is part of. */
  size_t    max_num_h2_leafs_per_row;

  /** @brief OpenCL buffer objects for the number of
    *        @f$\mathcal{H}^2@f$-matrices a row cluster is part of. */
  cl_mem    *buf_num_h2_leafs_per_row;

  /** @brief Workload of a row cluster, basically the sum of all products of
    *        the dimensions of the farfield matrices the row cluster is a part
    *        of. */
  uint      *workload_per_row;

  /** @brief Index offset to acces nondimensionalized arrays. E.g.
    *        @p host_xtoffs[idx_off[i] + j] is equivalent to xtoffs[i][j]. */
  uint      *idx_off;

  /** @brief OpenCL buffer objects for the index offsets. */
  cl_mem    *buf_idx_off;

  /** @brief The @f$\mathcal{H}^2@f$-matrices a row cluster is part of. */
  ph2matrix **h2_leafs_per_row;

  /** @brief The names of the column cluster of the
    *        @f$\mathcal{H}^2@f$-matrices a row cluster is part of. */
  uint      **col_names_per_row;

  /***/
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

  cl_event  *xt_events;

  uint      *ytoffs;

  cl_mem    *buf_ytoffs;

  cl_mem    *buf_yt;

  cl_event  *yt_events;

  uint roff;

  uint coff;
};

HEADER_PREFIX void
init_gcopencl(pgcopencl gcocl);

HEADER_PREFIX void
uninit_gcopencl(pgcopencl gcocl);

HEADER_PREFIX void
del_gcopencl(pgcopencl gcocls);

HEADER_PREFIX void
del_gcopencls(const uint num_gcocls, pgcopencl *gcocls);

HEADER_PREFIX size_t
getsize_gcopencl(pcgcopencl gcocl);

HEADER_PREFIX void
get_ocl_tree_informations_gcocl(ph2matrix H2,
                                uint      *num_levels,
                                pgcopencl **gcocls);

#endif // GCOPENCL_H