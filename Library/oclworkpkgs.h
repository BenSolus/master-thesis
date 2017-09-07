/* ------------------------------------------------------------
 * This is the file "oclworkpkgs.h" of this master thesis.
 * All rights reserved, Bennet Carstensen 2017
 * ------------------------------------------------------------ */

/**
 * @file      oclworkpkgs.h
 * @author    Bennet Carstensen
 * @date      2017
 * @copyright All rights reserved, Bennet Carstensen 2017
 */

#ifndef OCLWORKPKGS_H
#define OCLWORKPKGS_H

#include "basic.h"
#include "gcopencl.h"

#ifdef __APPLE__
#include "OpenCL/opencl.h"
#else
#include "CL/cl.h"
#endif

/** @defgroup oclworkpkgs oclworkpkgs
 *  @brief
 *
 * @{*/

/** @brief @ref gcopencl is just an abbreviation for the struct @ref
                _gcopencl. */
typedef struct _oclworkpgs oclworkpgs;

/** @brief Pointer to a @ref gcopencl object. */
typedef oclworkpgs *poclworkpgs;

/** @brief Pointer to a constant @ref gcopencl object. */
typedef const oclworkpgs *pcoclworkpgs;

struct _oclworkpgs
{
  uint   num_wrk_pkgs;

  uint   *wrk_per_pkg;

  uint   *first_idx_of_pkgs;

  uint   *last_idx_of_pkg;

  uint   *num_rows_per_pkg;

  uint   **rows_per_pkg;

  cl_mem *buf_rows_this_device;
};

HEADER_PREFIX void
init_oclwork(poclworkpgs oclwrk);

HEADER_PREFIX void
uninit_oclwork(poclworkpgs oclwrk);

HEADER_PREFIX poclworkpgs
new_equidistant_distributed_oclwork(pcgcopencl gcocl,
                                    ph2matrix  H2,
                                    const uint num_packages);

HEADER_PREFIX void
del_oclwrk(poclworkpgs oclwrk);

/** @} */

#endif // OCLWORKPKGS_H
