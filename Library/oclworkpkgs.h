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

  uint   *num_rows_per_pkg;

  uint   **rows_per_pkg;

  cl_mem *buf_rows_this_device;

  uint   offset;

  uint   size;
};

void
init_oclwork(poclworkpgs oclwrk);

void
uninit_oclwork(poclworkpgs oclwrk);

/** @} */

#endif // OCLWORKPKGS_H
