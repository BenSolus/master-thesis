/* ------------------------------------------------------------
 * This is the file "ocl_system.h" of this master thesis.
 * All rights reserved, Sven Christophersen 2015
 * ------------------------------------------------------------ */

/**
 * @file      ocl_system.h
 * @author    Bennet Carstensen
 * @date      2017
 * @copyright All rights reserved, Sven Christophersen 2015
 */

#ifndef OCL_SYSTEM_H
#define OCL_SYSTEM_H

#include "basic.h"
#include "opencl.h"

/** @defgroup ocl_system ocl_system
 *  @brief Additional functions and fixes for the
 *  <a href="http://www.h2lib.org/doc/d5/de5/group__opencl.html">opencl module of the H2Lib</a>.
 *  @{ */

/**
 * @brief Reads the source code specified by (multiple) strings @p src_strs and
 *        compiles all OpenCL kernels given by the array @p kernel_names into
 *        OpenCL kernels @p kernels.
 *
 * @param n            Number of source code strings.
 * @param src_strs     Source code strings.
 * @param num_kernels  Number of kernels that should be compiled given by @p
 *                     kernel_names.
 * @param kernel_names Array of function names that should be compiled as OpenCL
 *                     kernels.
 * @param kernels      Resulting array of OpenCL kernels.
 */
HEADER_PREFIX void
setup_kernels_fix(const uint n,
                  const char **src_strs,
                  const uint num_kernels,
                  const char **kernel_names,
                  cl_kernel **kernels);

/** @} */

#endif // OCL_SYSTEM_H
