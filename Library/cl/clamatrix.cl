/* ------------------------------------------------------------
 * This is the file "clamatrix.cl" of this master thesis.
 * All rights reserved, Steffen Boerm 2009
 * ------------------------------------------------------------ */

/**
 * @file   clamatrix.cl
 * @author Bennet Carstensen
 * @date   2017
 */

#ifndef CLAMATRIX_CL
#define CLAMATRIX_CL

#include "clbasic.cl"

/** @defgroup amatrix amatrix
 *  @brief Representation of a matrix as an array in column-major order.
 *
 *  The @ref amatrix class is used to handle standard linear algebra
 *  operations like matrix multiplication, factorization, solving
 *  linear systems or eigenvalue problems.
 *  @{ */

/** @brief Representation of a matrix as a column-order array. */
typedef struct _amatrix amatrix;

/** @brief Pointer to @ref amatrix object. */
typedef amatrix *pamatrix;

/** @brief Pointer to constant @ref amatrix object. */
typedef const amatrix *pcamatrix;

/** @brief Representation of a matrix as an array in column-major order. */
struct _amatrix {
  /** @brief Matrix coefficients in column-major order, i.e., @f$a_{ij}@f$
             corresponds to `a[i+j*ld]`.  */
  global field *a;

  /** @brief Leading dimension, i.e., increment used to switch from one column
             to the next.  */
  uint ld;

  /** @brief Number of rows. */
  uint rows;
  /** @brief Number of columns.  */
  uint cols;
};

/** @brief Create a new @ref amatrix object.
 *
 *  Sets up the pointer to the matrix coefficients and the other components.
 *
 *  @remark @p a must be (sizof(field) * rows * cols) bytes allocated global
 *          memory.
 *
 *  @param rows Number of rows.
 *  @param cols Number of columns.
 *  @param a    Pointer to the matrix coefficients.
 *  @returns New @ref amatrix object. */
pamatrix
new_amatrix(global field *a, uint rows, uint cols)
{
  amatrix m;

  m.a    = a;
  m.ld   = rows;
  m.rows = rows;
  m.cols = cols;

  return &m;
}

/** @} */

#endif // CLAMATRIX_CL
