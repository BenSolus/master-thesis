/* ------------------------------------------------------------
 * This is the file "clamatrix.cl" of this master thesis.
 * All rights reserved, Steffen Boerm 2009
 * ------------------------------------------------------------ */

/**
 * @file      cl/clamatrix.cl
 * @author    Bennet Carstensen
 * @date      2017
 * @copyright All rights reserved, Steffen Boerm 2009
 */

#ifndef CLAMATRIX_CL
#define CLAMATRIX_CL

#include "clbasic.cl"
#include "clavector.cl"

/** @defgroup amatrix amatrix
 *  @brief Representation of a matrix as an array in column-major order in
 *         OpenCL.
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
    *        corresponds to `a[i+j*ld]`.  */
  global field *a;

  /** @brief Leading dimension, i.e., increment used to switch from one column
    *        to the next.  */
  size_t ld;

  /** @brief Number of rows. */
  size_t rows;
  /** @brief Number of columns.  */
  size_t cols;
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
new_amatrix(global field *a, size_t rows, size_t cols)
{
  amatrix m;

  m.a    = a;
  m.ld   = rows;
  m.rows = rows;
  m.cols = cols;

  return &m;
}

/** @brief Print a matrix.
 *
 *  @param a Matrix object. */
void
print_amatrix(pcamatrix a)
{
  size_t rows = a->rows;
  size_t cols = a->cols;
  size_t lda = a->ld;

  if((get_local_id(0) == 0) && (get_local_id(1) == 0) && (get_local_id(2) == 0))
  {
    (void) printf("amatrix(%u,%u,%u)\n", rows, cols, a->ld);
    if (rows == 0 || cols == 0)
      return;

    for(size_t i = 0; i < rows; i++) {
      (void) printf("  (%.5e", a->a[i]);
      for(size_t j = 1; j < cols; j++)
        (void) printf(" | %.5e", a->a[i + j * lda]);
      (void) printf(")\n");
    }
  }
}

/** @brief Multiply a matrix @f$A@f$ by a vector @f$x@f$,
 *  @f$y \gets y + \alpha A x@f$.
 *
 *  The matrix @f$A@f$ is multiplied by the source vector @f$x@f$,
 *  the result is scaled by @f$\alpha@f$ and added to the
 *  target vector @f$y@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param a Matrix @f$A@f$.
 *  @param src Source vector @f$x@f$.
 *  @param trg Target vector @f$y@f$. */
void
addeval_amatrix_avector(field alpha, pcamatrix a, pcavector src, pavector trg)
{
  const size_t lcl_size = get_local_size(0);
  const size_t lcl_id   = get_local_id(0);

  const size_t rows = a->rows;
  const size_t cols = a->cols;

  /* Number of iterations, where ALL threads perform. */
  const size_t n    = (size_t) floor((real) rows / (real) lcl_size);

  /* Calculate lcl_size entries of the result vector at once. */
  for(size_t k = 0; k < n; ++k)
  {
    const size_t i = lcl_size * k + lcl_id; // Row index

    field y = f_zero;

    for(size_t j = 0; j < cols; ++j)
      y += a->a[i + j * rows] * src->v[j];

    trg->v[i] += alpha * y;
  }

  /* Calculate the remaining entries of the result vector. */
  if(lcl_id < (rows - n * lcl_size))
  {
    const size_t i = lcl_size * n + lcl_id;

    field y = f_zero;

    for(size_t j = 0; j < cols; ++j)
      y += a->a[i + j * rows] * src->v[j];

    trg->v[i] += alpha * y;
  }
}

/** @brief Multiply the adjoint of a matrix @f$A@f$ by a vector @f$x@f$,
 *  @f$y \gets y + \alpha A^* x@f$.
 *
 *  The adjoint @f$A^*@f$ is multiplied by the source vector @f$x@f$,
 *  the result is scaled by @f$\alpha@f$ and added to the
 *  target vector @f$y@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param a Matrix @f$A@f$.
 *  @param src Source vector @f$x@f$.
 *  @param trg Target vector @f$y@f$. */
void
addevaltrans_amatrix_avector(field     alpha,
                             pcamatrix a,
                             pcavector src,
                             pavector trg)
{
  const size_t lcl_size = get_local_size(0);
  const size_t lcl_id   = get_local_id(0);

  const size_t rows = a->rows;
  const size_t cols = a->cols;

  /* Number of iterations, where ALL threads perform. */
  const size_t n    = (size_t) floor((real) cols / (real) lcl_size);

  /* Calculate lcl_size entries of the result vector at once. */
  for(size_t k = 0; k < n; ++k)
  {
    const size_t j = lcl_size * k + lcl_id; // Row index

    field y = f_zero;

    for(size_t i = 0; i < rows; ++i)
      y += a->a[i + j * rows] * src->v[i];

    trg->v[j] += alpha * y;
  }

  /* Calculate the remaining entries of the result vector. */
  if(lcl_id < (cols - n * lcl_size))
  {
    const size_t j = lcl_size * n + lcl_id;

    field y = f_zero;

    for(size_t i = 0; i < rows; ++i)
      y += a->a[i + j * rows] * src->v[i];

    trg->v[j] += alpha * y;
  }
}

/** @brief Multiply a matrix @f$A@f$ or its adjoint @f$A^*@f$ by a
 *  vector, @f$y \gets y + \alpha A x@f$ or @f$y \gets y + \alpha A^* x@f$.
 *
 *  The matrix or its adjoint is multiplied by the source vector @f$x@f$,
 *  the result is scaled by @f$\alpha@f$ and added to the target vector
 *  @f$y@f$.
 *
 *  @param alpha Scaling factor @f$\alpha@f$.
 *  @param atrans Set if @f$A^*@f$ is to be used instead of @f$A@f$.
 *  @param a Matrix @f$A@f$.
 *  @param src Source vector @f$x@f$.
 *  @param trg Target vector @f$y@f$. */
void
mvm_amatrix_avector(field     alpha,
                    bool      atrans,
                    pcamatrix a,
                    pcavector src,
                    pavector  trg)
{
  if(atrans)
    addevaltrans_amatrix_avector(alpha, a, src, trg);
  else
    addeval_amatrix_avector(alpha, a, src, trg);
}

// kernel void
// mvm_amatrices_avectors(global field *alphas,
//                        global bool  *atrans,
//                        global field *as,
//                        global uint  *as_off,
//                        global uint  *n_rows,
//                        global uint  *n_cols,
//                        global field *srcs,
//                        global uint  *srcs_off,
//                        global field *trgs,
//                        global uint  *trgs_off)
// {
//   const size_t grp_id = get_group_id(0);
//
//   const size_t rows   = n_rows[grp_id];
//   const size_t cols   = n_cols[grp_id];
//
//   const bool   trans  = atrans[grp_id];
//
//   pamatrix a   = new_amatrix(as + as_off[grp_id], rows, cols);
//   pavector src = new_avector(srcs + srcs_off[grp_id], trans ? rows : cols);
//   pavector trg = new_avector(trgs + trgs_off[grp_id], trans ? cols : rows);
//
//   if(get_group_id(0) == 0)
//     printf("%u\n", trgs_off[grp_id]);
//
//   mvm_amatrix_avector(1.0, trans, a, src, trg);
// }

/** @} */

#endif // CLAMATRIX_CL
