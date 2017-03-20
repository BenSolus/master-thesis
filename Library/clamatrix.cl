/* ------------------------------------------------------------
 * This is the file "clamatrix.cl" of this master thesis.
 * All rights reserved, Steffen Boerm 2009
 * ------------------------------------------------------------ */

/**
 * @file greencross.cl
 * @author Bennet Carstensen
 * @date 2017
 */

#ifndef CLAMATRIX_CL
#define CLAMATRIX_CL

static const char clamatrix_src[] =
{
  "\n"
  "/** @defgroup amatrix amatrix\n"
  " *  @brief Representation of a matrix as an array in column-major order.\n"
  " *\n"
  " *  The @ref amatrix class is used to handle standard linear algebra\n"
  " *  operations like matrix multiplication, factorization, solving\n"
  " *  linear systems or eigenvalue problems.\n"
  " *  @{ */\n"
  "\n"
  "/** @brief Representation of a matrix as a column-order array. */\n"
  "typedef struct _amatrix amatrix;\n"
  "\n"
  "/** @brief Pointer to @ref amatrix object. */\n"
  "typedef amatrix *pamatrix;\n"
  "\n"
  "/** @brief Pointer to constant @ref amatrix object. */\n"
  "typedef const amatrix *pcamatrix;\n"
  "\n"
  "/** @brief Representation of a matrix as an array in column-major order. */\n"
  "struct _amatrix {\n"
  "  /** @brief Matrix coefficients in column-major order, i.e., @f$a_{ij}@f$\n"
  "             corresponds to `a[i+j*ld]`.  */\n"
  "  global field *a;\n"
  "\n"
  "  /** @brief Leading dimension, i.e., increment used to switch from one column\n"
  "             to the next.  */\n"
  "  uint ld;\n"
  "\n"
  "  /** @brief Number of rows. */\n"
  "  uint rows;\n"
  "  /** @brief Number of columns.  */\n"
  "  uint cols;\n"
  "};\n"
  "\n"
  "/** @brief Create a new @ref amatrix object.\n"
  " *\n"
  " *  Sets up the pointer to the matrix coefficients and the other components.\n"
  " *\n"
  " *  @remark @p a must be (sizof(field) * rows * cols) bytes allocated global\n"
  " *          memory.\n"
  " *\n"
  " *  @param rows Number of rows.\n"
  " *  @param cols Number of columns.\n"
  " *  @param a    Pointer to the matrix coefficients.\n"
  " *  @returns New @ref amatrix object. */\n"
  "pamatrix\n"
  "new_amatrix(global field *a, uint rows, uint cols)\n"
  "{\n"
  "  amatrix m;\n"
  "\n"
  "  m.a    = a;\n"
  "  m.ld   = rows;\n"
  "  m.rows = rows;\n"
  "  m.cols = cols;\n"
  "\n"
  "  return &m;\n"
  "}\n"
};

#endif // CLAMATRIX_CL
