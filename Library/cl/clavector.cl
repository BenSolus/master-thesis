/* ------------------------------------------------------------
 This is the file "clavector.cl" of this master thesis.
 All rights reserved, Steffen Boerm 2009
 ------------------------------------------------------------ */

 /**
  * @file      cl/clavector.cl
  * @author    Bennet Carstensen
  * @date      2017
  * @copyright All rights reserved, Steffen Boerm 2009
  */

#ifndef CLAVECTOR_H
#define CLAVECTOR_H

/** @defgroup avector avector
 *  @brief Representation of a vector as an array in OpenCL.
 *
 *  The @ref avector class is used to handle standard linear algebra
 *  operations like adding, scaling and multiplying vectors.
 *  @{ */

/** Representation of a vector as an array. */
typedef struct _avector avector;

/** Pointer to a @ref avector object. */
typedef avector *pavector;

/** Pointer to a constant @ref avector object. */
typedef const avector *pcavector;

#include "clbasic.cl"

/** Representation of a vector as an array. */
struct _avector {
  /** @brief Vector coefficients. */
  global field *v;

  /** @brief Vector dimension. */
  size_t dim;
};

/** @brief Create a new @ref avector object.
 *
 *  Sets up the pointer to the vector coefficients and the other components.
 *
 *  @param v   Coefficients of the new vector.
 *  @param dim Dimension of the new vector.
 *  @returns New @ref avector object. */
pavector
new_avector(global field *v, size_t dim)
{
  avector a;

  a.v   = v;
  a.dim = dim;

  return &a;
}

/** @brief Print a vector.
 *
 *  @param v Vector @f$v@f$. */
void
print_avector(pcavector v)
{
  const size_t dim = v->dim;

  if((get_local_id(0) == 0) && (get_local_id(1) == 0) && (get_local_id(2) == 0))
  {
    (void) printf("Group (%u, %u, %u): avector(%u)\n", get_group_id(0),
                                                       get_group_id(1),
                                                       get_group_id(2),
                                                       dim);
    if(dim == 0)
      return;

    (void) printf("  (%.5e", v->v[0]);
    for(size_t i = 1; i < dim; i++)
      (void) printf(" %.5e", v->v[i]);
    (void) printf(")\n");
  }
}

/** @} */

#endif // CLAVECTOR_H
