/* ------------------------------------------------------------
 * This is the file "objects.h" of this master thesis.
 * All rights reserved, Bennet Carstensen 2017
 * ------------------------------------------------------------ */

/** @file objects.h
 *  @author Bennet Carstensen
 */

#ifndef OBJECTS_H
#define OBJECTS_H

/** @defgroup objects objects
 *
 *  @brief Representation of objects in d-dimensional space.
 *
 * The @ref objects class represents objects with a location in d-dimensional
 * space and a corresponding mass. Basing on this object the clustergeometry
 * object and the H2-matrix are build.
 *
 * @{ */

/** @brief Representation of objects in d-dimensional space. */
typedef struct _objects objects;

/** @brief Pointer to @ref objects. */
typedef objects* pobjects;

/** @brief Pointer to constant @ref objects.*/
typedef const objects* pcobjects;

#include "basic.h"

/** @brief Representation of objects in d-dimensional space.
 *
 * The @ref objects class represents objects with a location in d-dimensional
 * space and a corresponding mass. Basing on this object the clustergeometry
 * object and the H2-matrix are build.
 */
struct _objects {
  /** @brief Location of the objects. */
  real **x;

  /** @brief Number of objects. */
  uint n;

  /** @brief Dimension of the space containing the objects. */
  uint dim;

  /** @brief Internal field containing all informations about the objects */
  real *mem;
};

/** @brief Initialize the @ref objects.
 *
 * Allocates storage for the objects and sets up pointers.
 *
 * @param objs Objects to initialize.
 */
HEADER_PREFIX void
init_objects(pobjects objs);

/** @brief Uninitialize the @ref objects.
 *
 * Release the storage corresponding to the objects if allocated previously.
 *
 * @param objs Objects to uninitialize.
 */
HEADER_PREFIX void
uninit_objects(pobjects objs);

/** @brief Create new @ref objects.
 *
 *  Allocates storage for the objects and sets up its components.
 *
 *  @remark Should always be matched by a call to @ref del_objects.
 *
 *  @param n   Number of objects
 *  @param dim Dimension of the space containing the objects
 *
 *  @return Returns the newly created @ref objects.
 */
HEADER_PREFIX pobjects
new_objects(uint n, uint dim);

/** @brief Delete the @ref objects.
 *
 *  Releases the storage corresponding to the objects.
 *
 *  @param objs Objects to be deleted.
 */
HEADER_PREFIX void
del_objects(pobjects objs);

/** @} */

#endif // OBJECTS_H
