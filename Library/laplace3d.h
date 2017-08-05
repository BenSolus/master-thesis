/* ------------------------------------------------------------
 * This is the file "laplace3d.h" of this master thesis.
 * All rights reserved, Bennet Carstensen 2017
 * ------------------------------------------------------------ */

/**
 * @file      laplace3d.h
 * @author    Bennet Carstensen
 * @date      2017
 * @copyright All rights reserved, Bennet Carstensen 2017
 */

#ifndef LAPLACE3D_H
#define LAPLACE3D_H

#include "greencross.h"

#include "surface3d.h"

/** @defgroup laplace2d laplace2d
 *  @brief Algorithms and functions to perform the Green cross approximation
 *         method on the fundamental solution of the negative Laplace operator
 *         in 3D.
 *
 * @{*/

/** @brief Prepares a @ref _greencross "greencross" object for approximating the
 *         integral equation
 *
 * Given a subspace @f$\Omega \subset \mathbb{R}^3 @f$ by @p s3d, an object
 * to approximate the integral equation
 * is given with an approximation order @p m, a quadrature order @p q,
 * necessary data for the admissibility condition in the hierarchical
 * clustering @p eta and a maximal size of leaf clusters @p res.
 *
 * @param[in] s3d Geometry describing the subspace
 *                @f$\Omega \subset \mathbb{R}^3 @f$.
 * @param[in] res Maximal size of leaf clusters.
 * @param[in] q   Quadratur order.
 * @param[in] m   Approximation order.
 * @result        A @ref _greencross "greencross" object ready for
 *                approximating the integral equation
 */
HEADER_PREFIX pgreencross
new_greencross_laplace3d(psurface3d s3d, uint res, uint q, uint m, real accur);

/** @} */

#endif // LAPLACE3D_H
