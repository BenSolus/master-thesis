/* ------------------------------------------------------------
 * This is the file "laplace2d.h" of this master thesis.
 * All rights reserved, Bennet Carstensen 2017
 * ------------------------------------------------------------ */

/**
 * @file      laplace2d.h
 * @author    Bennet Carstensen
 * @date      2017
 * @copyright All rights reserved, Bennet Carstensen 2017
 */

#ifndef LAPLACE2D_H
#define LAPLACE2D_H

#include "greencross.h"

#include "curve2d.h"

/** @defgroup laplace2d laplace2d
 *  @brief Algorithms and functions to perform the Green cross approximation
 *         method on the fundamental solution of the negative Laplace operator
 *         in 2D.
 *
 * @{*/

/** @brief Prepares a @ref _greencross "greencross" object for approximating the
 *         integral equation @f$ \int\limits_{\Omega} -\frac{1}{2 \pi}
 *         \log{\left \Vert x - y \right \Vert_2} u(y) dy = f(x) @f$, where
 *         @f$\Omega \subset \mathbb{R}^2 @f$.
 *
 * Given a subspace @f$\Omega \subset \mathbb{R}^2 @f$ by @p c2d, an object
 * to approximate the integral equation @f$ \int\limits_{\Omega}
 * -\frac{1}{2 \pi} \log{\left \Vert x - y \right \Vert_2} u(y) dy = f(x) @f$
 * is given with an approximation order @p m, a quadrature order @p q,
 * necessary data for the admissibility condition in the hierarchical
 * clustering @p eta and a maximal size of leaf clusters @p res.
 *
 * @param[in] c2d Geometry describing the subspace
 *                @f$\Omega \subset \mathbb{R}^2 @f$.
 * @param[in] res Maximal size of leaf clusters.
 * @param[in] q   Quadratur order.
 * @param[in] m   Approximation order.
 * @result        A @ref _greencross "greencross" object ready for
 *                approximating the integral equation
 *                @f$ \int\limits_{\Omega} -\frac{1}{2 \pi}
 *                \log{\left \Vert x - y \right \Vert_2} u(y) dy = f(x) @f$. */
HEADER_PREFIX pgreencross
new_greencross_laplace2d(pcurve2d c2d, uint res, uint q, uint m);

/** @} */

#endif // LAPLACE2D_H
