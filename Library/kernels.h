/* ------------------------------------------------------------
 * This is the file "kernels.h" of this master thesis.
 * All rights reserved, Bennet Carstensen 2017
 * ------------------------------------------------------------ */

/**
 * @file      kernels.h
 * @author    Bennet Carstensen
 * @date      2017
 * @copyright All rights reserved, Bennet Carstensen 2017
 */

#ifndef KERNELS_H
#define KERNELS_H

#include "basic.h"

/** @defgroup kernels kernels
 *  @brief Set of kernel functions and their partial derivative usable in the
 *         Green cross approximation method.
 *
 * @{*/

/** @brief @f$- \frac{1}{2 \pi}@f$ */
static const real r_minus_two_pi = -0.159154943091895336;

/** @brief @f$\frac{1}{4 \pi}@f$ */
static const real r_four_pi = 0.0795774715459476679;

/** @brief @f$- \frac{1}{\pi}@f$ */
static const real r_minus_pi = -0.318309886183790671;

/** @brief Implementation of the fundamental solution of the negative Laplace
 *         operator @f$- \Delta@f$ in 2D and 3D.
 *
 * Calculates the fundamental solution of the negative Laplace operator:
 *
 * \f[
 *   g(x, y) :=
 *     \begin{cases}
 *      -\frac{1}{2 \pi} \log{\left \Vert x - y \right \Vert_2} &
 *      \text{if } x \neq y, \ dim = 2 \\
 *      \\
 *      \frac{1}{4 \pi \left \Vert x - y \right \Vert_2} &
 *      \text{if } x \neq y, \ dim = 3 \\
 *      \\
 *      0 & \text{otherwise}
 *    \end{cases}
 * \f]
 *
 *  @param[in] x First operand as an array of @p dim coefficients.
 *  @param[in] y Second operand as an array of @p dim coefficients.
 *  @param[in] dim Dimension in which to evaluate the negative Laplace operator.
 *  @return        The result of the @p dim dimensional negative Laplace
 *                 operator in the parameters @p x and @p y.
 */
HEADER_PREFIX real
laplace_kernel(const field* x, const field* y, const uint dim);

/** @brief Calculates the partial derivative of the negative Laplace operator
 *         @f$- \Delta@f$ in the \p i th component of the first operand in 2D
 *         and 3D.
 *
 * Calculates the partial derivative of the fundamental solution of the negative
 * Laplace operatorin in the \p i th component of the first operand:
 *
 * \f[
 *   \frac{\delta g}{\delta x_{i}}(x, y) :=
 *     \begin{cases}
 *      - \frac{x_i - y_i}{\pi \left \Vert x - y \right \Vert_{2}^{2}} &
 *      \text{if } x \neq y, \ dim = 2 \\
 *      \\
 *      - \frac{x_i - y_i}{4 \pi \left \Vert x - y \right \Vert_{2}^{3}} &
 *      \text{if } x \neq y, \ dim = 3 \\
 *      \\
 *      0 & \text{otherwise}
 *    \end{cases}
 * \f]
 *
 *  @param[in] x   First operand as an array of @p dim coefficients. Parital
 *                 derivating the @p i th component of this operand.
 *  @param[in] y   Second operand as an array of @p dim coefficients.
 *  @param[in] dim Dimension in which to evaluate the partial derivative of the
 *                 negative Laplace operator
 *  @param[in] i   Index of the component of @p x in which to parital derive.
 *  @return        The result of the @p dim dimensional partial derivative of
 *                 the fundamental solution of the negative Laplace operator. */
HEADER_PREFIX real
pdx_laplace_kernel(const field *x,
                   const field *y,
                   const uint dim,
                   const uint i);

/** @brief Calculates the partial derivative of the negative Laplace operator
 *         @f$- \Delta@f$ in the \p i th component of the second operand in 2D
 *         and 3D.
 *
 * Calculates the partial derivative of the fundamental solution of the negative
 * Laplace operatorin in the \p i th component of the second operand:
 *
 * \f[
 *   \frac{\delta g}{\delta y_{i}}(x, y) :=
 *     \begin{cases}
 *      \frac{x_i - y_i}{\pi \left \Vert x - y \right \Vert_{2}^{2}} &
 *      \text{if } x \neq y, \ dim = 2 \\
 *      \\
 *      \frac{x_i - y_i}{4 \pi \left \Vert x - y \right \Vert_{2}^{3}} &
 *      \text{if } x \neq y, \ dim = 3 \\
 *      \\
 *      0 & \text{otherwise}
 *    \end{cases}
 * \f]
 *
 *  @param[in] x   First operand as an array of @p dim coefficients.
 *  @param[in] y   Second operand as an array of @p dim coefficients. Parital
 *                 derivating the @p i th component of this operand.
 *  @param[in] dim Dimension in which to evaluate the partial derivative of the
 *                 negative Laplace operator
 *  @param[in] i   Index of the component of @p y in which to parital derive.
 *  @return        The result of the @p dim dimensional partial derivative of
 *                 the fundamental solution of the negative Laplace operator. */
HEADER_PREFIX real
pdy_laplace_kernel(const field *x,
                   const field *y,
                   const uint   dim,
                   const uint   i);

/** @} */

#endif // GREENCROSS_H
