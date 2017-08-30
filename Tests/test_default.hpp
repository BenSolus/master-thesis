/* ------------------------------------------------------------
 * This is the file "test_default.hpp" of this master thesis.
 * All rights reserved, Bennet Carstensen 2017
 * ------------------------------------------------------------ */

/**
 * @file      test_default.hpp
 * @author    Bennet Carstensen
 * @date      2017
 * @copyright All rights reserved, Bennet Carstensen 2017
 */

#ifndef TEST_DEFAULT_H
#define TEST_DEFAULT_H

#include "gtest/gtest.h"

#include "basic.h"

static real eta   = 1.0;   // Parameter for the accuracy of hierarchical
                           // clustering
static uint n     = 2048;  // Problem size
static uint m     = 8;     // Approximation order
static uint q     = 2;     // Quadratur order
static uint res   = m * m; // Cluster resolution
static real aca   = 10e-4; // ACA resolution

static real accur = 1e-4;
static uint tests = 100;   // Number of random/benchmark tests performed

#endif // TEST_DEFAULT
