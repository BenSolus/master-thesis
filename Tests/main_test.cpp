/* ------------------------------------------------------------
 * This is the file "main_test.hpp" of this master thesis.
 * All rights reserved, Bennet Carstensen 2017
 * ------------------------------------------------------------ */

/**
 * @file      main_test.hpp
 * @author    Bennet Carstensen
 * @date      2017
 * @copyright All rights reserved, Bennet Carstensen 2017
 */

#ifndef MAIN_TEST_HPP
#define MAIN_TEST_HPP

#include "test_default.hpp"

#include "test_laplace3d.hpp"

int
main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}

#endif // MAIN_TEST_HPP
