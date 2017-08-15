/* ------------------------------------------------------------
 * This is the file global_fixture.h" of this master thesis.
 * All rights reserved, Bennet Carstensen 2017
 * ------------------------------------------------------------ */

/** @file      global_fixture.h
 *  @author    Bennet Carstensen
 *  @date      2017
 *  @copyright All rights reserved, Bennet Carstensen 2017
 */

#ifndef GLOBAL_FIXTURE_H
#define GLOBAL_FIXTURE_H

#include <boost/test/unit_test.hpp>

#include "amatrix.h"
#include "basic.h"
#include "clusterbasis.h"
#include "opencl.h"

#include <iostream>
#include <limits>

using namespace boost::unit_test::framework;

/** @defgroup global_fixture global_fixture
 *  @brief General initializations and cleanup routines for all test moduls.
 *
 *  @{ */

static const real eps = std::numeric_limits<real>::epsilon();

static real eta   = 1.0;     // Parameter for the accuracy of hierarchical
                             // clustering
static uint n     = 8192;    // Problem size
static uint m     = 16;       // Approximation order
static uint q     = 2;       // Quadratur order
static uint res   = 8 * m * m;     // Cluster resolution
static real aca   = 10e-4;   // ACA resolution

/** @brief Global fixture for the
  *        <a href="https://boost.org/libs/test">Boost Test Library</a>
  *
  * Performes initialization routines before any test case of an module and
  * cleanup after all test cases. */
struct global_fixture
{
  global_fixture()
  {
    init_h2lib(&master_test_suite().argc, &master_test_suite().argv);

    std::cout << "\nMachine epsilon: " << eps << "\n\n";

    cl_device_id *devices;
    cl_uint      ndevices;



    get_opencl_devices(&devices, &ndevices);

    set_opencl_devices(devices, ndevices, 1);
  }

  ~global_fixture()
  {
    for(uint i = 0; i < ocl_system.num_devices; ++i)
      CL_CHECK(clReleaseContext(ocl_system.contexts[i]));

    printf("  %u vectors, %u matrices and %u clusterbasis\n",
      getactives_avector(),
      getactives_amatrix(),
      getactives_clusterbasis());

    uninit_h2lib();
  }
};

/** @} */

#endif // GLOBAL_FIXTURE_H
