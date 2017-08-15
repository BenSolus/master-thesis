#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE laplace3d

#include <boost/test/unit_test.hpp>

#include "global_fixture.h"

#include "laplace3d.h"

#include "macrosurface3d.h"

#include <chrono>

using namespace boost::unit_test::framework;

BOOST_GLOBAL_FIXTURE(global_fixture);

static pgreencross gc;
static ph2matrix H2;
static pavector x;
static pavector y_ref;

BOOST_AUTO_TEST_CASE(reference)
{
  for(uint i = 0; i < ocl_system.num_devices; ++i)
  {
    cl_ulong                 cache_size;
    cl_device_mem_cache_type cache_type;

    CL_CHECK(clGetDeviceInfo(ocl_system.devices[i],
                             CL_DEVICE_GLOBAL_MEM_CACHE_SIZE,
                             sizeof(cl_ulong),
                             &cache_size,
                             NULL));

    CL_CHECK(clGetDeviceInfo(ocl_system.devices[i],
                             CL_DEVICE_GLOBAL_MEM_CACHE_TYPE,
                             sizeof(cl_device_mem_cache_type),
                             &cache_type,
                             NULL));

    printf("Device %u: Global memory cache of type %s with %lu bytes in size\n",
           i,
           cache_type == CL_NONE
             ? "CL_NONE"
             : cache_type == CL_READ_ONLY_CACHE
               ? "CL_READ_ONLY_CACHE"
               : cache_type == CL_READ_WRITE_CACHE
                 ? "CL_READ_WRITE_CACHE"
                 : "UNKNOWN",
           cache_size);
  }
  pmacrosurface3d mg;
  psurface3d      gr;

  mg      = new_sphere_macrosurface3d();

  gr      = build_from_macrosurface3d_surface3d(mg, REAL_SQRT((real) n * 0.125));

  gc      = new_greencross_laplace3d(gr, res, q, m, aca);

  printf("Cluster resolution: %u\n", res);

  H2      = build_green_cross_h2matrix_greencross(gc, (void *) &eta);

  printf("H2-Matrix): %u blocks\n"
         "            %.5f GB\n",
         H2->desc,
         getsize_h2matrix(H2) / (1024.0 * 1024.0 * 1024.0));

  x       = new_avector(n);
  random_avector(x);

  y_ref   = new_avector(n);
  clear_avector(y_ref);

  auto begin(std::chrono::high_resolution_clock::now());

  mvm_h2matrix_avector(1.0, false, H2, x, y_ref);

  auto end(std::chrono::high_resolution_clock::now());

  std::cout << "Reference H2-MVM:  "
            << std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count()
            << " ns\n";

  del_macrosurface3d(mg);
}

BOOST_AUTO_TEST_CASE(greencross1)
{
  pavector y = new_avector(n);
  clear_avector(y);

  auto begin(std::chrono::high_resolution_clock::now());

  mvm_h2matrix_avector_greencross(gc, 1.0, false, H2, x, y);

  auto end(std::chrono::high_resolution_clock::now());

  add_avector(r_minusone, y_ref, y);

  std::cout << "Greencross H2-MVM: "
            << std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count()
            << " ns, rel. error: " << std::scientific
            << norm2_avector(y) / norm2_avector(y_ref)
            << "\n";

  del_avector(y);
}

BOOST_AUTO_TEST_CASE(cleanup)
{
  del_avector(y_ref);
  del_avector(x);
  del_h2matrix(H2);
  del_greencross(gc);
}