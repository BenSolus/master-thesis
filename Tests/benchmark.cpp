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
  pmacrosurface3d mg;
  psurface3d      gr;

  mg      = new_sphere_macrosurface3d();

  gr      = build_from_macrosurface3d_surface3d(mg, REAL_SQRT((real) n * 0.125));

  gc      = new_greencross_laplace3d(gr, res, q, m, 10e-15);

  H2      = build_green_cross_h2matrix_greencross(gc, (void *) &eta);

  x       = new_avector(n);
  random_avector(x);

  y_ref   = new_avector(n);
  clear_avector(y_ref);

  auto begin(std::chrono::high_resolution_clock::now());

  mvm_h2matrix_avector(1.0, false, H2, x, y_ref);

  auto end(std::chrono::high_resolution_clock::now());

  std::cout << "Reference H2-MVM: "
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