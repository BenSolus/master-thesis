#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE laplace3d

#include <boost/test/unit_test.hpp>

#include "global_fixture.h"

#include "laplace3d.h"

#include "laplacebem3d.h"
#include "matrixnorms.h"

using namespace boost::unit_test::framework;

BOOST_GLOBAL_FIXTURE(global_fixture);

BOOST_AUTO_TEST_CASE(test_full)
{
  pamatrix        Gbem, Ggc;
  pbem3d          bem;
  pgreencross     gc;
  pmacrosurface3d mg;
  psurface3d      gr;

  real        rel_err;

  mg      = new_sphere_macrosurface3d();

  gr      = build_from_macrosurface3d_surface3d(mg, REAL_SQRT(n * 0.125));

  gc      = new_greencross_laplace3d(gr, res, q, m);

  bem     = new_slp_laplace_bem3d(gr,
                                  q,
                                  q + 2,
                                  BASIS_CONSTANT_BEM3D,
                                  BASIS_CONSTANT_BEM3D);

  Ggc     = build_amatrix_greencross(gc, n, NULL, n, NULL, false);

  Gbem    = new_amatrix(n, n);
  bem->nearfield(NULL, NULL, bem, false, Gbem);

  rel_err = norm2diff_amatrix(Ggc, Gbem) / norm2_amatrix(Gbem);

  std::cout << "\n-----------------------------------------------------------\n"
            << "\nTesting constructing full matrix of the fundamental "
            << "solution of the 3D Laplace equation...\n"
            << "\nRelative error: "
            << rel_err;

  if(rel_err < eps)
    std::cout << " (Test succeeded)\n\n";
  else
    std::cout << "\n\n";

  BOOST_CHECK_EQUAL(rel_err < eps, true);

  del_amatrix(Gbem);
  del_amatrix(Ggc);
  del_bem3d(bem);
  del_clusterbasis(gc->rb);
  del_clusterbasis(gc->cb);
  del_greencross(gc);
}

BOOST_AUTO_TEST_CASE(test_rk_laplace2d)
{
  pamatrix        G;
  pgreencross     gc;
  pmacrosurface3d mg;
  prkmatrix       RK;
  psurface3d      gr;

  real        rel_err;

  mg      = new_sphere_macrosurface3d();

  gr      = build_from_macrosurface3d_surface3d(mg, REAL_SQRT(n * 0.125));

  gc      = new_greencross_laplace3d(gr, res, q, m);

  RK      = build_green_rkmatrix_greencross(gc, gc->rc, gc->cc);

  G       = build_amatrix_greencross(gc, n, NULL, n, NULL, false);

  rel_err = norm2diff_amatrix_rkmatrix(RK, G) / norm2_amatrix(G);

  std::cout << "\n-----------------------------------------------------------\n"
            << "\nTesting constructing low-rank approximation of the "
            << "fundamental solution of the 3D Laplace equation...\n"
            << "\n|Relative error - 1|: "
            << ABS(rel_err - r_one);

  if(ABS(rel_err - r_one) < accur)
    std::cout << " (Test succeeded)\n\n";
  else
    std::cout << "\n\n";

  if(sizeof(real) == 8) // real = double
    BOOST_CHECK_EQUAL
      (ABS(rel_err - r_one) < std::numeric_limits<float>::epsilon(), true);

  // TODO: Error tolerance for real = float

  del_amatrix(G);
  del_rkmatrix(RK);
  del_clusterbasis(gc->rb);
  del_clusterbasis(gc->cb);
  del_greencross(gc);
}

BOOST_AUTO_TEST_CASE(test_green_laplace2d)
{
  pamatrix        G;
  pgreencross     gc;
  phmatrix        H;
  pmacrosurface3d mg;
  psurface3d      gr;

  real        rel_err;

  mg      = new_sphere_macrosurface3d();

  gr      = build_from_macrosurface3d_surface3d(mg, REAL_SQRT(n * 0.125));

  gc      = new_greencross_laplace3d(gr, res, q, m);

  H       = build_green_hmatrix_greencross(gc, (void *) &eta);

  G       = build_amatrix_greencross(gc, n, NULL, n, NULL, false);

  rel_err = norm2diff_amatrix_hmatrix(H, G) / norm2_amatrix(G);

  std::cout << "\n-----------------------------------------------------------\n"
            << "\nTesting constructing hierarchical matrix of the "
            << "fundamental solution of the 3D Laplace equation...\n"
            << "\nRelative error: "
            << rel_err;

  if(rel_err < std::numeric_limits<float>::epsilon())
    std::cout << " (Test succeeded)\n\n";
  else
    std::cout << "\n\n";

  if(sizeof(real) == 8) // real = double
    BOOST_CHECK_EQUAL(rel_err < accur, true);

  // TODO: Error tolerance for real = float

  del_amatrix(G);
  del_hmatrix(H);
  del_clusterbasis(gc->rb);
  del_clusterbasis(gc->cb);
  del_greencross(gc);
}

BOOST_AUTO_TEST_CASE(test_green_cross)
{
  pamatrix        G;
  pgreencross     gc;
  ph2matrix       H2;
  pmacrosurface3d mg;
  psurface3d      gr;

  real        rel_err;

  mg      = new_sphere_macrosurface3d();

  gr      = build_from_macrosurface3d_surface3d(mg, REAL_SQRT(n * 0.125));

  gc      = new_greencross_laplace3d(gr, res, q, m);

  H2      = build_green_cross_h2matrix_greencross(gc, (void *) &eta);

  G       = build_amatrix_greencross(gc, n, NULL, n, NULL, false);

  rel_err = norm2diff_amatrix_h2matrix(H2, G) / norm2_amatrix(G);

  std::cout << "\n-----------------------------------------------------------\n"
            << "\nTesting constructing H^2 matrix of the "
            << "fundamental solution of the 3D Laplace equation...\n"
            << "\nRelative error: "
            << rel_err;

  if(rel_err < std::numeric_limits<float>::epsilon())
    std::cout << " (Test succeeded)\n\n";
  else
    std::cout << "\n\n";

  if(sizeof(real) == 8) // real = double
    BOOST_CHECK_EQUAL(rel_err < std::numeric_limits<float>::epsilon(), true);

  // TODO: Error tolerance for real = float

  del_amatrix(G);
  del_h2matrix(H2);
  del_greencross(gc);
}
