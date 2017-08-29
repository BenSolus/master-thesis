#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE laplace2d

#include <boost/test/unit_test.hpp>

#include "global_fixture.h"

#include "laplace2d.h"

#include "laplacebem2d.h"
#include "matrixnorms.h"

using namespace boost::unit_test::framework;

BOOST_GLOBAL_FIXTURE(global_fixture);

static const real accur = 10e-04;

BOOST_AUTO_TEST_CASE(test_full)
{
  pamatrix    Gbem, Ggc;
  pbem2d      b2d;
  pcurve2d    c2d;
  pgreencross gc;

  real        rel_err;

  c2d  = new_circle_curve2d(n, 0.333);

  gc   = new_greencross_laplace2d(c2d, res, q, m);

  b2d  = new_slp_laplace_bem2d(c2d, q, BASIS_CONSTANT_BEM2D);

  Ggc  = build_amatrix_greencross(gc, n, NULL, n, NULL, false);

  Gbem = new_amatrix(n, n);
  b2d->nearfield(NULL, NULL, b2d, false, Gbem);

  rel_err = norm2diff_amatrix(Ggc, Gbem) / norm2_amatrix(Gbem);

  std::cout << "\n-----------------------------------------------------------\n"
            << "\nTesting constructing full matrix of the fundamental "
            << "solution of the 2D Laplace equation...\n"
            << "\nRelative error: "
            << rel_err;

  if(rel_err < eps)
    std::cout << " (Test succeeded)\n\n";
  else
    std::cout << "\n\n";

  BOOST_CHECK_EQUAL(rel_err < eps, true);

  del_amatrix(Gbem);
  del_amatrix(Ggc);
  del_bem2d(b2d);
  del_clusterbasis(gc->rb);
  del_clusterbasis(gc->cb);
  del_greencross(gc);
}

BOOST_AUTO_TEST_CASE(test_rk)
{
  pamatrix    G;
  pcurve2d    c2d;
  pgreencross gc;
  prkmatrix   RK;

  real        rel_err;

  c2d     = new_circle_curve2d(n, 0.333);

  gc      = new_greencross_laplace2d(c2d, res, q, m);

  RK      = build_green_rkmatrix_greencross(gc, gc->rc, gc->cc);

  G       = build_amatrix_greencross(gc, n, NULL, n, NULL, false);

  rel_err = norm2diff_amatrix_rkmatrix(RK, G) / norm2_amatrix(G);

  std::cout << "\n-----------------------------------------------------------\n"
            << "\nTesting constructing low-rank approximation of the "
            << "fundamental solution of the 2D Laplace equation...\n"
            << "\n|Relative error - 1|: "
            << ABS(rel_err - r_one);

  if(ABS(rel_err - r_one) < accur)
    std::cout << " (Test succeeded)\n\n";
  else
    std::cout << "\n\n";

  if(sizeof(real) == 8) // real = double
    BOOST_CHECK_EQUAL
      (ABS(rel_err - r_one) < accur, true);

  // TODO: Error tolerance for real = float

  del_amatrix(G);
  del_rkmatrix(RK);
  del_clusterbasis(gc->rb);
  del_clusterbasis(gc->cb);
  del_greencross(gc);
}

BOOST_AUTO_TEST_CASE(test_green)
{
  pamatrix    G;
  pcurve2d    c2d;
  pgreencross gc;
  phmatrix    H;

  real        rel_err;

  c2d  = new_circle_curve2d(n, 0.333);

  gc   = new_greencross_laplace2d(c2d, res, q, m);

  H    = build_green_hmatrix_greencross(gc, (void *) &eta);

  G    = build_amatrix_greencross(gc, n, NULL, n, NULL, false);

  rel_err = norm2diff_amatrix_hmatrix(H, G) / norm2_amatrix(G);

  std::cout << "\n-----------------------------------------------------------\n"
            << "\nTesting constructing hierarchical matrix of the "
            << "fundamental solution of the 2D Laplace equation...\n"
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
  pamatrix    G;
  pcurve2d    gr;
  pgreencross gc;
  ph2matrix   H2;

  real        rel_err;

  gr = new_circle_curve2d(n, r_one / (real) 3.0);

  gc = new_greencross_laplace2d(gr, res, q, m);

  H2 = build_green_cross_h2matrix_greencross(gc, (void *) &eta);

  G  = build_amatrix_greencross(gc, n, NULL, n, NULL, false);

  rel_err = norm2diff_amatrix_h2matrix(H2, G) / norm2_amatrix(G);

  std::cout << "\n-----------------------------------------------------------\n"
            << "\nTesting constructing H^2 matrix of the "
            << "fundamental solution of the 2D Laplace equation...\n"
            << "\nRelative error: "
            << rel_err;

  if(rel_err < accur)
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

// BOOST_AUTO_TEST_CASE(test_green_cross_mvm_laplace2d)
// {
//   pavector    x, y_gc, y_ref;
//   pcurve2d    gr;
//   pgreencross gc;
//   ph2matrix   H2;
//
//   real        rel_err;
//
//   gr = new_circle_curve2d(n, r_one / (real) 3.0);
//
//   gc = new_greencross_laplace2d(gr, res, q, m);
//
//   H2 = build_green_cross_h2matrix_greencross(gc, (void *) &eta);
//
//   x  = new_avector(n);
//   random_avector(x);
//
//   y_gc  = new_avector(n);
//   clear_avector(y_gc);
//   y_ref = new_avector(n);
//   clear_avector(y_ref);
//
//   std::cout << "\n-----------------------------------------------------------\n"
//             << "\nTesting calculating the Matrix-Vector-Product of an "
//             << "H^2-matrix resulting from a greencross \nobject describing the "
//             << "fundamental solution of the 2D Laplace equation...\n";
//
//   mvm_h2matrix_avector_greencross(1.0, false, H2, x, y_gc);
//   mvm_h2matrix_avector(1.0, false, H2, x, y_ref);
//
//   add_avector(r_minusone, y_ref, y_gc);
//
//   rel_err = norm2_avector(y_gc) / norm2_avector(y_ref);
//
//   std::cout << "\nRelative error: "
//             << rel_err;
//
//   del_avector(y_ref);
//   del_avector(y_gc);
//   del_avector(x);
//   del_h2matrix(H2);
//   del_greencross(gc);
// }
