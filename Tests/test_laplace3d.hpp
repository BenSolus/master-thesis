
#ifndef LAPLACE3D_H
#define LAPLACE3D_H

#include "test_default.hpp"

#include "laplacebem3d.h"
#include "laplace3dgca.h"
#include "macrosurface3d.h"
#include "matrixnorms.h"
#include "ocl_system.h"

static pamatrix        G;
static pgreencross     gc;
static pmacrosurface3d mg;

static real            norm_full;

TEST(laplace3d, full_matrix)
{
  pamatrix        Gbem;
  pbem3d          bem;
  psurface3d      gr;

  real        rel_err;

  mg      = new_sphere_macrosurface3d();

  gr      = build_from_macrosurface3d_surface3d(mg, REAL_SQRT(n * 0.125));

  gc      = new_greencross_laplace3d(gr, res, q, m, aca);

  bem     = new_slp_laplace_bem3d(gr,
                                  q,
                                  q + 2,
                                  BASIS_CONSTANT_BEM3D,
                                  BASIS_CONSTANT_BEM3D);

  G       = build_amatrix_greencross(gc, n, NULL, n, NULL, false);

  Gbem    = new_amatrix(n, n);
  bem->nearfield(NULL, NULL, bem, false, Gbem);

  rel_err = norm2diff_amatrix(G, Gbem) / norm2_amatrix(Gbem);

  std::cout << "\n-----------------------------------------------------------\n"
            << "\nTesting constructing full matrix of the fundamental "
            << "solution of the 3D Laplace equation...\n"
            << "\nRelative error: "
            << rel_err;

  if(rel_err < accur)
    std::cout << " (Test succeeded)\n\n";
  else
    std::cout << "\n\n";

  ASSERT_EQ(rel_err < accur, true);

  del_amatrix(Gbem);
}

TEST(laplace3d, low_rank_approximation)
{
  prkmatrix   RK;

  real        rel_err;

  RK        = build_green_rkmatrix_greencross(gc, gc->rc, gc->cc);

  norm_full = norm2_amatrix(G);

  rel_err   = norm2diff_amatrix_rkmatrix(RK, G) / norm_full;

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
  {
    ASSERT_EQ
      (ABS(rel_err - r_one) < std::numeric_limits<float>::epsilon(), true);
  }

   // TODO: Error tolerance for real = float

  del_rkmatrix(RK);
}
//
// BOOST_AUTO_TEST_CASE(test_green_laplace2d)
// {
//   pamatrix        G;
//   pgreencross     gc;
//   phmatrix        H;
//   pmacrosurface3d mg;
//   psurface3d      gr;
//
//   real        rel_err;
//
//   mg      = new_sphere_macrosurface3d();
//
//   gr      = build_from_macrosurface3d_surface3d(mg, REAL_SQRT(n * 0.125));
//
//   gc      = new_greencross_laplace3d(gr, res, q, m);
//
//   H       = build_green_hmatrix_greencross(gc, (void *) &eta);
//
//   G       = build_amatrix_greencross(gc, n, NULL, n, NULL, false);
//
//   rel_err = norm2diff_amatrix_hmatrix(H, G) / norm2_amatrix(G);
//
//   std::cout << "\n-----------------------------------------------------------\n"
//             << "\nTesting constructing hierarchical matrix of the "
//             << "fundamental solution of the 3D Laplace equation...\n"
//             << "\nRelative error: "
//             << rel_err;
//
//   if(rel_err < std::numeric_limits<float>::epsilon())
//     std::cout << " (Test succeeded)\n\n";
//   else
//     std::cout << "\n\n";
//
//   if(sizeof(real) == 8) // real = double
//     BOOST_CHECK_EQUAL(rel_err < accur, true);
//
//   // TODO: Error tolerance for real = float
//
//   del_amatrix(G);
//   del_hmatrix(H);
//   del_clusterbasis(gc->rb);
//   del_clusterbasis(gc->cb);
//   del_greencross(gc);
// }

static ph2matrix   H2;

TEST(laplace3d, building_h2_matrix)
{
  (void) tests;

  real rel_err;

  H2      = build_green_cross_h2matrix_greencross(gc, (void *) &eta);

  rel_err = norm2diff_amatrix_h2matrix(H2, G) / norm_full;

  std::cout << "\n-----------------------------------------------------------\n"
            << "\nTesting constructing H^2 matrix of the "
            << "fundamental solution of the 3D Laplace equation...\n"
            << "\nRelative error: "
            << std::scientific << rel_err;

  if (rel_err < accur)
    std::cout << " (Test succeeded)\n\n";
  else
    std::cout << "\n\n";

  if(sizeof(real) == 8) // real = double
  {
    ASSERT_EQ(rel_err < accur, true);
  }
  // TODO: Error tolerance for real = float
}

TEST(laplace3d, h2_mvm)
{
  pavector  x, y_gc, y_ref;

  real      rel_err;

  x       = new_avector(n);
  random_avector(x);

  y_gc    = new_avector(n);
  clear_avector(y_gc);
  y_ref   = new_avector(n);
  clear_avector(y_ref);

  std::cout << "\n-----------------------------------------------------------\n"
            << "\nTesting calculating the Matrix-Vector-Product of an "
            << "H^2-matrix resulting from a greencross \nobject describing the "
            << "fundamental solution of the 3D Laplace equation...\n";

  mvm_h2matrix_avector_greencross(gc, 1.0, false, H2, x, y_gc, 0);
  mvm_h2matrix_avector(1.0, false, H2, x, y_ref);

  add_avector(r_minusone, y_ref, y_gc);

  rel_err = norm2_avector(y_gc) / norm2_avector(y_ref);

  std::cout << "\nRelative error: " << std::scientific
            << rel_err;

  if(rel_err < std::numeric_limits<float>::epsilon())
    std::cout << " (Test succeeded)\n\n";
  else
    std::cout << "\n\n";

  if(sizeof(real) == 8) // real = double
  {
    ASSERT_EQ(rel_err < std::numeric_limits<float>::epsilon(), true);
  }

  del_avector(y_ref);
  del_avector(y_gc);
  del_avector(x);
 }

TEST(laplace3d, cleanup)
{
  del_h2matrix(H2);
  del_amatrix(G);
  del_greencross(gc);
  del_macrosurface3d(mg);

  for(uint i = 0; i < ocl_system.num_devices; ++i)
    CL_CHECK(clReleaseContext(ocl_system.contexts[i]));

  std::cout << "[----------] "
            << getactives_avector() << " vectors, "
            << getactives_amatrix() << " matrices and "
            << getactives_clusterbasis() << " clusterbases still reachable\n";

  ASSERT_EQ(getactives_avector(), 0);
  ASSERT_EQ(getactives_amatrix(), 0);
  ASSERT_EQ(getactives_clusterbasis(), 0);

  uninit_h2lib();
}
#endif // LAPLACE3D