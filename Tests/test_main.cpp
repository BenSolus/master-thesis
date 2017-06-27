#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE main

#include <boost/test/unit_test.hpp>

#include "global_fixture.h"

#include "greencross.h"

#include "laplacebem2d.h"
#include "matrixnorms.h"

#include <iostream>
#include <limits>
#include <iomanip>
using namespace boost::unit_test::framework;

BOOST_GLOBAL_FIXTURE(global_fixture);

/* ------------------------------------------------------------
 * "Header" for static functions
 * ------------------------------------------------------------ */

INLINE_PREFIX void
norm2diff_h2matrix_leafs_amatrix(ph2matrix h2, pamatrix a);

BOOST_AUTO_TEST_CASE(test_full_laplace2d)
{
  pamatrix    Gbem, Ggc;
  pbem2d      b2d;
  pcurve2d    c2d;
  pgreencross gc;

  real        rel_err;

  c2d  = new_circle_curve2d(n, 0.333);

  gc   = new_laplace2d_greencross(c2d, res, q, m);

  b2d  = new_slp_laplace_bem2d(c2d, q, BASIS_CONSTANT_BEM2D);

  Ggc  = new_amatrix(n, n);
  nearfield_greencross(gc, n, NULL, n, NULL, Ggc);

  Gbem = new_amatrix(n, n);
  b2d->nearfield(NULL, NULL, b2d, false, Gbem);

  rel_err = norm2diff_amatrix(Ggc, Gbem) / norm2_amatrix(Gbem);

  std::cout << "\n-----------------------------------------------------------\n"
            << "\nTesting constructing full matrix of the fundamental "
            << "solution of the 2D Laplace equation...\n"
            << "\nRelative error: "
            << rel_err
            << "\n\nTest "
            << ((rel_err < eps) ? "succeeded"
                                : ("failed (Relative Error >= " +
                                   std::to_string(eps) +
                                   " (Machine epsilon))"))
            << ".\n";

  BOOST_REQUIRE_EQUAL(rel_err < std::numeric_limits<real>::epsilon(), true);

  del_amatrix(Gbem);
  del_amatrix(Ggc);
  del_bem2d(b2d);
  del_greencross(gc);
}

BOOST_AUTO_TEST_CASE(test_rk_laplace2d)
{
  pamatrix    G, Gbem;
  pbem2d      b2d;
  pblock      broot;
  pcluster    root;
  pcurve2d    c2d;
  pgreencross gc;
  prkmatrix   RKbem, RKgc;

  real        rel_err;

  c2d   = new_circle_curve2d(n, 0.333);

  gc    = new_laplace2d_greencross(c2d, res, q, m);

  b2d   = new_slp_laplace_bem2d(c2d, q, BASIS_CONSTANT_BEM2D);

  root  = build_bem2d_cluster(b2d, res, BASIS_CONSTANT_BEM2D);

  broot = build_nonstrict_block(root, root, &eta, admissible_max_cluster);

  setup_hmatrix_aprx_green_row_bem2d(b2d,
                                     root,
                                     root,
                                     broot,
                                     m,
                                     1,
                                     0.5,
                                     build_bem2d_rect_quadpoints);

  RKbem = build_bem2d_rkmatrix(root, root, (void *) b2d);

  Gbem = new_amatrix(n, n);
  b2d->nearfield(NULL, NULL, b2d, false, Gbem);

  // print_amatrix(&RKbem->A);
  // print_amatrix(&RKbem->B);

  RKgc  = build_green_rkmatrix_greencross(gc, gc->rc, gc->cc);

  G     = new_amatrix(n, n);
  nearfield_greencross(gc, n, NULL, n, NULL, G);

  rel_err = norm2diff_amatrix_rkmatrix(RKgc, G) / norm2_amatrix(G);

  std::cout << "\n-----------------------------------------------------------\n"
            << "\nTesting constructing low-rank approximation of the "
            << "fundamental solution of the 2D Laplace equation...\n"
            << "\nRelative error - 1: "
            << rel_err - r_one
            << "\n\nTest "
            << ((rel_err - r_one <= 10e-07) ? "succeeded"
                                            : ("failed (Relative Error - 1 > " +
                                              std::to_string(10e-07)))
            << ".\n";

  BOOST_REQUIRE_EQUAL(rel_err - r_one <= 10e-07, true);

  del_amatrix(G);
  del_rkmatrix(RKbem);
  del_rkmatrix(RKgc);
  del_block(broot);
  del_cluster(root);
  del_bem2d(b2d);
  del_greencross(gc);
}

// BOOST_AUTO_TEST_CASE(test_green)
// {
//   pamatrix    G, Gbem;
//   pblock      broot;
//   pbem2d      b2d;
//   pcluster    root;
//   pcurve2d    c2d;
//   pgreencross gc;
//   phmatrix    H, Hbem;
//
//   c2d  = new_circle_curve2d(n, 0.333);
//
//   b2d  = new_slp_laplace_bem2d(c2d, q, BASIS_CONSTANT_BEM2D);
//
//     root  = build_bem2d_cluster(b2d, res, BASIS_CONSTANT_BEM2D);
//
//     broot = build_nonstrict_block(root, root, &eta, admissible_max_cluster);
//
//     setup_hmatrix_aprx_green_row_bem2d(b2d,
//                                        root,
//                                        root,
//                                        broot,
//                                        m,
//                                        1,
//                                        0.5,
//                                        build_bem2d_rect_quadpoints);
//
//   Hbem = build_from_block_hmatrix(broot, 0);
//   assemble_bem2d_hmatrix(b2d, broot, Hbem);
//
//   Gbem = new_amatrix(n, n);
//   b2d->nearfield(NULL, NULL, b2d, false, Gbem);
//
//   std::cout << norm2diff_amatrix_hmatrix(Hbem, Gbem) / norm2_amatrix(Gbem) << "\n";
//   gc   = new_laplace2d_greencross(c2d, res, q, m);
//
//   H    = fill_green_hmatrix_greencross(gc, (void *) &eta);
//
//   G    = new_amatrix(gc->rc->size, gc->cc->size);
//
//   nearfield_greencross(gc,
//                        gc->rc->size,
//                        NULL,
//                        gc->cc->size,
//                        NULL,
//                        G);
//
//   std::cout << "Relative error: "
//             << norm2diff_amatrix_hmatrix(H, G) / norm2_amatrix(G)
//             << "\n";
//
//   del_amatrix(G);
//   del_hmatrix(H);
//   del_greencross(gc);
// }
//
// BOOST_AUTO_TEST_CASE(test_leafs)
// {
//   pamatrix    G;
//   pcurve2d    gr;
//   pgreencross gc;
//   ph2matrix   h2;
//
//   gr = new_circle_curve2d(n, r_one / (real) 3.0);
//
//   gc = new_laplace2d_greencross(gr, res, q, m);
//
//   h2 = green_cross_approximation(gc, (void *) &eta);
//
//   G  = new_amatrix(gc->rc->size, gc->cc->size);
//
//   nearfield_greencross(gc,
//                        gc->rc->size,
//                        gc->rc->idx,
//                        gc->cc->size,
//                        gc->cc->idx,
//                        G);
//
//   norm2diff_h2matrix_leafs_amatrix(h2, G);
//
//   del_amatrix(G);
//   del_h2matrix(h2);
//   del_greencross(gc);
// }

void
norm2diff_h2matrix_leafs_amatrix(ph2matrix h2, pamatrix a)
{
  if(h2->son)
  {
    const uint rsons = h2->rsons;
    const uint csons = h2->csons;

    uint coff = 0;

    amatrix  tmp;
    pamatrix asub;

    for(uint j = 0; j < csons; ++j)
    {
      const uint cols = h2->cb->son ? h2->cb->son[j]->t->size : h2->cb->t->size;

      uint roff = 0;

      for(uint i = 0; i < rsons; ++i)
      {
        const uint rows = h2->rb->son ? h2->rb->son[i]->t->size
                                      : h2->rb->t->size;

        asub  = init_sub_amatrix(&tmp, a, rows, roff, cols, coff);

        norm2diff_h2matrix_leafs_amatrix(h2->son[i + j * rsons], asub);

        roff += rows;
      }

      coff += cols;
    }
  }
  else
  {
    pamatrix G, tmp;

    real     err;

    if(h2->u)
    {
      tmp = new_zero_amatrix(h2->rb->t->size, h2->u->S.cols);

      addmul_amatrix(1.0, false, &h2->rb->V, false, &h2->u->S, tmp);

      G   = new_zero_amatrix(h2->rb->t->size, h2->cb->t->size);

      addmul_amatrix(1.0, false, tmp, true, &h2->cb->V, G);

      err = norm2diff_amatrix(G, a);

      //if(err >= accur)
        std::cout << std::scientific
                  << "Error was: "
                  << err
                  << ". Maximum fault tolerance: "
                  << accur
                  << "\n";

      BOOST_REQUIRE_EQUAL(err < accur, true);

      del_amatrix(G);
      del_amatrix(tmp);
    }
    else if(h2->f)
    {
      err = norm2diff_amatrix(h2->f, a);

      BOOST_REQUIRE_EQUAL(err < std::numeric_limits<real>::epsilon(), true);
    }
  }
}

/** @brief Test case for testing the whole library in an usual use case.
 *
 * Tests the whole library for segmentation faults, double frees or
 * corruptions in an usual use case. */
// BOOST_AUTO_TEST_CASE(test_main)
// {
//   pamatrix    G;
//   pgreencross gc;
//   pcurve2d    gr;
//
//   real        err;
//
//   gr = new_circle_curve2d(n, (real) 1.0 / 3.0);
//
//   gc = new_laplace2d_greencross(gr, res, (void *) &eta, q, m);
//
//   assemble_leafs_h2matrix_greencross(gc, gc->h2);
//
//   G  = new_amatrix(gc->h2->rb->t->size, gc->h2->cb->t->size);
//
//   nearfield_greencross(gc, gc->h2->rb, gc->h2->cb, G);
//
//   print_amatrix(G);
//
//   err = norm2diff_h2matrix_leafs_amatrix(gc->h2, G);
//
//   printf("Error: %.5e\n", err);
//
//   del_amatrix(G);
//   del_greencross(gc);
// }
