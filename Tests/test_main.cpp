#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE main

#include <boost/test/unit_test.hpp>

#include "global_fixture.h"

#include "greencross.h"

#include "matrixnorms.h"

#include <iostream>
#include <limits>

using namespace boost::unit_test::framework;

BOOST_GLOBAL_FIXTURE(global_fixture);

/* ------------------------------------------------------------
 * "Header" for static functions
 * ------------------------------------------------------------ */

INLINE_PREFIX void
norm2diff_h2matrix_leafs_amatrix(ph2matrix h2, pamatrix a);

BOOST_AUTO_TEST_CASE(test_green)
{
  pamatrix    G;
  pcurve2d    c2d;
  pgreencross gc;
  phmatrix    H;

  c2d = new_circle_curve2d(n, r_one / 3.0);

  gc  = new_laplace2d_greencross(c2d, res, q, m);

  H   = fill_green_hmatrix_greencross(gc, (void *) &eta);

  G   = new_amatrix(gc->rc->size, gc->cc->size);

  nearfield_greencross(gc,
                       gc->rc->size,
                       gc->rc->idx,
                       gc->cc->size,
                       gc->cc->idx,
                       G);

  std::cout << norm2diff_amatrix_hmatrix(H, G) / norm2_amatrix(G) << "\n";

  del_amatrix(G);
  del_hmatrix(H);
  del_greencross(gc);
}

BOOST_AUTO_TEST_CASE(test_leafs)
{
  pamatrix    G;
  pcurve2d    gr;
  pgreencross gc;
  ph2matrix   h2;

  gr = new_circle_curve2d(n, r_one / (real) 3.0);

  gc = new_laplace2d_greencross(gr, res, q, m);

  h2 = green_cross_approximation(gc, (void *) &eta);

  G  = new_amatrix(gc->rc->size, gc->cc->size);

  nearfield_greencross(gc,
                       gc->rc->size,
                       gc->rc->idx,
                       gc->cc->size,
                       gc->cc->idx,
                       G);

  norm2diff_h2matrix_leafs_amatrix(h2, G);

  del_amatrix(G);
  del_h2matrix(h2);
  del_greencross(gc);
}

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
    pamatrix G;

    real     err;

    if(h2->u)
    {
      G = new_zero_amatrix(h2->rb->t->size, h2->u->S.cols);

      addmul_amatrix(1.0, false, &h2->rb->V, false, &h2->u->S, G);

      err = norm2diff_amatrix(G, a);

      if(err >= accur)
        std::cout << std::scientific
                  << "Error was: "
                  << err
                  << ". Maximum fault tolerance: "
                  << accur
                  << "\n";

      BOOST_REQUIRE_EQUAL(err < accur, true);

      del_amatrix(G);

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
