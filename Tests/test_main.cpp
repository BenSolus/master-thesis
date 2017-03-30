#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE main

/** @defgroup test_module
 *  @brief Test modules for testing this library.
 *
 * @{ */

#include <boost/test/unit_test.hpp>

#include "basic.h"
#include "ie1d.h"
#include "opencl.h"

#include <time.h>

using namespace boost::unit_test::framework;

static const real eta = 1.0;   // Parameter for the accuracy
static const uint k   = 8;     // Approximation order
static const uint n   = 1024; // Problem size
static const uint res = 4 * k; // Cluster resolution

// void
// forward_leaf(pclusterbasis cb, pcavector x, pavector xt)
// {
//   const uint sons = cb->sons;
//
//   avector  locxc;
//   pavector xc   = init_sub_avector(&locxc, xt, cb->k, 0);
//
//   clear_avector(xc);
//
//   if(sons)
//   {
//     uint xtoff = cb->k;
//
//     for(uint i = 0; i < sons; ++i)
//     {
//       avector  locxt1;
//       pavector xt1 = init_sub_avector(&locxt1, xt, cb->son[i]->ktree, xtoff);
//
//       forward_leaf(cb->son[i], x, xt1);
//
//       uninit_avector(xt1);
//
//       xtoff += cb->son[i]->ktree;
//     }
//
//     assert(xtoff == cb->ktree);
//   }
//   else
//   {
//     avector locxp;
//     /* Permuted entries of x */
//     pavector xp = init_sub_avector(&locxp, xt, cb->t->size, cb->k);
//
//     /* Find and copy entries */
//     for(size_t i = 0; i < cb->t->size; i++)
//       xp->v[i] = x->v[cb->t->idx[i]];
//
//     /* Multiply by leaf matrix */
//     mvm_amatrix_avector(1.0, true, &cb->V, xp, xc);
//
//     uninit_avector(xp);
//   }
//
//   uninit_avector(xc);
// }

/** @brief Test case for testing the whole library in an usual use case.
 *
 * Tests the whole library for segment segmentation faults, double frees or
 * corruptions in an usual use case. */
BOOST_AUTO_TEST_CASE(test_main)
{
  avector       locx;
  pavector      x;
  pavector      xt_cpu;
  pavector      xt_gpu;
  pblock        broot; // Root of the block tree
  pie1d         ie;    // Problems
  pcluster      root;  // Root of the cluster tree
  pclusterbasis rb;    // Row clusterbasis
  pclusterbasis cb;    // Column clusterbasis
  ph2matrix     H2;    // H2-matrix
  pstopwatch    sw;

  size_t        sz;    // Storage requirement

  srand(time(NULL));   // Initialize RNG

  omp_set_num_threads(2);

  init_h2lib(&master_test_suite().argc, &master_test_suite().argv);

  /* Create ie1d (problem description) */
  ie = new_ie1d(n, 3);

  /* Create cluster tree */
  printf("Creating cluster tree\n");
  root = build_ie1d_cluster(ie, res);
  printf("  %u clusters\n", root->desc);

  /* Prepare ie1d objects */
  printf("Prepare for interpolation\n");
  setup_aprx_interpolation_ie1d(ie, root, k, 0.0);

  /* Create block tree */
  printf("Creating block tree\n");
  broot = build_nonstrict_block(root,
                                root,
                                (void *) &eta,
                                admissible_2_cluster);
  printf("  %u blocks\n",broot->desc);

  /* Create cluster basis*/
  printf("Creating cluster basis\n");
  rb = build_from_cluster_clusterbasis(root);
  cb = build_from_cluster_clusterbasis(root);

  /* Filling cluster basis*/
  printf("Filling cluster basis\n");
  fill_clusterbasis_ie1d(ie, rb);
  sz = getsize_clusterbasis(rb);
  printf("  Row cluster basis %.1f KB (%.1f KB/DoF)\n", sz / 1024.0,
                                                        sz / 1024.0 / n);
  fill_clusterbasis_ie1d(ie, cb);
  sz = getsize_clusterbasis(cb);
  printf("  Column cluster basis %.1f KB (%.1f KB/DoF)\n", sz / 1024.0,
                                                           sz / 1024.0 / n);

  /* Create H2-matrix*/
  printf("Creating H2-matrix\n");
  H2 = build_from_block_h2matrix(broot, rb, cb);

  /* Fill H2-matrix */
  printf("Filling H2-matrix\n");
  fill_h2matrix_ie1d(ie, H2);
  sz = getsize_h2matrix(H2);
  printf("  H2-matrix %.1f KB (%.1f KB/DoF)\n", sz / 1024.0, sz / 1024.0 / n);

  sw     = new_stopwatch();

  x      = init_avector(&locx, H2->cb->t->size);
  xt_cpu = new_coeffs_clusterbasis_avector(H2->cb);
  xt_gpu = new_coeffs_clusterbasis_avector(H2->cb);

  random_avector(x);

  del_stopwatch(sw);
  del_avector(xt_gpu);
  del_avector(xt_cpu);
  uninit_avector(x);
  del_h2matrix(H2);
  del_block(broot);
  freemem(root->idx);
  del_cluster(root);
  del_ie1d(ie);

  printf("  %u vectors, %u matrices and\n"
	       "  %u cluster basis still allocated\n", getactives_avector(),
                                                 getactives_amatrix(),
                                                 getactives_clusterbasis());

  uninit_h2lib();
}

/** @} */
