#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE main

#include <boost/test/unit_test.hpp>

#include "greencross.h"
#include "basic.h"
#include "ie1d.h"
#include "opencl.h"

using namespace boost::unit_test::framework;

static const real eta = 1.0;   // Parameter for the accuracy
static const uint k   = 4;     // Approximation order
static const uint n   = 1024;  // Problem size
static const uint res = 4 * k; // Cluster resolution

BOOST_AUTO_TEST_CASE(test_main)
{
  pavector      x;
  pavector      xt;
  pblock        broot; // Root of the block tree
  pie1d         ie;    // Problems
  pcluster      root;  // Root of the cluster tree
  pclusterbasis rb;    // Row clusterbasis
  pclusterbasis cb;    // Column clusterbasis
  ph2matrix     H2;    // H2-matrix

  size_t        sz;    // Storage requirement

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

  init_greencross_ocl();

  x  = new_avector(H2->cb->t->size);
  xt = new_avector(H2->cb->t->size);

  forward_greencross(H2->cb, x, xt);

  uninit_greencross_ocl();

  del_avector(xt);
  del_avector(x);
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
