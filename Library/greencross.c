/* ------------------------------------------------------------
 * This is the file "greencross.c" of this master thesis.
 * All rights reserved, Bennet Carstensen 2017
 * ------------------------------------------------------------ */

/**
 * @file      greencross.c
 * @author    Bennet Carstensen
 * @date      2017
 * @copyright All rights reserved, Bennet Carstensen 2017
 */

#include "greencross.h"

#include "aca.h"
#include "bem2d.h"
#include "bem3d.h"
#include "gaussquad.h"
#include "opencl.h"
#include "singquad1d.h"
#include "singquad2d.h"

#include <string.h>

/** @addtogroup greencross
 *
 * @{*/

struct _gcopencl
{
  ph2matrix *H2n;

  uint num_row_leafs;

  uint *names_row_leafs;

  uint *num_h2matrix_leafs_per_cluster;

  ph2matrix **h2matrix_leafs_per_cluster;

  uint      **buf_ridx_off;

  cl_mem    *buf_ridx;

  uint      **buf_cidx_off;

  cl_mem    *buf_cidx;

  uint      *xtoffs;

  uint      **ytoffs;

  pcluster  *row_clusters;

  uint      *current_h2matrix_leafs;

  uint roff;
  uint coff;
};

static void
iterate_h2matrix_serial(ph2matrix G,
                        uint      mname,
                        uint      rname,
                        uint      cname,
		                    uint      pardepth,
		                    void      (*pre)  (ph2matrix G,
                                           uint      mname,
                                           uint      rname,
                                           uint      cname,
			                                     uint      pardepth,
                                           void      *data),
		                    void      (*post) (ph2matrix G,
                                           uint      mname,
                                           uint      rname,
			                                     uint      cname,
                                           uint      pardepth,
                                           void      *data),
		                    void      *data);

/** @} */

/* ------------------------------------------------------------
 * Source of static functions.
 * ------------------------------------------------------------ */

void
iterate_h2matrix_serial(ph2matrix G,
                        uint      mname,
                        uint      rname,
                        uint      cname,
		                    uint      pardepth,
		                    void      (*pre)  (ph2matrix G,
                                           uint      mname,
                                           uint      rname,
                                           uint      cname,
			                                     uint      pardepth,
                                           void      *data),
		                    void      (*post) (ph2matrix G,
                                           uint      mname,
                                           uint      rname,
			                                     uint      cname,
                                           uint      pardepth,
                                           void      *data),
		                    void      *data)
{
  pccluster rc, cc;
  uint      rsons, csons;
  ph2matrix *Gn;
  uint     *mnamen, *rnamen, *cnamen;
  uint      mname1, rname1, cname1;

  uint      i, j, k, l;

  /* Call a priori callback function */
  if (pre)
    pre(G, mname, rname, cname, pardepth, data);

  /* Handle sons */
  if (G->son) {
    rc = G->rb->t;
    cc = G->cb->t;

    rsons = G->rsons;
    csons = G->csons;

    if (G->rb == G->son[0]->rb) {
      if (G->cb == G->son[0]->cb) {
	/* Just one son */
	iterate_h2matrix(G->son[0], mname + 1, rname, cname,
			 (pardepth > 0 ? pardepth - 1 : 0), pre, post, data);
      }
      else {
	/* No row son, but column sons */
	mname1 = mname + 1;
	cname1 = cname + 1;
	for (j = 0; j < csons; j++) {
	  iterate_h2matrix(G->son[j * rsons], mname1, rname, cname1,
			   (pardepth > 0 ? pardepth - 1 : 0), pre, post,
			   data);
	  mname1 += G->son[j * rsons]->desc;
	  cname1 += cc->son[j]->desc;
	}
	assert(mname1 == mname + G->desc);
	assert(cname1 == cname + cc->desc);
      }
    }
    else {
      if (G->cb == G->son[0]->cb) {
	/* No column son, but row sons */
	mname1 = mname + 1;
	rname1 = rname + 1;
	for (i = 0; i < rsons; i++) {
	  iterate_h2matrix(G->son[i], mname1, rname1, cname,
			   (pardepth > 0 ? pardepth - 1 : 0), pre, post,
			   data);
	  mname1 += G->son[i]->desc;
	  rname1 += rc->son[i]->desc;
	}
	assert(mname1 == mname + G->desc);
	assert(rname1 == rname + rc->desc);
      }
      else {
	/* Row and column sons, so we can parallelize something */
	if (rsons >= csons) {
	  /* Allocate arrays for submatrices and names */
	  Gn = (ph2matrix *) allocmem(sizeof(ph2matrix) * csons);
	  mnamen = (uint *) allocmem(sizeof(uint) * csons);
	  rnamen = (uint *) allocmem(sizeof(uint) * csons);
	  cnamen = (uint *) allocmem(sizeof(uint) * csons);

	  /* Consider all parallelizable subsets, in this case we use
	     row-wise shifted cyclic diagonals */
	  for (k = 0; k < rsons; k++) {
	    cname1 = cname + 1;

	    for (j = 0; j < csons; j++) {
	      i = (k + j) % rsons;

	      Gn[j] = G->son[i + j * rsons];

	      /* Compute submatrix number */
	      mname1 = mname + 1;
	      for (l = 0; l < i + j * rsons; l++)
		mname1 += G->son[l]->desc;

	      /* Compute row number */
	      rname1 = rname;
	      if (G->rb != G->son[0]->rb) {
		rname1 = rname + 1;
		for (l = 0; l < i; l++)
		  rname1 += rc->son[l]->desc;
	      }

	      mnamen[j] = mname1;
	      rnamen[j] = rname1;
	      cnamen[j] = cname1;

	      cname1 += cc->son[j]->desc;
	    }
	    assert(cname1 == cname + cc->desc);

	    /* Recursive call */
	    for (j = 0; j < csons; j++)
	      iterate_h2matrix(Gn[j], mnamen[j], rnamen[j], cnamen[j],
			       (pardepth > 0 ? pardepth - 1 : 0), pre, post,
			       data);
	  }

	  /* Clean up */
	  freemem(cnamen);
	  freemem(rnamen);
	  freemem(mnamen);
	  freemem(Gn);
	}
	else {
	  /* Allocate arrays for submatrices and names */
	  Gn = (ph2matrix *) allocmem(sizeof(ph2matrix) * rsons);
	  mnamen = (uint *) allocmem(sizeof(uint) * rsons);
	  rnamen = (uint *) allocmem(sizeof(uint) * rsons);
	  cnamen = (uint *) allocmem(sizeof(uint) * rsons);

	  /* Consider all parallelizable subsets, in this case we use
	     column-wise shifted cyclic diagonals */
	  for (k = 0; k < csons; k++) {
	    rname1 = rname + 1;

	    for (i = 0; i < rsons; i++) {
	      j = (k + i) % csons;

	      Gn[i] = G->son[i + j * rsons];

	      /* Compute submatrix number */
	      mname1 = mname + 1;
	      for (l = 0; l < i + j * rsons; l++)
		mname1 += G->son[l]->desc;

	      /* Compute colum number */
	      cname1 = cname;
	      if (G->cb != G->son[0]->cb) {
		cname1 = cname + 1;
		for (l = 0; l < j; l++)
		  cname1 += cc->son[l]->desc;
	      }

	      mnamen[i] = mname1;
	      rnamen[i] = rname1;
	      cnamen[i] = cname1;

	      rname1 += rc->son[i]->desc;
	    }
	    assert(rname1 == rname + rc->desc);

	    /* Recursive call */
	    for (j = 0; j < rsons; j++)
	      iterate_h2matrix(Gn[j], mnamen[j], rnamen[j], cnamen[j],
			       (pardepth > 0 ? pardepth - 1 : 0), pre, post,
			       data);
	  }

	  /* Clean up */
	  freemem(cnamen);
	  freemem(rnamen);
	  freemem(mnamen);
	  freemem(Gn);
	}
      }
    }
  }

  /* Call a posteriori callback function */
  if (post)
    post(G, mname, rname, cname, pardepth, data);
}

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

void
init_gcopencl(pgcopencl gcocl)
{
  assert(gcocl != NULL);

  gcocl->num_row_leafs = 0;
  gcocl->roff          = 0;
  gcocl->coff          = 0;

  gcocl->H2n                            = NULL;
  gcocl->names_row_leafs                = NULL;
  gcocl->num_h2matrix_leafs_per_cluster = NULL;
  gcocl->h2matrix_leafs_per_cluster     = NULL;
  gcocl->buf_ridx_off                   = NULL;
  gcocl->buf_ridx                       =
    (cl_mem *) allocmem(ocl_system.num_devices * sizeof(cl_mem));
  gcocl->buf_cidx_off                   = NULL;
  gcocl->buf_cidx                       =
    (cl_mem *) allocmem(ocl_system.num_devices * sizeof(cl_mem));
  gcocl->xtoffs                         = NULL;
  gcocl->ytoffs                         = NULL;
  gcocl->row_clusters                   = NULL;
  gcocl->current_h2matrix_leafs         = NULL;
}

void
uninit_gcopencl(pgcopencl gcocl)
{
  if(gcocl->H2n != NULL)
  {
    freemem(gcocl->H2n);
    gcocl->H2n = NULL;
  }

  if(gcocl->names_row_leafs != NULL)
  {
    freemem(gcocl->names_row_leafs);
    gcocl->names_row_leafs = NULL;
  }

  if(gcocl->num_h2matrix_leafs_per_cluster != NULL)
  {
    freemem(gcocl->num_h2matrix_leafs_per_cluster);
    gcocl->num_h2matrix_leafs_per_cluster = NULL;
  }

  if(gcocl->h2matrix_leafs_per_cluster != NULL)
  {
    for(uint i = 0; i < gcocl->num_row_leafs; ++i)
      if(gcocl->h2matrix_leafs_per_cluster[i] != NULL)
      {
        freemem(gcocl->h2matrix_leafs_per_cluster[i]);
        gcocl->h2matrix_leafs_per_cluster[i] = NULL;
      }

    freemem(gcocl->h2matrix_leafs_per_cluster);
    gcocl->h2matrix_leafs_per_cluster = NULL;
  }

  if(gcocl->buf_ridx_off != NULL)
  {
    for(uint i = 0; i < gcocl->num_row_leafs; ++i)
      if(gcocl->buf_ridx_off[i] != NULL)
      {
        freemem(gcocl->buf_ridx_off[i]);
        gcocl->buf_ridx_off[i] = NULL;
      }

    freemem(gcocl->buf_ridx_off);
    gcocl->buf_ridx_off = NULL;
  }

  if(gcocl->buf_ridx != NULL)
  {
    freemem(gcocl->buf_ridx);
    gcocl->buf_ridx = NULL;
  }

  if(gcocl->buf_cidx_off != NULL)
  {
    for(uint i = 0; i < gcocl->num_row_leafs; ++i)
      if(gcocl->buf_cidx_off[i] != NULL)
      {
        freemem(gcocl->buf_cidx_off[i]);
        gcocl->buf_cidx_off[i] = NULL;
      }

    freemem(gcocl->buf_cidx_off);
    gcocl->buf_cidx_off = NULL;
  }

  if(gcocl->buf_cidx != NULL)
  {
    freemem(gcocl->buf_cidx);
    gcocl->buf_cidx = NULL;
  }

  if(gcocl->xtoffs != NULL)
  {
    freemem(gcocl->xtoffs);
    gcocl->xtoffs = NULL;
  }

  if(gcocl->ytoffs != NULL)
  {
    for(uint i = 0; i < gcocl->num_row_leafs; ++i)
      if(gcocl->ytoffs[i] != NULL)
      {
        freemem(gcocl->ytoffs[i]);
        gcocl->ytoffs[i] = NULL;
      }

    freemem(gcocl->ytoffs);
    gcocl->ytoffs = NULL;
  }

  if(gcocl->row_clusters != NULL)
  {
    freemem(gcocl->row_clusters);
    gcocl->row_clusters = NULL;
  }

  if(gcocl->current_h2matrix_leafs != NULL)
  {
    freemem(gcocl->current_h2matrix_leafs);
    gcocl->current_h2matrix_leafs = NULL;
  }
}

void
init_greencross(pgreencross gc, uint dim)
{
  assert(gc != NULL);
  assert((dim >= greencross_min_dim) && (greencross_max_dim <= 3));

  /* Scalar components */
  gc->dim                 = dim;
  gc->n                   = 0;
  gc->m                   = 0;
  gc->K                   = 0;

  gc->gcocl               = (pgcopencl) allocmem(sizeof(gcopencl));

  init_gcopencl(gc->gcocl);

  /* Pointer components */
  gc->geom                = NULL;
  gc->bem                 = NULL;
  gc->idx                 = NULL;
  gc->rb                  = NULL;
  gc->cb                  = NULL;
  gc->rc                  = NULL;
  gc->cc                  = NULL;
}

void
uninit_greencross(pgreencross gc)
{
  assert(gc != NULL);

  if(gc->geom != NULL)
  {
    gc->dim == 2 ? del_curve2d((pcurve2d) gc->geom)
                 : del_surface3d((psurface3d) gc->geom);
    gc->geom = NULL;
  }

  if(gc->bem != NULL)
  {
    gc->dim == 2 ? del_bem2d((pbem2d) gc->bem)
                 : del_bem3d((pbem3d) gc->bem);
    gc->bem = NULL;
  }

  if(gc->idx != NULL)
  {
    freemem(gc->idx);
    gc->idx = NULL;
  }

  if(gc->rc == gc->cc)
    gc->cc = NULL; // Only delete cluster once if row cluster = column cluster.

  if(gc->rc != NULL)
  {
    del_cluster(gc->rc);
    gc->rc = NULL;
  }

  if(gc->cc != NULL)
  {
    del_cluster(gc->cc);
    gc->cc = NULL;
  }

  uninit_gcopencl(gc->gcocl);
}

void
del_greencross(pgreencross gc)
{
  assert(gc != NULL);

  uninit_greencross(gc);

  freemem(gc);

  gc = NULL;
}

pclustergeometry
build_clustergeometry_greencross(const void *data, const uint dim, uint **idx)
{
  pclustergeometry cg;

  real (*x)[dim];
  uint (*p)[dim];
  uint n = 0;

  if(dim == 2)
  {
    pccurve2d gr = (pccurve2d) data;

    x  = gr->x;
    p  = gr->e;

    n  = gr->edges;
  }
  else if(dim == 3)
  {
    pcsurface3d gr = (pcsurface3d) data;

    x = gr->x;
    p = gr->t;

    n = gr->triangles;
  }
  else
  {
    printf("Unknown dimension for triangular mesh: %u\n", dim);
    exit(1);
  }

  cg   = new_clustergeometry(dim, n);
  *idx = allocuint(n);

  for(uint i = 0; i < n; ++i)
  {
    (*idx)[i] = i;

    for(uint d = 0; d < dim; ++d)
    {
      /* Center of gravity as characteristic point */
      cg->x[i][d] = r_zero;
      for(uint j = 0; j < dim; ++j)
        cg->x[i][d] += x[p[i][j]][d];
      cg->x[i][d] *= 0.5;

      /* Lower front left corner of bounding box */
      cg->smin[i][d] =
        dim == 2 ? REAL_MIN(x[p[i][0]][d], x[p[i][1]][d])
                 : REAL_MIN3(x[p[i][0]][d], x[p[i][1]][d], x[p[i][2]][d]);

      /* Upper back right corner of bounding box */
      cg->smax[i][d] =
        dim == 2 ? REAL_MAX(x[p[i][0]][d], x[p[i][1]][d])
                 : REAL_MAX3(x[p[i][0]][d], x[p[i][1]][d], x[p[i][2]][d]);
    }
  }

  return cg;
}

void
assemble_amatrix_greencross(pcgreencross gc,
                            const uint   *ridx,
                            const uint   *cidx,
                            bool         ntrans,
                            pamatrix     G)
{
  if(gc->dim == 2)
    ((pbem2d) gc->bem)->nearfield(ridx, cidx, (pbem2d) gc->bem, ntrans, G);
  else
    ((pbem3d) gc->bem)->nearfield(ridx, cidx, (pbem3d) gc->bem, ntrans, G);
}

pamatrix
build_amatrix_greencross(pcgreencross gc,
                         const uint rows,
                         const uint *ridx,
                         const uint cols,
                         const uint *cidx,
                         bool       ntrans)
{
  pamatrix G = new_amatrix(rows, cols);

  assemble_amatrix_greencross(gc, ridx, cidx, ntrans, G);

  return G;
}

prkmatrix
build_green_rkmatrix_greencross(pcgreencross gc, pccluster row, pccluster col)
{
  prkmatrix AB;

  if(gc->dim == 2)
  {
    pbem2d b2d = (pbem2d) gc->bem;

    setup_hmatrix_aprx_green_row_bem2d(b2d,
                                       NULL,
                                       NULL,
                                       NULL,
                                       gc->m,
                                       1,
                                       0.5,
    				                           build_bem2d_rect_quadpoints);

    AB = build_bem2d_rkmatrix(row, col, gc->bem);

    setup_h2matrix_aprx_greenhybrid_bem2d(b2d,
                                          gc->rb,
                                          gc->cb,
                                          NULL,
                                          gc->m,
                                          1,
                                          0.5,
                                          10e-07,
                                          build_bem2d_rect_quadpoints);
  }
  else
  {
    pbem3d b3d = (pbem3d) gc->bem;

    setup_hmatrix_aprx_green_row_bem3d(b3d,
                                       NULL,
                                       NULL,
                                       NULL,
                                       gc->m,
                                       1,
                                       0.5,
                                       build_bem3d_cube_quadpoints);

    AB = build_bem3d_rkmatrix(row, col, gc->bem);

    setup_h2matrix_aprx_greenhybrid_bem3d(b3d,
                                          gc->rb,
                                          gc->cb,
                                          NULL,
                                          gc->m,
                                          1,
                                          0.5,
                                          10e-07,
                                          build_bem3d_cube_quadpoints);
  }

  return AB;
}

phmatrix
build_green_hmatrix_greencross(pgreencross gc, void *eta)
{
  pblock   broot;
  phmatrix H;

  broot = build_nonstrict_block(gc->rc, gc->cc, eta, admissible_max_cluster);

  H     = build_from_block_hmatrix(broot, gc->K);

  if(gc->dim == 2)
  {
    pbem2d b2d = (pbem2d) gc->bem;

    setup_hmatrix_aprx_green_row_bem2d(b2d,
                                       NULL,
                                       NULL,
    				                           NULL,
                                       gc->m,
                                       1,
                                       0.5,
    				                           build_bem2d_rect_quadpoints);

    assemble_bem2d_hmatrix(b2d, broot, H);

    setup_h2matrix_aprx_greenhybrid_bem2d(b2d,
                                          gc->rb,
                                          gc->cb,
                                          NULL,
                                          gc->m,
                                          1,
                                          0.5,
                                          10e-07,
                                          build_bem2d_rect_quadpoints);
  }
  else
  {
    pbem3d b3d = (pbem3d) gc->bem;

    setup_hmatrix_aprx_green_row_bem3d(b3d,
                                       NULL,
                                       NULL,
                                       NULL,
                                       gc->m,
                                       1,
                                       0.5,
                                       build_bem3d_cube_quadpoints);

    assemble_bem3d_hmatrix(b3d, broot, H);

    setup_h2matrix_aprx_greenhybrid_bem3d(b3d,
                                          gc->rb,
                                          gc->cb,
                                          NULL,
                                          gc->m,
                                          1,
                                          0.5,
                                          10e-07,
                                          build_bem3d_cube_quadpoints);
  }

  del_block(broot);

  return H;
}

static void
get_leaf_block_informations(ph2matrix H2,
                            uint      mname,
                            uint      rname,
                            uint      cname,
                            uint      pardepth,
                            void      *data)
{
  (void) mname;
  (void) cname;
  (void) pardepth;

  pgcopencl gcocl = (pgcopencl) data;

  if(H2->u || H2->f)
  {
    bool row_already_stored = false;

    for(uint i = 0; i < gcocl->num_row_leafs; ++i)
      if(gcocl->names_row_leafs[i] == rname)
      {
        row_already_stored = true;

        // Update H2-matrix leaf counter for this cluster
        gcocl->num_h2matrix_leafs_per_cluster[i] += 1;

        // Add new H2-matrix leafs entry for this cluster
        gcocl->h2matrix_leafs_per_cluster[i] =
          realloc(gcocl->h2matrix_leafs_per_cluster[i],
                  gcocl->num_h2matrix_leafs_per_cluster[i] * sizeof(ph2matrix));

        assert(gcocl->h2matrix_leafs_per_cluster[i] != NULL);

        // Set H2-matrix leafs entry for this cluster
        gcocl->h2matrix_leafs_per_cluster[i][gcocl->num_h2matrix_leafs_per_cluster[i] - 1] = H2;

        gcocl->buf_ridx_off[i] =
          realloc(gcocl->buf_ridx_off[i],
                  gcocl->num_h2matrix_leafs_per_cluster[i] * sizeof(uint));

        assert(gcocl->buf_ridx_off[i] != NULL);

        gcocl->buf_ridx_off[i][gcocl->num_h2matrix_leafs_per_cluster[i] - 1]
          = gcocl->roff;

        gcocl->roff += H2->rb->t->size;

        break;
      }

    if(!row_already_stored)
    {
      cl_int res;

      // Update number of leafs
      gcocl->num_row_leafs += 1;

      // Add new name entry for this cluster
      gcocl->names_row_leafs = realloc(gcocl->names_row_leafs,
                                       gcocl->num_row_leafs * sizeof(uint));

      assert(gcocl->names_row_leafs != NULL);

      // Set new name entry for this cluster
      gcocl->names_row_leafs[gcocl->num_row_leafs - 1] = rname;

      // Add new H2-matrix leaf counter entry for this cluster
      gcocl->num_h2matrix_leafs_per_cluster =
        realloc(gcocl->num_h2matrix_leafs_per_cluster,
                gcocl->num_row_leafs * sizeof(uint));

      assert(gcocl->num_h2matrix_leafs_per_cluster != NULL);

      // Set new H2-matrix leaf counter entry for this cluster
      gcocl->num_h2matrix_leafs_per_cluster[gcocl->num_row_leafs - 1] = 1;

      // Add new pointer to H2-matrix leafs array for this cluster
      gcocl->h2matrix_leafs_per_cluster =
        realloc(gcocl->h2matrix_leafs_per_cluster,
                gcocl->num_row_leafs * sizeof(ph2matrix*));

      assert(gcocl->h2matrix_leafs_per_cluster != NULL);

      // Add new H2-matrix leafs array for this cluster with one entry
      gcocl->h2matrix_leafs_per_cluster[gcocl->num_row_leafs - 1] =
        (ph2matrix*) allocmem(sizeof(ph2matrix));

      // Set H2-matrix leafs entry
      gcocl->h2matrix_leafs_per_cluster[gcocl->num_row_leafs - 1][0] = H2;

      // Add new pointer to row offset array for this H2-matrx
      gcocl->buf_ridx_off = realloc(gcocl->buf_ridx_off,
                                    gcocl->num_row_leafs * sizeof(uint*));

      assert(gcocl->buf_ridx_off != NULL);

      // Add new row offset array for this matrx with one entry
      gcocl->buf_ridx_off[gcocl->num_row_leafs - 1] =
        (uint*) allocmem(sizeof(uint));

      // Set row offset entry
      gcocl->buf_ridx_off[gcocl->num_row_leafs - 1][0] = gcocl->roff;

      res = clEnqueueWriteBuffer(ocl_system.queues[0],
                                 gcocl->buf_ridx[0],
                                 false,
                                 gcocl->roff,
                                 H2->rb->t->size,
                                 H2->rb->t->idx,
                                 0,
                                 NULL,
                                 NULL);

      CL_CHECK(res);

      // Update row offset for next leaf H2-matrix
      gcocl->roff += H2->rb->t->size;

      // Add new pointer to column offset array for this H2-matrx
      gcocl->buf_cidx_off = realloc(gcocl->buf_cidx_off,
                                    gcocl->num_row_leafs * sizeof(uint*));

      assert(gcocl->buf_cidx_off != NULL);

      // Add new column offset array for this matrx with one entry
      gcocl->buf_cidx_off[gcocl->num_row_leafs - 1] =
        (uint*) allocmem(sizeof(uint));

      // Set column offset entry
      gcocl->buf_cidx_off[gcocl->num_row_leafs - 1][0] = gcocl->coff;

      res = clEnqueueWriteBuffer(ocl_system.queues[0],
                                 gcocl->buf_cidx[0],
                                 false,
                                 gcocl->coff,
                                 H2->cb->t->size,
                                 H2->cb->t->idx,
                                 0,
                                 NULL,
                                 NULL);

      CL_CHECK(res);

      // Update column offset for next leaf H2-matrix
      gcocl->coff += H2->cb->t->size;
    }
  }
}

static bool
are_equal_clusters(pccluster a, pccluster b)
{
  if(a->size != b->size)
    return false;

  for(uint i = 0; i < a->size; ++i)
    if(a->idx[i] != b->idx[i])
      return false;

  return true;
}

static void
iterate_recursively_h2matrix(ph2matrix H2,
                             uint      xtoff,
                             uint      ytoff,
                             pgcopencl gcocl)
{
  if(H2->u || H2->f)
  {
    int  i = -1;
    uint j;

    for(uint k = 0; k < gcocl->num_row_leafs; ++k)
      if(are_equal_clusters(H2->rb->t, gcocl->row_clusters[gcocl->names_row_leafs[k]]))
      {
        i = k;
        break;
      }

    assert(i >= 0);

    j = gcocl->current_h2matrix_leafs[i];

    gcocl->xtoffs[i]    = xtoff;
    gcocl->ytoffs[i][j] = ytoff; // TODO: Change to be similar to cidxs

    gcocl->current_h2matrix_leafs[i] += 1;
  }
  else
  {
    pcclusterbasis rb = H2->rb;
    pcclusterbasis cb = H2->cb;

    const uint rsons  = H2->rsons;
    const uint csons  = H2->csons;

    uint xtoff1 = cb->k;

    for(uint j = 0; j < csons; ++j)
    {
      assert(csons == 1 || cb->sons > 0);

      uint ytoff1 = rb->k;

      for(uint i = 0; i < rsons; ++i)
      {
        assert(rsons == 1 || rb->sons > 0);

        iterate_recursively_h2matrix(H2->son[i + j * rsons],
                                     xtoff + xtoff1,
                                     ytoff + ytoff1,
                                     gcocl);

       ytoff1 += (rb->sons > 0 ? rb->son[i]->ktree : rb->t->size);
      }

      assert(ytoff1 == rb->ktree);

      xtoff1 += (cb->sons > 0 ? cb->son[j]->ktree : cb->t->size);
    }

    assert(xtoff1 == cb->ktree);
  }
}

static void
get_ocl_informations_gcocl(ph2matrix H2, pgcopencl gcocl)
{
  cl_int res;

  gcocl->num_row_leafs                  = 0;
  gcocl->roff                           = 0;
  gcocl->coff                           = 0;
  gcocl->names_row_leafs                = (uint *) calloc(0, sizeof(uint));
  gcocl->num_h2matrix_leafs_per_cluster = (uint *) calloc(0, sizeof(uint));
  gcocl->h2matrix_leafs_per_cluster     = (ph2matrix**) calloc
                                                          (0,
                                                           sizeof(ph2matrix*));
  gcocl->buf_ridx_off                   = (uint **) calloc(0, sizeof(uint*));
  gcocl->buf_ridx[0]                    = clCreateBuffer
                                            (ocl_system.contexts[0],
                                             CL_MEM_READ_WRITE,
                                             ocl_system.max_package_size,
                                             NULL,
                                             &res);
  gcocl->buf_cidx_off                   = (uint **) calloc(0, sizeof(uint*));
  gcocl->buf_cidx[0]                    = clCreateBuffer
                                            (ocl_system.contexts[0],
                                             CL_MEM_READ_WRITE,
                                             ocl_system.max_package_size,
                                             NULL,
                                             &res);

  CL_CHECK(res);

  iterate_h2matrix_serial(H2,
                          0,
                          0,
                          0,
                          0,
                          get_leaf_block_informations,
                          NULL,
                          (void*) gcocl);

  gcocl->row_clusters           = enumerate_cluster((pcluster) H2->rb->t);

  gcocl->current_h2matrix_leafs =
    (uint *) calloc(gcocl->num_row_leafs, sizeof(uint));

  gcocl->xtoffs                 = (uint *)  calloc(gcocl->num_row_leafs,
                                                   sizeof(uint));
  gcocl->ytoffs                 = (uint **) calloc(gcocl->num_row_leafs,
                                                   sizeof(uint*));

  for(uint i = 0; i < gcocl->num_row_leafs; ++i)
    gcocl->ytoffs[i]            =
      (uint *)  calloc(gcocl->num_h2matrix_leafs_per_cluster[i], sizeof(uint));

  iterate_recursively_h2matrix(H2, 0, 0, gcocl);

  freemem(gcocl->current_h2matrix_leafs);
  gcocl->current_h2matrix_leafs = NULL;

  freemem(gcocl->row_clusters);
  gcocl->row_clusters = NULL;
}

ph2matrix
build_green_cross_h2matrix_greencross(pgreencross gc, void *eta)
{
  pblock    broot;
  ph2matrix H2;

  broot = build_strict_block(gc->rc, gc->cc, eta, admissible_2_cluster);

  H2    = build_from_block_h2matrix(broot,
                                    gc->rb,
                                    gc->cb);

  if(gc->dim == 2)
  {
    pbem2d bem = (pbem2d) gc->bem;

    assemble_bem2d_h2matrix_row_clusterbasis(bem, H2->rb);
    assemble_bem2d_h2matrix_col_clusterbasis(bem, H2->cb);
    assemble_bem2d_h2matrix(bem, broot, H2);
  }
  else
  {
    pbem3d bem = (pbem3d) gc->bem;

    assemble_bem3d_h2matrix_row_clusterbasis(bem, H2->rb);
    assemble_bem3d_h2matrix_col_clusterbasis(bem, H2->cb);
    assemble_bem3d_h2matrix(bem, H2);
  }

  get_ocl_informations_gcocl(H2, gc->gcocl);

  del_block(broot);

  return H2;
}

void
fastaddeval_h2matrix_avector_greencross(field      alpha,
                                        pch2matrix H2,
                                        pavector   xt,
			                                  pavector   yt)
{
  (void) alpha;
  (void) H2;
  (void) xt;
  (void) yt;
}

void
addeval_h2matrix_avector_greencross(field      alpha,
                                    pch2matrix H2,
                                    pcavector  x,
                                    pavector   y)
{
  pavector  xt, yt;

  xt = new_coeffs_clusterbasis_avector(H2->cb);
  yt = new_coeffs_clusterbasis_avector(H2->rb);

  clear_avector(yt);

  forward_clusterbasis_avector(H2->cb, x, xt);

  fastaddeval_h2matrix_avector(alpha, H2, xt, yt);

  backward_clusterbasis_avector(H2->rb, yt, y);

  del_avector(yt);
  del_avector(xt);
}

void
addevaltrans_h2matrix_avector_greencross(field      alpha,
                                         pch2matrix H2,
                                         pcavector  x,
			                                   pavector   y)
{
  pavector  xt, yt;

  xt = new_coeffs_clusterbasis_avector(H2->rb);
  yt = new_coeffs_clusterbasis_avector(H2->cb);

  clear_avector(yt);

  forward_clusterbasis_avector(H2->rb, x, xt);

  fastaddevaltrans_h2matrix_avector(alpha, H2, xt, yt);

  backward_clusterbasis_avector(H2->cb, yt, y);

  del_avector(yt);
  del_avector(xt);
}

void
mvm_h2matrix_avector_greencross(field      alpha,
                                bool       h2trans,
                                pch2matrix H2,
                                pcavector  x,
		                            pavector   y)
{
  if (h2trans)
    addevaltrans_h2matrix_avector(alpha, H2, x, y);
  else
    addeval_h2matrix_avector(alpha, H2, x, y);
}
