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

#include "bem2d.h"
#include "ocl_system.h"

#include "clfastaddevalgca.cl"
#include "clgcidxinfo.cl"
#include "clgeom.cl"
#include "clgreencross.cl"
#include "clkernels.cl"
#include "clintegralinfo.cl"
#include "clsingquad.cl"

#include <string.h>

const char *src_code_strs[] = { clgcidxinfo_src,
                                clgeom_src,
                                clkernels_src,
                                clsingquad_src,
                                clfastaddevalgca_src,
                                clintegralinfo_src,
                                clgreencross_src };

const char *kernel_names[]  = { "fastaddeval_h2matrix_avector_0",
                                "fastaddeval_nf_common",
                                "fastaddeval_nf_min_vert",
                                "fastaddeval_nf_min_edge"
                              };

/** @addtogroup greencross
 *
 * @{*/

/*
 * Just an abbreviation for the struct _greencluster3d .
 */
typedef struct _greencluster3d greencluster3d;

/*
 * Pointer to a @ref greencluster3d object.
 */
typedef greencluster3d *pgreencluster3d;

/*
 * Just an abbreviation for the struct _greencluster3d .
 */
typedef struct _greenclusterbasis3d greenclusterbasis3d;

/*
 * Pointer to a @ref greenclusterbasis3d object.
 */
typedef greenclusterbasis3d *pgreenclusterbasis3d;

struct _greenclusterbasis3d {
  uint     *xi;
  /** local indices of pivot elements */
  uint     *xihat;
  /** global indices of pivot elements */
  pamatrix  Qinv;
  /** Triangular factor of QR-decomposition */
  pcclusterbasis cb; /** corresponding clusterbasis */
  uint      sons;
  /** number of sons for current clusterbasis */
  uint      m;
};

struct _parbem3d {
  /*
   * special members for different H- and H2-matrix approximations
   */

  phmatrix *hn;			/* temporary enumerated list of hmatrices */
  ph2matrix *h2n;		/* temporary enumerated list of h2matrices */
  pdh2matrix *dh2n;		/* temporary enumerated list of h2matrices */
  pclusterbasis *rbn;
  pclusterbasis *cbn;
  pdclusterbasis *drbn;
  pdclusterbasis *dcbn;
  pclusteroperator *rwn;	/* temporary enumerated list of clusteroperators for row-cluster */
  pclusteroperator *cwn;	/* temporary enumerated list of clusteroperators for col-cluster */
  pdclusteroperator *ron;	/* temporary enumerated list of dclusteroperators for row-cluster */
  pdclusteroperator *con;	/* temporary enumerated list of dclusteroperators for col-cluster */
  uint     *leveln;		/* temporary list of levelnumber for each block in blocktree. */
  pgreencluster3d *grcn;
  uint      grcnn;
  pgreencluster3d *gccn;
  uint      gccnn;
  pgreenclusterbasis3d *grbn;
  uint      grbnn;
  pgreenclusterbasis3d *gcbn;
  uint      gcbnn;
};

/** @} */

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

void
init_greencross(pgreencross gc, uint dim)
{
  assert(gc != NULL);
  assert((dim >= greencross_min_dim) && (dim <= greencross_max_dim));

  /* Scalar components */
  gc->dim        = dim;
  gc->n          = 0;
  gc->m          = 0;
  gc->K          = 0;
  gc->aca_accur  = 0;

  gc->ocl_info_nf = NULL;

  /* Pointer components */
  gc->geom       = NULL;
  gc->buf_x      =
    (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));
  gc->buf_p      =
    (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));
  gc->event_p    =
    (cl_event *) calloc(ocl_system.num_devices, sizeof(cl_event));
  gc->buf_g      =
    (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));
  gc->bem        = NULL;
  gc->sq_gca     = NULL;
  gc->sq_partial_min_vert = NULL;
  gc->sq_partial_min_edge = NULL;
  gc->feval      = NULL;
  gc->kernel_3d  = NULL;
#ifdef USE_SIMD
  gc->kernel_simd_3d = NULL;
#endif
  gc->idx             = NULL;
  gc->rb              = NULL;
  gc->cb              = NULL;
  gc->rc              = NULL;
  gc->cc              = NULL;
  gc->ocl_kernels     = NULL;
  gc->kernel_events   = NULL;
  gc->gcocl           = NULL;
  gc->ocl_info_nf     = NULL;
  gc->iinfos_min_edge = NULL;
  gc->oclwrk          = NULL;
  gc->ocl_wrk_nf      = NULL;
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

  if(gc->buf_x != NULL)
  {
    for(uint i = 0; i < ocl_system.num_devices; ++i)
      if(gc->buf_x[i] != NULL)
      {
        CL_CHECK(clReleaseMemObject(gc->buf_x[i]));
        gc->buf_x[i] = NULL;
      }

    freemem(gc->buf_x);
    gc->buf_x = NULL;
  }

  if(gc->buf_p != NULL)
  {
    for(uint i = 0; i < ocl_system.num_devices; ++i)
      if(gc->buf_p[i] != NULL)
      {
        CL_CHECK(clReleaseMemObject(gc->buf_p[i]));
        gc->buf_p[i] = NULL;
      }

    freemem(gc->buf_p);
    gc->buf_p = NULL;
  }

  if(gc->event_p != NULL)
  {
    for(uint i = 0; i < ocl_system.num_devices; ++i)
      if(gc->event_p[i] != NULL)
      {
        CL_CHECK(clReleaseEvent(gc->event_p[i]));
        gc->event_p[i] = NULL;
      }

    freemem(gc->event_p);
  }

  if(gc->buf_g != NULL)
  {
    for(uint i = 0; i < ocl_system.num_devices; ++i)
      if(gc->buf_g[i] != NULL)
      {
        CL_CHECK(clReleaseMemObject(gc->buf_g[i]));
        gc->buf_g[i] = NULL;
      }

    freemem(gc->buf_g);
    gc->buf_g = NULL;
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

  if(gc->sq_gca != NULL)
    del_singquadgca(gc->sq_gca);

  if(gc->sq_partial_min_vert != NULL)
    del_singquadgca(gc->sq_partial_min_vert);

  if(gc->sq_partial_min_edge != NULL)
    del_singquadgca(gc->sq_partial_min_edge);

  if(gc->feval != NULL)
    del_fastaddevalgca(gc->feval);

/*  if(delete_kernels(num_kernels, &gc->kernels) != NULL)
    fprintf(stderr, "warning: failed to delete kernels");*/

  freemem(gc->ocl_kernels);

  if(gc->kernel_events != NULL)
  {
    for(uint i = 0; i < num_kernels; ++i)
      clReleaseEvent(gc->kernel_events[i]);

    freemem(gc->kernel_events);
  }
  del_gcopencl(gc->gcocl);
  gc->gcocl = NULL;

  del_gcopencl(gc->ocl_info_nf);

//  uninit_gcopencl(gc->ocl_gca_nf);
//  freemem(gc->ocl_gca_nf);
//  gc->ocl_gca_nf = NULL;

  //del_oclwrk(gc->oclwrk);
  gc->oclwrk = NULL;
}

void
del_greencross(pgreencross gc)
{
  assert(gc != NULL);

  uninit_greencross(gc);

  freemem(gc);
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
                                          gc->aca_accur,
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
                                          gc->aca_accur,
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
                                          gc->aca_accur,
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
  (void) pardepth;

  pgreencross gc = (pgreencross) data;

  pgcopencl gcocl = gc->gcocl;

  if(H2->u)
  {
    bool row_already_stored = false;

    for(uint i = 0; i < gcocl->num_row_leafs; ++i)
      if(gcocl->names_row_leafs[i] == rname)
      {
        row_already_stored = true;

        bool column_already_used = false;

        // Update H2-matrix leaf counter for this cluster
        gcocl->num_h2_leafs_per_row[i] += 1;

        /************************* h2_leafs_per_row ***************************/

        gcocl->h2_leafs_per_row[i] =
          realloc(gcocl->h2_leafs_per_row[i],
                  gcocl->num_h2_leafs_per_row[i] * sizeof(ph2matrix*));

        gcocl->h2_leafs_per_row[i][gcocl->num_h2_leafs_per_row[i] - 1] = H2;

        /************************* col_names_per_row **************************/

        gcocl->col_names_per_row[i] =
          realloc(gcocl->col_names_per_row[i],
                  gcocl->num_h2_leafs_per_row[i] * sizeof(ph2matrix*));

        gcocl->col_names_per_row[i][gcocl->num_h2_leafs_per_row[i] - 1] = cname;

        for(uint j = 0; j < gcocl->num_row_leafs; ++j)
          for(uint k = 0; k < gcocl->num_h2_leafs_per_row[j]; ++k)
            if(i != j || ((gcocl->num_h2_leafs_per_row[i] - 1) != k))
              if(gcocl->col_names_per_row[j][k] == cname)
                {
                  column_already_used = true;

                  /********************** buf_cidx_off ************************/

                  gcocl->cidx_off[i] =
                    realloc(gcocl->cidx_off[i],
                            gcocl->num_h2_leafs_per_row[i] * sizeof(ph2matrix*));

                  gcocl->cidx_off[i][gcocl->num_h2_leafs_per_row[i] - 1] =
                    gcocl->cidx_off[j][k];

                  break;
                }

        if(!column_already_used)
        {
          /************************** buf_cidx_off ****************************/

          gcocl->cidx_off[i] =
            realloc(gcocl->cidx_off[i],
                    gcocl->num_h2_leafs_per_row[i] * sizeof(ph2matrix*));

          gcocl->cidx_off[i][gcocl->num_h2_leafs_per_row[i] - 1] = gcocl->coff;

          /**************************** host_cidx *****************************/

          gcocl->host_cidx = realloc(gcocl->host_cidx,
                                     (gcocl->coff + H2->cb->k) * sizeof(uint));

          assert(gcocl->host_cidx != NULL);

          memcpy(gcocl->host_cidx + gcocl->coff,
                 ((pbem3d) gc->bem)->par->gcbn[cname]->xihat,
                 H2->cb->k * sizeof(uint));

          gcocl->coff += H2->cb->k;
        }

        break;
      }

    if(!row_already_stored)
    {
      bool column_already_used = false;

      // Update number of leafs
      gcocl->num_row_leafs += 1;

      /**************************** names_row_leafs ***************************/

      gcocl->names_row_leafs = realloc(gcocl->names_row_leafs,
                                       gcocl->num_row_leafs * sizeof(uint));

      assert(gcocl->names_row_leafs != NULL);

      gcocl->names_row_leafs[gcocl->num_row_leafs - 1] = rname;

      /************************* num_h2_leafs_per_row *************************/

      gcocl->num_h2_leafs_per_row = realloc(gcocl->num_h2_leafs_per_row,
                                            gcocl->num_row_leafs * sizeof(uint));

      assert(gcocl->num_h2_leafs_per_row != NULL);

      gcocl->num_h2_leafs_per_row[gcocl->num_row_leafs - 1] = 1;

      /*************************** h2_leafs_per_row ***************************/

      gcocl->h2_leafs_per_row = realloc(gcocl->h2_leafs_per_row,
                                        gcocl->num_row_leafs * sizeof(ph2matrix*));

      assert(gcocl->h2_leafs_per_row != NULL);

      gcocl->h2_leafs_per_row[gcocl->num_row_leafs - 1] =
        calloc(1, sizeof(ph2matrix*));

      gcocl->h2_leafs_per_row[gcocl->num_row_leafs - 1][0] = H2;


      /*************************** col_names_per_row **************************/

      gcocl->col_names_per_row = realloc(gcocl->col_names_per_row,
                                         gcocl->num_row_leafs * sizeof(uint*));

      assert(gcocl->col_names_per_row != NULL);

      gcocl->col_names_per_row[gcocl->num_row_leafs - 1] =
        calloc(1, sizeof(uint*));

      gcocl->col_names_per_row[gcocl->num_row_leafs - 1][0] = cname;

      /******************************* ridx_off *******************************/

      gcocl->ridx_off = realloc(gcocl->ridx_off,
                                gcocl->num_row_leafs * sizeof(uint));

      assert(gcocl->ridx_off != NULL);

      gcocl->ridx_off[gcocl->num_row_leafs - 1] = gcocl->roff;

      /****************************** host_ridx *******************************/

      gcocl->host_ridx = realloc(gcocl->host_ridx,
                                 (gcocl->roff + H2->rb->k) * sizeof(uint));

      assert(gcocl->host_ridx != NULL);

      memcpy(gcocl->host_ridx + gcocl->roff,
             ((pbem3d) gc->bem)->par->grbn[rname]->xihat,
             H2->rb->k * sizeof(uint));

      gcocl->roff += H2->rb->k;

      for(int i = 0; i < (int) gcocl->num_row_leafs - 2; ++i)
        for(uint j = 0; j < gcocl->num_h2_leafs_per_row[i]; ++j)
          if(gcocl->col_names_per_row[i][j] == cname)
          {
            column_already_used = true;

            /**************************** cidx_off ****************************/

            gcocl->cidx_off = realloc(gcocl->cidx_off,
                                      gcocl->num_row_leafs * sizeof(uint*));

            assert(gcocl->cidx_off != NULL);

            gcocl->cidx_off[gcocl->num_row_leafs - 1] =
              calloc(1, sizeof(uint*));

            gcocl->cidx_off[gcocl->num_row_leafs - 1][0] =
              gcocl->cidx_off[i][j];

            break;
          }

      if(!column_already_used)
      {
        /**************************** cidx_off ****************************/

        gcocl->cidx_off = realloc(gcocl->cidx_off,
                                  gcocl->num_row_leafs * sizeof(uint*));

        assert(gcocl->cidx_off != NULL);

        gcocl->cidx_off[gcocl->num_row_leafs - 1] =
          calloc(1, sizeof(uint*));

        gcocl->cidx_off[gcocl->num_row_leafs - 1][0] = gcocl->coff;

        /**************************** host_cidx ****************************/

        gcocl->host_cidx = realloc(gcocl->host_cidx,
                                   (gcocl->coff + H2->cb->k) * sizeof(uint));

        assert(gcocl->host_cidx != NULL);

        memcpy(gcocl->host_cidx + gcocl->coff,
               ((pbem3d) gc->bem)->par->gcbn[cname]->xihat,
               H2->cb->k * sizeof(uint));

        gcocl->coff += H2->cb->k;
      }
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
  if(H2->u)
  {
//    if(H2->f && !gcocl->is_farfield)
//    {
//      xtoff += H2->cb->k;
//      ytoff += H2->rb->k;
//    }

    int  i = -1;
    uint j;

    for(uint k = 0; k < gcocl->num_row_leafs; ++k)
      if(are_equal_clusters(H2->rb->t, gcocl->h2_leafs_per_row[k][0]->rb->t))
      {
        i = k;
        break;
      }

    assert(i >= 0);

    gcocl->ytoffs[i] = ytoff;

    for(j = 0; j < gcocl->num_h2_leafs_per_row[i]; ++j)
      if(are_equal_clusters(H2->cb->t, gcocl->h2_leafs_per_row[i][j]->cb->t))
      {
        gcocl->xtoffs[i][j] = xtoff;
        break;
      }

    assert(j < gcocl->num_h2_leafs_per_row[i]);
  }
  else if(H2->son)
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
                                     xtoff + xtoff1 - (cb->sons > 0 ? 0 : cb->k),
                                     ytoff + ytoff1 - (rb->sons > 0 ? 0 : rb->k),
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
get_ocl_informations_gcocl(ph2matrix H2, pgreencross gc)
{
  pgcopencl gcocl = gc->gcocl;

  gcocl->num_row_leafs            = 0;
  gcocl->max_num_h2_leafs_per_row = 0;
  gcocl->size_ridx                = 0;
  gcocl->size_cidx                = 0;
  gcocl->roff                     = 0;
  gcocl->coff                     = 0;

  gcocl->names_row_leafs          = (uint *) calloc(0, sizeof(uint));
  gcocl->num_h2_leafs_per_row     = (uint *) calloc(0, sizeof(uint));
  gcocl->h2_leafs_per_row         = (ph2matrix **) calloc(0, sizeof(ph2matrix *));
  gcocl->col_names_per_row        = (uint **) calloc(0, sizeof(uint *));
  gcocl->ridx_off                 = (uint *) calloc(0, sizeof(uint));
  gcocl->host_ridx                = (uint *) calloc(0, sizeof(uint));
  gcocl->cidx_off                 = (uint **) calloc(0, sizeof(uint *));
  gcocl->host_cidx                = (uint *) calloc(0, sizeof(uint));

  iterate_h2matrix(H2,
                   0,
                   0,
                   0,
                   0,
                   get_leaf_block_informations,
                   NULL,
                   (void *) gc);

  gcocl->size_ridx  = gcocl->roff;
  gcocl->size_cidx  = gcocl->coff;
  gcocl->workload_per_row = (uint *) calloc(gcocl->num_row_leafs, sizeof(uint));
  gcocl->idx_off    = (uint *) calloc(gcocl->num_row_leafs, sizeof(uint));
  gcocl->ridx_sizes = (uint *) calloc(gcocl->num_row_leafs, sizeof(uint));

  for(uint i = 0; i < gcocl->num_row_leafs; ++i)
  {
    /* Get maximum number of H2-matrix leafs which share the same writing
     * cluster. */
    if(gcocl->num_h2_leafs_per_row[i] > gcocl->max_num_h2_leafs_per_row)
      gcocl->max_num_h2_leafs_per_row = gcocl->num_h2_leafs_per_row[i];

    /* Get workload of the current writing cluster. */
    for(uint j = 0; j < gcocl->num_h2_leafs_per_row[i]; ++j)
      gcocl->workload_per_row[i] +=
        gcocl->h2_leafs_per_row[i][j]->u
          ? gcocl->h2_leafs_per_row[i][j]->u->S.rows *
            gcocl->h2_leafs_per_row[i][j]->u->S.cols
          : gcocl->h2_leafs_per_row[i][j]->f->rows *
            gcocl->h2_leafs_per_row[i][j]->f->cols;

    /* Get the index offset for the current writing cluster needed by the
     * devices. */
    gcocl->idx_off[i] =
      (i == 0 ? 0 : gcocl->idx_off[i - 1] + gcocl->num_h2_leafs_per_row[i - 1]);

    /* Sum up the number of indices needed by the devices. */
    gcocl->ridx_sizes[i] = gcocl->h2_leafs_per_row[i][0]->rb->k;
  }

  gcocl->cidx_sizes =
    (uint *) calloc(gcocl->idx_off[gcocl->num_row_leafs - 1] +
                    gcocl->num_h2_leafs_per_row[gcocl->num_row_leafs - 1],
                    sizeof(uint));

  gcocl->host_cidx_off =
    (uint *) calloc(gcocl->idx_off[gcocl->num_row_leafs - 1] +
                    gcocl->num_h2_leafs_per_row[gcocl->num_row_leafs - 1],
                    sizeof(uint));

  /* Get number of indices for all H2-matrix leafs, ready to be written to
   * OpenCL devices. */
  for (uint i = 0; i < gcocl->num_row_leafs; ++i)
  {
    for(uint j = 0; j < gcocl->num_h2_leafs_per_row[i]; ++j)
      gcocl->cidx_sizes[gcocl->idx_off[i] + j] =
        gcocl->h2_leafs_per_row[i][j]->cb->k;

    memcpy(gcocl->host_cidx_off + gcocl->idx_off[i],
           gcocl->cidx_off[i],
           gcocl->num_h2_leafs_per_row[i] * sizeof(uint));
  }

  /* Write H2-matrix informations to the devices. */
  for(uint i = 0; i < ocl_system.num_devices; ++i)
  {
    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           gcocl->num_row_leafs,
                           sizeof(uint),
                           gcocl->num_h2_leafs_per_row,
                           NULL,
                           &gcocl->buf_num_h2_leafs_per_row[i]);

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           gcocl->num_row_leafs,
                           sizeof(uint),
                           gcocl->idx_off,
                           NULL,
                           &gcocl->buf_idx_off[i]);

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           gcocl->num_row_leafs,
                           sizeof(uint),
                           gcocl->ridx_sizes,
                           NULL,
                           &gcocl->buf_ridx_sizes[i]);

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           gcocl->num_row_leafs,
                           sizeof(uint),
                           gcocl->ridx_off,
                           NULL,
                           &gcocl->buf_ridx_off[i]);

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           gcocl->size_ridx,
                           sizeof(uint),
                           gcocl->host_ridx,
                           NULL,
                           &gcocl->buf_ridx[i]);

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           gcocl->idx_off[gcocl->num_row_leafs - 1] +
                           gcocl->num_h2_leafs_per_row[gcocl->num_row_leafs - 1],
                           sizeof(uint),
                           gcocl->cidx_sizes,
                           NULL,
                           &gcocl->buf_cidx_sizes[i]);

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           gcocl->idx_off[gcocl->num_row_leafs - 1] +
                           gcocl->num_h2_leafs_per_row[gcocl->num_row_leafs - 1],
                           sizeof(uint),
                           gcocl->host_cidx_off,
                           NULL,
                           &gcocl->buf_cidx_off[i]);

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           gcocl->size_cidx,
                           sizeof(uint),
                           gcocl->host_cidx,
                           NULL,
                           &gcocl->buf_cidx[i]);
  }

  /* Get offsets of individual writing clusters relative to the whole
   * H^2-matrix. */

  gcocl->xtoffs             = (uint **) calloc(gcocl->num_row_leafs,
                                               sizeof(uint*));

  gcocl->ytoffs             = (uint *) calloc(gcocl->num_row_leafs,
                                              sizeof(uint));

  for(uint i = 0; i < gcocl->num_row_leafs; ++i)
    gcocl->xtoffs[i]       =
      (uint *) calloc(gcocl->num_h2_leafs_per_row[i], sizeof(uint));

  iterate_recursively_h2matrix(H2, 0, 0, gcocl);

  gcocl->host_xtoffs   =
    (uint *) calloc(gcocl->idx_off[gcocl->num_row_leafs - 1] +
                    gcocl->num_h2_leafs_per_row[gcocl->num_row_leafs - 1],
                    sizeof(uint));

  /* Make the offsets ready to be written to OpenCL devices. */
  for(uint i = 0; i < gcocl->num_row_leafs; ++i)
  {
    memcpy(gcocl->host_xtoffs + gcocl->idx_off[i],
           gcocl->xtoffs[i],
           gcocl->num_h2_leafs_per_row[i] * sizeof(uint));
  }

  /* Write offsets to the devices. */
  for(uint i = 0; i < ocl_system.num_devices; ++i)
  {
    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           gcocl->idx_off[gcocl->num_row_leafs - 1] +
                           gcocl->num_h2_leafs_per_row[gcocl->num_row_leafs - 1],
                           sizeof(uint),
                           gcocl->host_xtoffs,
                           NULL,
                           &gcocl->buf_xtoffs[i]);

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           gcocl->num_row_leafs,
                           sizeof(uint),
                           gcocl->ytoffs,
                           NULL,
                           &gcocl->buf_ytoffs[i]);
  }
}

ph2matrix
build_green_cross_h2matrix_greencross(pgreencross gc, void *eta)
{
  pblock    broot;
  ph2matrix H2;

  broot = build_strict_block(gc->rc, gc->cc, eta, admissible_max_cluster);

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

  /* Get farfield leaf informations. */

  gc->gcocl = new_gcopencl(FARFIELD_GCA);

  get_ocl_informations_gcocl(H2, gc);

  /* Get the informations an OpenCL device needs to perform the MVM-part of
   * nearfield matrices. */
  if(gc->dim == 2)
  {
    pcurve2d gr = (pcurve2d) gc->geom;

    gc->ocl_info_nf = new_nearfield_gcopencl(H2, gc->dim, gr->e);

    gc->iinfos_min_edge = build_from_idxinfos_integralinfos(gc->ocl_info_nf,
                                                            H2,
                                                            gc->dim,
                                                            (void*) gr->e,
                                                            1);
  }
  else
  {
    psurface3d gr = (psurface3d) gc->geom;

    gc->ocl_info_nf = new_nearfield_gcopencl(H2, gc->dim, gr->t);

    gc->iinfos_min_edge = build_from_idxinfos_integralinfos(gc->ocl_info_nf,
                                                            H2,
                                                            gc->dim,
                                                            (void*) gr->t,
                                                            1);
  }

  gc->oclwrk      = new_equidistant_distributed_oclwork(gc->gcocl,
                                                        H2,
                                                        1);

  /* TODO: This must be synchronized with the above work packages. */
  gc->ocl_wrk_nf  = new_equidistant_distributed_oclwork(gc->ocl_info_nf,
                                                        H2,
                                                        1);

  /* Set up buffer for fastaddeval. */
  gc->feval = new_fastaddevalgca(H2->rb, H2->cb, 1);

  gc->kernel_events = (cl_event*) allocmem(num_kernels * sizeof(cl_event));

  for(uint i = 0; i < num_kernels; ++i)
    gc->kernel_events[i] = CL_SUCCESS;

  if(delete_kernels(num_kernels, &gc->ocl_kernels) != NULL)
  {
    fprintf(stderr, "error: can't create new kernels\n");
    exit(1);
  }

  const uint max_num_h2_leafs =
    UINT_MAX((uint) gc->gcocl->max_num_h2_leafs_per_row,
             (uint) gc->ocl_info_nf->max_num_h2_leafs_per_row);

  char add_flags[100];

  int n = sprintf(add_flags,
                  "-cl-nv-verbose -DQUADRATUR_ORDER=%u -DSIZE=%u "
                  "-DMAX_WRK_GRP_SIZE=%u",
                  ((pbem3d) gc->bem)->sq->n_dist,
                  16,
                  max_num_h2_leafs);

  if(n <= 0)
  {
    fprintf(stderr, "error: can't define compiler flags for OpenCL!");
    exit(1);
  }

  setup_kernels_fix(7,
                    src_code_strs,
                    add_flags,
                    num_kernels,
                    kernel_names,
                    &gc->ocl_kernels);


  for(uint i = 0; i < gc->oclwrk->num_wrk_pkgs; ++i)
  {
      for(uint k = 0; k < ocl_system.queues_per_device; ++k)
      {
        for(uint j = 0; j < num_kernels; ++j)
        {
          cl_kernel kernel =
            gc->ocl_kernels[k + i * ocl_system.queues_per_device +
                            j * ocl_system.num_devices *
                            ocl_system.queues_per_device];

          CL_CHECK(clSetKernelArg(kernel, 0, sizeof(uint), &gc->dim));
          CL_CHECK(clSetKernelArg(kernel, 1, sizeof(uint), &gc->n));
          CL_CHECK(clSetKernelArg(kernel, 2, sizeof(cl_mem), &gc->buf_x[i]));
          CL_CHECK(clSetKernelArg(kernel, 3, sizeof(cl_mem), &gc->buf_p[i]));
          CL_CHECK(clSetKernelArg(kernel, 4, sizeof(cl_mem), &gc->buf_g[i]));
        }

        cl_kernel kernel = gc->ocl_kernels[k + i * ocl_system.queues_per_device];

        CL_CHECK(clSetKernelArg(kernel, 5, sizeof(uint),   &gc->sq_gca->nq));
        CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), &gc->sq_gca->buf_xqs[i]));
        CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), &gc->sq_gca->buf_yqs[i]));
        CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), &gc->sq_gca->buf_wqs[i]));
        CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), &gc->sq_gca->buf_bases[i]));

        CL_CHECK(clSetKernelArg(kernel, 10, sizeof(uint),   &gc->oclwrk->num_rows_per_pkg[i]));
        CL_CHECK(clSetKernelArg(kernel, 11, sizeof(cl_mem), &gc->oclwrk->buf_rows_this_device[i]));
        CL_CHECK(clSetKernelArg(kernel, 12, sizeof(cl_mem), &gc->gcocl->buf_num_h2_leafs_per_row[i]));
        CL_CHECK(clSetKernelArg(kernel, 13, sizeof(cl_mem), &gc->gcocl->buf_idx_off[i]));
        CL_CHECK(clSetKernelArg(kernel, 14, sizeof(cl_mem), &gc->gcocl->buf_ridx_sizes[i]));
        CL_CHECK(clSetKernelArg(kernel, 15, sizeof(cl_mem), &gc->gcocl->buf_cidx_sizes[i]));
        CL_CHECK(clSetKernelArg(kernel, 16, sizeof(cl_mem), &gc->gcocl->buf_ridx_off[i]));
        CL_CHECK(clSetKernelArg(kernel, 17, sizeof(cl_mem), &gc->gcocl->buf_cidx_off[i]));
        CL_CHECK(clSetKernelArg(kernel, 18, sizeof(cl_mem), &gc->gcocl->buf_ridx[i]));
        CL_CHECK(clSetKernelArg(kernel, 19, sizeof(cl_mem), &gc->gcocl->buf_cidx[i]));
        CL_CHECK(clSetKernelArg(kernel, 20, sizeof(cl_mem), &gc->gcocl->buf_xtoffs[i]));
        CL_CHECK(clSetKernelArg(kernel, 21, sizeof(cl_mem), &gc->gcocl->buf_ytoffs[i]));

        if(gc->dim == 2)
        {
          CL_CHECK(clSetKernelArg(kernel, 22, sizeof(real), &((pbem2d) gc->bem)->alpha));
        }
        else
        {
          CL_CHECK(clSetKernelArg(kernel, 22, sizeof(real), &((pbem3d) gc->bem)->alpha));
        }

        CL_CHECK(clSetKernelArg(kernel, 24, sizeof(cl_mem), &gc->feval->buf_xt[i]));
//#ifndef USE_OPENMP
        CL_CHECK(clSetKernelArg(kernel, 25, sizeof(cl_mem), &gc->feval->buf_yt[i]));
//#else
//        CL_CHECK(clSetKernelArg(kernel, 25, sizeof(cl_mem), &gc->feval->buf_yt_ff[i]));
//#endif

        kernel = gc->ocl_kernels[k + i * ocl_system.queues_per_device +
                                 ocl_system.num_devices *
                                 ocl_system.queues_per_device];

        CL_CHECK(clSetKernelArg(kernel, 5, sizeof(uint),   &gc->sq_gca->nq));
        CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), &gc->sq_gca->buf_xqs[i]));
        CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), &gc->sq_gca->buf_yqs[i]));
        CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), &gc->sq_gca->buf_wqs[i]));
        CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), &gc->sq_gca->buf_bases[i]));

        CL_CHECK(clSetKernelArg(kernel, 10, sizeof(uint),   &gc->ocl_wrk_nf->num_rows_per_pkg[i]));
        CL_CHECK(clSetKernelArg(kernel, 11, sizeof(cl_mem), &gc->ocl_wrk_nf->buf_rows_this_device[i]));
        CL_CHECK(clSetKernelArg(kernel, 12, sizeof(cl_mem), &gc->ocl_info_nf->buf_num_h2_leafs_per_row[i]));
        CL_CHECK(clSetKernelArg(kernel, 13, sizeof(cl_mem), &gc->ocl_info_nf->buf_idx_off[i]));
        CL_CHECK(clSetKernelArg(kernel, 14, sizeof(cl_mem), &gc->ocl_info_nf->buf_ridx_sizes[i]));
        CL_CHECK(clSetKernelArg(kernel, 15, sizeof(cl_mem), &gc->ocl_info_nf->buf_cidx_sizes[i]));
        CL_CHECK(clSetKernelArg(kernel, 16, sizeof(cl_mem), &gc->ocl_info_nf->buf_ridx_off[i]));
        CL_CHECK(clSetKernelArg(kernel, 17, sizeof(cl_mem), &gc->ocl_info_nf->buf_cidx_off[i]));
        CL_CHECK(clSetKernelArg(kernel, 18, sizeof(cl_mem), &gc->ocl_info_nf->buf_ridx[i]));
        CL_CHECK(clSetKernelArg(kernel, 19, sizeof(cl_mem), &gc->ocl_info_nf->buf_cidx[i]));
        CL_CHECK(clSetKernelArg(kernel, 20, sizeof(cl_mem), &gc->ocl_info_nf->buf_xtoffs[i]));
        CL_CHECK(clSetKernelArg(kernel, 21, sizeof(cl_mem), &gc->ocl_info_nf->buf_ytoffs[i]));

        if(gc->dim == 2)
        {
          CL_CHECK(clSetKernelArg(kernel, 22, sizeof(real), &((pbem2d) gc->bem)->alpha));
        }
        else
        {
          CL_CHECK(clSetKernelArg(kernel, 22, sizeof(real), &((pbem3d) gc->bem)->alpha));
        }

        CL_CHECK(clSetKernelArg(kernel, 24, sizeof(cl_mem), &gc->feval->buf_xt[i]));
//#ifndef USE_OPENMP
        CL_CHECK(clSetKernelArg(kernel, 25, sizeof(cl_mem), &gc->feval->buf_yt[i]));
//#else
//        CL_CHECK(clSetKernelArg(kernel, 25, sizeof(cl_mem), &gc->feval->buf_yt_nf_common[i]));
//#endif

        kernel = gc->ocl_kernels[k + i * ocl_system.queues_per_device +
                                 2 * ocl_system.num_devices *
                                 ocl_system.queues_per_device];

        CL_CHECK(clSetKernelArg(kernel, 5, sizeof(uint),   &gc->sq_partial_min_vert->nq));
        CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), &gc->sq_partial_min_vert->buf_xqs[i]));
        CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), &gc->sq_partial_min_vert->buf_yqs[i]));
        CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), &gc->sq_partial_min_vert->buf_wqs[i]));
        CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), &gc->sq_partial_min_vert->buf_bases[i]));

        CL_CHECK(clSetKernelArg(kernel, 10, sizeof(uint),   &gc->ocl_wrk_nf->num_rows_per_pkg[i]));
        CL_CHECK(clSetKernelArg(kernel, 11, sizeof(cl_mem), &gc->ocl_wrk_nf->buf_rows_this_device[i]));
        CL_CHECK(clSetKernelArg(kernel, 12, sizeof(cl_mem), &gc->ocl_info_nf->buf_num_h2_leafs_per_row[i]));
        CL_CHECK(clSetKernelArg(kernel, 13, sizeof(cl_mem), &gc->ocl_info_nf->buf_idx_off[i]));
        CL_CHECK(clSetKernelArg(kernel, 14, sizeof(cl_mem), &gc->ocl_info_nf->buf_ridx_sizes[i]));
        CL_CHECK(clSetKernelArg(kernel, 15, sizeof(cl_mem), &gc->ocl_info_nf->buf_cidx_sizes[i]));
        CL_CHECK(clSetKernelArg(kernel, 16, sizeof(cl_mem), &gc->ocl_info_nf->buf_ridx_off[i]));
        CL_CHECK(clSetKernelArg(kernel, 17, sizeof(cl_mem), &gc->ocl_info_nf->buf_cidx_off[i]));
        CL_CHECK(clSetKernelArg(kernel, 18, sizeof(cl_mem), &gc->ocl_info_nf->buf_ridx[i]));
        CL_CHECK(clSetKernelArg(kernel, 19, sizeof(cl_mem), &gc->ocl_info_nf->buf_cidx[i]));
        CL_CHECK(clSetKernelArg(kernel, 20, sizeof(cl_mem), &gc->ocl_info_nf->buf_xtoffs[i]));
        CL_CHECK(clSetKernelArg(kernel, 21, sizeof(cl_mem), &gc->ocl_info_nf->buf_ytoffs[i]));

        CL_CHECK(clSetKernelArg(kernel, 22, sizeof(cl_mem), &gc->ocl_info_nf->buf_num_any_id[i]));
        CL_CHECK(clSetKernelArg(kernel, 23, sizeof(cl_mem), &gc->ocl_info_nf->buf_idx_off_any_id[i]));
        CL_CHECK(clSetKernelArg(kernel, 24, sizeof(cl_mem), &gc->ocl_info_nf->buf_rows_any_id[i]));
        CL_CHECK(clSetKernelArg(kernel, 25, sizeof(cl_mem), &gc->ocl_info_nf->buf_cols_any_id[i]));
        CL_CHECK(clSetKernelArg(kernel, 26, sizeof(cl_mem), &gc->ocl_info_nf->buf_cidx_any_id[i]));

        if(gc->dim == 2)
        {
          CL_CHECK(clSetKernelArg(kernel, 27, sizeof(real), &((pbem2d) gc->bem)->alpha));
        }
        else
        {
          CL_CHECK(clSetKernelArg(kernel, 27, sizeof(real), &((pbem3d) gc->bem)->alpha));
        }

        CL_CHECK(clSetKernelArg(kernel, 29, sizeof(cl_mem), &gc->feval->buf_xt[i]));
//#ifndef USE_OPENMP
        CL_CHECK(clSetKernelArg(kernel, 30, sizeof(cl_mem), &gc->feval->buf_yt[i]));
//#else
//        CL_CHECK(clSetKernelArg(kernel, 30, sizeof(cl_mem), &gc->feval->buf_yt_nf_min_vert[i]));
//#endif

        kernel = gc->ocl_kernels[k + i * ocl_system.queues_per_device +
                                 3 * ocl_system.num_devices *
                                 ocl_system.queues_per_device];

        CL_CHECK(clSetKernelArg(kernel, 5, sizeof(uint),   &gc->sq_partial_min_edge->nq));
        CL_CHECK(clSetKernelArg(kernel, 6, sizeof(cl_mem), &gc->sq_partial_min_edge->buf_xqs[i]));
        CL_CHECK(clSetKernelArg(kernel, 7, sizeof(cl_mem), &gc->sq_partial_min_edge->buf_yqs[i]));
        CL_CHECK(clSetKernelArg(kernel, 8, sizeof(cl_mem), &gc->sq_partial_min_edge->buf_wqs[i]));
        CL_CHECK(clSetKernelArg(kernel, 9, sizeof(cl_mem), &gc->sq_partial_min_edge->buf_bases[i]));

        CL_CHECK(clSetKernelArg(kernel, 10, sizeof(uint),   &gc->ocl_wrk_nf->num_rows_per_pkg[i]));
        CL_CHECK(clSetKernelArg(kernel, 11, sizeof(cl_mem), &gc->ocl_wrk_nf->buf_rows_this_device[i]));
        CL_CHECK(clSetKernelArg(kernel, 12, sizeof(cl_mem), &gc->ocl_info_nf->buf_num_h2_leafs_per_row[i]));
        CL_CHECK(clSetKernelArg(kernel, 13, sizeof(cl_mem), &gc->ocl_info_nf->buf_idx_off[i]));
        CL_CHECK(clSetKernelArg(kernel, 14, sizeof(cl_mem), &gc->ocl_info_nf->buf_ridx_sizes[i]));
        CL_CHECK(clSetKernelArg(kernel, 15, sizeof(cl_mem), &gc->ocl_info_nf->buf_cidx_sizes[i]));
        CL_CHECK(clSetKernelArg(kernel, 16, sizeof(cl_mem), &gc->ocl_info_nf->buf_ridx_off[i]));
        CL_CHECK(clSetKernelArg(kernel, 17, sizeof(cl_mem), &gc->ocl_info_nf->buf_cidx_off[i]));
        CL_CHECK(clSetKernelArg(kernel, 18, sizeof(cl_mem), &gc->ocl_info_nf->buf_ridx[i]));
        CL_CHECK(clSetKernelArg(kernel, 19, sizeof(cl_mem), &gc->ocl_info_nf->buf_cidx[i]));
        CL_CHECK(clSetKernelArg(kernel, 20, sizeof(cl_mem), &gc->ocl_info_nf->buf_xtoffs[i]));
        CL_CHECK(clSetKernelArg(kernel, 21, sizeof(cl_mem), &gc->ocl_info_nf->buf_ytoffs[i]));

        CL_CHECK(clSetKernelArg(kernel, 22, sizeof(cl_mem), &gc->iinfos_min_edge->buf_num_integrals[i]));
        CL_CHECK(clSetKernelArg(kernel, 23, sizeof(cl_mem), &gc->iinfos_min_edge->buf_idx_off[i]));
        CL_CHECK(clSetKernelArg(kernel, 24, sizeof(cl_mem), &gc->iinfos_min_edge->buf_rows[i]));
        CL_CHECK(clSetKernelArg(kernel, 25, sizeof(cl_mem), &gc->iinfos_min_edge->buf_cols[i]));
        CL_CHECK(clSetKernelArg(kernel, 26, sizeof(cl_mem), &gc->iinfos_min_edge->buf_cidx[i]));

        if(gc->dim == 2)
        {
          CL_CHECK(clSetKernelArg(kernel, 27, sizeof(real), &((pbem2d) gc->bem)->alpha));
        }
        else
        {
          CL_CHECK(clSetKernelArg(kernel, 27, sizeof(real), &((pbem3d) gc->bem)->alpha));
        }

        CL_CHECK(clSetKernelArg(kernel, 29, sizeof(cl_mem), &gc->feval->buf_xt[i]));
//#ifndef USE_OPENMP
        CL_CHECK(clSetKernelArg(kernel, 30, sizeof(cl_mem), &gc->feval->buf_yt[i]));
//#else
//        CL_CHECK(clSetKernelArg(kernel, 30, sizeof(cl_mem), &gc->feval->buf_yt_nf_min_edge[i]));
//#endif
      }
  }

  del_block(broot);

  return H2;
}

void
fastaddeval_farfield_h2matrix_avector(field      alpha,
                                      pch2matrix h2,
                                      pavector   xt,
                                      pavector   yt)
{
  avector        loc1, loc2;
  pavector       xt1, yt1;
  pcclusterbasis rb = h2->rb;
  pcclusterbasis cb = h2->cb;
  uint      rsons = h2->rsons;
  uint      csons = h2->csons;
  uint      xtoff, ytoff;
  uint      i, j;

  if(h2->u)
    addeval_amatrix_avector(alpha, &h2->u->S, xt, yt);
  else if(h2->son)
  {
    xtoff = cb->k;
    for (j = 0; j < csons; j++) {
      assert(csons == 1 || cb->sons > 0);
      xt1 =
        (cb->sons >
         0 ? init_sub_avector(&loc1, xt, cb->son[j]->ktree,
                              xtoff) : init_sub_avector(&loc1, xt, cb->ktree,
                                                        0));

      ytoff = rb->k;
      for (i = 0; i < rsons; i++) {
        assert(rsons == 1 || rb->sons > 0);
        yt1 = (rb->sons > 0 ?
               init_sub_avector(&loc2, yt, rb->son[i]->ktree, ytoff) :
               init_sub_avector(&loc2, yt, rb->ktree, 0));

        fastaddeval_farfield_h2matrix_avector(alpha, h2->son[i + j * rsons], xt1, yt1);

        uninit_avector(yt1);

        ytoff += (rb->sons > 0 ? rb->son[i]->ktree : rb->t->size);
      }
      assert(ytoff == rb->ktree);

      uninit_avector(xt1);

      xtoff += (cb->sons > 0 ? cb->son[j]->ktree : cb->t->size);
    }
    assert(xtoff == cb->ktree);
  }
}

void
fastaddeval_nearfield_h2matrix_avector(field      alpha,
                                       pch2matrix h2,
                                       pavector   xt,
                                       pavector   yt)
{
  avector   loc1, loc2;
  pavector  xp, yp, xt1, yt1;
  pcclusterbasis rb = h2->rb;
  pcclusterbasis cb = h2->cb;
  uint      rsons = h2->rsons;
  uint      csons = h2->csons;
  uint      xtoff, ytoff;
  uint      i, j;

  if (h2->f) {
    xp = init_sub_avector(&loc1, xt, cb->t->size, cb->k);
    yp = init_sub_avector(&loc2, yt, rb->t->size, rb->k);

    addeval_amatrix_avector(alpha, h2->f, xp, yp);

    uninit_avector(yp);
    uninit_avector(xp);
  }
  else if (h2->son) {
    xtoff = cb->k;
    for (j = 0; j < csons; j++) {
      assert(csons == 1 || cb->sons > 0);
      xt1 =
        (cb->sons >
         0 ? init_sub_avector(&loc1, xt, cb->son[j]->ktree,
                              xtoff) : init_sub_avector(&loc1, xt, cb->ktree,
                                                        0));

      ytoff = rb->k;
      for (i = 0; i < rsons; i++) {
        assert(rsons == 1 || rb->sons > 0);
        yt1 = (rb->sons > 0 ?
               init_sub_avector(&loc2, yt, rb->son[i]->ktree, ytoff) :
               init_sub_avector(&loc2, yt, rb->ktree, 0));

        fastaddeval_nearfield_h2matrix_avector(alpha, h2->son[i + j * rsons], xt1, yt1);

        uninit_avector(yt1);

        ytoff += (rb->sons > 0 ? rb->son[i]->ktree : rb->t->size);
      }
      assert(ytoff == rb->ktree);

      uninit_avector(xt1);

      xtoff += (cb->sons > 0 ? cb->son[j]->ktree : cb->t->size);
    }
    assert(xtoff == cb->ktree);
  }
}

void
fastaddeval_farfield_h2matrix_avector_gca(pgreencross gca
//#ifndef USE_OPENMP
)
//#else
//                                        , pavector    yt)
//#endif
{
  const size_t max_num_h2_leafs         = gca->gcocl->max_num_h2_leafs_per_row;

  const size_t max_num_writing_clusters = gca->gcocl->num_row_leafs;

  const size_t num_global_work_items    = max_num_h2_leafs *
                                          max_num_writing_clusters;

  pfastaddevalgca feval  = gca->feval;
  poclworkpgs     oclwrk = gca->oclwrk;

  cl_command_queue *queues = ocl_system.queues;

  for(uint i = 0; i < oclwrk->num_wrk_pkgs; ++i)
  {
//#ifndef USE_OPENMP
    cl_event events[2] = { feval->events_xt[i], feval->events_yt[i] };
//#else
//    cl_event events[2] = { feval->events_xt[i], feval->events_yt[4 * i] };
//#endif

    clWaitForEvents(2, events);

    /* Start kernel */
    CL_CHECK(clEnqueueNDRangeKernel
               (queues[i * ocl_system.queues_per_device],
                gca->ocl_kernels[i * ocl_system.queues_per_device],
                1,
                NULL,
                &num_global_work_items,
                &max_num_h2_leafs,
                0,
                NULL,
                &gca->kernel_events[0]));
  }

//#ifdef USE_OPENMP
//
//  cl_event events[oclwrk->num_wrk_pkgs];
//
//  for(uint i = 0; i < oclwrk->num_wrk_pkgs; ++i)
//  {
//    const uint offset   = oclwrk->first_idx_of_pkgs[i];
//    const uint num_comp = oclwrk->last_idx_of_pkg[i] - offset;
//
//    avector  tmp;
//
//    pavector sub_yt = init_sub_avector(&tmp, yt, num_comp, offset);
//
//    clWaitForEvents(1, &gca->kernel_events[0]);
//
//    CL_CHECK(clEnqueueReadBuffer (queues[i * ocl_system.queues_per_device],
//                                  feval->buf_yt_ff[i],
//                                  CL_FALSE,
//                                  offset * sizeof(real),
//                                  num_comp * sizeof(real),
//                                  sub_yt->v,
//                                  0,
//                                  NULL,
//                                  &feval->events_yt[4 * i]));
//
//    events[i] = feval->events_yt[4 * i];
//
//    uninit_avector(sub_yt);
//  }
//
//  clWaitForEvents(oclwrk->num_wrk_pkgs, events);
//#endif
}

void
fastaddeval_nearfield_common_h2matrix_avector_gca(pgreencross gc
//#ifndef USE_OPENMP
  )
//#else
//  , pavector    yt)
//#endif
{
  const size_t max_num_h2_leafs      = gc->ocl_info_nf->max_num_h2_leafs_per_row;
  const size_t num_writing_clusters  = gc->ocl_info_nf->num_row_leafs;
  const size_t num_global_work_items = max_num_h2_leafs * num_writing_clusters;

  pcoclworkpgs     oclwrk  = gc->oclwrk;

//#ifdef USE_OPENMP
//  pcfastaddevalgca feval   = gc->feval;
//#endif

  cl_command_queue *queues = ocl_system.queues;

  for(uint i = 0; i < oclwrk->num_wrk_pkgs; ++i)
  {
//#ifndef USE_OPENMP
//    clFinish(queues[i * ocl_system.queues_per_device]);
    clWaitForEvents(1, &gc->kernel_events[0]);
//#else
//    cl_event events[2] = { feval->events_xt[i], feval->events_yt[4 * i + 1] };
//
//    clWaitForEvents(2, events);
//#endif

    /* Start kernel */
    CL_CHECK(clEnqueueNDRangeKernel
               (queues[i * ocl_system.queues_per_device],
                gc->ocl_kernels[i * ocl_system.queues_per_device +
                                ocl_system.num_devices *
                                ocl_system.queues_per_device],
                1,
                NULL,
                &num_global_work_items,
                &max_num_h2_leafs,
                0,
                NULL,
                &gc->kernel_events[1]));
  }

//#ifdef USE_OPENMP
//
//  cl_event events[oclwrk->num_wrk_pkgs];
//
//  for(uint i = 0; i < oclwrk->num_wrk_pkgs; ++i)
//  {
//    const uint offset   = oclwrk->first_idx_of_pkgs[i];
//    const uint num_comp = oclwrk->last_idx_of_pkg[i] - offset;
//
//    avector  tmp;
//
//    pavector sub_yt = init_sub_avector(&tmp, yt, num_comp, offset);
//
//    clWaitForEvents(1, &gc->kernel_events[1]);
//
//    CL_CHECK(clEnqueueReadBuffer (queues[i * ocl_system.queues_per_device],
//                                  feval->buf_yt_nf_common[i],
//                                  CL_FALSE,
//                                  offset * sizeof(real),
//                                  num_comp * sizeof(real),
//                                  sub_yt->v,
//                                  0,
//                                  NULL,
//                                  &feval->events_yt[4 * i + 1]));
//
//    events[i] = feval->events_yt[4 * i + 1];
//
//    uninit_avector(sub_yt);
//  }
//
//  clWaitForEvents(oclwrk->num_wrk_pkgs, events);
//#endif
}

void
fastaddeval_nearfield_min_vert_h2matrix_avector_gca(pgreencross gc
//#ifndef USE_OPENMP
  )
//#else
//  , pavector    yt)
//#endif
{
  const size_t max_num_h2_leafs      = gc->ocl_info_nf->max_num_h2_leafs_per_row;
  const size_t num_writing_clusters  = gc->ocl_info_nf->num_row_leafs;
  const size_t num_global_work_items = max_num_h2_leafs * num_writing_clusters;

  pcoclworkpgs     oclwrk  = gc->oclwrk;

//#ifdef USE_OPENMP
//  pcfastaddevalgca feval   = gc->feval;
//#endif

  cl_command_queue *queues = ocl_system.queues;

  for(uint i = 0; i < oclwrk->num_wrk_pkgs; ++i)
  {
//#ifndef USE_OPENMP
//    clFinish(queues[i * ocl_system.queues_per_device]);
    clWaitForEvents(1, &gc->kernel_events[1]);
//#else
//    cl_event events[2] = { feval->events_xt[i], feval->events_yt[4 * i + 2] };
//
//    clWaitForEvents(2, events);
//#endif
    /* Start kernel */
    CL_CHECK(clEnqueueNDRangeKernel
               (queues[i * ocl_system.queues_per_device],
                gc->ocl_kernels[i * ocl_system.queues_per_device +
                                2 * ocl_system.num_devices *
                                ocl_system.queues_per_device],
                1,
                NULL,
                &num_global_work_items,
                &max_num_h2_leafs,
                0,
                NULL,
                &gc->kernel_events[2]));
  }

//#ifdef USE_OPENMP
//
//  cl_event events[oclwrk->num_wrk_pkgs];
//
//  for(uint i = 0; i < oclwrk->num_wrk_pkgs; ++i)
//  {
//    const uint offset   = oclwrk->first_idx_of_pkgs[i];
//    const uint num_comp = oclwrk->last_idx_of_pkg[i] - offset;
//
//    avector  tmp;
//
//    pavector sub_yt = init_sub_avector(&tmp, yt, num_comp, offset);
//
//    clWaitForEvents(1, &gc->kernel_events[2]);
//
//    CL_CHECK(clEnqueueReadBuffer (queues[i * ocl_system.queues_per_device],
//                                  feval->buf_yt_nf_min_vert[i],
//                                  CL_FALSE,
//                                  offset * sizeof(real),
//                                  num_comp * sizeof(real),
//                                  sub_yt->v,
//                                  0,
//                                  NULL,
//                                  &feval->events_yt[4 * i + 2]));
//
//    events[i] = feval->events_yt[4 * i + 2];
//
//    uninit_avector(sub_yt);
//  }
//
//  clWaitForEvents(oclwrk->num_wrk_pkgs, events);
//#endif
}

void
fastaddeval_nearfield_min_edge_h2matrix_avector_gca(pgreencross gc
//#ifndef USE_OPENMP
  )
//#else
//  , pavector    yt)
//#endif
{
  const size_t max_num_h2_leafs      = gc->ocl_info_nf->max_num_h2_leafs_per_row;
  const size_t num_groups            = gc->iinfos_min_edge->num_integral_grps;
  const size_t num_global_work_items = max_num_h2_leafs * num_groups;

  pcoclworkpgs     oclwrk            = gc->oclwrk;

//#ifdef USE_OPENMP
//  pcfastaddevalgca feval   = gc->feval;
//#endif

  cl_command_queue *queues           = ocl_system.queues;

  for(uint i = 0; i < oclwrk->num_wrk_pkgs; ++i)
  {
//#ifndef USE_OPENMP
//    clFinish(queues[i * ocl_system.queues_per_device]);
    clWaitForEvents(1, &gc->kernel_events[2]);
//#else
//    cl_event events[2] = { feval->events_xt[i], feval->events_yt[4 * i + 3] };
//
//    clWaitForEvents(2, events);
//#endif

    /* Start kernel */
    CL_CHECK(clEnqueueNDRangeKernel
               (queues[i * ocl_system.queues_per_device],
                gc->ocl_kernels[i * ocl_system.queues_per_device +
                                3 * ocl_system.num_devices *
                                ocl_system.queues_per_device],
                1,
                NULL,
                &num_global_work_items,
                &max_num_h2_leafs,
                0,
                NULL,
                &gc->kernel_events[3]));
  }

//#ifdef USE_OPENMP
//
//  cl_event events[oclwrk->num_wrk_pkgs];
//
//  for(uint i = 0; i < oclwrk->num_wrk_pkgs; ++i)
//  {
//    const uint offset   = oclwrk->first_idx_of_pkgs[i];
//    const uint num_comp = oclwrk->last_idx_of_pkg[i] - offset;
//
//    avector  tmp;
//
//    pavector sub_yt = init_sub_avector(&tmp, yt, num_comp, offset);
//
//    clWaitForEvents(1, &gc->kernel_events[3]);
//
//    CL_CHECK(clEnqueueReadBuffer (queues[i * ocl_system.queues_per_device],
//                                  feval->buf_yt_nf_min_edge[i],
//                                  CL_FALSE,
//                                  offset * sizeof(real),
//                                  num_comp * sizeof(real),
//                                  sub_yt->v,
//                                  0,
//                                  NULL,
//                                  &feval->events_yt[4 * i + 3]));
//
//    events[i] = feval->events_yt[4 * i + 3];
//
//    uninit_avector(sub_yt);
//  }
//
//  clWaitForEvents(oclwrk->num_wrk_pkgs, events);
//#endif
}

static void
nearfield_3d_partial_gca(pcgreencross gc,
                         const uint   rows,
                         const uint   cols,
                         const uint   *ridx,
                         const uint   *cidx,
                         const field  alpha,
                         pcavector    xt,
                         pavector     yt)
{
  if(gc->dim != 3)
  {
    fprintf(stderr, "error: trying to call nearfield_3d_partial_gca for "
                    "problems other than in 3D.");
    exit(1);
  }

  kernel_func3d kernel = gc->kernel_3d;

  pcbem3d    bem     = (pbem3d) gc->bem;

  const uint n_dist = bem->sq->n_dist;

  const real (*v)[3] = (const real(*)[3]) ((psurface3d) gc->geom)->x;
  const uint (*p)[3] = (const uint(*)[3]) ((psurface3d) gc->geom)->t;
  const real *g      = ((psurface3d) gc->geom)->g;

  for(uint i = 0; i < rows; ++i)
  {
    const uint ii  = (ridx == NULL ? i : ridx[i]);
    const real tmp = g[ii] * bem->kernel_const;

    real result = r_zero;

    for(uint j = 0; j < cols; ++j)
    {
      const uint jj     = (cidx == NULL ? j : cidx[j]);
      const real factor = tmp * g[jj];

      real base;
      uint nq;
      uint px[3], py[3];
      real *xq, *yq, *wq;

      /* Choose quadrature rule, ensuring symmetry of the matrix. */

      select_quadrature_singquad2d(bem->sq,
                                   p[ii],
                                   p[jj],
                                   px,
                                   py,
                                   &xq,
                                   &yq,
                                   &wq,
                                   &nq,
                                   &base);

      /* Only compute matrix entries where the corresponding triangulars have an
       * identical vertex, edge or are identical in themselfe. */
      if((nq - n_dist) > 0)
      {
        const uint vnq = ROUNDUP(nq, VREAL);

        wq += vnq * 9;

        /* Copy permuted vertex numbers. */

        px[0] = p[ii][px[0]];
        px[1] = p[ii][px[1]];
        px[2] = p[ii][px[2]];

        py[0] = p[jj][py[0]];
        py[1] = p[jj][py[1]];
        py[2] = p[jj][py[2]];

        /* Copy permuted vertices */

        const real x[3][3] =
          { { v[px[0]][0], v[px[0]][1], v[px[0]][2] },
            { v[px[1]][0], v[px[1]][1], v[px[1]][2] },
            { v[px[2]][0], v[px[2]][1], v[px[2]][2] } };

        const real y[3][3] =
          { { v[py[0]][0], v[py[0]][1], v[py[0]][2] },
            { v[py[1]][0], v[py[1]][1], v[py[1]][2] },
            { v[py[2]][0], v[py[2]][1], v[py[2]][2] } };

        real sum = base;

        uint q   = 0;

#ifdef __AVX__
#ifdef USE_FLOAT
        for(; (q + 8) <= nq; q += 8)
        {
          const __m256 one = _mm256_set1_ps(r_one);
          const __m256 w   = _mm256_load_ps(wq + q);

          __m256 a1 = _mm256_load_ps(xq + q);
          __m256 a2 = _mm256_load_ps(xq + vnq + q);

          __m256 xx[3] = { _mm256_set1_ps(r_zero),
                           _mm256_set1_ps(r_zero),
                           _mm256_set1_ps(r_zero) };

          __m256 scalar;

          scalar = _mm256_set1_ps(x[0][0]);
          xx[0]  = _mm256_mul_ps(_mm256_sub_ps(one, a1), scalar);
          scalar = _mm256_set1_ps(x[1][0]);
          xx[0]  = _mm256_add_ps(_mm256_mul_ps(_mm256_sub_ps(a1, a2), scalar), xx[0]);
          scalar = _mm256_set1_ps(x[2][0]);
          xx[0]  = _mm256_add_ps(_mm256_mul_ps(a2, scalar), xx[0]);

          scalar = _mm256_set1_ps(x[0][1]);
          xx[1]  = _mm256_mul_ps(_mm256_sub_ps(one, a1), scalar);
          scalar = _mm256_set1_ps(x[1][1]);
          xx[1]  = _mm256_add_ps(_mm256_mul_ps(_mm256_sub_ps(a1, a2), scalar), xx[1]);
          scalar = _mm256_set1_ps(x[2][1]);
          xx[1]  = _mm256_add_ps(_mm256_mul_ps(a2, scalar), xx[1]);

          scalar = _mm256_set1_ps(x[0][2]);
          xx[2]  = _mm256_mul_ps(_mm256_sub_ps(one, a1), scalar);
          scalar = _mm256_set1_ps(x[1][2]);
          xx[2]  = _mm256_add_ps(_mm256_mul_ps(_mm256_sub_ps(a1, a2), scalar), xx[2]);
          scalar = _mm256_set1_ps(x[2][2]);
          xx[2]  = _mm256_add_ps(_mm256_mul_ps(a2, scalar), xx[2]);

          a1 = _mm256_load_ps(yq + q);
          a2 = _mm256_load_ps(yq + vnq + q);

          __m256 yy[3] = { _mm256_set1_ps(r_zero),
                           _mm256_set1_ps(r_zero),
                           _mm256_set1_ps(r_zero) };

          scalar = _mm256_set1_ps(y[0][0]);
          yy[0]  = _mm256_mul_ps(_mm256_sub_ps(one, a1), scalar);
          scalar = _mm256_set1_ps(y[1][0]);
          yy[0]  = _mm256_add_ps(_mm256_mul_ps(_mm256_sub_ps(a1, a2), scalar), yy[0]);
          scalar = _mm256_set1_ps(y[2][0]);
          yy[0]  = _mm256_add_ps(_mm256_mul_ps(a2, scalar), yy[0]);

          scalar = _mm256_set1_ps(y[0][1]);
          yy[1]  = _mm256_mul_ps(_mm256_sub_ps(one, a1), scalar);
          scalar = _mm256_set1_ps(y[1][1]);
          yy[1]  = _mm256_add_ps(_mm256_mul_ps(_mm256_sub_ps(a1, a2), scalar), yy[1]);
          scalar = _mm256_set1_ps(y[2][1]);
          yy[1]  = _mm256_add_ps(_mm256_mul_ps(a2, scalar), yy[1]);

          scalar = _mm256_set1_ps(y[0][2]);
          yy[2]  = _mm256_mul_ps(_mm256_sub_ps(one, a1), scalar);
          scalar = _mm256_set1_ps(y[1][2]);
          yy[2]  = _mm256_add_ps(_mm256_mul_ps(_mm256_sub_ps(a1, a2), scalar), yy[2]);
          scalar = _mm256_set1_ps(y[2][2]);
          yy[2]  = _mm256_add_ps(_mm256_mul_ps(a2, scalar), yy[2]);

          __m256 kernel_result;

          gc->kernel_simd_3d(xx, yy, NULL, NULL, NULL, &kernel_result, NULL);

          kernel_result = _mm256_mul_ps(w, kernel_result);

          real results[8];

          _mm256_store_ps(results, kernel_result);

          sum += results[0] + results[1] + results[2] + results[3] +
                 results[4] + results[5] + results[6] + results[7];
        }
#endif // USE_FLOAT
#endif // __AVX__
        
        for(; q < (nq - n_dist); ++q)
        {
          real a1 = xq[q];
          real a2 = xq[q + nq];

          const real xx[3] =
            { (r_one - a1) * x[0][0] + (a1 - a2) * x[1][0] + a2 * x[2][0],
              (r_one - a1) * x[0][1] + (a1 - a2) * x[1][1] + a2 * x[2][1],
              (r_one - a1) * x[0][2] + (a1 - a2) * x[1][2] + a2 * x[2][2] };

          a1 = yq[q];
          a2 = yq[q + nq];

          const real yy[3] =
            { (r_one - a1) * y[0][0] + (a1 - a2) * y[1][0] + a2 * y[2][0],
              (r_one - a1) * y[0][1] + (a1 - a2) * y[1][1] + a2 * y[2][1],
              (r_one - a1) * y[0][2] + (a1 - a2) * y[1][2] + a2 * y[2][2] };

          sum += wq[q] * kernel(xx, yy, NULL, NULL, NULL);
        }

        const real factor2 =
          ((bem->alpha != 0.0) && (ii == jj))
            ? 0.5 * bem->alpha * g[ii]
            : r_zero;

        result += (sum * factor + factor2) * getentry_avector(xt, j);
      }
    }

    yt->v[i] += alpha * result;
  }
}

static void
nearfield_3d_partial_min_id_edge(pcgreencross gc,
                                 const uint   rows,
                                 const uint   cols,
                                 const uint   *ridx,
                                 const uint   *cidx,
                                 const field  alpha,
                                 pcavector    xt,
                                 pavector     yt)
{
  if(gc->dim != 3)
  {
    fprintf(stderr, "error: trying to call nearfield_3d_partial_gca for "
      "problems other than in 3D.");
    exit(1);
  }

  kernel_func3d kernel = gc->kernel_3d;

  pcbem3d    bem     = (pbem3d) gc->bem;

//  const uint n_vert = bem->sq->n_vert;

  const real (*v)[3] = (const real(*)[3]) ((psurface3d) gc->geom)->x;
  const real *g      = ((psurface3d) gc->geom)->g;

  uint (*p)[3] = ((psurface3d) gc->geom)->t;

  for(uint i = 0; i < rows; ++i)
  {
    const uint ii  = (ridx == NULL ? i : ridx[i]);
    const real tmp = g[ii] * bem->kernel_const;

    real result = r_zero;

    for(uint j = 0; j < cols; ++j)
    {
      const uint jj = (cidx == NULL ? j : cidx[j]);

      uint id_verts = 0;

      if(gc->dim == 2)
        fprintf(stderr, "error: nearfield_3d_partial_min_id_edge not yet "
          "implemented for 2D problems!\n");
      else
        id_verts = fast_select_quadrature(p, ii, jj);

      if(id_verts > 1)
      {
        const real factor = tmp * g[jj];

        real sum;
        uint px[3], py[3];
        real *xq, *yq, *wq;

        select_quadrature_singquadgca(gc->sq_partial_min_edge,
                                      p[ii],
                                      p[jj],
                                      px,
                                      py,
                                      &xq,
                                      &yq,
                                      &wq,
                                      &sum);

        /* Copy permuted vertex numbers. */

        px[0] = p[ii][px[0]];
        px[1] = p[ii][px[1]];
        px[2] = p[ii][px[2]];

        py[0] = p[jj][py[0]];
        py[1] = p[jj][py[1]];
        py[2] = p[jj][py[2]];

        /* Copy permuted vertices */

        const real x[3][3] =
          { { v[px[0]][0], v[px[0]][1], v[px[0]][2] },
            { v[px[1]][0], v[px[1]][1], v[px[1]][2] },
            { v[px[2]][0], v[px[2]][1], v[px[2]][2] } };

        const real y[3][3] =
          { { v[py[0]][0], v[py[0]][1], v[py[0]][2] },
            { v[py[1]][0], v[py[1]][1], v[py[1]][2] },
            { v[py[2]][0], v[py[2]][1], v[py[2]][2] } };

        const uint nq = gc->sq_partial_min_edge->nq;

        uint q = 0;

#ifdef __AVX__
#ifdef USE_FLOAT
        for(; (q + 8) <= nq; q += 8)
        {
          const __m256 one = _mm256_set1_ps(r_one);
          const __m256 w   = _mm256_load_ps(wq + q);

          __m256 a1 = _mm256_load_ps(xq + q);
          __m256 a2 = _mm256_load_ps(xq + nq + q);

          __m256 xx[3] = { _mm256_set1_ps(r_zero),
                           _mm256_set1_ps(r_zero),
                           _mm256_set1_ps(r_zero) };

          __m256 scalar;

          scalar = _mm256_set1_ps(x[0][0]);
          xx[0]  = _mm256_mul_ps(_mm256_sub_ps(one, a1), scalar);
          scalar = _mm256_set1_ps(x[1][0]);
          xx[0]  = _mm256_add_ps(_mm256_mul_ps(_mm256_sub_ps(a1, a2), scalar), xx[0]);
          scalar = _mm256_set1_ps(x[2][0]);
          xx[0]  = _mm256_add_ps(_mm256_mul_ps(a2, scalar), xx[0]);

          scalar = _mm256_set1_ps(x[0][1]);
          xx[1]  = _mm256_mul_ps(_mm256_sub_ps(one, a1), scalar);
          scalar = _mm256_set1_ps(x[1][1]);
          xx[1]  = _mm256_add_ps(_mm256_mul_ps(_mm256_sub_ps(a1, a2), scalar), xx[1]);
          scalar = _mm256_set1_ps(x[2][1]);
          xx[1]  = _mm256_add_ps(_mm256_mul_ps(a2, scalar), xx[1]);

          scalar = _mm256_set1_ps(x[0][2]);
          xx[2]  = _mm256_mul_ps(_mm256_sub_ps(one, a1), scalar);
          scalar = _mm256_set1_ps(x[1][2]);
          xx[2]  = _mm256_add_ps(_mm256_mul_ps(_mm256_sub_ps(a1, a2), scalar), xx[2]);
          scalar = _mm256_set1_ps(x[2][2]);
          xx[2]  = _mm256_add_ps(_mm256_mul_ps(a2, scalar), xx[2]);

          a1 = _mm256_load_ps(yq + q);
          a2 = _mm256_load_ps(yq + nq + q);

          __m256 yy[3] = { _mm256_set1_ps(r_zero),
                           _mm256_set1_ps(r_zero),
                           _mm256_set1_ps(r_zero) };

          scalar = _mm256_set1_ps(y[0][0]);
          yy[0]  = _mm256_mul_ps(_mm256_sub_ps(one, a1), scalar);
          scalar = _mm256_set1_ps(y[1][0]);
          yy[0]  = _mm256_add_ps(_mm256_mul_ps(_mm256_sub_ps(a1, a2), scalar), yy[0]);
          scalar = _mm256_set1_ps(y[2][0]);
          yy[0]  = _mm256_add_ps(_mm256_mul_ps(a2, scalar), yy[0]);

          scalar = _mm256_set1_ps(y[0][1]);
          yy[1]  = _mm256_mul_ps(_mm256_sub_ps(one, a1), scalar);
          scalar = _mm256_set1_ps(y[1][1]);
          yy[1]  = _mm256_add_ps(_mm256_mul_ps(_mm256_sub_ps(a1, a2), scalar), yy[1]);
          scalar = _mm256_set1_ps(y[2][1]);
          yy[1]  = _mm256_add_ps(_mm256_mul_ps(a2, scalar), yy[1]);

          scalar = _mm256_set1_ps(y[0][2]);
          yy[2]  = _mm256_mul_ps(_mm256_sub_ps(one, a1), scalar);
          scalar = _mm256_set1_ps(y[1][2]);
          yy[2]  = _mm256_add_ps(_mm256_mul_ps(_mm256_sub_ps(a1, a2), scalar), yy[2]);
          scalar = _mm256_set1_ps(y[2][2]);
          yy[2]  = _mm256_add_ps(_mm256_mul_ps(a2, scalar), yy[2]);

          __m256 kernel_result;

          gc->kernel_simd_3d(xx, yy, NULL, NULL, NULL, &kernel_result, NULL);

          kernel_result = _mm256_mul_ps(w, kernel_result);

          real results[8];

          _mm256_store_ps(results, kernel_result);

          sum += results[0] + results[1] + results[2] + results[3] +
                 results[4] + results[5] + results[6] + results[7];
        }
#endif // USE_FLOAT
#endif // __AVX__

        for(; q < nq; ++q)
        {
          real a1 = xq[q];
          real a2 = xq[q + nq];

          const real xx[3] =
            { (r_one - a1) * x[0][0] + (a1 - a2) * x[1][0] + a2 * x[2][0],
              (r_one - a1) * x[0][1] + (a1 - a2) * x[1][1] + a2 * x[2][1],
              (r_one - a1) * x[0][2] + (a1 - a2) * x[1][2] + a2 * x[2][2] };

          a1 = yq[q];
          a2 = yq[q + nq];

          const real yy[3] =
            { (r_one - a1) * y[0][0] + (a1 - a2) * y[1][0] + a2 * y[2][0],
              (r_one - a1) * y[0][1] + (a1 - a2) * y[1][1] + a2 * y[2][1],
              (r_one - a1) * y[0][2] + (a1 - a2) * y[1][2] + a2 * y[2][2] };

          sum += wq[q] * kernel(xx, yy, NULL, NULL, NULL);
        }

        const real factor2 =
          ((bem->alpha != 0.0) && (ii == jj))
          ? 0.5 * bem->alpha * g[ii]
          : r_zero;

        result += (sum * factor + factor2) * getentry_avector(xt, j);
      }

//      const real factor = tmp * g[jj];
//
//      real base;
//      uint nq;
//      uint px[3], py[3];
//      real *xq, *yq, *wq;
//
//       Choose quadrature rule, ensuring symmetry of the matrix.
//
//      select_quadrature_singquad2d(bem->sq,
//                                   p[ii],
//                                   p[jj],
//                                   px,
//                                   py,
//                                   &xq,
//                                   &yq,
//                                   &wq,
//                                   &nq,
//                                   &base);
//
//      const int num_quad_points = nq - n_vert;
//
//       Only compute matrix entries where the corresponding triangulars have an
//       * identical edge or are identical in themselves.
//      if(num_quad_points > 0)
//      {
//        const uint vnq = ROUNDUP(nq, VREAL);
//
//        wq += vnq * 9;
//
//         Copy permuted vertex numbers.
//
//        px[0] = p[ii][px[0]];
//        px[1] = p[ii][px[1]];
//        px[2] = p[ii][px[2]];
//
//        py[0] = p[jj][py[0]];
//        py[1] = p[jj][py[1]];
//        py[2] = p[jj][py[2]];
//
//         Copy permuted vertices
//
//        const real x[3][3] =
//          { { v[px[0]][0], v[px[0]][1], v[px[0]][2] },
//            { v[px[1]][0], v[px[1]][1], v[px[1]][2] },
//            { v[px[2]][0], v[px[2]][1], v[px[2]][2] } };
//
//        const real y[3][3] =
//          { { v[py[0]][0], v[py[0]][1], v[py[0]][2] },
//            { v[py[1]][0], v[py[1]][1], v[py[1]][2] },
//            { v[py[2]][0], v[py[2]][1], v[py[2]][2] } };
//
//        real sum = base;
//
//        uint q   = 0;

//#ifdef USE_SIMD
//        for(q = 0; (i + q) < num_quad_points; q += VREAL)
//        {
//          const vreal one = vset1(r_one);
//          const vreal w   = vload(wq + q);
//
//          vreal a1 = vload(xq + q);
//          vreal a2 = vload(xq + nq + q);
//
//          vreal xx[3] = { vset1(r_zero), vset1(r_zero), vset1(r_zero) };
//
//          vreal scalar;
//
//          scalar = vset1(x[0][0]); xx[0] = vfmadd(vsub(one, a1), scalar, xx[0]);
//          scalar = vset1(x[1][0]); xx[0] = vfmadd(vsub(a1,  a2), scalar, xx[0]);
//          scalar = vset1(x[2][0]); xx[0] = vfmadd(a2,            scalar, xx[0]);
//
//          scalar = vset1(x[0][1]); xx[1] = vfmadd(vsub(one, a1), scalar, xx[1]);
//          scalar = vset1(x[1][1]); xx[1] = vfmadd(vsub(a1,  a2), scalar, xx[1]);
//          scalar = vset1(x[2][1]); xx[1] = vfmadd(a2,            scalar, xx[1]);
//
//          scalar = vset1(x[0][2]); xx[2] = vfmadd(vsub(one, a1), scalar, xx[2]);
//          scalar = vset1(x[1][2]); xx[2] = vfmadd(vsub(a1,  a2), scalar, xx[2]);
//          scalar = vset1(x[2][2]); xx[2] = vfmadd(a2,            scalar, xx[2]);
//
//          a1 = vload(yq + q);
//          a2 = vload(yq + nq + q);
//
//          vreal yy[3] = { vset1(r_zero), vset1(r_zero), vset1(r_zero) };
//
//          scalar = vset1(y[0][0]); yy[0] = vfmadd(vsub(one, a1), scalar, yy[0]);
//          scalar = vset1(y[1][0]); yy[0] = vfmadd(vsub(a1,  a2), scalar, yy[0]);
//          scalar = vset1(y[2][0]); yy[0] = vfmadd(a2,            scalar, yy[0]);
//
//          scalar = vset1(y[0][1]); yy[1] = vfmadd(vsub(one, a1), scalar, yy[1]);
//          scalar = vset1(y[1][1]); yy[1] = vfmadd(vsub(a1,  a2), scalar, yy[1]);
//          scalar = vset1(y[2][1]); yy[1] = vfmadd(a2,            scalar, yy[1]);
//
//          scalar = vset1(y[0][2]); yy[2] = vfmadd(vsub(one, a1), scalar, yy[2]);
//          scalar = vset1(y[1][2]); yy[2] = vfmadd(vsub(a1,  a2), scalar, yy[2]);
//          scalar = vset1(y[2][2]); yy[2] = vfmadd(a2,            scalar, yy[2]);
//
//          vreal kernel_result;
//
//          gc->kernel_simd_3d(xx, yy, NULL, NULL, NULL, &kernel_result, NULL);
//
//          kernel_result = vmul(w, kernel_result);
//
//          real results[VREAL];
//
//          vstore(results, kernel_result);
//
//          for(uint k = (VREAL >> 1); k > 0; k >>= 1)
//            for(uint l = 0; l < k; ++l)
//              results[l] += results[k + l];
//
//          sum += results[0];
//        }
//#endif

//        for(; (int) q < num_quad_points; ++q)
//        {
//          real a1 = xq[q];
//          real a2 = xq[q + nq];
//
//          const real xx[3] =
//            { (r_one - a1) * x[0][0] + (a1 - a2) * x[1][0] + a2 * x[2][0],
//              (r_one - a1) * x[0][1] + (a1 - a2) * x[1][1] + a2 * x[2][1],
//              (r_one - a1) * x[0][2] + (a1 - a2) * x[1][2] + a2 * x[2][2] };
//
//          a1 = yq[q];
//          a2 = yq[q + nq];
//
//          const real yy[3] =
//            { (r_one - a1) * y[0][0] + (a1 - a2) * y[1][0] + a2 * y[2][0],
//              (r_one - a1) * y[0][1] + (a1 - a2) * y[1][1] + a2 * y[2][1],
//              (r_one - a1) * y[0][2] + (a1 - a2) * y[1][2] + a2 * y[2][2] };
//
//          sum += wq[q] * kernel(xx, yy, NULL, NULL, NULL);
//        }
//
//        const real factor2 =
//          ((bem->alpha != 0.0) && (ii == jj))
//          ? 0.5 * bem->alpha * g[ii]
//          : r_zero;
//
//        result += (sum * factor + factor2) * getentry_avector(xt, j);
//      }
//    }
//
//    yt->v[i] += alpha * result;
    }

    yt->v[i] += alpha * result;
  }
}

static void
nearfield_3d_cpu_gca(pcgreencross gc,
                     const uint   rows,
                     const uint   cols,
                     const uint   *ridx,
                     const uint   *cidx,
                     const field  alpha,
                     pcavector    xt,
                     pavector     yt)
{
  if(gc->dim != 3)
  {
    fprintf(stderr, "error: trying to call nearfield_3d_cpu_gca for "
      "problems other than in 3D.");
    exit(1);
  }

  kernel_func3d kernel = gc->kernel_3d;

  pcbem3d    bem     = (pbem3d) gc->bem;

  const uint vdim  = VREAL;
//  const uint nq      = gc->sq_gca->nq;

  const real (*v)[3] = (const real(*)[3]) ((psurface3d) gc->geom)->x;
  const uint (*p)[3] = (const uint(*)[3]) ((psurface3d) gc->geom)->t;
  const real *g      = ((psurface3d) gc->geom)->g;

  for(uint i = 0; i < rows; ++i)
  {
    const uint ii  = (ridx == NULL ? i : ridx[i]);
    const real tmp = g[ii] * bem->kernel_const;

    real result = r_zero;

    for(uint j = 0; j < cols; ++j)
    {
      const uint jj     = (cidx == NULL ? j : cidx[j]);
      const real factor = tmp * g[jj];

      real sum;
      uint nq;
      uint px[3], py[3];
      real *xq, *yq, *wq;

      /* Choose quadrature rule. */

      select_quadrature_singquad2d(bem->sq,
                                   p[ii],
                                   p[jj],
                                   px,
                                   py,
                                   &xq,
                                   &yq,
                                   &wq,
                                   &nq,
                                   &sum);

      const uint vnq = ROUNDUP(nq, vdim);
      wq += 9 * vnq;
//      select_quadrature_singquadgca(gc->sq_gca,
//                                    p[ii],
//                                    p[jj],
//                                    px,
//                                    py,
//                                    &xq,
//                                    &yq,
//                                    &wq,
//                                    &sum);

      /* Copy permuted vertex numbers. */

      px[0] = p[ii][px[0]];
      px[1] = p[ii][px[1]];
      px[2] = p[ii][px[2]];

      py[0] = p[jj][py[0]];
      py[1] = p[jj][py[1]];
      py[2] = p[jj][py[2]];

      /* Copy permuted vertices */

      const real x[3][3] =
        { { v[px[0]][0], v[px[0]][1], v[px[0]][2] },
          { v[px[1]][0], v[px[1]][1], v[px[1]][2] },
          { v[px[2]][0], v[px[2]][1], v[px[2]][2] } };

      const real y[3][3] =
        { { v[py[0]][0], v[py[0]][1], v[py[0]][2] },
          { v[py[1]][0], v[py[1]][1], v[py[1]][2] },
          { v[py[2]][0], v[py[2]][1], v[py[2]][2] } };

        uint q = 0;

#ifdef __AVX__
#ifdef USE_FLOAT
      for(; (q + 8) <= nq; q += 8)
      {
        const __m256 one = _mm256_set1_ps(r_one);
        const __m256 w   = _mm256_load_ps(wq + q);

        __m256 a1 = _mm256_load_ps(xq + q);
        __m256 a2 = _mm256_load_ps(xq + vnq + q);

        __m256 xx[3] = { _mm256_set1_ps(r_zero),
                         _mm256_set1_ps(r_zero),
                         _mm256_set1_ps(r_zero) };

        __m256 scalar;

        scalar = _mm256_set1_ps(x[0][0]);
        xx[0]  = _mm256_mul_ps(_mm256_sub_ps(one, a1), scalar);
        scalar = _mm256_set1_ps(x[1][0]);
        xx[0]  = _mm256_add_ps(_mm256_mul_ps(_mm256_sub_ps(a1, a2), scalar), xx[0]);
        scalar = _mm256_set1_ps(x[2][0]);
        xx[0]  = _mm256_add_ps(_mm256_mul_ps(a2, scalar), xx[0]);

        scalar = _mm256_set1_ps(x[0][1]);
        xx[1]  = _mm256_mul_ps(_mm256_sub_ps(one, a1), scalar);
        scalar = _mm256_set1_ps(x[1][1]);
        xx[1]  = _mm256_add_ps(_mm256_mul_ps(_mm256_sub_ps(a1, a2), scalar), xx[1]);
        scalar = _mm256_set1_ps(x[2][1]);
        xx[1]  = _mm256_add_ps(_mm256_mul_ps(a2, scalar), xx[1]);

        scalar = _mm256_set1_ps(x[0][2]);
        xx[2]  = _mm256_mul_ps(_mm256_sub_ps(one, a1), scalar);
        scalar = _mm256_set1_ps(x[1][2]);
        xx[2]  = _mm256_add_ps(_mm256_mul_ps(_mm256_sub_ps(a1, a2), scalar), xx[2]);
        scalar = _mm256_set1_ps(x[2][2]);
        xx[2]  = _mm256_add_ps(_mm256_mul_ps(a2, scalar), xx[2]);

        a1 = _mm256_load_ps(yq + q);
        a2 = _mm256_load_ps(yq + vnq + q);

        __m256 yy[3] = { _mm256_set1_ps(r_zero),
                         _mm256_set1_ps(r_zero),
                         _mm256_set1_ps(r_zero) };

        scalar = _mm256_set1_ps(y[0][0]);
        yy[0]  = _mm256_mul_ps(_mm256_sub_ps(one, a1), scalar);
        scalar = _mm256_set1_ps(y[1][0]);
        yy[0]  = _mm256_add_ps(_mm256_mul_ps(_mm256_sub_ps(a1, a2), scalar), yy[0]);
        scalar = _mm256_set1_ps(y[2][0]);
        yy[0]  = _mm256_add_ps(_mm256_mul_ps(a2, scalar), yy[0]);

        scalar = _mm256_set1_ps(y[0][1]);
        yy[1]  = _mm256_mul_ps(_mm256_sub_ps(one, a1), scalar);
        scalar = _mm256_set1_ps(y[1][1]);
        yy[1]  = _mm256_add_ps(_mm256_mul_ps(_mm256_sub_ps(a1, a2), scalar), yy[1]);
        scalar = _mm256_set1_ps(y[2][1]);
        yy[1]  = _mm256_add_ps(_mm256_mul_ps(a2, scalar), yy[1]);

        scalar = _mm256_set1_ps(y[0][2]);
        yy[2]  = _mm256_mul_ps(_mm256_sub_ps(one, a1), scalar);
        scalar = _mm256_set1_ps(y[1][2]);
        yy[2]  = _mm256_add_ps(_mm256_mul_ps(_mm256_sub_ps(a1, a2), scalar), yy[2]);
        scalar = _mm256_set1_ps(y[2][2]);
        yy[2]  = _mm256_add_ps(_mm256_mul_ps(a2, scalar), yy[2]);

        __m256 kernel_result;

        gc->kernel_simd_3d(xx, yy, NULL, NULL, NULL, &kernel_result, NULL);

        kernel_result = _mm256_mul_ps(w, kernel_result);

        real results[8];

        _mm256_store_ps(results, kernel_result);

        sum += results[0] + results[1] + results[2] + results[3] +
               results[4] + results[5] + results[6] + results[7];
      }
#endif // USE_FLOAT
#endif // __AVX__

      q = 0; sum = r_zero;

      /* Perform remaining quadrature. */
      for(; q < nq; ++q)
      {
        real a1 = xq[q];
        real a2 = xq[q + vnq];

        const real xx[3] =
          { (r_one - a1) * x[0][0] + (a1 - a2) * x[1][0] + a2 * x[2][0],
            (r_one - a1) * x[0][1] + (a1 - a2) * x[1][1] + a2 * x[2][1],
            (r_one - a1) * x[0][2] + (a1 - a2) * x[1][2] + a2 * x[2][2] };

        a1 = yq[q];
        a2 = yq[q + vnq];

        const real yy[3] =
          { (r_one - a1) * y[0][0] + (a1 - a2) * y[1][0] + a2 * y[2][0],
            (r_one - a1) * y[0][1] + (a1 - a2) * y[1][1] + a2 * y[2][1],
            (r_one - a1) * y[0][2] + (a1 - a2) * y[1][2] + a2 * y[2][2] };

        sum += wq[q] * kernel(xx, yy, NULL, NULL, NULL);
      }

      const real factor2 =
        ((bem->alpha != 0.0) && (ii == jj))
          ? 0.5 * bem->alpha * g[ii]
          : r_zero;

      result += (sum * factor + factor2) * getentry_avector(xt, j);
    }

    yt->v[i] += alpha * result;
  }
}

void
fastaddeval_nearfield_partial_h2matrix_avector_gca(pcgreencross gc,
                                                   field        alpha,
                                                   pavector     xt,
                                                   pavector     yt)
{
  pgcopencl ocl_info = gc->ocl_info_nf;

  avector tmp1, tmp2;

//  #pragma omp parallel for private(tmp1, tmp2)
  for(uint i = 0; i < ocl_info->num_row_leafs; ++i)
  {
    pavector yt1 = init_sub_avector(&tmp1,
                                    yt,
                                    ocl_info->ridx_sizes[i],
                                    ocl_info->ytoffs[i]);

    const uint idx_off = ocl_info->idx_off[i];
    const uint *ridx   = ocl_info->host_ridx + ocl_info->ridx_off[i];

    for(uint j = 0; j < ocl_info->num_h2_leafs_per_row[i]; ++j)
    {
      pavector xt1 = init_sub_avector(&tmp2,
                                      xt,
                                      ocl_info->cidx_sizes[idx_off + j],
                                      ocl_info->xtoffs[i][j]);

      nearfield_3d_partial_gca(gc,
                               yt1->dim,
                               xt1->dim,
                               ridx,
                               ocl_info->host_cidx + ocl_info->cidx_off[i][j],
                               alpha,
                               xt1,
                               yt1);

      uninit_avector(xt1);
    }

    uninit_avector(yt1);
  }
}

void
fastaddeval_nearfield_partial_min_id_edge_gca_ref(pcgreencross gc,
                                                  field        alpha,
                                                  pavector     xt,
                                                  pavector     yt)
{
  pgcopencl ocl_info = gc->ocl_info_nf;

  avector tmp1, tmp2;

//  #pragma omp parallel for private(tmp1, tmp2)
  for(uint i = 0; i < ocl_info->num_row_leafs; ++i)
  {
    pavector yt1 = init_sub_avector(&tmp1,
                                    yt,
                                    ocl_info->ridx_sizes[i],
                                    ocl_info->ytoffs[i]);

    const uint idx_off = ocl_info->idx_off[i];
    const uint *ridx   = ocl_info->host_ridx + ocl_info->ridx_off[i];

    for(uint j = 0; j < ocl_info->num_h2_leafs_per_row[i]; ++j)
    {
      pavector xt1 = init_sub_avector(&tmp2,
                                      xt,
                                      ocl_info->cidx_sizes[idx_off + j],
                                      ocl_info->xtoffs[i][j]);

      nearfield_3d_partial_min_id_edge
        (gc,
         yt1->dim,
         xt1->dim,
         ridx,
         ocl_info->host_cidx + ocl_info->cidx_off[i][j],
         alpha,
         xt1,
         yt1);

      uninit_avector(xt1);
    }

    uninit_avector(yt1);
  }
}

void
fastaddeval_nearfield_partial_min_id_edge_gca(pcgreencross gc,
                                              field        alpha,
                                              pavector     xt,
                                              pavector     yt)
{
  const real kernel_const = ((pbem3d) gc->bem)->kernel_const;
  const real bem_alpha    = ((pbem3d) gc->bem)->alpha;

  const real (*v)[3]      = (const real(*)[3]) ((psurface3d) gc->geom)->x;
  const uint (*p)[3]      = (const uint(*)[3]) ((psurface3d) gc->geom)->t;
  const real *g           = ((psurface3d) gc->geom)->g;

  kernel_func3d kernel    = gc->kernel_3d;

  pgcopencl      ocl_info = gc->ocl_info_nf;
  pintegralinfos iinfos   = gc->iinfos_min_edge;
  psingquadgca   sq       = gc->sq_partial_min_edge;

  for(uint i = 0; i < iinfos->num_integral_grps; ++i)
  {
    const uint ytoff = ocl_info->ytoffs[i];

    const uint *ridx = ocl_info->host_ridx + ocl_info->ridx_off[i];

    const uint *rows = iinfos->rows[i];
    const uint *cidx = iinfos->cidx[i];
    const uint *cols = iinfos->cols[i];

    for(uint j = 0; j < iinfos->num_integrals[i]; ++j)
    {
      const uint ii = ridx[rows[j]];
      const uint jj = cidx[j];

      const real factor = g[ii] * g[jj] * kernel_const;

      real sum;
      uint px[3], py[3];
      real *xq, *yq, *wq;

      select_quadrature_singquadgca(sq,
                                    p[ii],
                                    p[jj],
                                    px,
                                    py,
                                    &xq,
                                    &yq,
                                    &wq,
                                    &sum);

      px[0] = p[ii][px[0]];
      px[1] = p[ii][px[1]];
      px[2] = p[ii][px[2]];

      py[0] = p[jj][py[0]];
      py[1] = p[jj][py[1]];
      py[2] = p[jj][py[2]];

      /* Copy permuted vertices */

      const real x[3][3] =
        { { v[px[0]][0], v[px[0]][1], v[px[0]][2] },
          { v[px[1]][0], v[px[1]][1], v[px[1]][2] },
          { v[px[2]][0], v[px[2]][1], v[px[2]][2] } };

      const real y[3][3] =
        { { v[py[0]][0], v[py[0]][1], v[py[0]][2] },
          { v[py[1]][0], v[py[1]][1], v[py[1]][2] },
          { v[py[2]][0], v[py[2]][1], v[py[2]][2] } };

      for(uint q = 0; q < sq->nq; ++q)
      {
        real a1 = xq[q];
        real a2 = xq[q + sq->nq];

        const real xx[3] =
          { (r_one - a1) * x[0][0] + (a1 - a2) * x[1][0] + a2 * x[2][0],
            (r_one - a1) * x[0][1] + (a1 - a2) * x[1][1] + a2 * x[2][1],
            (r_one - a1) * x[0][2] + (a1 - a2) * x[1][2] + a2 * x[2][2] };

        a1 = yq[q];
        a2 = yq[q + sq->nq];

        const real yy[3] =
          { (r_one - a1) * y[0][0] + (a1 - a2) * y[1][0] + a2 * y[2][0],
            (r_one - a1) * y[0][1] + (a1 - a2) * y[1][1] + a2 * y[2][1],
            (r_one - a1) * y[0][2] + (a1 - a2) * y[1][2] + a2 * y[2][2] };

        sum += wq[q] * kernel(xx, yy, NULL, NULL, NULL);
      }

      const real factor2 =
        ((bem_alpha != 0.0) && (ii == jj))
        ? 0.5 * bem_alpha * g[ii]
        : r_zero;

      const uint yt_idx = ytoff + rows[j];

      setentry_avector(yt,
                       yt_idx,
                       getentry_avector(yt, yt_idx)     +
                       alpha * (sum * factor + factor2) *
                       getentry_avector(xt, cols[j]));
    }
  }
}

void
fastaddeval_nearfield_partial_min_id_vert_gca(pcgreencross gc,
                                              field        alpha,
                                              pavector     xt,
                                              pavector     yt)
{
  const real kernel_const = ((pbem3d) gc->bem)->kernel_const;
  const real bem_alpha    = ((pbem3d) gc->bem)->alpha;

  const real (*v)[3]      = (const real(*)[3]) ((psurface3d) gc->geom)->x;
  const uint (*p)[3]      = (const uint(*)[3]) ((psurface3d) gc->geom)->t;
  const real *g           = ((psurface3d) gc->geom)->g;

  kernel_func3d kernel    = gc->kernel_3d;

  pgcopencl    ocl_info   = gc->ocl_info_nf;
  psingquadgca sq         = gc->sq_partial_min_vert;

  for(uint i = 0; i < ocl_info->num_row_leafs; ++i)
  {
    const uint ytoff        = ocl_info->ytoffs[i];

    const uint *ridx        = ocl_info->host_ridx + ocl_info->ridx_off[i];
    const uint *rows_any_id = ocl_info->rows_any_id[i];
    const uint *cidx_any_id = ocl_info->cidx_any_id[i];
    const uint *cols_any_id = ocl_info->cols_any_id[i];

    for(uint j = 0; j < ocl_info->num_any_id[i]; ++j)
    {
      const uint ii = ridx[rows_any_id[j]];
      const uint jj = cidx_any_id[j];

      const real factor = g[ii] * g[jj] * kernel_const;

      real sum;
      uint px[3], py[3];
      real *xq, *yq, *wq;

      select_quadrature_singquadgca(sq,
                                    p[ii],
                                    p[jj],
                                    px,
                                    py,
                                    &xq,
                                    &yq,
                                    &wq,
                                    &sum);

      px[0] = p[ii][px[0]];
      px[1] = p[ii][px[1]];
      px[2] = p[ii][px[2]];

      py[0] = p[jj][py[0]];
      py[1] = p[jj][py[1]];
      py[2] = p[jj][py[2]];

      /* Copy permuted vertices */

      const real x[3][3] =
        { { v[px[0]][0], v[px[0]][1], v[px[0]][2] },
          { v[px[1]][0], v[px[1]][1], v[px[1]][2] },
          { v[px[2]][0], v[px[2]][1], v[px[2]][2] } };

      const real y[3][3] =
        { { v[py[0]][0], v[py[0]][1], v[py[0]][2] },
          { v[py[1]][0], v[py[1]][1], v[py[1]][2] },
          { v[py[2]][0], v[py[2]][1], v[py[2]][2] } };

      uint q = 0;

#ifdef __AVX__
#ifdef USE_FLOAT
      for(; (q + 8) <= sq->nq; q += 8)
      {
        const __m256 one = _mm256_set1_ps(r_one);
        const __m256 w   = _mm256_load_ps(wq + q);

        __m256 a1 = _mm256_load_ps(xq + q);
        __m256 a2 = _mm256_load_ps(xq + sq->nq + q);

        __m256 xx[3] = { _mm256_set1_ps(r_zero),
                         _mm256_set1_ps(r_zero),
                         _mm256_set1_ps(r_zero) };

        __m256 scalar;

        scalar = _mm256_set1_ps(x[0][0]);
        xx[0]  = _mm256_mul_ps(_mm256_sub_ps(one, a1), scalar);
        scalar = _mm256_set1_ps(x[1][0]);
        xx[0]  = _mm256_add_ps(_mm256_mul_ps(_mm256_sub_ps(a1, a2), scalar), xx[0]);
        scalar = _mm256_set1_ps(x[2][0]);
        xx[0]  = _mm256_add_ps(_mm256_mul_ps(a2, scalar), xx[0]);

        scalar = _mm256_set1_ps(x[0][1]);
        xx[1]  = _mm256_mul_ps(_mm256_sub_ps(one, a1), scalar);
        scalar = _mm256_set1_ps(x[1][1]);
        xx[1]  = _mm256_add_ps(_mm256_mul_ps(_mm256_sub_ps(a1, a2), scalar), xx[1]);
        scalar = _mm256_set1_ps(x[2][1]);
        xx[1]  = _mm256_add_ps(_mm256_mul_ps(a2, scalar), xx[1]);

        scalar = _mm256_set1_ps(x[0][2]);
        xx[2]  = _mm256_mul_ps(_mm256_sub_ps(one, a1), scalar);
        scalar = _mm256_set1_ps(x[1][2]);
        xx[2]  = _mm256_add_ps(_mm256_mul_ps(_mm256_sub_ps(a1, a2), scalar), xx[2]);
        scalar = _mm256_set1_ps(x[2][2]);
        xx[2]  = _mm256_add_ps(_mm256_mul_ps(a2, scalar), xx[2]);

        a1 = _mm256_load_ps(yq + q);
        a2 = _mm256_load_ps(yq + sq->nq + q);

        __m256 yy[3] = { _mm256_set1_ps(r_zero),
                         _mm256_set1_ps(r_zero),
                         _mm256_set1_ps(r_zero) };

        scalar = _mm256_set1_ps(y[0][0]);
        yy[0]  = _mm256_mul_ps(_mm256_sub_ps(one, a1), scalar);
        scalar = _mm256_set1_ps(y[1][0]);
        yy[0]  = _mm256_add_ps(_mm256_mul_ps(_mm256_sub_ps(a1, a2), scalar), yy[0]);
        scalar = _mm256_set1_ps(y[2][0]);
        yy[0]  = _mm256_add_ps(_mm256_mul_ps(a2, scalar), yy[0]);

        scalar = _mm256_set1_ps(y[0][1]);
        yy[1]  = _mm256_mul_ps(_mm256_sub_ps(one, a1), scalar);
        scalar = _mm256_set1_ps(y[1][1]);
        yy[1]  = _mm256_add_ps(_mm256_mul_ps(_mm256_sub_ps(a1, a2), scalar), yy[1]);
        scalar = _mm256_set1_ps(y[2][1]);
        yy[1]  = _mm256_add_ps(_mm256_mul_ps(a2, scalar), yy[1]);

        scalar = _mm256_set1_ps(y[0][2]);
        yy[2]  = _mm256_mul_ps(_mm256_sub_ps(one, a1), scalar);
        scalar = _mm256_set1_ps(y[1][2]);
        yy[2]  = _mm256_add_ps(_mm256_mul_ps(_mm256_sub_ps(a1, a2), scalar), yy[2]);
        scalar = _mm256_set1_ps(y[2][2]);
        yy[2]  = _mm256_add_ps(_mm256_mul_ps(a2, scalar), yy[2]);

        __m256 kernel_result;

        gc->kernel_simd_3d(xx, yy, NULL, NULL, NULL, &kernel_result, NULL);

        kernel_result = _mm256_mul_ps(w, kernel_result);

        real results[8];

        _mm256_store_ps(results, kernel_result);

        sum += results[0] + results[1] + results[2] + results[3] +
               results[4] + results[5] + results[6] + results[7];
      }
#endif // USE_FLOAT
#endif // __AVX__

      for(; q < sq->nq; ++q)
      {
        real a1 = xq[q];
        real a2 = xq[q + sq->nq];

        const real xx[3] =
          { (r_one - a1) * x[0][0] + (a1 - a2) * x[1][0] + a2 * x[2][0],
            (r_one - a1) * x[0][1] + (a1 - a2) * x[1][1] + a2 * x[2][1],
            (r_one - a1) * x[0][2] + (a1 - a2) * x[1][2] + a2 * x[2][2] };

        a1 = yq[q];
        a2 = yq[q + sq->nq];

        const real yy[3] =
          { (r_one - a1) * y[0][0] + (a1 - a2) * y[1][0] + a2 * y[2][0],
            (r_one - a1) * y[0][1] + (a1 - a2) * y[1][1] + a2 * y[2][1],
            (r_one - a1) * y[0][2] + (a1 - a2) * y[1][2] + a2 * y[2][2] };

        sum += wq[q] * kernel(xx, yy, NULL, NULL, NULL);
      }

      const real factor2 =
        ((bem_alpha != 0.0) && (ii == jj))
        ? 0.5 * bem_alpha * g[ii]
        : r_zero;

      const uint yt_idx = ytoff + rows_any_id[j];

      setentry_avector(yt,
                       yt_idx,
                       getentry_avector(yt, yt_idx)     +
                       alpha * (sum * factor + factor2) *
                       getentry_avector(xt, cols_any_id[j]));
    }
  }
}

void
fastaddeval_nearfield_cpu_h2matrix_avectors_gca(pgreencross gc,
                                                field       alpha,
                                                pavector    xt,
                                                pavector    yt)
{
  pgcopencl ocl_info = gc->ocl_info_nf;

  avector tmp1, tmp2;

  #pragma omp parallel for num_threads(8) private(tmp1, tmp2)
  for(uint i = 0; i < ocl_info->num_row_leafs; ++i)
  {
    pavector yt1 = init_sub_avector(&tmp1,
                                    yt,
                                    ocl_info->ridx_sizes[i],
                                    ocl_info->ytoffs[i]);

    const uint idx_off = ocl_info->idx_off[i];
    const uint *ridx   = ocl_info->host_ridx + ocl_info->ridx_off[i];

    for(uint j = 0; j < ocl_info->num_h2_leafs_per_row[i]; ++j)
    {
      pavector xt1 = init_sub_avector(&tmp2,
                                      xt,
                                      ocl_info->cidx_sizes[idx_off + j],
                                      ocl_info->xtoffs[i][j]);

      nearfield_3d_cpu_gca(gc,
                           yt1->dim,
                           xt1->dim,
                           ridx,
                           ocl_info->host_cidx + ocl_info->cidx_off[i][j],
                           alpha,
                           xt1,
                           yt1);

      uninit_avector(xt1);
    }

    uninit_avector(yt1);
  }
}

void
fastaddeval_farfield_cpu_h2matrix_avectors_gca(pcgreencross gca,
                                               const field  alpha,
                                               pavector     xt,
                                               pavector     yt)
{
  pgcopencl gcocl = gca->gcocl;

  avector   tmp1, tmp2;

//  #pragma omp parallel for num_threads(4) private(tmp1, tmp2)
  for(uint i = 0; i < gcocl->num_row_leafs; ++i)
  {
    pavector yt1 = init_sub_avector(&tmp1,
                                    yt,
                                    gcocl->ridx_sizes[i],
                                    gcocl->ytoffs[i]);

    const uint idx_off = gcocl->idx_off[i];

    for (uint j = 0; j < gcocl->num_h2_leafs_per_row[i]; ++j)
    {
      const uint cidx_size = gcocl->cidx_sizes[idx_off + j];

      pavector xt1 = init_sub_avector(&tmp2,
                                      xt,
                                      cidx_size,
                                      gcocl->xtoffs[i][j]);

      nearfield_3d_cpu_gca(gca,
                           gca->gcocl->ridx_sizes[i],
                           cidx_size,
                           gca->gcocl->host_ridx + gca->gcocl->ridx_off[i],
                           gca->gcocl->host_cidx + gca->gcocl->cidx_off[i][j],
                           alpha,
                           xt1,
                           yt1);

      uninit_avector(xt1);
    }

    uninit_avector(yt1);
  }
}

void
fastaddeval_h2matrix_avector_greencross(pgreencross gca,
                                        field       alpha,
                                        pch2matrix  H2,
                                        pavector    xt,
			                                  pavector    yt)
{
  (void) H2;

  pcfastaddevalgca feval   = gca->feval;
  pcoclworkpgs     oclwrk  = gca->oclwrk;

  cl_command_queue *queues = ocl_system.queues;

  for(uint i = 0; i < oclwrk->num_wrk_pkgs; ++i)
  {
    CL_CHECK(clEnqueueWriteBuffer
               (queues[i * ocl_system.queues_per_device],
                feval->buf_xt[i],
                CL_FALSE,
                0,
                xt->dim * sizeof(real),
                xt->v,
                0,
                NULL,
                &feval->events_xt[i]));

    for(uint k = 0; k < ocl_system.queues_per_device; ++k)
    {
      cl_kernel kernel = gca->ocl_kernels[k + i * ocl_system.queues_per_device];

      CL_CHECK(clSetKernelArg(kernel, 23, sizeof(real), &alpha));

      kernel = gca->ocl_kernels[k + i * ocl_system.queues_per_device +
                                ocl_system.num_devices *
                                ocl_system.queues_per_device];

      CL_CHECK(clSetKernelArg(kernel, 23, sizeof(real), &alpha));

      kernel = gca->ocl_kernels[k + i * ocl_system.queues_per_device +
                                2 * ocl_system.num_devices *
                                ocl_system.queues_per_device];

      CL_CHECK(clSetKernelArg(kernel, 28, sizeof(real), &alpha));

      kernel = gca->ocl_kernels[k + i * ocl_system.queues_per_device +
                                3 * ocl_system.num_devices *
                                ocl_system.queues_per_device];

      CL_CHECK(clSetKernelArg(kernel, 28, sizeof(real), &alpha));
    }
  }

//#ifndef USE_OPENMP
  fastaddeval_farfield_h2matrix_avector_gca(gca);

  fastaddeval_nearfield_common_h2matrix_avector_gca(gca);

  fastaddeval_nearfield_min_vert_h2matrix_avector_gca(gca);

  fastaddeval_nearfield_min_edge_h2matrix_avector_gca(gca);

  printf("yt->dim: %u\n", yt->dim);

  for(uint i = 0; i < oclwrk->num_wrk_pkgs; ++i)
  {
    const uint offset   = oclwrk->first_idx_of_pkgs[i];
    const uint num_comp = oclwrk->last_idx_of_pkg[i] - offset;

    avector  tmp;

    pavector sub_yt = init_sub_avector(&tmp, yt, num_comp, offset);

    clFinish(queues[i * ocl_system.queues_per_device]);

    CL_CHECK(clEnqueueReadBuffer(queues[i * ocl_system.queues_per_device],
                                 feval->buf_yt[i],
                                 CL_FALSE,
                                 offset * sizeof(real),
                                 num_comp * sizeof(real),
                                 sub_yt->v,
                                 0,
                                 NULL,
                                 &feval->events_yt[i]));

    uninit_avector(sub_yt);
  }

  clWaitForEvents(oclwrk->num_wrk_pkgs, feval->events_yt);
//#else
//
//  pavector yt_1 = new_zero_avector(yt->dim);
//  pavector yt_2 = new_zero_avector(yt->dim);
//  pavector yt_3 = new_zero_avector(yt->dim);
//  pavector yt_4 = new_zero_avector(yt->dim);
//
//  fastaddeval_farfield_h2matrix_avector_gca(gca, yt_1);
//
//  fastaddeval_nearfield_common_h2matrix_avector_gca(gca, yt_2);
//
//  fastaddeval_nearfield_min_vert_h2matrix_avector_gca(gca, yt_3);
//
//  fastaddeval_nearfield_min_edge_h2matrix_avector_gca(gca, yt_4);
//
//  add_avector(f_one, yt_3, yt_1);
//  add_avector(f_one, yt_4, yt_2);
//
//  add_avector(f_one, yt_2, yt_1);
//
//  add_avector(f_one, yt_1, yt);
//
//  del_avector(yt_4);
//  del_avector(yt_3);
//  del_avector(yt_2);
//  del_avector(yt_1);
//
//#endif

//  fastaddeval_farfield_cpu_h2matrix_avectors_gca(gca, alpha, xt, yt);
//  fastaddeval_nearfield_cpu_h2matrix_avectors_gca(gca, alpha, xt, yt);
}

void
addeval_h2matrix_avector_greencross(pgreencross gc,
                                    field       alpha,
                                    pch2matrix  H2,
                                    pcavector   x,
                                    pavector    y)
{
  pavector  xt, yt;

  xt = new_coeffs_clusterbasis_avector(H2->cb);
  yt = new_coeffs_clusterbasis_avector(H2->rb);

  clear_avector(yt);

  for(uint i = 0; i < gc->oclwrk->num_wrk_pkgs; ++i)
  {
//#ifndef USE_OPENMP
    CL_CHECK(clEnqueueWriteBuffer
               (ocl_system.queues[i * ocl_system.queues_per_device],
                gc->feval->buf_yt[i],
                CL_FALSE,
                0,
                yt->dim * sizeof(real),
                yt->v,
                0,
                NULL,
                &gc->feval->events_yt[i]));
//#else
//    CL_CHECK(clEnqueueWriteBuffer
//               (ocl_system.queues[i * ocl_system.queues_per_device],
//                gc->feval->buf_yt_ff[i],
//                CL_FALSE,
//                0,
//                yt->dim * sizeof(real),
//                yt->v,
//                0,
//                NULL,
//                &gc->feval->events_yt[4 * i]));
//
//    CL_CHECK(clEnqueueWriteBuffer
//               (ocl_system.queues[i * ocl_system.queues_per_device],
//                gc->feval->buf_yt_nf_common[i],
//                CL_FALSE,
//                0,
//                yt->dim * sizeof(real),
//                yt->v,
//                0,
//                NULL,
//                &gc->feval->events_yt[4 * i + 1]));
//
//    CL_CHECK(clEnqueueWriteBuffer
//               (ocl_system.queues[i * ocl_system.queues_per_device],
//                gc->feval->buf_yt_nf_min_vert[i],
//                CL_FALSE,
//                0,
//                yt->dim * sizeof(real),
//                yt->v,
//                0,
//                NULL,
//                &gc->feval->events_yt[4 * i + 2]));
//    CL_CHECK(clEnqueueWriteBuffer
//               (ocl_system.queues[i * ocl_system.queues_per_device],
//                gc->feval->buf_yt_nf_min_edge[i],
//                CL_FALSE,
//                0,
//                yt->dim * sizeof(real),
//                yt->v,
//                0,
//                NULL,
//                &gc->feval->events_yt[4 * i + 3]));
//#endif
  }

  forward_clusterbasis_avector(H2->cb, x, xt);

  fastaddeval_h2matrix_avector_greencross(gc, alpha, H2, xt, yt);

  backward_clusterbasis_avector(H2->rb, yt, y);

  del_avector(yt);
  del_avector(xt);
}

void
addevaltrans_h2matrix_avector_greencross(pgreencross gc,
                                         field       alpha,
                                         pch2matrix  H2,
                                         pcavector   x,
			                                   pavector    y)
{
  (void) gc;

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
mvm_h2matrix_avector_greencross(pgreencross  gc,
                                field        alpha,
                                bool         h2trans,
                                pch2matrix   H2,
                                pcavector    x,
		                            pavector     y)
{
  if (h2trans)
    addevaltrans_h2matrix_avector_greencross(gc, alpha, H2, x, y);
  else
    addeval_h2matrix_avector_greencross(gc, alpha, H2, x, y);
}
