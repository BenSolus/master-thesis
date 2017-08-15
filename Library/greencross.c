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
#include "bem3d.h"
#include "ocl_system.h"

#include <string.h>

const size_t OCL_LOCAL_WORK_GROUP_SIZE = 64;

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

struct _oclwork
{
  uint num_wrk_pkgs;

  uint *wrk_per_pkg;

  uint *num_rows_per_pkg;

  uint **rows_per_pkg;
};

/** @} */

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

void
init_oclwork(poclwork oclwrk)
{
  oclwrk->num_wrk_pkgs     = 0;

  oclwrk->wrk_per_pkg      = NULL;

  oclwrk->num_rows_per_pkg = NULL;

  oclwrk->rows_per_pkg     = NULL;
}

void
uninit_oclwork(poclwork oclwrk)
{
  if(oclwrk->wrk_per_pkg != NULL)
  {
    freemem(oclwrk->wrk_per_pkg);
    oclwrk->wrk_per_pkg = NULL;
  }

  if(oclwrk->num_rows_per_pkg)
  {
    freemem(oclwrk->num_rows_per_pkg);
    oclwrk->num_rows_per_pkg = NULL;
  }

  if(oclwrk->rows_per_pkg != NULL)
  {
    for(uint i = 0; i < oclwrk->num_wrk_pkgs; ++i)
      if(oclwrk->rows_per_pkg[i] != NULL)
      {
        freemem(oclwrk->rows_per_pkg[i]);
        oclwrk->rows_per_pkg[i] = NULL;
      }

    freemem(oclwrk->rows_per_pkg);
    oclwrk->rows_per_pkg = NULL;
  }
}

void
init_greencross(pgreencross gc, uint dim)
{
  assert(gc != NULL);
  assert((dim >= greencross_min_dim) && (dim <= greencross_max_dim));

  /* Scalar components */
  gc->dim     = dim;
  gc->n       = 0;
  gc->m       = 0;
  gc->K       = 0;
  gc->accur   = 0;

  gc->gcocl   = (pgcopencl) allocmem(sizeof(gcopencl));

  init_gcopencl(gc->gcocl);

  gc->oclwrk  = (poclwork) allocmem(sizeof(oclwork));

  init_oclwork(gc->oclwrk);

  /* Pointer components */
  gc->geom    = NULL;
  gc->buf_x   =
    (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));
  gc->buf_p   =
    (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));
  gc->event_p =
    (cl_event *) calloc(ocl_system.num_devices, sizeof(cl_event));
  gc->buf_g   =
    (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));
  gc->bem     = NULL;
  gc->idx     = NULL;
  gc->rb      = NULL;
  gc->cb      = NULL;
  gc->rc      = NULL;
  gc->cc      = NULL;
  gc->buf_qx  =
    (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));
  gc->buf_qy  =
    (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));
  gc->buf_w   =
    (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));
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

  if(gc->buf_qx != NULL)
  {
    for(uint i = 0; i < ocl_system.num_devices; ++i)
      if(gc->buf_qx != NULL)
      {
        CL_CHECK(clReleaseMemObject(gc->buf_qx[i]));
        gc->buf_qx[i] = NULL;
      }

    freemem(gc->buf_qx);
    gc->buf_qx = NULL;
  }

  if(gc->buf_qy != NULL)
  {
    for(uint i = 0; i < ocl_system.num_devices; ++i)
      if(gc->buf_qy != NULL)
      {
        CL_CHECK(clReleaseMemObject(gc->buf_qy[i]));
        gc->buf_qy[i] = NULL;
      }

    freemem(gc->buf_qy);
    gc->buf_qy = NULL;
  }

  if(gc->buf_w != NULL)
  {
    for(uint i = 0; i < ocl_system.num_devices; ++i)
      if(gc->buf_w != NULL)
      {
        CL_CHECK(clReleaseMemObject(gc->buf_w[i]));
        gc->buf_w[i] = NULL;
      }

    freemem(gc->buf_w);
    gc->buf_w = NULL;
  }

  uninit_gcopencl(gc->gcocl);

  uninit_oclwork(gc->oclwrk);
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
                                          gc->accur,
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
                                          gc->accur,
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
                                          gc->accur,
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
      {//if(H2->u->S.rows == 115) printf("%u: %u\n", i, rname);
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
    {//if(H2->u->S.rows == 115) printf("%u: %u\n", gcocl->num_row_leafs, rname);
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
    int  i = -1;
    uint j;

    for(uint k = 0; k < gcocl->num_row_leafs; ++k)
      if(are_equal_clusters(H2->rb->t, gcocl->row_clusters[gcocl->names_row_leafs[k]]))
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

  gcocl->num_row_leafs = 0;
  gcocl->size_ridx = 0;
  gcocl->size_cidx = 0;
  gcocl->roff = 0;
  gcocl->coff = 0;
  gcocl->names_row_leafs = (uint *) calloc(0, sizeof(uint));
  gcocl->num_h2_leafs_per_row = (uint *) calloc(0, sizeof(uint));
  gcocl->h2_leafs_per_row = (ph2matrix **) calloc(0, sizeof(ph2matrix *));
  gcocl->col_names_per_row = (uint **) calloc(0, sizeof(uint *));
  gcocl->ridx_off = (uint *) calloc(0, sizeof(uint));
  gcocl->host_ridx = (uint *) calloc(0, sizeof(uint));
  gcocl->cidx_off = (uint **) calloc(0, sizeof(uint *));
  gcocl->host_cidx = (uint *) calloc(0, sizeof(uint));

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
    for(uint j = 0; j < gcocl->num_h2_leafs_per_row[i]; ++j)
      gcocl->workload_per_row[i] += gcocl->h2_leafs_per_row[i][j]->u->S.rows *
                                    gcocl->h2_leafs_per_row[i][j]->u->S.cols;
    gcocl->idx_off[i] =
      (i == 0 ? 0 : gcocl->idx_off[i - 1] + gcocl->num_h2_leafs_per_row[i - 1]);
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

  for (uint i = 0; i < gcocl->num_row_leafs; ++i)
  {
    for(uint j = 0; j < gcocl->num_h2_leafs_per_row[i]; ++j)
      gcocl->cidx_sizes[gcocl->idx_off[i] + j] =
        gcocl->h2_leafs_per_row[i][j]->cb->k;

    memcpy(gcocl->host_cidx_off + gcocl->idx_off[i],
           gcocl->cidx_off[i],
           gcocl->num_h2_leafs_per_row[i] * sizeof(uint));
  }

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

    CL_CHECK(clEnqueueWriteBuffer
               (ocl_system.queues[i * ocl_system.queues_per_device], // Currently using only one of everything)
                gcocl->buf_ridx[i],
                false,
                0,
                gcocl->size_ridx * sizeof(uint),
                gcocl->host_ridx,
                0,
                NULL,
                NULL));

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

    CL_CHECK(clEnqueueWriteBuffer
               (ocl_system.queues[0], // Currently using only one of everything)
                gcocl->buf_cidx[i * ocl_system.queues_per_device],
                false,
                0,
                gcocl->size_cidx * sizeof(uint),
                gcocl->host_cidx,
                0,
                NULL,
                NULL));
  }

  gcocl->row_clusters       = enumerate_cluster((pcluster) H2->rb->t);

  gcocl->cur_h2matrix_leafs =
    (uint *) calloc(gcocl->num_row_leafs, sizeof(uint));

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

  for(uint i = 0; i < gcocl->num_row_leafs; ++i)
  {
    memcpy(gcocl->host_xtoffs + gcocl->idx_off[i],
           gcocl->xtoffs[i],
           gcocl->num_h2_leafs_per_row[i] * sizeof(uint));
  }

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

  freemem(gcocl->cur_h2matrix_leafs);
  gcocl->cur_h2matrix_leafs = NULL;

  freemem(gcocl->row_clusters);
  gcocl->row_clusters = NULL;
}

static void
distribute_opencl_work(pcgcopencl gcocl, poclwork oclwrk)
{
  oclwrk->num_wrk_pkgs = ocl_system.num_devices;

  oclwrk->wrk_per_pkg      =
    (uint *) calloc(oclwrk->num_wrk_pkgs, sizeof(uint));

  oclwrk->num_rows_per_pkg =
    (uint *) calloc(oclwrk->num_wrk_pkgs, sizeof(uint));

  oclwrk->rows_per_pkg     =
    (uint **) calloc(oclwrk->num_wrk_pkgs, sizeof(uint*));

  for(uint i = 0; i < oclwrk->num_wrk_pkgs; ++i)
    oclwrk->rows_per_pkg[i] = (uint *) calloc(0, sizeof(uint));

  for(uint i = 0; i < gcocl->num_row_leafs; ++i)
  {
    uint min_wrk     = oclwrk->wrk_per_pkg[0];
    uint min_wrk_pkg = 0;

    for(uint j = 1; j < oclwrk->num_wrk_pkgs; ++j)
      if(oclwrk->wrk_per_pkg[j] < min_wrk)
      {
        min_wrk     = oclwrk->wrk_per_pkg[j];
        min_wrk_pkg = j;
      }

    oclwrk->wrk_per_pkg[min_wrk_pkg]      += gcocl->workload_per_row[i];

    oclwrk->num_rows_per_pkg[min_wrk_pkg] += 1;

    oclwrk->rows_per_pkg[min_wrk_pkg] =
      realloc(oclwrk->rows_per_pkg[min_wrk_pkg],
              oclwrk->num_rows_per_pkg[min_wrk_pkg] * sizeof(uint));

    oclwrk->rows_per_pkg[min_wrk_pkg][oclwrk->num_rows_per_pkg[min_wrk_pkg] - 1]
      = i;
  }
}

ph2matrix
build_green_cross_h2matrix_greencross(pgreencross gc, void *eta)
{
  pblock    broot;
  ph2matrix H2;

  broot = build_strict_block(gc->rc, gc->cc, eta, admissible_max_cluster);
//  view_block(broot);
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

  get_ocl_informations_gcocl(H2, gc);

  for(uint i = 0; i < ocl_system.num_devices; ++i)
  {
    for(uint j = 0; i < num_kernels; ++i)
    {
      for(uint k = 0; j < ocl_system.queues_per_device; ++j)
      {
        cl_kernel kernel = gc->gcocl->kernels[k +
                                              i * ocl_system.queues_per_device +
                                              j * ocl_system.num_devices *
                                              ocl_system.queues_per_device];

        CL_CHECK(clSetKernelArg(kernel, 10, sizeof(cl_mem), &gc->gcocl->buf_num_h2_leafs_per_row[i]));
        CL_CHECK(clSetKernelArg(kernel, 11, sizeof(cl_mem), &gc->gcocl->buf_idx_off[i]));
        CL_CHECK(clSetKernelArg(kernel, 12, sizeof(cl_mem), &gc->gcocl->buf_ridx_sizes[i]));
        CL_CHECK(clSetKernelArg(kernel, 13, sizeof(cl_mem), &gc->gcocl->buf_cidx_sizes[i]));
        CL_CHECK(clSetKernelArg(kernel, 14, sizeof(cl_mem), &gc->gcocl->buf_ridx_off[i]));
        CL_CHECK(clSetKernelArg(kernel, 15, sizeof(cl_mem), &gc->gcocl->buf_cidx_off[i]));
        CL_CHECK(clSetKernelArg(kernel, 16, sizeof(cl_mem), &gc->gcocl->buf_ridx[i]));
        CL_CHECK(clSetKernelArg(kernel, 17, sizeof(cl_mem), &gc->gcocl->buf_cidx[i]));

        CL_CHECK(clSetKernelArg(kernel, 19, sizeof(cl_mem), &gc->gcocl->buf_xtoffs[i]));
        CL_CHECK(clSetKernelArg(kernel, 20, sizeof(cl_mem), &gc->gcocl->buf_ytoffs[i]));
        CL_CHECK(clSetKernelArg(kernel, 21, sizeof(cl_mem), &gc->gcocl->buf_xt[i]));
        CL_CHECK(clSetKernelArg(kernel, 22, sizeof(cl_mem), &gc->gcocl->buf_yt[i]));
      }
    }
  }

  distribute_opencl_work(gc->gcocl, gc->oclwrk);

//  for(uint i = 0; i < gc->oclwrk->num_wrk_pkgs; ++i)
//  {
//    printf("Group %u (Work load: %u):\n", i, gc->oclwrk->wrk_per_pkg[i]);
//
//    for(uint j = 0; j < gc->oclwrk->num_rows_per_pkg[i]; ++j)
//      printf("  %u\n", gc->oclwrk->rows_per_pkg[i][j]);
//  }

  del_block(broot);

  return H2;
}

void
fastaddeval_nearfield_h2matrix_avector_greencross(field      alpha,
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

        fastaddeval_nearfield_h2matrix_avector_greencross(alpha, h2->son[i + j * rsons], xt1, yt1);

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
fastaddeval_farfield_cpu_h2matrix_avector_greencross(pgreencross gc,
                                                     field      alpha,
                                                     pch2matrix h2,
                                                     pavector   xt,
                                                     pavector   yt)
{
  avector   loc1, loc2;
  pavector  xt1, yt1;
  pcclusterbasis rb = h2->rb;
  pcclusterbasis cb = h2->cb;
  uint      rsons = h2->rsons;
  uint      csons = h2->csons;
  uint      xtoff, ytoff;
  uint      i, j;

  if(h2->u)
  {
//    if(are_equal_clusters(gc->gcocl->h2_leafs_per_row[0][0]->rb->t, h2->rb->t))
//    {
      addeval_amatrix_avector(alpha, &h2->u->S, xt, yt);
//    }
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

        fastaddeval_farfield_cpu_h2matrix_avector_greencross(gc, alpha, h2->son[i + j * rsons], xt1, yt1);

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
fastaddeval_farfield_h2matrix_avector_greencross(pgreencross  gc,
                                                 field        alpha,
                                                 pavector     xt,
                                                 pavector     yt)
{
  const size_t num_global_work_items = OCL_LOCAL_WORK_GROUP_SIZE *
                                       gc->gcocl->num_row_leafs;

  pgcopencl gcocl  = gc->gcocl;
  poclwork  oclwrk = gc->oclwrk;

  cl_event  evnt;

  cl_kernel kernel = gcocl->kernels[0];

  CL_CHECK(clEnqueueWriteBuffer(ocl_system.queues[0],
                                gcocl->buf_xt[0],
                                false,
                                0,
                                xt->dim * sizeof(real),
                                xt->v,
                                0,
                                NULL,
                                &evnt));

  CL_CHECK(clWaitForEvents(1, &gc->event_p[0]));
  CL_CHECK(clReleaseEvent(gc->event_p[0]));
  gc->event_p[0] = NULL;

  CL_CHECK(clWaitForEvents(1, &evnt));
  CL_CHECK(clReleaseEvent(evnt));

  CL_CHECK(clSetKernelArg(kernel, 9, sizeof(uint), &oclwrk->num_rows_per_pkg[0]));
  CL_CHECK(clSetKernelArg(kernel, 18, sizeof(real), &alpha));


  CL_CHECK(clEnqueueNDRangeKernel(ocl_system.queues[0],
                                  kernel,
                                  1,
                                  NULL,
                                  &num_global_work_items,
                                  &OCL_LOCAL_WORK_GROUP_SIZE,
                                  0,
                                  NULL,
                                  &evnt));

  CL_CHECK(clWaitForEvents(1, &evnt));
  CL_CHECK(clReleaseEvent(evnt));

  CL_CHECK(clEnqueueReadBuffer(ocl_system.queues[0],
                               gcocl->buf_yt[0],
                               true,
                               0,
                               yt->dim * sizeof(real),
                               yt->v,
                               0,
                               NULL,
                               NULL));

//  for(uint i = 0; i < oclwrk->num_rows_per_pkg[0]; ++i)
//  {
//    CL_CHECK(clEnqueueReadBuffer(ocl_system.queues[0],
//                                 gcocl->buf_yt[0],
//                                 true,
//                                 gcocl->ytoffs[oclwrk->rows_per_pkg[0][i]] *
//                                 sizeof(real),
//                                 gcocl->ridx_sizes[oclwrk->rows_per_pkg[0][i]] *
//                                 sizeof(real),
//                                 yt->v +
//                                 gcocl->ytoffs[oclwrk->rows_per_pkg[0][i]],
//                                 0,
//                                 NULL,
//                                 NULL));
//  }
}

void
fastaddeval_h2matrix_avector_greencross(pgreencross gc,
                                        field       alpha,
                                        pch2matrix  H2,
                                        pavector    xt,
			                                  pavector    yt)
{
//  fastaddeval_farfield_cpu_h2matrix_avector_greencross(gc, alpha, H2, xt, yt);
//
//  print_avector(yt);
//  clear_avector(yt);

  fastaddeval_farfield_h2matrix_avector_greencross(gc, alpha, xt, yt);

//  print_avector(yt);

  fastaddeval_nearfield_h2matrix_avector_greencross(alpha, H2, xt, yt);
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
