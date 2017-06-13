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

#include "kernels.h"

#include "aca.h"
#include "gaussquad.h"
#include "singquad1d.h"
#include "singquad2d.h"

/* ------------------------------------------------------------
 * "Header" for static functions
 * ------------------------------------------------------------ */

/** @addtogroup greencross
 *
 * @{*/

/** @brief Retrieve the maximum level of a cluster basis.
 *
 *  @attention @p level should be 0 since it's only used insid the function to
 *             keep track of the level of the current viewed cluster basis.
 *
 *  @param cb    The cluster basis of which we search the maximum level.
 *  @param level The level of the current viewed cluster basis. */
INLINE_PREFIX uint
get_max_lvl_clusterbasis(pclusterbasis cb, const uint level);

/** @brief Retrieve the number of cluster basis for each level from a given root
 *         cluster basis.
 *
 *  After execution, each entry of @p on_level with index i will hold the
 *  number of cluster bases whose level are i.
 *
 *  @attention @p level should be 0 since it's only used insid the function to
 *             keep track of the level of the current viewed cluster basis.
 *  @attention @p on_level should be an zero initialized array of length
 *             @ref get_max_lvl_clusterbasis
 *             "get_max_lvl_clusterbasis(cb, level)" + 1.
 *
 *  @param cb       The root cluster basis.
 *  @param level    The level of the current viewed cluster basis.
 *  @param on_level An array whose entries will hold the number of cluster bases
 *                  for each level. */
INLINE_PREFIX void
set_on_lvl_info_clusterbasis(pclusterbasis cb, const uint level, uint *on_level);

/** @brief Retrieve the cluster bases of a given root cluster basis divided into
 *         the level they are on.
 *
 *  After execution, each entry of @p cb_per_level will hold the pointers to all
 *  cluster bases whose level is i, where i is the index of the entry inside
 *  cb_per_level.
 *
 *  @attention @p level should be 0 since it's only used insid the function to
 *             keep track of the level of the current viewed cluster basis.
 *  @attention @p lvl_idx should be an zero initialized array of length
 *             @ref get_max_lvl_clusterbasis
 *             "get_max_lvl_clusterbasis(cb, level)" + 1.
 *  @attention @p cb_per_level should be an initialized array of length
 *             get_max_lvl_clusterbasis(cb, level) + 1, where each entry with
 *             index should be an initialized array of length @p on_level[i],
 *             where @p on_level is retrieved by
 *             @ref set_on_lvl_info_clusterbasis
 *             "set_on_lvl_info_clusterbasis(cb, 0, on_level)".
 *
 *  @param cb           The root cluster basis.
 *  @param level        The level of the current viewed cluster basis.
 *  @param lvl_idx      The number of cluster bases already found on each level
 *                      respectively.
 *  @param cb_per_level An arrays whose entries will hold pointers to all
 *                      cluster bases of the corresponding level. */
INLINE_PREFIX void
set_per_level_clusterbasis(pclusterbasis cb,
                           const uint level,
                           uint *lvl_idx,
                           pclusterbasis **cb_per_level);

/** @brief Creats a <a href="http://www.h2lib.org/doc/d8/da5/struct__clustergeometry.html">
 *         clustergeometry</a> object from a triangular mesh.
 *
 *  The geometry is taken from @p data interpreted as different objects
 *  depending on the dimension of the problem @p dim.
 *
 *  @attention For @p dim = 2, @p data is expected to be an object of type
 *             <a href="http://www.h2lib.org/doc/d8/df9/struct__curve2d.html">
 *             curve2d</a>.
 *
 *  @param[in]  data Pointer to an object containing the informations to
 *                   construct the cluster geometry.
 *  @param[in]  dim  Dimension of the cluster geometry.
 *  @param[out] idx  Index set which will be allocated by this method to index
 *                   each triangle in data. The entries of idx will be
 *                   0, 1, ..., N-1, where N is equal to the number of
 *                   triangles.
 *  @return          A valid cluster geometry object which can be used to
 *                   construct a cluster tree along with the index set @p idx.*/
INLINE_PREFIX pclustergeometry
build_clustergeometry_greencross(const void *data, const uint dim, uint **idx);

/** @} */

/* ------------------------------------------------------------
 * Source of static functions.
 * ------------------------------------------------------------ */

uint
get_max_lvl_clusterbasis(pclusterbasis cb, const uint level)
{
  const uint sons = cb->sons;

  uint max_level = 0;

  if(sons)
    for(uint i = 0; i < sons; ++i)
    {
      const uint max_level_son = get_max_lvl_clusterbasis(cb->son[i],
                                                          level + 1);
      /* Get the max level of all sons. */
      max_level = max_level < max_level_son ? max_level_son : max_level;
    }
  else
    max_level = level; // Max level on this path found.

  return max_level;
}

void
set_on_lvl_info_clusterbasis(pclusterbasis cb, const uint level, uint *on_level)
{
  const uint sons = cb->sons;

  on_level[level] += 1; // Increment the corresponding level counter.

  /* Iterate over all sons. */
  if(sons)
    for(uint i = 0; i < sons; ++i)
      set_on_lvl_info_clusterbasis(cb->son[i], level + 1, on_level);
}

void
set_per_level_clusterbasis(pclusterbasis cb,
                           const uint level,
                           uint *lvl_idx,
                           pclusterbasis **cb_per_level)
{
  const uint sons = cb->sons;

  cb_per_level[level][lvl_idx[level]] = cb; // Add this cluster basis to this
                                            // level.

  lvl_idx[level] += 1; // Increment the index for the next cluster basis on
                       // this level.

  /* Iterate over all sons. */
  if(sons)
    for(uint i = 0; i < sons; ++i)
      set_per_level_clusterbasis(cb->son[i], level + 1, lvl_idx, cb_per_level);
}

pclustergeometry
build_clustergeometry_greencross(const void *data, const uint dim, uint **idx)
{
  pclustergeometry cg;
  pcurve2d         gr;

  real (*x)[dim];
  uint (*e)[2];
  uint (*t)[3];
  uint n;

  if(dim == 2)
  {
    gr = (pcurve2d) data;
    x  = gr->x;
    e  = gr->e;
    (void) t;
    n  = gr->edges;
  }
  else if(dim == 3)
  {
    // TODO: Implement 3-dimensional version
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

    for(uint j = 0; j < dim; ++j)
    {
      /* Center of gravity as characteristic point */
      cg->x[i][j] = (x[e[i][0]][j] + x[e[i][1]][j]) * 0.5;

      /* Lower front left corner of bounding box */
      cg->smin[i][j] = REAL_MIN(x[e[i][0]][j], x[e[i][1]][j]);

      /* Upper back right corner of bounding box */
      cg->smax[i][j] = REAL_MAX(x[e[i][0]][j], x[e[i][1]][j]);
    }
  }

  return cg;
}

uint uipow(uint base, uint exp)
{
  uint result = 1;
  while (exp)
  {
    if(exp & 1)
      result *= base;
    exp >>= 1;
    base *= base;
  }

  return result;
}

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

void
init_greencross(pgreencross gc, uint dim)
{
  assert(gc != NULL);
  assert((dim >= greencross_min_dim) && (greencross_max_dim <= 3));

  /* Scalar components */
  gc->n                   = 0;
  gc->nq                  = 0;
  gc->ng                  = 0;
  gc->K                   = 0;
  gc->dim                 = dim;

  /* Pointer components */
  gc->idx                 = NULL;
  gc->zg                  = NULL;
  gc->wg                  = NULL;
  gc->h2                  = NULL;
  gc->kernel_function     = NULL;
  gc->pdx_kernel_function = NULL;
  gc->pdy_kernel_function = NULL;
}

void
uninit_greencross(pgreencross gc)
{
  assert(gc != NULL);

  pcluster rc, cc;

  rc = cc = NULL;

  if(gc->idx != NULL)
  {
    freemem(gc->idx);
    gc->idx = NULL;
  }

  if(gc->geom != NULL)
  {
    gc->dim == 2 ? del_curve2d((pcurve2d) gc->geom)
                 : del_surface3d((psurface3d) gc->geom);
    gc->geom = NULL;
  }

  if(gc->zq != NULL)
  {
    freemem(gc->zq);
    gc->zq = NULL;
  }

  if(gc->wq != NULL)
  {
    freemem(gc->wq);
    gc->wq = NULL;
  }

  if(gc->sq != NULL)
  {
    gc->dim == 2 ? del_singquad1d((psingquad1d) gc->sq)
                 : del_singquad2d((psingquad2d) gc->sq);
    gc->sq = NULL;
  }

  if(gc->zg != NULL)
  {
    freemem(gc->zg);
    gc->zg = NULL;
  }

  if(gc->wg != NULL)
  {
    freemem(gc->wg);
    gc->wg = NULL;
  }

  if(gc->h2 != NULL)
  {
    rc = (pcluster) gc->h2->rb->t;
    cc = (pcluster) gc->h2->cb->t;

    del_h2matrix(gc->h2);
    gc->h2 = NULL;
  }

  if(rc == cc)
    cc = NULL; // Only delete cluster once if row cluster = column cluster.

  if(rc != NULL)
  {
    del_cluster(rc);
    rc = NULL;
  }

  if(cc != NULL)
  {
    del_cluster(cc);
    cc = NULL;
  }

  gc->kernel_function     = NULL;
  gc->pdx_kernel_function = NULL;
  gc->pdy_kernel_function = NULL;
}

pgreencross
new_laplace2d_greencross(pcurve2d gr, uint res, void *eta, uint q, uint m)
{
  pblock           broot;
  pcluster         root;
  pclusterbasis    cb, rb;
  pclustergeometry cg;
  pgreencross      gc;

  /* Allocate and init object */

  gc = (pgreencross) allocmem(sizeof(greencross));

  init_greencross(gc, 2);

  gc->geom = (void *) gr;

  /* Set up basis data structures for H2-matrix approximations */

  gc->n  = gr->edges;

  cg     = build_clustergeometry_greencross((const void *) gr, 2, &gc->idx);

  root   = build_adaptive_cluster(cg, gr->edges, gc->idx, res);

  broot  = build_strict_block(root, root, eta, admissible_2_cluster);

  cb     = build_from_cluster_clusterbasis(root);
  rb     = build_from_cluster_clusterbasis(root);

  gc->h2 = build_from_block_h2matrix(broot, cb, rb);

  gc->nq = q;
  gc->zq = allocreal(q);
  gc->wq = allocreal(q);
  assemble_gauss(q, gc->zq, gc->wq);

  /* TODO: 3d and DLP case implementation */
  gc->sq = (void *) build_log_singquad1d(q, gc->zq, gc->wq); //

  gc->ng = m;
  gc->zg = allocreal(m);
  gc->wg = allocreal(m);
  assemble_gauss(m, gc->zg, gc->wg);

  gc->K                   = 4 * m;
  gc->kernel_function     = laplace_kernel;
  gc->pdx_kernel_function = pdx_laplace_kernel;
  gc->pdy_kernel_function = pdy_laplace_kernel;

  del_clustergeometry(cg);

  return gc;
}

void
del_greencross(pgreencross gc)
{
  assert(gc != NULL);

  uninit_greencross(gc);

  freemem(gc);

  gc = NULL;
}

void
nearfield_greencross(pcgreencross gc,
                     const uint   rsize,
                     const uint   *ridx,
                     const uint   csize,
                     const uint   *cidx,
                     pamatrix     G)
{
  const uint dim = gc->dim;

  assert((dim == 2) || (dim == 3));

  pcurve2d    c2d;
  psingquad1d sq1d;

  real base;
  uint nq;

  uint pe0[greencross_max_dim], pe1[greencross_max_dim];
  real e0[greencross_max_dim], e1[greencross_max_dim];

  preal e0q, e1q, g, wq;

  real  (*x)[dim];
  uint  (*e)[2];

  if(dim == 2)
  {
    c2d  = (pcurve2d) gc->geom;
    x    = c2d->x;
    e    = c2d->e;
    g    = c2d->g;

    sq1d = (psingquad1d) gc->sq;
  }
  else
  {
    // TODO: 3d case
  }

  resize_amatrix(G, rsize, csize);

  for(uint j = 0; j < csize; ++j)
  {
    const uint jj = cidx ? cidx[j] : j;

    for(uint i = 0; i < rsize; ++i)
    {
      const uint ii = ridx ? ridx[i] : i;

      real sum;

      if(ii < jj)
        select_quadrature_singquad1d(sq1d,
                                     e[ii],
                                     e[jj],
                                     pe0,
                                     pe1,
                                     &e0q,
                                     &e1q,
                                     &wq,
                                     &nq,
                                     &base);
      else
        select_quadrature_singquad1d(sq1d,
                                     e[jj],
                                     e[ii],
                                     pe1,
                                     pe0,
                                     &e1q,
                                     &e0q,
                                     &wq,
                                     &nq,
                                     &base);

      for(uint d = 0; d < dim; ++d)
      {
        pe0[d] = e[ii][pe0[d]];
        pe1[d] = e[jj][pe1[d]];
      }

      sum = r_minus_two_pi * base;

      for(uint q = 0; q < nq; ++q)
      {
        real a0, a1;

        a1 = e0q[q];
        a0 = r_one - a1;

        for(uint d = 0; d < dim; ++d)
          e0[d] = a0 * x[pe0[0]][d] + a1 * x[pe1[1]][d];

        a1 = e1q[q];
        a0 = r_one - a1;

        for(uint d = 0; d < dim; ++d)
          e1[d] = a0 * x[pe0[0]][d] + a1 * x[pe1[1]][d];

        sum += wq[q] * gc->kernel_function(e0, e1, dim);
      }

      setentry_amatrix(G, i, j, g[ii] * g[jj] * sum);
    }
  }
}

void
fill_green_left_greencross(pcgreencross gc, pccluster c, pamatrix A)
{
  const uint dim  = gc->dim;
  const uint K    = gc->K;
  const uint ng   = gc->ng;
  const uint nq   = gc->nq;
  const uint size = c->size;

  pcreal wg = gc->wg;
  pcreal wq = gc->wq;
  pcreal zg = gc->zg;
  pcreal zq = gc->zq;

  amatrix tmp1, tmp2;

  pamatrix A1, A2;
  pcurve2d c2d;

  real delta;

  field bpa[greencross_max_dim], bma[greencross_max_dim];
  field xx[greencross_max_dim];
  field z[greencross_max_dim];

  uint M[greencross_max_dim - 1];

  real *g;

  real (*x)[greencross_max_dim];
  uint (*e)[2];

  if(dim == 2)
  {
    c2d  = (pcurve2d) gc->geom;
    x    = c2d->x;
    e    = c2d->e;
    g    = c2d->g;
  }
  else
  {
    // TODO: 3d case
  }

  delta = r_zero;

  for(uint d = 0; d < dim; ++d)
  {
    const real diam = REAL_ABS(c->bmax[d] - c->bmin[d]);

    delta  = delta < diam ? diam : delta;

    bpa[d] = (c->bmax[d] + c->bmax[d]) * 0.5;
    bma[d] = diam;
  }

  delta *= 0.5;

  for(uint d = 0; d < dim; ++d)
    bma[d] = (bma[d] + 2.0 * delta) * 0.5;

  resize_amatrix(A, size, 2 * K);

  A1 = init_sub_amatrix(&tmp1, A, size, 0, K, 0);
  A2 = init_sub_amatrix(&tmp2, A, size, 0, K, K);

  for(uint ny = 0; ny < K; ++ny)
  {
    const uint iota   = ny % (2 * dim);
    const real z_iota = iota & 1 ? r_one : r_minusone;

    uint l, tmp;
    real det, w;

    tmp = ny / (2 * dim);

    for(uint d = 0; d < (dim - 1); ++d)
    {
      M[d] = tmp % ng;
      tmp /= ng;
    }

    l   = iota >> 1;

    det = w = r_one;

    for(uint d = 0; d < dim; ++d)
    {
      z[d] = bpa[d] + bma[d] * (d == l ? z_iota : zg[M[d - (d > l ? 1 : 0)]]);

      w   *= d == l ? r_one : wg[M[d - (d > l ? 1 : 0)]];

      det *= d == l ? r_one : bma[d] * bma[d];
    }

    w *= REAL_SQRT(det);
    w  = REAL_SQRT(w);

    for(uint i = 0; i < size; ++i)
    {
      const uint ii = c->idx[i];

      real int_kernel, int_pdx_kernel;

      int_kernel = int_pdx_kernel = r_zero;

      for(uint q = 0; q < nq; ++q)
      {
        const real a1 = zq[q];
        const real a0 = r_one - a1;

        for(uint d = 0; d < dim; ++d)
          // TODO: 3d case
          xx[d] = a0 * x[e[ii][0]][d] + a1 * x[e[ii][1]][d];

        int_kernel     += wq[q] * gc->kernel_function(xx, z, dim);
        int_pdx_kernel += wq[q] * gc->pdx_kernel_function(xx, z, dim, l);
      }

      setentry_amatrix(A1, i, ny, w * g[ii] * int_kernel);
      setentry_amatrix(A2, i, ny, delta * w * g[ii] * z_iota * int_pdx_kernel);
    }
  }

  uninit_amatrix(A1);
  uninit_amatrix(A2);
}

void
fill_green_right_greencross(pcgreencross gc,
                            pccluster    t,
                            pccluster    s,
                            pamatrix     B)
{
  const uint dim  = gc->dim;
  const uint K    = gc->K;
  const uint ng   = gc->ng;
  const uint nq   = gc->nq;
  const uint size = s->size;

  pcreal wg = gc->wg;
  pcreal wq = gc->wq;
  pcreal zg = gc->zg;
  pcreal zq = gc->zq;

  amatrix  tmp1, tmp2;

  pamatrix B1, B2;
  pcurve2d c2d;

  real delta;

  field bpa[greencross_max_dim], bma[greencross_max_dim];
  field yy[greencross_max_dim];
  field z[greencross_max_dim];

  uint M[greencross_max_dim - 1];

  real *g;

  real (*x)[greencross_max_dim];
  uint (*e)[2];

  if(dim == 2)
  {
    c2d  = (pcurve2d) gc->geom;
    x    = c2d->x;
    e    = c2d->e;
    g    = c2d->g;
  }
  else
  {
    // TODO: 3d case
  }

  delta = r_zero;

  for(uint d = 0; d < dim; ++d)
  {
    const real diam = t->bmax[d] - t->bmin[d];

    delta = delta < diam ? diam : delta;

    bpa[d] = (t->bmax[d] + t->bmax[d]) * 0.5;
    bma[d] = diam;
  }

  delta *= 0.5;

  for(uint d = 0; d < dim; ++d)
    bma[d] = (bma[d] + 2.0 * delta) * 0.5;

  resize_amatrix(B, size, 2 * K);

  B1 = init_sub_amatrix(&tmp1, B, size, 0, K, 0);
  B2 = init_sub_amatrix(&tmp2, B, size, 0, K, K);

  for(uint ny = 0; ny < K; ++ny)
  {
    const uint iota   = ny % (2 * dim);
    const real z_iota = iota & 1 ? r_one : r_minusone;

    uint l, tmp;
    real det, w;

    tmp = ny / (2 * dim);

    for(uint d = 0; d < (dim - 1); ++d)
    {
      M[d] = tmp % ng;
      tmp /= ng;
    }

    l   = iota >> 1;

    det = w = r_one;

    for(uint d = 0; d < dim; ++d)
    {
      z[d] = bpa[d] + bma[d] * (d == l ? z_iota : zg[M[d - (d > l ? 1 : 0)]]);

      w   *= d == l ? r_one : wg[M[d - (d > l ? 1 : 0)]];

      det *= d == l ? r_one : bma[d] * bma[d];
    }

    w *= REAL_SQRT(det);
    w  = REAL_SQRT(w);

    for(uint i = 0; i < size; ++i)
    {
      const uint ii = s->idx[i];

      real int_kernel, int_pdy_kernel;

      int_kernel = int_pdy_kernel = r_zero;

      for(uint q = 0; q < nq; ++q)
      {
        const real a1 = zq[q];
        const real a0 = r_one - a1;

        for(uint d = 0; d < dim; ++d)
          // TODO: 3d case
          yy[d] = a0 * x[e[ii][0]][d] + a1 * x[e[ii][1]][d];

        int_kernel     += wq[q] * gc->kernel_function(z, yy, dim);
        int_pdy_kernel += wq[q] * gc->pdy_kernel_function(z, yy, dim, l);
      }

      setentry_amatrix(B1, i, ny, w * g[ii] * z_iota * int_pdy_kernel);
      setentry_amatrix(B2, i, ny, - w * g[ii] * int_kernel / delta);
    }
  }

  uninit_amatrix(B1);
  uninit_amatrix(B2);
}

void
green_cross_approximation(pcgreencross gc, ph2matrix h2)
{
  if(h2->son)
  {
    const uint rsons = h2->rsons;
    const uint csons = h2->csons;

    for(uint j = 0; j < csons; ++j)
      for(uint i = 0; i < rsons; ++i)
        green_cross_approximation(gc, h2->son[i + j * rsons]);
  }
  else
    if(h2->f)
      nearfield_greencross(gc,
                           h2->rb->t->size,
                           h2->rb->t->idx,
                           h2->cb->t->size,
                           h2->cb->t->idx,
                           h2->f);
    else if(h2->u)
    {
      avector  tmpb;

      pamatrix A, Ai;
      pavector b;
      prkmatrix RK;

      uint csize, k;

      uint *rpidx;

      /* Get the matrix from greens formula we want to approximate. */
      A  = new_amatrix(0, 0);

      fill_green_left_greencross(gc, h2->rb->t, A);

      /* Get a rank factorization and the corresponding row pivot indices of *
       * the previously calculated matrix.                                   */
      RK  = new_rkmatrix(0, 0, 0);

      decomp_fullaca_rkmatrix(A, 0.0, &rpidx, NULL, RK);

      /* Permute the left factor of the rank factorization by the row pivot *
       * indices. This matrix is lower triangular with diagonal elements    *
       * equal to one.                                                      */
      k   = RK->k;

      resize_amatrix(A, k, k);

      for(uint j = 0; j < k; ++j)
        for(uint i = 0; i < k; ++i)
          setentry_amatrix(A, i, j, getentry_amatrix(&RK->A, rpidx[i], j));

      /* Invert the permuted factor of the rank factorization. */
      Ai = new_identity_amatrix(k, k);

      for(uint i = 0; i < k; ++i)
      {
        b = init_column_avector(&tmpb, Ai, i);

        triangularsolve_amatrix_avector(true, true, false, A, b);

        uninit_avector(b);
      }

      /* Multiply the left factor of the rank factorization with its permuted *
       * inverse and use this as the leaf matrix of the current row           *
       * clusterbasis.                                                        */
      resize_amatrix(&h2->rb->V, h2->rb->t->size, k);
      clear_amatrix(&h2->rb->V);

      addmul_amatrix(f_one, false, &RK->A, false, Ai, &h2->rb->V);

      nearfield_greencross(gc,
                           h2->rb->t->size,
                           h2->rb->t->idx,
                           h2->cb->t->size,
                           h2->cb->t->idx,
                           A);

      csize = h2->cb->t->size;

      resize_amatrix(&h2->u->S, k, csize);

      for(uint j = 0; j < csize; ++j)
        for(uint i = 0; i < k; ++i)
          setentry_amatrix(&h2->u->S,
                           i,
                           j,
                           getentry_amatrix(A, rpidx[i], j));

      freemem(rpidx);

      del_rkmatrix(RK);
      del_amatrix(Ai);
      del_amatrix(A);
    }
}
