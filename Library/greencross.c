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

#include <string.h>

/* ------------------------------------------------------------
 * "Header" for static functions
 * ------------------------------------------------------------ */

/** @addtogroup greencross
 *
 * @{*/

INLINE_PREFIX void
fill_green_recursive_hmatrix_greencross(pcgreencross gc, phmatrix h);

INLINE_PREFIX void
green_cross_approximation_recursive_greencross(pcgreencross gc, ph2matrix h2);

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

INLINE_PREFIX void
fill_green_recursive_hmatrix_greencross(pcgreencross gc, phmatrix h)
{
  if(h->son)
  {
    const uint rsons = h->rsons;
    const uint csons = h->csons;

    /* Approximate sons. */
    for(uint j = 0; j < csons; ++j)
      for(uint i = 0; i < rsons; ++i)
        fill_green_recursive_hmatrix_greencross(gc, h->son[i + j * rsons]);
  }
  else if(h->f) // Can't approximate matrix, use full matrix.
    nearfield_greencross(gc,
                         h->rc->size,
                         h->rc->idx,
                         h->cc->size,
                         h->cc->idx,
                         h->f);
  else if(h->r) // Approximate with greens degenerated approximation.
  {
    setrank_rkmatrix(h->r, gc->K);
    fill_green_left_greencross(gc, h->rc, &h->r->A);
    fill_green_right_greencross(gc, h->rc, h->cc, &h->r->B);
  }
}

uint
fill_green_cross_leaf_cluster_greencross(bool         row,
                                         pcgreencross gc,
                                         pccluster    c,
                                         uint         **pidx,
                                         pamatrix     V)
{
  avector   tmpb;

  pamatrix  A, Ai;
  pavector  b;
  prkmatrix RK;

  uint l;

  /* Get the matrix from greens formula we want to approximate. */
  A    = new_amatrix(0, 0);

  fill_green_left_greencross(gc, c, A);

  /* Get a rank factorization and the corresponding row pivot indices of *
   * the previously calculated matrix.                                   */
  RK   = new_rkmatrix(0, 0, 0);

  decomp_fullaca_rkmatrix(A, 0.0, row ? pidx : NULL, row ? NULL : pidx, RK);

  /* Permute the left factor of the rank factorization by the row pivot *
   * indices. This matrix is lower triangular with diagonal elements    *
   * equal to one.                                                      */
  l = RK->k;

  resize_amatrix(A, l, l);

  for(uint j = 0; j < l; ++j)
    for(uint i = 0; i < l; ++i)
      setentry_amatrix(A, i, j, getentry_amatrix(&RK->A, (*pidx)[i], j));

  /* Invert the permuted factor of the rank factorization. */
  Ai = new_identity_amatrix(l, l);

  for(uint i = 0; i < l; ++i)
  {
    b = init_column_avector(&tmpb, Ai, i);

    triangularsolve_amatrix_avector(true, true, false, A, b);

    uninit_avector(b);
  }

  /* Multiply the left factor of the rank factorization with its permuted *
   * inverse and use this as the leaf matrix of the current row           *
   * clusterbasis.                                                        */
  resize_amatrix(V, c->size, l);
  clear_amatrix(V);

  addmul_amatrix(f_one, false, &RK->A, false, Ai, V);

  del_rkmatrix(RK);
  del_amatrix(Ai);
  del_amatrix(A);

  return l;
}

void
green_cross_approximation_recursive_greencross(pcgreencross gc, ph2matrix h2)
{
  if(h2->son)
  {
    const uint rsons = h2->rsons;
    const uint csons = h2->csons;

    for(uint j = 0; j < csons; ++j)
      for(uint i = 0; i < rsons; ++i)
        green_cross_approximation_recursive_greencross(gc,
                                                       h2->son[i + j * rsons]);
  }
  else
  {
    if(h2->f)
      nearfield_greencross(gc,
                           h2->rb->t->size,
                           h2->rb->t->idx,
                           h2->cb->t->size,
                           h2->cb->t->idx,
                           h2->f);
    else if(h2->u)
    {
      pamatrix  A;

      uint csize, lr, lc;
      uint *rpidx, *cpidx;

      lr = fill_green_cross_leaf_cluster_greencross(true,
                                                    gc,
                                                    h2->rb->t,
                                                    &rpidx,
                                                    &h2->rb->V);

      lc = fill_green_cross_leaf_cluster_greencross(false,
                                                    gc,
                                                    h2->cb->t,
                                                    &cpidx,
                                                    &h2->cb->V);

      A = new_amatrix(0, 0);

      nearfield_greencross(gc,
                           h2->rb->t->size,
                           h2->rb->t->idx,
                           h2->cb->t->size,
                           h2->cb->t->idx,
                           A);

      resize_amatrix(&h2->u->S, lr, lc);
      nearfield_greencross(gc, lr, rpidx, lc, cpidx, &h2->u->S);

      freemem(cpidx);
      freemem(rpidx);

      del_amatrix(A);
    }
  }
}

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

    for(uint d = 0; d < dim; ++d)
    {
      /* Center of gravity as characteristic point */
      cg->x[i][d] = (x[e[i][0]][d] + x[e[i][1]][d]) * 0.5;

      /* Lower front left corner of bounding box */
      cg->smin[i][d] = REAL_MIN(x[e[i][0]][d], x[e[i][1]][d]);

      /* Upper back right corner of bounding box */
      cg->smax[i][d] = REAL_MAX(x[e[i][0]][d], x[e[i][1]][d]);
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
  gc->rc                  = NULL;
  gc->cc                  = NULL;
  gc->zg                  = NULL;
  gc->wg                  = NULL;
  gc->kernel_function     = NULL;
  gc->pdx_kernel_function = NULL;
  gc->pdy_kernel_function = NULL;
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

  gc->kernel_function     = NULL;
  gc->pdx_kernel_function = NULL;
  gc->pdy_kernel_function = NULL;
}

pgreencross
new_laplace2d_greencross(pcurve2d gr, uint res, uint q, uint m)
{
  pclustergeometry cg;
  pgreencross      gc;

  /* Allocate and init object. */

  gc = (pgreencross) allocmem(sizeof(greencross));

  init_greencross(gc, 2);

  /* Set and get cluster from geometry. */

  gc->geom = (void *) gr;

  gc->n  = gr->edges;

  cg     = build_clustergeometry_greencross((const void *) gr, 2, &gc->idx);

  gc->rc = gc->cc = build_adaptive_cluster(cg, gc->n, gc->idx, res);

  /* Get regular quadrature points. */

  gc->nq = q;
  gc->zq = allocreal(q);
  gc->wq = allocreal(q);
  assemble_gauss(q, gc->zq, gc->wq);

  /* Get quadrature for singular integrals. */

  /* TODO: 3d and DLP case implementation */
  gc->sq = (void *) build_log_singquad1d(q, gc->zq, gc->wq); //

  /* Get green quadrature points. */

  gc->ng = m;
  gc->zg = allocreal(m);
  gc->wg = allocreal(m);
  assemble_gauss(m, gc->zg, gc->wg);

  /* Local rank of greens degenerated approximation. */

  gc->K                   = 4 * m;

  /* Kernel function and partial derivatives. */

  gc->kernel_function     = laplace_kernel;
  gc->pdx_kernel_function = pdx_laplace_kernel;
  gc->pdy_kernel_function = pdy_laplace_kernel;

  /* Clean up */

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
    const uint jj = cidx != NULL ? cidx[j] : j;

    for(uint i = 0; i < rsize; ++i)
    {
      const uint ii = ridx != NULL ? ridx[i] : i;
      const uint c  = select_quadrature_singquad1d(sq1d,
                                                   e[ii],
                                                   e[jj],
                                                   pe0,
                                                   pe1,
                                                   &e0q,
                                                   &e1q,
                                                   &wq,
                                                   &nq,
                                                   &base);

      real sum;

      switch(c)
      {
        case 0: case 2:
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
              e0[d] = a0 * x[pe0[0]][d] + a1 * x[pe0[1]][d];

            a1 = e1q[q];
            a0 = r_one - a1;

            for(uint d = 0; d < dim; ++d)
              e1[d] = a0 * x[pe1[0]][d] + a1 * x[pe1[1]][d];

            sum += wq[q] * gc->kernel_function(e0, e1, dim);
          }
          setentry_amatrix(G, i, j, g[ii] * g[jj] * sum);
          break;
        case 1:
          // TODO: 3D-Case
          if(e[ii][0] == e[jj][0])
          {
            pe0[0] = e[ii][1];
            pe0[1] = e[ii][0];
            pe1[0] = e[jj][1];
          }
          else if(e[ii][0] == e[jj][1])
          {
            pe0[0] = e[ii][1];
            pe0[1] = e[ii][0];
            pe1[0] = e[jj][0];
          }
          else if(e[ii][1] == e[jj][0])
          {
            pe0[0] = e[ii][0];
            pe0[1] = e[ii][1];
            pe1[0] = e[jj][1];
          }
          else if(e[ii][1] == e[jj][1])
          {
            pe0[0] = e[ii][0];
            pe0[1] = e[ii][1];
            pe1[0] = e[jj][0];
          }
          else
          {
            printf("ERROR!\n");
            exit(0);
          }

          sum = 0.0;

          for(uint q = 0; q < nq; ++q) {

            real dx, dy, tx, ty;

            tx = e0q[q];
            ty = e1q[q];

            dx = x[pe0[0]][0] * (-tx) + x[pe0[1]][0] * (tx + ty) + x[pe1[0]][0] * (-ty);
            dy = x[pe0[0]][1] * (-tx) + x[pe0[1]][1] * (tx + ty) + x[pe1[0]][1] * (-ty);

            sum += wq[q] * REAL_LOG(dx * dx + dy * dy);
          }
          setentry_amatrix(G, i, j, (2.0 * base + sum) * r_minus_two_pi * g[ii] * g[jj] * 0.5);
          break;

        default:
          break;
      }
    }
  }
}

void
fill_green_left_greencross(pcgreencross gc, pccluster t, pamatrix A)
{
  const uint dim  = gc->dim; // Dimension of the problem
  const uint dim2 = 2 * dim;
  const uint K    = gc->K;   // Local rank of greens low-rank approximation.
  const uint m    = gc->ng;  // Approximation order in 1st dimension.
  const uint m1   = dim == 2 ? 1 : m; // Approximation order in 2nd dimension
  const uint nq   = gc->nq;  // Number of regular quadrature points/weigts.
  const uint size = t->size; // Number of base functions.

  pcreal wg = gc->wg; // Green quadrature weights.
  pcreal wq = gc->wq; // Regular quadrature weights.
  pcreal zg = gc->zg; // Green quadrature points in [-1, 1].
  pcreal zq = gc->zq; // Regular quadrature points in [-1, 1].

  amatrix tmp1, tmp2;

  pamatrix A1, A2;
  pcurve2d c2d;

  real delta;

  /* Affine mapping factors */
  field bpa[dim], bma[dim];
  field xx[dim];   // Transformed regular quadrature points
  field z[dim];    // Transformed green quadrature points

  uint M[dim - 1]; // Index set for green quadrature points.

  real *g;

  real (*x)[dim];
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

  /* Get mapping factors for boundary domain. */
  delta = r_zero;

  for(uint d = 0; d < dim; ++d)
  {
    const real diam = t->bmax[d] - t->bmin[d];

    delta  = delta < diam ? diam : delta;

    bpa[d] = (t->bmax[d] + t->bmin[d]) * 0.5;
    bma[d] = diam;
  }

  delta *= 0.5;

  for(uint d = 0; d < dim; ++d)
    bma[d] = bma[d] * 0.5 + delta;

  resize_amatrix(A, size, 2 * K);


  A1 = init_sub_amatrix(&tmp1, A, size, 0, K, 0);
  A2 = init_sub_amatrix(&tmp2, A, size, 0, K, K);

  for(uint iota = 0; iota < dim2; ++iota)
  {
    /* Quadrature point we induce for greens representation formula. */
    const real z_iota = iota & 1 ? r_one : r_minusone;
    /* Position in which we induce the quadrature point. */
    const uint l      = iota >> 1;

    for(uint mu0 = 0; mu0 < m; ++mu0)
    {
      M[0] = mu0;

      for(uint mu1 = 0; mu1 < m1; ++mu1)
      {
        real det, w;

        M[1] = mu1;

        det  = w = r_one;

        /* Get green quadrature points and weight */
        for(uint d = 0; d < dim; ++d)
        {
          z[d] = bpa[d] +
                 bma[d] * (d == l ? z_iota : zg[M[(d > l ? d - 1 : d)]]);

          w   *= d == l ? r_one : wg[M[(d > l ? d - 1 : d)]];

          det *= d == l ? r_one : bma[d] * bma[d];
        }

        w *= REAL_SQRT(det);
        w  = REAL_SQRT(w);

        for(uint i = 0; i < size; ++i)
        {
          const uint ii = t->idx[i];

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

          setentry_amatrix(A1,
                           i,
                           iota + dim2 * (mu0 + m * mu1),
                           w * g[ii] * int_kernel);
          /*
           * Since our faces are axially parallel, the outer normal vectors are
           * the (negative) l-th unit vectors and thus we can reduce the
           * directionl derivative to the l-th partial derivative * z_iota.
           */
          setentry_amatrix(A2,
                           i,
                           iota + dim2 * (mu0 + m * mu1),
                           w * g[ii] * z_iota * int_pdx_kernel);
        }
      }
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
  const uint dim  = gc->dim; // Dimension of the problem
  const uint dim2 = 2 * dim;
  const uint K    = gc->K;   // Local rank of greens low-rank approximation.
  const uint m    = gc->ng;  // Approximation order in 1st dimension.
  const uint m1   = dim == 2 ? 1 : gc->ng; // Approximation order in 2nd dimension
  const uint nq   = gc->nq;  // Number of regular quadrature points/weigts.
  const uint size = s->size; // Number of base functions.

  pcreal wg = gc->wg; // Green quadrature weights.
  pcreal wq = gc->wq; // Regular quadrature weights.
  pcreal zg = gc->zg; // Green quadrature points in [-1, 1].
  pcreal zq = gc->zq; // Regular quadrature points in [-1, 1].

  amatrix tmp1, tmp2;

  pamatrix B1, B2;
  pcurve2d c2d;

  real delta;

  /* Affine mapping factors. */
  field bpa[dim], bma[dim];
  field yy[dim]; // Transformed regular quadrature points.
  field z[dim];  // Transformed green quadrature points.

  uint M[dim - 1]; // Index set for green quadrature points.

  real *g;

  real (*x)[dim];
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

  /* Get mapping factors for boundary domain. */
  delta = r_zero;

  for(uint d = 0; d < dim; ++d)
  {
    const real diam = REAL_ABS(t->bmax[d] - t->bmin[d]);

    delta  = delta < diam ? diam : delta;

    bpa[d] = (t->bmax[d] + t->bmin[d]) * 0.5;
    bma[d] = diam;
  }

  delta *= 0.5;

  for(uint d = 0; d < dim; ++d)
    bma[d] = (bma[d] + 2.0 * delta) * 0.5;

  resize_amatrix(B, size, 2 * K);

  B1 = init_sub_amatrix(&tmp1, B, size, 0, K, 0);
  B2 = init_sub_amatrix(&tmp2, B, size, 0, K, K);

  for(uint iota = 0; iota < dim2; ++iota)
  {
    /* Quadrature point we induce for greens representation formula. */
    const real z_iota = iota & 1 ? r_one : r_minusone;
    /* Position in which we induce the quadrature point. */
    const uint l      = iota >> 1;

    for(uint mu0 = 0; mu0 < m; ++mu0)
    {
      M[0] = mu0;

      for(uint mu1 = 0; mu1 < m1; ++mu1)
      {
        real det, w;

        M[1] = mu1;

        det  = w = r_one;

        /* Get green quadrature points and weight */
        for(uint d = 0; d < dim; ++d)
        {
          z[d] = bpa[d] +
                 bma[d] * (d == l ? z_iota : zg[M[(d > l ? d - 1 : d)]]);

          w   *= d == l ? r_one : wg[M[(d > l ? d - 1 : d)]];

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

            int_pdy_kernel += wq[q] * gc->pdy_kernel_function(z, yy, dim, l);
            int_kernel     += wq[q] * gc->kernel_function(z, yy, dim);
          }

          /*
           * Since our faces are axially parallel, the outer normal vectors are
           * the (negative) l-th unit vectors and thus we can reduce the
           * directionl derivative to the l-th partial derivative * z_iota.
           */
          setentry_amatrix(B1,
                           i,
                           iota + dim2 * (mu0 + m * mu1),
                           w * g[ii] * z_iota * int_pdy_kernel);
          setentry_amatrix(B2,
                           i,
                           iota + dim2 * (mu0 + m * mu1),
                           - w * g[ii] * int_kernel);
        }
      }
    }
  }

  uninit_amatrix(B1);
  uninit_amatrix(B2);
}

void
assemble_green_rkmatrix_greencross(pcgreencross gc,
                                   pccluster    row,
                                   pccluster    col,
                                   prkmatrix    AB)
{
  setrank_rkmatrix(AB, 2 * gc->K);
  fill_green_left_greencross(gc, row, &AB->A);
  fill_green_right_greencross(gc, row, col, &AB->B);
}

prkmatrix
build_green_rkmatrix_greencross(pcgreencross gc, pccluster row, pccluster col)
{
  prkmatrix AB;

  AB = new_rkmatrix(row->size, col->size, 2 * gc->K);

  assemble_green_rkmatrix_greencross(gc, row, col, AB);

  return AB;
}

phmatrix
fill_green_hmatrix_greencross(pcgreencross gc, void *eta)
{
  pblock   broot;
  phmatrix h;

  broot = build_strict_block(gc->rc, gc->cc, eta, admissible_2_cluster);

  h     = build_from_block_hmatrix(broot, gc->K);

  fill_green_recursive_hmatrix_greencross(gc, h);

  return h;
}

ph2matrix
green_cross_approximation(pcgreencross gc, void *eta)
{
  pblock        broot;
  pclusterbasis rb, cb;
  ph2matrix     h2;

  broot = build_strict_block(gc->rc, gc->cc, eta, admissible_2_cluster);

  rb    = build_from_cluster_clusterbasis(gc->rc);
  cb    = build_from_cluster_clusterbasis(gc->cc);

  h2    = build_from_block_h2matrix(broot, rb, cb);

  green_cross_approximation_recursive_greencross(gc, h2);

  return h2;
}
