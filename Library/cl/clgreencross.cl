/* ------------------------------------------------------------
 * This is the file "clgreencross.cl" of this master thesis.
 * All rights reserved, Steffen Boerm 2009
 * ------------------------------------------------------------ */

/**
 * @file      cl/clgreencross.cl
 * @author    Bennet Carstensen
 * @date      2017
 * @copyright All rights reserved, Bennet Carstensen 2017
 */

#ifndef GREENCROSS_CL
#define GREENCROSS_CL

#include "clbasic.cl"
#include "clgcidxinfo.cl"
#include "clgeom.cl"
#include "clquad.cl"

#ifdef USE_FLOAT
typedef float2  real2;
typedef float4  real4;
typedef float8  real8;
typedef float16 real16;
#else
typedef double2  real2;
typedef double4  real4;
typedef double8  real8;
typedef double16 real16;
#endif

/** @brief @f$\frac{1}{4 \pi}@f$ */
constant static real r_four_pi = 0.0795774715459476;

real
laplace3d(private real *x, private real *y)
{
  real norm;

  const real norm2 = (x[0] - y[0]) * (x[0] - y[0]) +
                     (x[1] - y[1]) * (x[1] - y[1]) +
                     (x[2] - y[2]) * (x[2] - y[2]);

#ifndef USE_FLOAT
  norm = convert_double(native_rsqrt(convert_float(norm2)));
#else
  norm = native_rsqrt(norm2);
#endif

  return r_four_pi * norm;
}

void
mvm_sub_row(pgeom      sur,
            pquad      quad,
            pgcidxinfo gcii,
                    const uint i,
            global  const uint *cidx_sizes,
            global  const uint *cidx_offs,
            global  const uint *cidxs,
            global  const uint *xt_offs,
            global  const real *xts,
            private       real x[3][3],
            private       real y[3][3],
            private       real xx[3],
            private       real yy[3])
{
  const size_t lid0    = get_local_id(0);
  const size_t lsize0  = get_local_size(0);

  real gr;

  if((i + lid0) < gcii->ridx_size)
  {
    const uint ii = gcii->ridx[i + lid0];
    gr = sur->g[ii];

    for(uint d1 = 0; d1 < sur->dim; ++d1)
      for(uint d2 = 0; d2 < sur->dim; ++d2)
        x[d1][d2] = sur->v[d2 + sur->dim * sur->p[sur->dim * ii + d1]];
  }

  real res = r_zero;

  for(uint k = 0; k < gcii->num_h2_leafs; ++k)
  {
    set_column_info_gcidxinfo(gcii,
                              cidx_sizes[gcii->idx_off + k],
                              cidx_offs[gcii->idx_off + k],
                              cidxs,
                              xt_offs[gcii->idx_off + k],
                              xts);

    for(uint j = 0; j < gcii->cidx_size; ++j)
    {
      const uint jj = gcii->cidx[j]; // Bank conflict
      const real gc = sur->g[jj];

      for(uint d1 = 0; d1 < sur->dim; ++d1)
        for(uint d2 = 0; d2 < sur->dim; ++d2)
          y[d1][d2] = sur->v[d2 + sur->dim * sur->p[sur->dim * jj + d1]];

      if((i + lid0) < gcii->ridx_size)
      {
        real sum = r_zero;

        for(uint q = 0; q < quad->nq; ++q)
        {
          real t = quad->qx[q];
          real s = quad->qx[q + quad->nq];

          xx[0] = (r_one - t) * x[0][0] + (t - s) * x[1][0] + s * x[2][0];
          xx[1] = (r_one - t) * x[0][1] + (t - s) * x[1][1] + s * x[2][1];
          xx[2] = (r_one - t) * x[0][2] + (t - s) * x[1][2] + s * x[2][2];

          t = quad->qy[q];
          s = quad->qy[q + quad->nq];

          yy[0] = (r_one - t) * y[0][0] + (t - s) * y[1][0] + s * y[2][0];
          yy[1] = (r_one - t) * y[0][1] + (t - s) * y[1][1] + s * y[2][1];
          yy[2] = (r_one - t) * y[0][2] + (t - s) * y[1][2] + s * y[2][2];

          sum += quad->w[q] * laplace3d(xx, yy);
        }

        res += gcii->xt[j] * gc * gr * sum;
      }
    }
  }

  if((i + lid0) < gcii->ridx_size)
    gcii->yt[i + lid0] += res;
}

kernel void
fastaddeval_h2matrix_avector(       const uint dim,
                                    const uint n,
                             global const real *vs,
                             global const uint *p,
                             global const real *g,
                                    const uint nq,
                             constant     real *qx,
                             constant     real *qy,
                             constant     real *w,
                                    const uint num_row_leafs,
                             global const uint *num_h2_leafs_per_cluster,
                             global const uint *idx_offs,
                             global const uint *ridx_sizes,
                             global const uint *cidx_sizes,
                             global const uint *ridx_offs,
                             global const uint *cidx_offs,
                             global const uint *ridxs,
                             global const uint *cidxs,
                                    const real alpha,
                             global const uint *xt_offs,
                             global const uint *yt_offs,
                             global const real *xts,
                             global       real *yt)
{
  const size_t grpid0 = get_group_id(0);

  if(grpid0 >= num_row_leafs)
    return;
  else
  {
    local real xt_tmp[WRK_GRP_SIZE0];
    local uint cidx_tmp[WRK_GRP_SIZE0];

    pgeom sur  = new_geom(dim, n, vs, p, g);
    pquad quad = new_quad(nq, qx, qy, w);
    pgcidxinfo gcii = new_row_gcidxinfo(num_h2_leafs_per_cluster[grpid0],
                                        idx_offs[grpid0],
                                        ridx_sizes[grpid0],
                                        ridx_offs[grpid0],
                                        ridxs,
                                        yt_offs[grpid0],
                                        yt,
                                        cidx_tmp,
                                        xt_tmp);

    const size_t lid0    = get_local_id(0);
    const size_t lsize0  = get_local_size(0);

    private real x[3][3], y[3][3], xx[3], yy[3];

    uint i = 0;

    for(i = 0; (i + lsize0) < gcii->ridx_size; i += lsize0)
      mvm_sub_row(sur,
                  quad,
                  gcii,
                  i,
                  cidx_sizes,
                  cidx_offs,
                  cidxs,
                  xt_offs,
                  xts,
                  x,
                  y,
                  xx,
                  yy);
    // {
    //   const uint ii = gcii->ridx[i + lid0];
    //   const real gr = sur->g[ii];
    //
    //   for(uint d1 = 0; d1 < dim; ++d1)
    //     for(uint d2 = 0; d2 < dim; ++d2)
    //       x[d1][d2] = sur->v[d2 + dim * p[dim * ii + d1]];
    //
    //   real res = r_zero;
    //
    //   for(uint k = 0; k < gcii->num_h2_leafs; ++k)
    //   {
    //     set_column_info_gcidxinfo(gcii,
    //                               cidx_sizes[gcii->idx_off + k],
    //                               cidx_offs[gcii->idx_off + k],
    //                               cidxs,
    //                               xt_offs[gcii->idx_off + k],
    //                               xt);
    //
    //     // const uint cidx_size = cidx_sizes[gcii->idx_off + k];
    //     // const uint cidx_off  = cidx_offs[gcii->idx_off + k];
    //     // const uint xt_off    = xt_offs[gcii->idx_off + k];
    //
    //     for(uint j = 0; j < gcii->cidx_size; ++j)
    //     {
    //       const uint jj = gcii->cidx[j]; // Bank conflict
    //       const real gc = sur->g[jj];
    //
    //       for(uint d1 = 0; d1 < dim; ++d1)
    //         for(uint d2 = 0; d2 < dim; ++d2)
    //           y[d1][d2] = sur->v[d2 + dim * p[dim * jj + d1]];
    //
    //       real sum = r_zero;
    //
    //       for(uint q = 0; q < quad->nq; ++q)
    //       {
    //         real t = quad->qx[q];
    //         real s = quad->qx[q + nq];
    //
    //         xx[0] = (r_one - t) * x[0][0] + (t - s) * x[1][0] + s * x[2][0];
    //         xx[1] = (r_one - t) * x[0][1] + (t - s) * x[1][1] + s * x[2][1];
    //         xx[2] = (r_one - t) * x[0][2] + (t - s) * x[1][2] + s * x[2][2];
    //
    //         t = quad->qy[q];
    //         s = quad->qy[q + nq];
    //
    //         yy[0] = (r_one - t) * y[0][0] + (t - s) * y[1][0] + s * y[2][0];
    //         yy[1] = (r_one - t) * y[0][1] + (t - s) * y[1][1] + s * y[2][1];
    //         yy[2] = (r_one - t) * y[0][2] + (t - s) * y[1][2] + s * y[2][2];
    //
    //         sum += quad->w[q] * laplace3d(xx, yy);
    //       }
    //
    //       // res += xt[xt_off + j] * gc * gr * sum;
    //       res += gcii->xt[j] * gc * gr * sum;
    //     }
    //   }
    //
    //   gcii->yt[i + lid0] += res;
    // }

    mvm_sub_row(sur,
                quad,
                gcii,
                i,
                cidx_sizes,
                cidx_offs,
                cidxs,
                xt_offs,
                xts,
                x,
                y,
                xx,
                yy);
  //   {
  //     const uint ii = gcii->ridx[i + lid0];
  //     const real gr = g[ii];
  //
  //     for(uint d1 = 0; d1 < dim; ++d1)
  //       for(uint d2 = 0; d2 < dim; ++d2)
  //         x[d1][d2] = vs[d2 + dim * p[dim * ii + d1]];
  //
  //     real res = r_zero;
  //
  //     for(uint k = 0; k < gcii->num_h2_leafs; ++k)
  //     {
  //       set_column_info_gcidxinfo(gcii,
  //                                 cidx_sizes[gcii->idx_off + k],
  //                                 cidx_offs[gcii->idx_off + k],
  //                                 cidxs,
  //                                 xt_offs[gcii->idx_off + k],
  //                                 xts);
  //
  //       // const uint cidx_size = cidx_sizes[gcii->idx_off + k];
  //       // const uint cidx_off  = cidx_offs[gcii->idx_off + k];
  //       // const uint xt_off    = xt_offs[gcii->idx_off + k];
  //
  //       for(uint j = 0; j < gcii->cidx_size; ++j)
  //       {
  //         const uint jj = gcii->cidx[j]; // Bank conflict
  //         const real gc = g[jj];
  //
  //         for(uint d1 = 0; d1 < dim; ++d1)
  //           for(uint d2 = 0; d2 < dim; ++d2)
  //             y[d1][d2] = vs[d2 + dim * p[dim * jj + d1]];
  //
  //         real sum = r_zero;
  //
  //         for(uint q = 0; q < nq; ++q)
  //         {
  //           real t = qx[q];
  //           real s = qx[q + nq];
  //
  //           xx[0] = (r_one - t) * x[0][0] + (t - s) * x[1][0] + s * x[2][0];
  //           xx[1] = (r_one - t) * x[0][1] + (t - s) * x[1][1] + s * x[2][1];
  //           xx[2] = (r_one - t) * x[0][2] + (t - s) * x[1][2] + s * x[2][2];
  //
  //           t = qy[q];
  //           s = qy[q + nq];
  //
  //           yy[0] = (r_one - t) * y[0][0] + (t - s) * y[1][0] + s * y[2][0];
  //           yy[1] = (r_one - t) * y[0][1] + (t - s) * y[1][1] + s * y[2][1];
  //           yy[2] = (r_one - t) * y[0][2] + (t - s) * y[1][2] + s * y[2][2];
  //
  //           sum += w[q] * laplace3d(xx, yy);
  //         }
  //
  //         res += gcii->xt[j] * gc * gr * sum;
  //       }
  //     }
  //
  //     gcii->yt[i + lid0] += res;
  //   }
  }
}

#endif
