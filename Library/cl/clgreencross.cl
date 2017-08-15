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
typedef float4  real4;
typedef float8  real8;
typedef float16 real16;
#else
typedef double4  real4;
typedef double8  real8;
typedef double16 real16;
#endif

/** @brief @f$\frac{1}{4 \pi}@f$ */
constant static real r_four_pi = 0.0795774715459476;

real
laplace3d(private real3 x, private real3 y)
{
  private real3 xmy2 = x - y;
  xmy2 *= xmy2;


  const real norm2 = xmy2.x + xmy2.y + xmy2.z;

  real norm;

#ifndef USE_FLOAT
  norm = convert_double(native_rsqrt(convert_float(norm2)));
#else
  norm = native_rsqrt(norm2);
#endif

  return r_four_pi * norm;
}

void
mvm_sub_row(        pgeom      sur,
                    pquad      quad,
                    pgcidxinfo  gcii,
                    const uint  i,
            global  const uint  *cidx_sizes,
            global  const uint  *cidx_offs,
            global  const uint  *cidxs,
            global  const uint  *xt_offs,
            global  const real  *xts,
            private       real3 x[3],
            private       real3 y[3],
            private       real3 xx,
            private       real3 yy)
{
  const size_t lid0    = get_local_id(0);

  local real3 (*v_tmp)[3] = (local real3 (*)[3]) sur->v_tmp;

  real gr;

  if((i + lid0) < gcii->ridx_size)
  {
    const uint ii = gcii->ridx[i + lid0];
    gr = sur->g[ii];

    private uint3 p = vload3(ii, sur->p);

    x[0] = vload3(p.x, sur->v);
    x[1] = vload3(p.y, sur->v);
    x[2] = vload3(p.z, sur->v);
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


    for(uint l = 0; l < gcii->cidx_size; l += WRK_GRP_SIZE0)
    {
      const uint loops = ((l + WRK_GRP_SIZE0) < gcii->cidx_size)
                           ? WRK_GRP_SIZE0 : (gcii->cidx_size - WRK_GRP_SIZE0);

      barrier(CLK_LOCAL_MEM_FENCE);

      if(lid0 < loops)
      {
        const uint jj = gcii->cidx[l + lid0];
        const uint3 p = vload3(jj, sur->p);

        //local real3 (*v_tmp)[3] = (local real3 (*)[3]) sur->v_tmp;

        sur->g_tmp[lid0]   = sur->g[jj];

        v_tmp[lid0][0] = vload3(p.x, sur->v);
        v_tmp[lid0][1] = vload3(p.y, sur->v);
        v_tmp[lid0][2] = vload3(p.z, sur->v);

        gcii->xt_tmp[lid0] = gcii->xt[l + lid0];

      }

      barrier(CLK_LOCAL_MEM_FENCE);

      if((i + lid0) < gcii->ridx_size)
      {
        for(uint j = 0; j < loops; ++j)
        {
          const real gc = sur->g_tmp[j];

          y[0] = vload3(0, (local real*) &v_tmp[j][0]);
          y[1] = vload3(0, (local real*) &v_tmp[j][1]);
          y[2] = vload3(0, (local real*) &v_tmp[j][2]);
          // private uint3 p = vload3(jj, sur->p);
          //
          // y[0] = vload3(p.x, sur->v);
          // y[1] = vload3(p.y, sur->v);
          // y[2] = vload3(p.z, sur->v);

          real sum = r_zero;

          for(uint q = 0; q < quad->nq; ++q)
          {
            real t = quad->qx[q];
            real s = quad->qx[q + quad->nq];

            xx = (r_one - t) * x[0] + (t - s) * x[1] + s * x[2];

            t = quad->qy[q];
            s = quad->qy[q + quad->nq];

            yy = (r_one - t) * y[0] + (t - s) * y[1] + s * y[2];

            sum += quad->w[q] * laplace3d(xx, yy);
          }

          res += gcii->xt_tmp[j] * gc * gr * sum;
        }
      }
    }
  }

    // const uint loops = gcii->cidx_size / WRK_GRP_SIZE0;
    //
    // uint l = 0;
    //
    // for(l = 0; l < loops; ++l)
    // {
    //   gcii->xt_tmp[lid0] = gcii->xt[l * WRK_GRP_SIZE0 + lid0];
    //
    //   if((i + lid0) < gcii->ridx_size)
    //   {
    //     for(uint j = 0; j < WRK_GRP_SIZE0; ++j)
    //     {
    //       const uint jj = gcii->cidx[l * WRK_GRP_SIZE0 + j]; // Bank conflict
    //       const real gc = sur->g[jj];
    //
    //       private uint3 p = vload3(jj, sur->p);
    //
    //       y[0] = vload3(p.x, sur->v);
    //       y[1] = vload3(p.y, sur->v);
    //       y[2] = vload3(p.z, sur->v);
    //
    //       real sum = r_zero;
    //
    //       for(uint q = 0; q < quad->nq; ++q)
    //       {
    //         real t = quad->qx[q];
    //         real s = quad->qx[q + quad->nq];
    //
    //         xx = (r_one - t) * x[0] + (t - s) * x[1] + s * x[2];
    //
    //         t = quad->qy[q];
    //         s = quad->qy[q + quad->nq];
    //
    //         yy = (r_one - t) * y[0] + (t - s) * y[1] + s * y[2];
    //
    //         sum += quad->w[q] * laplace3d(xx, yy);
    //       }
    //
    //       res += gcii->xt[l * WRK_GRP_SIZE0 + j] * gc * gr * sum;
    //     }
    //   }
    // }
    //
    // if(lid0 < (gcii->cidx_size - l * WRK_GRP_SIZE0))
    // {
    //   gcii->xt_tmp[lid0] = gcii->xt[l * WRK_GRP_SIZE0 + lid0];
    // }
    //
    // barrier(CLK_LOCAL_MEM_FENCE);
    //
    // if((i + lid0) < gcii->ridx_size)
    // {
    //   for(uint j = 0; j < (gcii->cidx_size - l * WRK_GRP_SIZE0); ++j)
    //   {
    //     const uint jj = gcii->cidx[l * WRK_GRP_SIZE0 + j]; // Bank conflict
    //     const real gc = sur->g[jj];
    //
    //     private uint3 p = vload3(jj, sur->p);
    //
    //     y[0] = vload3(p.x, sur->v);
    //     y[1] = vload3(p.y, sur->v);
    //     y[2] = vload3(p.z, sur->v);
    //
    //     real sum = r_zero;
    //
    //     for(uint q = 0; q < quad->nq; ++q)
    //     {
    //       real t = quad->qx[q];
    //       real s = quad->qx[q + quad->nq];
    //
    //       xx = (r_one - t) * x[0] + (t - s) * x[1] + s * x[2];
    //
    //       t = quad->qy[q];
    //       s = quad->qy[q + quad->nq];
    //
    //       yy = (r_one - t) * y[0] + (t - s) * y[1] + s * y[2];
    //
    //       sum += quad->w[q] * laplace3d(xx, yy);
    //     }
    //
    //     res += gcii->xt[l * WRK_GRP_SIZE0 + j] * gc * gr * sum;
    //   }
    // }

  if((i + lid0) < gcii->ridx_size)
    gcii->yt[i + lid0] += res;
}

kernel void
fastaddeval_h2matrix_avector_1(       const uint dim,
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
    local real g_tmp[WRK_GRP_SIZE0], xt_tmp[WRK_GRP_SIZE0];
    local real3 v_tmp[WRK_GRP_SIZE0][3];

    gcidxinfo gcii;
    geom      sur;
    quad      q;

    init_geom(&sur, dim, n, vs, p, g, (local void *) v_tmp, g_tmp);
    init_quad(&q, nq, qx, qy, w);
    init_row_gcidxinfo(&gcii,
                       num_h2_leafs_per_cluster[grpid0],
                       idx_offs[grpid0],
                       ridx_sizes[grpid0],
                       ridx_offs[grpid0],
                       ridxs,
                       yt_offs[grpid0],
                       yt,
                       xt_tmp);

    private real3 x[3], y[3], xx, yy;

    uint i = 0;

    for(i = 0; i < gcii.ridx_size; i += WRK_GRP_SIZE0)
      mvm_sub_row(&sur,
                  &q,
                  &gcii,
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
  }
}

void
clPrefixSum(local real* x, const size_t n)
{
  const size_t lid0 = get_local_id(0);

  for(uint i = 2; i <= n; i <<= 1)
  {
    if((lid0 % i) == (i - 1))
      x[lid0] += x[lid0 - (i >> 1)];

    barrier(CLK_LOCAL_MEM_FENCE);
  }

  for(uint i = (n >> 1); i >= 2; i >>= 1)
  {
    if((lid0 >= i) && ((lid0 % i) == ((i >> 1) - 1)))
      x[lid0] += x[lid0 - (i >> 1)];

    barrier(CLK_LOCAL_MEM_FENCE);
  }
}

kernel void
fastaddeval_h2matrix_avector_2(       const uint dim,
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
    const uint   num_h2_leafs = num_h2_leafs_per_cluster[grpid0];
    const size_t lid0         = get_local_id(0);

    gcidxinfo gcii;
    geom      sur;
    quad      quad;

    local   real  sum[WRK_GRP_SIZE0], yt_tmp[WRK_GRP_SIZE0];
    private real3 x[3], y[3], xx, yy;

    event_t event = 0;

    init_geom(&sur, dim, n, vs, p, g, 0, 0);
    init_quad(&quad, nq, qx, qy, w);

    /* Same information for all threads in a work group */
    init_row_gcidxinfo(&gcii,
                       num_h2_leafs,
                       idx_offs[grpid0],
                       ridx_sizes[grpid0],
                       ridx_offs[grpid0],
                       ridxs,
                       yt_offs[grpid0],
                       yt,
                       0);

    /* Different information for all threads in a work group */
    if(lid0 < num_h2_leafs)
      set_column_info_gcidxinfo(&gcii,
                                cidx_sizes[gcii.idx_off + lid0],
                                cidx_offs[gcii.idx_off + lid0],
                                cidxs,
                                xt_offs[gcii.idx_off + lid0],
                                xts);

    for(uint k = 0; k < gcii.ridx_size; k += WRK_GRP_SIZE0)
    {
      const uint loops = (k + WRK_GRP_SIZE0) < gcii.ridx_size
                           ? WRK_GRP_SIZE0
                           : gcii.ridx_size - k;

      if(lid0 < num_h2_leafs)
      {
        for(uint i = 0; i < loops; ++i)
        {
          const uint ii = gcii.ridx[k + i];
          const real gr = sur.g[ii];

          private uint3 p = vload3(ii, sur.p);

          x[0] = vload3(p.x, sur.v);
          x[1] = vload3(p.y, sur.v);
          x[2] = vload3(p.z, sur.v);

          real res = r_zero;

          for(uint j = 0; j < gcii.cidx_size; ++j)
          {
            const uint jj = gcii.cidx[j];
            const real gc = sur.g[jj];

            private uint3 p = vload3(jj, sur.p);

            y[0] = vload3(p.x, sur.v);
            y[1] = vload3(p.y, sur.v);
            y[2] = vload3(p.z, sur.v);

            real sum = r_zero;

            for(uint q = 0; q < quad.nq; ++q)
            {
              real t = quad.qx[q];
              real s = quad.qx[q + quad.nq];

              xx = (r_one - t) * x[0] + (t - s) * x[1] + s * x[2];

              t = quad.qy[q];
              s = quad.qy[q + quad.nq];

              yy = (r_one - t) * y[0] + (t - s) * y[1] + s * y[2];

              sum += quad.w[q] * laplace3d(xx, yy);
            }

            res += gcii.xt[j] * gc * gr * sum;
          }

          /* Parallel (sum) reduction to get the entry of yt */
          sum[lid0] = res;

          uint num_leafs = gcii.num_h2_leafs;

          while(num_leafs > 1)
          {
            const uint corrector = num_leafs % 2;

            num_leafs >>= 1;

            barrier(CLK_LOCAL_MEM_FENCE);

            if((lid0 >= corrector) && (lid0 < (num_leafs + corrector)))
              sum[lid0] += sum[lid0 + num_leafs];

            num_leafs += corrector;
          }

          barrier(CLK_LOCAL_MEM_FENCE);

          if(lid0 == 0)
            yt_tmp[i] = alpha * sum[0];
          //  gcii.yt[k + i] += alpha * sum[0];
        }
      }

      wait_group_events(1, &event);
      event = async_work_group_copy(gcii.yt + k,
                                    yt_tmp,
                                    loops,
                                    0);
    }

    wait_group_events(1, &event);
    // if(lid0 < num_h2_leafs)
    // {
    //   for(uint i = 0; i < gcii.ridx_size; ++i)
    //   {
    //     const uint ii = gcii.ridx[i];
    //     const real gr = sur.g[ii];
    //
    //     private uint3 p = vload3(ii, sur.p);
    //
    //     x[0] = vload3(p.x, sur.v);
    //     x[1] = vload3(p.y, sur.v);
    //     x[2] = vload3(p.z, sur.v);
    //
    //     real res = r_zero;
    //
    //     for(uint j = 0; j < gcii.cidx_size; ++j)
    //     {
    //       const uint jj = gcii.cidx[j];
    //       const real gc = sur.g[jj];
    //
    //       private uint3 p = vload3(jj, sur.p);
    //
    //       y[0] = vload3(p.x, sur.v);
    //       y[1] = vload3(p.y, sur.v);
    //       y[2] = vload3(p.z, sur.v);
    //
    //       real sum = r_zero;
    //
    //       for(uint q = 0; q < quad.nq; ++q)
    //       {
    //         real t = quad.qx[q];
    //         real s = quad.qx[q + quad.nq];
    //
    //         xx = (r_one - t) * x[0] + (t - s) * x[1] + s * x[2];
    //
    //         t = quad.qy[q];
    //         s = quad.qy[q + quad.nq];
    //
    //         yy = (r_one - t) * y[0] + (t - s) * y[1] + s * y[2];
    //
    //         sum += quad.w[q] * laplace3d(xx, yy);
    //       }
    //
    //       res += gcii.xt[j] * gc * gr * sum;
    //     }
    //
    //     tmp[lid0] = res;
    //
    //     uint num_leafs = gcii.num_h2_leafs;
    //
    //     while(num_leafs > 1)
    //     {
    //       const uint corrector = num_leafs % 2;
    //
    //       num_leafs >>= 1;
    //
    //       barrier(CLK_LOCAL_MEM_FENCE);
    //
    //       if((lid0 >= corrector) && (lid0 < (num_leafs + corrector)))
    //         tmp[lid0] += tmp[lid0 + num_leafs];
    //
    //       num_leafs += corrector;
    //     }
    //
    //     barrier(CLK_LOCAL_MEM_FENCE);
    //
    //     if(lid0 == 0)
    //       gcii.yt[i] += alpha * tmp[0];
    //   }
    // }
  }
}

kernel void
fastaddeval_h2matrix_avector_3(       const uint dim,
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
    const uint   num_h2_leafs = num_h2_leafs_per_cluster[grpid0];
    const size_t lid0         = get_local_id(0);

    gcidxinfo gcii;
    geom      sur;
    quad      quad;

    local   real  sum[WRK_GRP_SIZE0], yt_tmp[WRK_GRP_SIZE0];
    private real3 x[3], y[3], xx, yy;

    event_t event = 0;

    init_geom(&sur, dim, n, vs, p, g, 0, 0);
    init_quad(&quad, nq, qx, qy, w);

    /* Same information for all threads in a work group */
    init_row_gcidxinfo(&gcii,
                       num_h2_leafs,
                       idx_offs[grpid0],
                       ridx_sizes[grpid0],
                       ridx_offs[grpid0],
                       ridxs,
                       yt_offs[grpid0],
                       yt,
                       0);

    /* Different information for all threads in a work group */
    if(lid0 < num_h2_leafs)
      set_column_info_gcidxinfo(&gcii,
                                cidx_sizes[gcii.idx_off + lid0],
                                cidx_offs[gcii.idx_off + lid0],
                                cidxs,
                                xt_offs[gcii.idx_off + lid0],
                                xts);

    for(uint k = 0; k < gcii.ridx_size; k += WRK_GRP_SIZE0)
    {
      const uint loops = (k + WRK_GRP_SIZE0) < gcii.ridx_size
                           ? WRK_GRP_SIZE0
                           : gcii.ridx_size - k;

      if(lid0 < num_h2_leafs)
      {
        for(uint i = 0; i < loops; ++i)
        {
          const uint ii = gcii.ridx[k + i];
          const real gr = sur.g[ii];

          private uint3 p = vload3(ii, sur.p);

          x[0] = vload3(p.x, sur.v);
          x[1] = vload3(p.y, sur.v);
          x[2] = vload3(p.z, sur.v);

          real res = r_zero;

          for(uint j = 0; j < gcii.cidx_size; ++j)
          {
            if(j != (gcii.cidx_size - 1))
            {
              const uint jj = gcii.cidx[j + 1];

              prefetch(sur.g + jj, 1);

              private uint3 p = vload3(jj, sur.p);

              prefetch(sur.v + sur.dim * p.x, 1);
              prefetch(sur.v + sur.dim * p.y, 1);
              prefetch(sur.v + sur.dim * p.z, 1);
            }

            const uint jj = gcii.cidx[j];
            const real gc = sur.g[jj];

            private uint3 p = vload3(jj, sur.p);

            y[0] = vload3(p.x, sur.v);
            y[1] = vload3(p.y, sur.v);
            y[2] = vload3(p.z, sur.v);

            real sum = r_zero;

            for(uint q = 0; q < quad.nq; ++q)
            {
              real t = quad.qx[q];
              real s = quad.qx[q + quad.nq];

              xx = (r_one - t) * x[0] + (t - s) * x[1] + s * x[2];

              t = quad.qy[q];
              s = quad.qy[q + quad.nq];

              yy = (r_one - t) * y[0] + (t - s) * y[1] + s * y[2];

              sum += quad.w[q] * laplace3d(xx, yy);
            }

            res += gcii.xt[j] * gc * gr * sum;
          }

          /* Parallel (sum) reduction to get the entry of yt */
          sum[lid0] = res;

          uint num_leafs = gcii.num_h2_leafs;

          while(num_leafs > 1)
          {
            const uint corrector = num_leafs % 2;

            num_leafs >>= 1;

            barrier(CLK_LOCAL_MEM_FENCE);

            if((lid0 >= corrector) && (lid0 < (num_leafs + corrector)))
              sum[lid0] += sum[lid0 + num_leafs];

            num_leafs += corrector;
          }

          barrier(CLK_LOCAL_MEM_FENCE);

          if(lid0 == 0)
            yt_tmp[i] = alpha * sum[0];
          //  gcii.yt[k + i] += alpha * sum[0];
        }
      }

      wait_group_events(1, &event);
      event = async_work_group_copy(gcii.yt + k,
                                    yt_tmp,
                                    loops,
                                    0);
    }

    wait_group_events(1, &event);
  }
}

#endif
