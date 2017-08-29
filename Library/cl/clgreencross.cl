/* ------------------------------------------------------------
 * This is the file "clgreencross.cl" of this master thesis.
 * All rights reserved, Bennet Carstensen 2017
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
laplace3d(private float3 x, private float3 y)
{
  private float3 xmy2 = x - y;
  xmy2 *= xmy2;


  const float norm2 = xmy2.x + xmy2.y + xmy2.z;

  float norm;

  norm = native_rsqrt(norm2);

  return r_four_pi * norm;
}

// kernel void
// fastaddeval_h2matrix_avector_1(       const uint dim,
//                                       const uint n,
//                                global const real *vs,
//                                global const uint *p,
//                                global const real *g,
//                                       const uint nq,
//                                constant     real *qx,
//                                constant     real *qy,
//                                constant     real *w,
//                                       const uint num_row_leafs,
//                                global const uint *rows_this_device,
//                                global const uint *num_h2_leafs_per_cluster,
//                                global const uint *idx_offs,
//                                global const uint *ridx_sizes,
//                                global const uint *cidx_sizes,
//                                global const uint *ridx_offs,
//                                global const uint *cidx_offs,
//                                global const uint *ridxs,
//                                global const uint *cidxs,
//                                       const real alpha,
//                                global const uint *xt_offs,
//                                global const uint *yt_offs,
//                                global const real *xts,
//                                global       real *yt)
// {
//   const size_t grpid0 = get_group_id(0);
//
//   if(grpid0 >= num_row_leafs)
//     return;
//   else
//   {
//     const uint   row_this_group = rows_this_device[grpid0];
//
//     const uint   num_h2_leafs   = num_h2_leafs_per_cluster[row_this_group];
//     const size_t lid0           = get_local_id(0);
//
//     gcidxinfo gcii;
//     geom      sur;
//     quad      quad;
//
//     local   real  yt_tmp[32], gr_tmp[32];
//     local   real3 x_tmp[32][3];
//     private real  gc0, gc1;
//     private real  *gc;
//     private real3 (*y)[3];
//     private real3 x[3], y0[3], y1[3], xx, yy;
//
//     uint    k = 0;
//
//     init_geom(&sur, dim, n, vs, p, g, 0, 0);
//     init_quad(&quad, nq, qx, qy, w);
//
//     /* Same information for all threads in a work group */
//     init_row_gcidxinfo(&gcii,
//                        num_h2_leafs,
//                        idx_offs[row_this_group],
//                        ridx_sizes[row_this_group],
//                        ridx_offs[row_this_group],
//                        ridxs,
//                        yt_offs[row_this_group],
//                        yt,
//                        0);
//
//     /* Different information for all threads in a work group */
//     if(lid0 < num_h2_leafs)
//       set_column_info_gcidxinfo(&gcii,
//                                 cidx_sizes[gcii.idx_off + lid0],
//                                 cidx_offs[gcii.idx_off + lid0],
//                                 cidxs,
//                                 xt_offs[gcii.idx_off + lid0],
//                                 xts);
//
//     //if(32 < gcii.ridx_size : 32 : gcii.ridx_size)
//     for(k = 0; k < gcii.ridx_size; k += 32)
//     {
//       const uint loops = (k + 32) < gcii.ridx_size ? 32 : gcii.ridx_size - k;
//
//       barrier(CLK_LOCAL_MEM_FENCE);
//
//       if(lid0 < loops)
//       {
//         const uint  ii = gcii.ridx[k + lid0];
//
//         gr_tmp[lid0]   = sur.g[ii];
//
//         const uint3 p  = vload3(ii, sur.p);
//
//         x_tmp[lid0][0] = vload3(p.x, sur.v);
//         x_tmp[lid0][1] = vload3(p.y, sur.v);
//         x_tmp[lid0][2] = vload3(p.z, sur.v);
//
//         yt_tmp[lid0]   = r_zero;
//       }
//
//       if(lid0 < num_h2_leafs)
//       {
//         for(uint i = 0; i < loops; ++i)
//         {
//           barrier(CLK_LOCAL_MEM_FENCE);
//
//           const uint ii = (lid0 + i) % loops;
//           const real gr = gr_tmp[ii];
//
//           x[0] = x_tmp[ii][0];
//           x[1] = x_tmp[ii][1];
//           x[2] = x_tmp[ii][2];
//
//           real res = r_zero;
//
//           uint  jj = gcii.cidx[0];
//           uint3 p  = vload3(jj, sur.p);
//
//           gc0 = sur.g[jj];
//
//           y0[0] = vload3(p.x, sur.v);
//           y0[1] = vload3(p.y, sur.v);
//           y0[2] = vload3(p.z, sur.v);
//
//           uchar prefetch_switch = 1;
//
//           for(uint j = 0; j < gcii.cidx_size; ++j)
//           {
//             if(j < (gcii.cidx_size - 1))
//             {
//               gc  = prefetch_switch == 0 ? &gc0 : &gc1;
//               y   = prefetch_switch == 0 ? &y0 : &y1;
//
//               jj  = gcii.cidx[j + 1];
//
//               *gc = sur.g[jj];
//
//               *y[0] = vload3(p.x, sur.v);
//               *y[1] = vload3(p.y, sur.v);
//               *y[2] = vload3(p.z, sur.v);
//
//               prefetch_switch = prefetch_switch == 0 ? 1 : 0;
//             }
//
//             gc  = prefetch_switch == 0 ? &gc0 : &gc1;
//             y   = prefetch_switch == 0 ? &y0 : &y1;
//
//             real sum = r_zero;
//
//             for(uint q = 0; q < quad.nq; ++q)
//             {
//               real t = quad.qx[q];
//               real s = quad.qx[q + quad.nq];
//
//               xx = (r_one - t) * x[0] + (t - s) * x[1] + s * x[2];
//
//               t = quad.qy[q];
//               s = quad.qy[q + quad.nq];
//
//               yy = (r_one - t) * *y[0] + (t - s) * *y[1] + s * *y[2];
//
//               sum += quad.w[q] * laplace3d((float3) xx, (float3) yy);
//             }
//
//             res += gcii.xt[j] * *gc * gr * sum;
//           }
//
//           for(uint j = 0; j < num_h2_leafs; j += loops)
//           {
//             const uint max = (j + loops) < num_h2_leafs
//                                ? j + loops
//                                : num_h2_leafs;
//
//
//             barrier(CLK_LOCAL_MEM_FENCE);
//
//             if((lid0 >= j) && (lid0 < max))
//               yt_tmp[ii] += res;
//           }
//         }
//       }
//
//       barrier(CLK_LOCAL_MEM_FENCE);
//
//       if(lid0 < loops)
//         gcii.yt[k + lid0] = alpha * yt_tmp[lid0];
//     }
//   }
// }

kernel void
fastaddeval_h2matrix_avector_0(       const uint dim,
                                      const uint n,
                               global const real *vs,
                               global const uint *p,
                               global const real *g,
                                      const uint nq,
                               constant     real *qx,
                               constant     real *qy,
                               constant     real *w,
                                      const uint num_row_leafs,
                               global const uint *rows_this_device,
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
    const uint   row_this_group = rows_this_device[grpid0];

    const uint   num_h2_leafs   = num_h2_leafs_per_cluster[row_this_group];
    const size_t lid0           = get_local_id(0);


    gcidxinfo gcii;
    geom      sur;
    quad      quad;

    local   real  yt_tmp[32], gr_tmp[32];
    local   real3 x_tmp[32][3];
    private float3 xx, yy;
    private real3 x[3], y[3];

    init_geom(&sur, dim, n, vs, p, g, 0, 0);
    init_quad(&quad, nq, qx, qy, w);

    /* Same information for all threads in a work group */
    init_row_gcidxinfo(&gcii,
                       num_h2_leafs,
                       idx_offs[row_this_group],
                       ridx_sizes[row_this_group],
                       ridx_offs[row_this_group],
                       ridxs,
                       yt_offs[row_this_group],
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

    for(uint k = 0; k < gcii.ridx_size; k += 32)
    {
      const uint loops = (k + 32) < gcii.ridx_size ? 32 : gcii.ridx_size - k;

      barrier(CLK_LOCAL_MEM_FENCE);

      if(lid0 < loops)
      {
        const uint  ii = gcii.ridx[k + lid0];

        gr_tmp[lid0]   = sur.g[ii];

        const uint3 p  = vload3(ii, sur.p);

        x_tmp[lid0][0] = vload3(p.x, sur.v);
        x_tmp[lid0][1] = vload3(p.y, sur.v);
        x_tmp[lid0][2] = vload3(p.z, sur.v);

        yt_tmp[lid0]   = r_zero;
      }

      if(lid0 < num_h2_leafs)
      {
        for(uint i = 0; i < loops; ++i)
        {
          barrier(CLK_LOCAL_MEM_FENCE);

          const uint ii = (lid0 + i) % loops;
          const real gr = gr_tmp[ii];

          x[0] = x_tmp[ii][0];
          x[1] = x_tmp[ii][1];
          x[2] = x_tmp[ii][2];

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
              float s, t;

#ifndef USE_FLOAT
              t  = convert_float(quad.qx[q]);
              s  = convert_float(quad.qx[q + quad.nq]);

              xx = (1.0f - t) * convert_float3(x[0]) +
                   (t - s)    * convert_float3(x[1]) +
                   s          * convert_float3(x[2]);

              t  = convert_float(quad.qy[q]);
              s  = convert_float(quad.qy[q + quad.nq]);

              yy = (1.0f - t) * convert_float3(y[0]) +
                   (t - s)    * convert_float3(y[1]) +
                   s          * convert_float3(y[2]);
#else
              t  = quad.qx[q];
              s  = quad.qx[q + quad.nq]

              xx = (r_one - t) * x[0] + (t - s) * x[1] + s * x[2];

              t = quad.qy[q];
              s = quad.qy[q + quad.nq];

              yy = (r_one - t) * y[0] + (t - s) * y[1] + s * y[2];
#endif

              sum += quad.w[q] * laplace3d(xx, yy);
            }

            res += gcii.xt[j] * gc * gr * sum;
          }

          for(uint j = 0; j < num_h2_leafs; j += loops)
          {
            const uint max = (j + loops) < num_h2_leafs
                               ? j + loops
                               : num_h2_leafs;


            barrier(CLK_LOCAL_MEM_FENCE);

            if((lid0 >= j) && (lid0 < max))
              yt_tmp[ii] += res;
          }
        }
      }

      barrier(CLK_LOCAL_MEM_FENCE);

      if(lid0 < loops)
        gcii.yt[k + lid0] = alpha * yt_tmp[lid0];
    }
  }
}

#endif
