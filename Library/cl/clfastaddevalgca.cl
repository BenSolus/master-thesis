/* ------------------------------------------------------------
 * This is the file "clfastaddeval.cl" of this master thesis.
 * All rights reserved, Bennet Carstensen 2017
 * ------------------------------------------------------------ */

/**
 * @file      cl/clfastaddeval.cl
 * @author    Bennet Carstensen
 * @date      2017
 * @copyright All rights reserved, Bennet Carstensen 2017
 */

#include "clgcidxinfo.cl"
#include "clgeom.cl"
#include "clkernels.cl"
#include "clsingquad.cl"

#ifndef FASTADDEVAL_CL
#define FASTADDEVAL_CL

void
fastaddeval_farfield(pgcidxinfo gcii, pgeom sur, psingquadl sq,
                     private const real bem_alpha,
                     private const real kernel_factor,
                     private const real alpha)
{
  const size_t lid0 = get_local_id(0);

  for(uint i = 0; i < gcii->ridx_size; ++i)
  {
    const uint ii = gcii->ridx[i];

    const real  gr = sur->g[ii];

    const uint3 q  = vload3(ii, sur->p);

    private real3 x[3];

    x[0] = vload3(q.x, sur->v);
    x[1] = vload3(q.y, sur->v);
    x[2] = vload3(q.z, sur->v);

    // const real yt_i = gcii->yt[i];

    real res = r_zero;

    for(uint j = 0; j < gcii->cidx_size; ++j)
    {
      const uint  jj = gcii->cidx[j];

      const real  gc = sur->g[jj];
      const uint3 p  = vload3(jj, sur->p);

      const real factor = (ii == jj) * 0.5 * bem_alpha * gr;

      private real3 y[3];

      y[0] = vload3(p.x, sur->v);
      y[1] = vload3(p.y, sur->v);
      y[2] = vload3(p.z, sur->v);

      real sum = sq->bases.x;

      uint q = 0;

      for(; (q + 4) <= sq->nq; q += 4)
      {
        const float4 w = vload4(q / 4, sq->wqs);

        float4 a1 = vload4(q / 4, sq->xqs);
        float4 a2 = vload4(q / 4, sq->xqs + sq->nq);

        float4 xx[3] = { 0, 0, 0 };

        xx[0] = (r_one - a1) * x[0].x + (a1 - a2) * x[1].x + a2 * x[2].x;
        xx[1] = (r_one - a1) * x[0].y + (a1 - a2) * x[1].y + a2 * x[2].y;
        xx[2] = (r_one - a1) * x[0].z + (a1 - a2) * x[1].z + a2 * x[2].z;

        float4 yy[3] = { 0, 0, 0 };

        yy[0] = (r_one - a1) * y[0].x + (a1 - a2) * y[1].x + a2 * y[2].x;
        yy[1] = (r_one - a1) * y[0].y + (a1 - a2) * y[1].y + a2 * y[2].y;
        yy[2] = (r_one - a1) * y[0].z + (a1 - a2) * y[1].z + a2 * y[2].z;

        const float4 quad = w * laplace3d4(xx, yy);

        sum += quad.x + quad.y + quad.z + quad.w;
      }

      for(; q < sq->nq; ++q)
      {
        float s, t;

#ifndef USE_FLOAT
        t  = convert_float(sq->xqs[q]);
        s  = convert_float(sq->xqs[q + sq->nq]);

        private float3 xx, yy;

        xx = (1.0f - t) * convert_float3(x[0]) +
             (t - s)    * convert_float3(x[1]) +
             s          * convert_float3(x[2]);

        t  = convert_float(sq->yqs[q]);
        s  = convert_float(sq->yqs[q + sq->nq]);

        yy = (1.0f - t) * convert_float3(y[0]) +
             (t - s)    * convert_float3(y[1]) +
             s          * convert_float3(y[2]);
#else
        t  = sq->xqs[q];
        s  = sq->xqs[q + sq->nq];

        // if(get_global_id(0) == 0 && i == 0 && j == 0)
        //   printf("%2.u: %.5e %.5e ", q, t, s);

        private float3 xx, yy;

        xx = (r_one - t) * x[0] + (t - s) * x[1] + s * x[2];

        t  = sq->yqs[q];
        s  = sq->yqs[q + sq->nq];

        yy = (r_one - t) * y[0] + (t - s) * y[1] + s * y[2];
#endif

        sum += sq->wqs[q] * laplace3d(xx, yy);
      }

      res += gcii->xt[j] * ((gc * gr * kernel_factor * sum) + factor);
    }

    for(uint j = 0; j < gcii->num_h2_leafs; ++j)
    {
      barrier(CLK_GLOBAL_MEM_FENCE);

      if(lid0 == j)
        // gcii->xt[j % gcii->cidx_size] *= r_one;
        gcii->yt[i] += alpha * res;
    }
  }
}

void
fastaddeval_farfield_back(pgcidxinfo gcii, pgeom sur, psingquadl sq,
                         private const real bem_alpha,
                         private const real kernel_factor,
                         private const real alpha)
{
  // prefetch(sq->xqs, 2 * QUADRATUR_ORDER);
  // prefetch(sq->yqs, 2 * QUADRATUR_ORDER);
  // prefetch(sq->wqs,     QUADRATUR_ORDER);

  for(uint i = 0; i < gcii->ridx_size; ++i)
  {
    const uint ii = gcii->ridx[i];

    const real  gr = sur->g[ii];

    const uint3 q  = vload3(ii, sur->p);

    private real3 x[3];

    x[0] = vload3(q.x, sur->v);
    x[1] = vload3(q.y, sur->v);
    x[2] = vload3(q.z, sur->v);

    real res = r_zero;

    for(uint j = 0; j < gcii->cidx_size; ++j)
    {
      const uint  jj = gcii->cidx[j];

      const real  gc = sur->g[jj];
      const uint3 p  = vload3(jj, sur->p);

      const real factor = (ii == jj) * 0.5 * bem_alpha * gr;

      private real3 y[3];

      y[0] = vload3(p.x, sur->v);
      y[1] = vload3(p.y, sur->v);
      y[2] = vload3(p.z, sur->v);

      real sum = sq->bases.x;

      uint q = 0;

#ifndef USE_FLOAT
      for(; (q + 2) < sq->nq; q += 2)
      {
        real2  tmp = vload2(q / 2, sq->xqs);
        float2 t   = { convert_float(tmp.x), convert_float(tmp.y) };

        tmp        = vload2(q / 2, sq->xqs + sq->nq);
        float2 s   = { convert_float(tmp.x), convert_float(tmp.y) };

        private float2 xx[3] = { { 0.0f, 0.0f },
                                 { 0.0f, 0.0f },
                                 { 0.0f, 0.0f } };

        float scalar = 0.0f;

        scalar = convert_float(x[0].x); xx[0] = (1.0f - t) * scalar;
        scalar = convert_float(x[0].y); xx[1] = (1.0f - t) * scalar;
        scalar = convert_float(x[0].z); xx[2] = (1.0f - t) * scalar;

        scalar = convert_float(x[1].x); xx[0] = mad(t - s, scalar, xx[0]);
        scalar = convert_float(x[1].y); xx[1] = mad(t - s, scalar, xx[1]);
        scalar = convert_float(x[1].z); xx[2] = mad(t - s, scalar, xx[2]);

        scalar = convert_float(x[2].x); xx[0] = mad(s,     scalar, xx[0]);
        scalar = convert_float(x[2].y); xx[1] = mad(s,     scalar, xx[1]);
        scalar = convert_float(x[2].z); xx[2] = mad(s,     scalar, xx[2]);

        tmp = vload2(q / 2, sq->yqs);
        t   = ( convert_float(tmp.x), convert_float(tmp.y) );

        tmp = vload2(q / 2, sq->yqs + sq->nq);
        s   = ( convert_float(tmp.x), convert_float(tmp.y) );

        private float2 yy[3] = { { 0.0f, 0.0f },
                                 { 0.0f, 0.0f },
                                 { 0.0f, 0.0f } };

        scalar = convert_float(y[0].x); yy[0] = (1.0f - t) * scalar;
        scalar = convert_float(y[0].y); yy[1] = (1.0f - t) * scalar;
        scalar = convert_float(y[0].z); yy[2] = (1.0f - t) * scalar;

        scalar = convert_float(y[1].x); yy[0] = mad(t - s, scalar, yy[0]);
        scalar = convert_float(y[1].y); yy[1] = mad(t - s, scalar, yy[1]);
        scalar = convert_float(y[1].z); yy[2] = mad(t - s, scalar, yy[2]);

        scalar = convert_float(y[2].x); yy[0] = mad(s,     scalar, yy[0]);
        scalar = convert_float(y[2].y); yy[1] = mad(s,     scalar, yy[1]);
        scalar = convert_float(y[2].z); yy[2] = mad(s,     scalar, yy[2]);

        tmp    = vload2(q / 2, sq->wqs);

        const float2 result = laplace3d2(xx, yy);



        sum += convert_float(tmp.x) * result.x +
               convert_float(tmp.y) * result.y;
      }
#else

#endif
      for(; q < sq->nq; ++q)
      {
        float s, t;

#ifndef USE_FLOAT
        t  = convert_float(sq->xqs[q]);
        s  = convert_float(sq->xqs[q + sq->nq]);

        private float3 xx, yy;

        xx = (1.0f - t) * convert_float3(x[0]) +
             (t - s)    * convert_float3(x[1]) +
             s          * convert_float3(x[2]);

        t  = convert_float(sq->yqs[q]);
        s  = convert_float(sq->yqs[q + sq->nq]);

        yy = (1.0f - t) * convert_float3(y[0]) +
             (t - s)    * convert_float3(y[1]) +
             s          * convert_float3(y[2]);
#else
        t  = sq->xqs[q];
        s  = sq->xqs[q + sq->nq];

        private float3 xx, yy;

        xx = (r_one - t) * x[0] + (t - s) * x[1] + s * x[2];

        t  = sq->yqs[q];
        s  = sq->yqs[q + sq->nq];

        yy = (r_one - t) * y[0] + (t - s) * y[1] + s * y[2];
#endif

        sum += sq->wqs[q] * laplace3d(xx, yy);
      }

      res += gcii->xt[j] * ((gc * gr * kernel_factor * sum) + factor);
    }

    for(uint j = 0; j < gcii->num_h2_leafs; ++j)
    {
      barrier(CLK_GLOBAL_MEM_FENCE);

      if(get_local_id(0) == j)
        gcii->yt[i] += alpha * res;
    }
  }
}

void
mvm_on_the_fly_gca(pgcidxinfo gcii, pgeom sur, psingquadg sq,
                   private const real bem_alpha,
                   private const real kernel_factor,
                   private const real alpha)
{
  const size_t lid0 = get_local_id(0);

  for(uint i = 0; i < gcii->ridx_size; ++i)
  {
    const uint ii = gcii->ridx[i];

    const real  gr = sur->g[ii];

    const uint3 q  = vload3(ii, sur->p);

    private real3 x[3];

    x[0] = vload3(q.x, sur->v);
    x[1] = vload3(q.y, sur->v);
    x[2] = vload3(q.z, sur->v);

    real res = r_zero;

    for(uint j = 0; j < gcii->cidx_size; ++j)
    {
      const uint  jj = gcii->cidx[j];

      const real  gc = sur->g[jj];
      const uint3 p  = vload3(jj, sur->p);

      const real factor = (ii == jj) * 0.5 * bem_alpha * gr;

      private real3 y[3];

      y[0] = vload3(p.x, sur->v);
      y[1] = vload3(p.y, sur->v);
      y[2] = vload3(p.z, sur->v);

      real sum;
      uint xp[3], yp[3];

      global const real *xq, *yq, *wq;

      select_quadratureg(sq, q, p, xp, yp, &xq, &yq, &wq, &sum);

  //     for(uint k = 0; k < sq->nq; k += SIZE)
  //     {
  //       const size_t size = (k + SIZE) < sq->nq ? SIZE : sq->nq - k;
  //
  //       barrier(CLK_LOCAL_MEM_FENCE);
  //
  //       if(lid0 < size)
  //       {
  //         sq->xql[lid0]        = xq[k + lid0];
  //         sq->xql[lid0 + SIZE] = xq[k + lid0 + sq->nq];
  //         sq->yql[lid0]        = yq[k + lid0];
  //         sq->yql[lid0 + SIZE] = yq[k + lid0 + sq->nq];
  //         sq->wql[lid0]        = wq[k + lid0];
  //       }
  //
  //       barrier(CLK_LOCAL_MEM_FENCE);
  //
  //       for(uint q = 0; q < size; ++q)
  //       {
  //         float s, t;
  //
  // #ifndef USE_FLOAT
  //         t  = convert_float(sq->xql[q]);
  //         s  = convert_float(sq->xql[q + SIZE]);
  //
  //         private float3 xx, yy;
  //
  //         xx = (1.0f - t) * convert_float3(x[xp[0]]) +
  //              (t - s)    * convert_float3(x[xp[1]]) +
  //              s          * convert_float3(x[xp[2]]);
  //
  //         t  = convert_float(sq->yql[q]);
  //         s  = convert_float(sq->yql[q + SIZE]);
  //
  //         yy = (1.0f - t) * convert_float3(y[yp[0]]) +
  //              (t - s)    * convert_float3(y[yp[1]]) +
  //              s          * convert_float3(y[yp[2]]);
  // #else
  //         t  = xq1[q];
  //         s  = xq2[q];
  //
  //         xx = (r_one - t) * x[xp[0]] + (t - s) * x[xp[1]] + s * x[xp[2]];
  //
  //         t  = yq1[q];
  //         s  = yq2[q];
  //
  //         yy = (r_one - t) * y[yp[0]] + (t - s) * y[yp[1]] + s * y[yp[2]];
  // #endif
  //
  //         sum += sq->wql[q] * laplace3d(xx, yy);
  //       }
  //     }

      for(uint q = 0; q < sq->nq; ++q)
      {
        float s, t;

#ifndef USE_FLOAT
        t  = convert_float(xq[q]);
        s  = convert_float(xq[q + sq->nq]);

        private float3 xx, yy;

        xx = (1.0f - t) * convert_float3(x[xp[0]]) +
             (t - s)    * convert_float3(x[xp[1]]) +
             s          * convert_float3(x[xp[2]]);

        t  = convert_float(yq[q]);
        s  = convert_float(yq[q + sq->nq]);

        yy = (1.0f - t) * convert_float3(y[yp[0]]) +
             (t - s)    * convert_float3(y[yp[1]]) +
             s          * convert_float3(y[yp[2]]);
#else
        t  = xq[q];
        s  = xq[q + sq->nq];

        private float3 xx, yy;

        xx = (r_one - t) * x[xp[0]] + (t - s) * x[xp[1]] + s * x[xp[2]];

        t  = yq[q];
        s  = yq[q + sq->nq];

        yy = (r_one - t) * y[yp[0]] + (t - s) * y[yp[1]] + s * y[yp[2]];
#endif

        sum += wq[q] * laplace3d(xx, yy);
      }

      res += gcii->xt[j] * ((gc * gr * kernel_factor * sum) + factor);
    }

    for(uint j = 0; j < gcii->num_h2_leafs; ++j)
    {
      barrier(CLK_GLOBAL_MEM_FENCE);

      if(get_local_id(0) == j)
        gcii->yt[i] += alpha * res;
    }
  }
}

#endif // FASTADDEVAL_CL
