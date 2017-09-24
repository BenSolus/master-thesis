/* ------------------------------------------------------------
 * This is the file "clintegralinfo.cl" of this master thesis.
 * All rights reserved, Bennet Carstensen 2017
 * ------------------------------------------------------------ */

/**
 * @file      cl/clintegralinfo.cl
 * @author    Bennet Carstensen
 * @date      2017
 * @copyright All rights reserved, Bennet Carstensen 2017
 */

#ifndef INTEGRALINFO_CL
#define INTEGRALINFO_CL

#include "clbasic.cl"
#include "clgcidxinfo.cl"
#include "clgeom.cl"
#include "clkernels.cl"
#include "clsingquad.cl"

/** @brief @ref intinfo is just an abbreviation for the struct @ref _intinfo. */
typedef struct _intinfo intinfo;

/** @brief Pointer to a @ref _intinfo "intinfo" object. */
typedef intinfo *pintinfo;

/** @brief Pointer to a constant @ref _intinfo "intinfo" object. */
typedef const intinfo *pcintinfo;

struct _intinfo
{
  uint num_integrals;

  global const uint *rows;

  global const uint *cols;

  global const uint *cidx;
};

void
init_intinfo(             pintinfo iinfo,
                    const uint     num_integrals,
                    const uint     idx_off,
             global const uint     *rows,
             global const uint     *cols,
             global const uint     *cidxs)
{
  iinfo->num_integrals = num_integrals;
  iinfo->rows          = rows  + idx_off;
  iinfo->cols          = cols  + idx_off;
  iinfo->cidx          = cidxs + idx_off;
}

void
eval_integrals(pgcidxinfo gcii, pcintinfo iinfo, pcgeom sur, psingquadc sq,
               private const real bem_alpha,
               private const real kernel_factor,
               private const real alpha,
               global  const real *xts)
{
  const size_t lid0 = get_local_id(0);

  for(uint k = 0; k < iinfo->num_integrals; ++k)
  {
    const uint i  = iinfo->rows[k];
    const uint j  = iinfo->cols[k];

    const uint ii = gcii->ridx[i];
    const uint jj = iinfo->cidx[k];

    const real gr     = sur->g[ii];

    const real factor = gr * sur->g[jj] * kernel_factor;

    const uint3 p = vload3(ii, sur->p);
    const uint3 q = vload3(jj, sur->p);

    private const real3 x[3] = { vload3(p.x, sur->v),
                                 vload3(p.y, sur->v),
                                 vload3(p.z, sur->v) };

    private const real3 y[3] = { vload3(q.x, sur->v),
                                 vload3(q.y, sur->v),
                                 vload3(q.z, sur->v) };

    real base;
    uint xp[3], yp[3];

    constant real *xq, *yq, *wq;

    select_quadraturec(sq, p, q, xp, yp, &xq, &yq, &wq, &base);

    const uint nq = sq->nq;

    real sum = r_zero;

    gcii->ytl[0] = r_zero;

    for(uint q = 0; q < nq; q += get_local_size(0))
    {
      if((q + lid0) < nq)
      {
        float s, t;

#ifndef USE_FLOAT
        t  = convert_float(xq[q + lid0]);
        s  = convert_float(xq[q + lid0 + nq]);

        private float3 xx, yy;

        xx = (1.0f - t) * convert_float3(x[xp[0]]) +
             (t - s)    * convert_float3(x[xp[1]]) +
             s          * convert_float3(x[xp[2]]);

        t  = convert_float(yq[q + lid0]);
        s  = convert_float(yq[q + lid0 + nq]);

        yy = (1.0f - t) * convert_float3(y[yp[0]]) +
             (t - s)    * convert_float3(y[yp[1]]) +
             s          * convert_float3(y[yp[2]]);
#else
        t  = xq[q + lid0];
        s  = xq[q + lid0 + nq];

        private float3 xx, yy;

        xx = (r_one - t) * x[xp[0]] + (t - s) * x[xp[1]] + s * x[xp[2]];

        t  = yq[q + lid0];
        s  = yq[q + lid0 + sq->nq];

        yy = (r_one - t) * y[yp[0]] + (t - s) * y[yp[1]] + s * y[yp[2]];
#endif

        sum += wq[q + lid0] * laplace3d(xx, yy);
      }
    }

    for(uint j = 0; j < get_local_size(0); ++j)
    {
      barrier(CLK_LOCAL_MEM_FENCE);

      if(lid0 == j)
        gcii->ytl[0] += factor * sum;
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    if(lid0 == 0)
      gcii->yt[i] +=
        alpha * xts[j] * (gcii->ytl[0] + (ii == jj) * 0.5 * bem_alpha * gr);
  }
}

#endif // INTEGRALINFO_CL
