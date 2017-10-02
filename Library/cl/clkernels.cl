/* ------------------------------------------------------------
 * This is the file "clkernels.cl" of this master thesis.
 * All rights reserved, Bennet Carstensen 2017
 * ------------------------------------------------------------ */

/**
 * @file      cl/clkernels.cl
 * @author    Bennet Carstensen
 * @date      2017
 * @copyright All rights reserved, Bennet Carstensen 2017
 */

#include "clbasic.cl"

#ifndef KERNELS_CL
#define KERNELS_CL

/** @brief @f$\frac{1}{4 \pi}@f$ */
constant static real r_four_pi = 0.0795774715459476;

float
laplace3d(private float3 x, private float3 y)
{
  private float3 xmy2 = x - y;
  xmy2 *= xmy2;


  const float norm2 = xmy2.x + xmy2.y + xmy2.z;

  float norm;

  norm = native_rsqrt(norm2);

  return norm;
}

float2
laplace3d2(private float2 x[3], private float2 y[3])
{
  private const float2 dist[3] = { x[0] - y[0], x[1] - y[1], x[2] - y[2] };

  return native_rsqrt(dist[0] * dist[0] +
                      dist[1] * dist[1] +
                      dist[2] * dist[2]);
}

float4
laplace3d4(private float4 x[3], private float4 y[3])
{
  private const float4 dist[3] = { x[0] - y[0], x[1] - y[1], x[2] - y[2] };

  return native_rsqrt(dist[0] * dist[0] +
                      dist[1] * dist[1] +
                      dist[2] * dist[2]);
}

float8
laplace3d8(private float8 x[3], private float8 y[3])
{
  private const float8 dist[3] = { x[0] - y[0], x[1] - y[1], x[2] - y[2] };

  return native_rsqrt(dist[0] * dist[0] +
                      dist[1] * dist[1] +
                      dist[2] * dist[2]);
}

#endif // KERNELS_CL
