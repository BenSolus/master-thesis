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

#ifndef CLKERNELS_CL
#define CLKERNELS_CL

/** @addtogroup kernels
 *  @{ */

/*  @brief The source code to performe operations on the module
 *         @ref kernels via OpenCL as a string. */
static const char clkernels_src[] =
{
  "\n"
  "#ifndef KERNELS_CL\n"
  "#define KERNELS_CL\n"
  "\n"
  "constant static real r_four_pi = 0.0795774715459476;\n"
  "\n"
  "float\n"
  "laplace3d(private float3 x, private float3 y)\n"
  "{\n"
  "  private float3 xmy2 = x - y;\n"
  "  xmy2 *= xmy2;\n"
  "\n"
  "\n"
  "  const float norm2 = xmy2.x + xmy2.y + xmy2.z;\n"
  "\n"
  "  float norm;\n"
  "\n"
  "  norm = native_rsqrt(norm2);\n"
  "\n"
  "  return norm;\n"
  "}\n"
  "\n"
  "float2\n"
  "laplace3d2(private float2 x[3], private float2 y[3])\n"
  "{\n"
  "  private const float2 dist[3] = { x[0] - y[0], x[1] - y[1], x[2] - y[2] };\n"
  "\n"
  "  return native_rsqrt(dist[0] * dist[0] +\n"
  "                      dist[1] * dist[1] +\n"
  "                      dist[2] * dist[2]);\n"
  "}\n"
  "\n"
  "float4\n"
  "laplace3d4(private float4 x[3], private float4 y[3])\n"
  "{\n"
  "  private const float4 dist[3] = { x[0] - y[0], x[1] - y[1], x[2] - y[2] };\n"
  "\n"
  "  return native_rsqrt(dist[0] * dist[0] +\n"
  "                      dist[1] * dist[1] +\n"
  "                      dist[2] * dist[2]);\n"
  "}\n"
  "\n"
  "float8\n"
  "laplace3d8(private float8 x[3], private float8 y[3])\n"
  "{\n"
  "  private const float8 dist[3] = { x[0] - y[0], x[1] - y[1], x[2] - y[2] };\n"
  "\n"
  "  return native_rsqrt(dist[0] * dist[0] +\n"
  "                      dist[1] * dist[1] +\n"
  "                      dist[2] * dist[2]);\n"
  "}\n"
  "\n"
  "#endif // KERNELS_CL\n"
};

/** @} */

#endif // CLKERNELS_CL