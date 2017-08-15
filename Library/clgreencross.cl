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

#ifndef CLGREENCROSS_CL
#define CLGREENCROSS_CL

/** @addtogroup greencross
 *  @{ */

/*  @brief The source code to performe operations on the module
 *         @ref greencross via OpenCL as a string. */
static const char clgreencross_src[] =
{
  "#ifndef GREENCROSS_CL\n"
  "#define GREENCROSS_CL\n"
  "\n"
  "\n"
  "#ifdef USE_FLOAT\n"
  "typedef float4  real4;\n"
  "typedef float8  real8;\n"
  "typedef float16 real16;\n"
  "#else\n"
  "typedef double4  real4;\n"
  "typedef double8  real8;\n"
  "typedef double16 real16;\n"
  "#endif\n"
  "\n"
  "constant static real r_four_pi = 0.0795774715459476;\n"
  "\n"
  "real\n"
  "laplace3d(private real3 x, private real3 y)\n"
  "{\n"
  "  private real3 xmy2 = x - y;\n"
  "  xmy2 *= xmy2;\n"
  "\n"
  "\n"
  "  const real norm2 = xmy2.x + xmy2.y + xmy2.z;\n"
  "\n"
  "  real norm;\n"
  "\n"
  "#ifndef USE_FLOAT\n"
  "  norm = convert_double(native_rsqrt(convert_float(norm2)));\n"
  "#else\n"
  "  norm = native_rsqrt(norm2);\n"
  "#endif\n"
  "\n"
  "  return r_four_pi * norm;\n"
  "}\n"
  "\n"
  "void\n"
  "mvm_sub_row(        pgeom      sur,\n"
  "                    pquad      quad,\n"
  "                    pgcidxinfo  gcii,\n"
  "                    const uint  i,\n"
  "            global  const uint  *cidx_sizes,\n"
  "            global  const uint  *cidx_offs,\n"
  "            global  const uint  *cidxs,\n"
  "            global  const uint  *xt_offs,\n"
  "            global  const real  *xts,\n"
  "            private       real3 x[3],\n"
  "            private       real3 y[3],\n"
  "            private       real3 xx,\n"
  "            private       real3 yy)\n"
  "{\n"
  "  const size_t lid0    = get_local_id(0);\n"
  "\n"
  "  local real3 (*v_tmp)[3] = (local real3 (*)[3]) sur->v_tmp;\n"
  "\n"
  "  real gr;\n"
  "\n"
  "  if((i + lid0) < gcii->ridx_size)\n"
  "  {\n"
  "    const uint ii = gcii->ridx[i + lid0];\n"
  "    gr = sur->g[ii];\n"
  "\n"
  "    private uint3 p = vload3(ii, sur->p);\n"
  "\n"
  "    x[0] = vload3(p.x, sur->v);\n"
  "    x[1] = vload3(p.y, sur->v);\n"
  "    x[2] = vload3(p.z, sur->v);\n"
  "  }\n"
  "\n"
  "  real res = r_zero;\n"
  "\n"
  "  for(uint k = 0; k < gcii->num_h2_leafs; ++k)\n"
  "  {\n"
  "    set_column_info_gcidxinfo(gcii,\n"
  "                              cidx_sizes[gcii->idx_off + k],\n"
  "                              cidx_offs[gcii->idx_off + k],\n"
  "                              cidxs,\n"
  "                              xt_offs[gcii->idx_off + k],\n"
  "                              xts);\n"
  "\n"
  "\n"
  "    for(uint l = 0; l < gcii->cidx_size; l += WRK_GRP_SIZE0)\n"
  "    {\n"
  "      const uint loops = ((l + WRK_GRP_SIZE0) < gcii->cidx_size)\n"
  "                           ? WRK_GRP_SIZE0 : (gcii->cidx_size - WRK_GRP_SIZE0);\n"
  "\n"
  "      barrier(CLK_LOCAL_MEM_FENCE);\n"
  "\n"
  "      if(lid0 < loops)\n"
  "      {\n"
  "        const uint jj = gcii->cidx[l + lid0];\n"
  "        const uint3 p = vload3(jj, sur->p);\n"
  "\n"
  "        //local real3 (*v_tmp)[3] = (local real3 (*)[3]) sur->v_tmp;\n"
  "\n"
  "        sur->g_tmp[lid0]   = sur->g[jj];\n"
  "\n"
  "        v_tmp[lid0][0] = vload3(p.x, sur->v);\n"
  "        v_tmp[lid0][1] = vload3(p.y, sur->v);\n"
  "        v_tmp[lid0][2] = vload3(p.z, sur->v);\n"
  "\n"
  "        gcii->xt_tmp[lid0] = gcii->xt[l + lid0];\n"
  "\n"
  "      }\n"
  "\n"
  "      barrier(CLK_LOCAL_MEM_FENCE);\n"
  "\n"
  "      if((i + lid0) < gcii->ridx_size)\n"
  "      {\n"
  "        for(uint j = 0; j < loops; ++j)\n"
  "        {\n"
  "          const real gc = sur->g_tmp[j];\n"
  "\n"
  "          y[0] = vload3(0, (local real*) &v_tmp[j][0]);\n"
  "          y[1] = vload3(0, (local real*) &v_tmp[j][1]);\n"
  "          y[2] = vload3(0, (local real*) &v_tmp[j][2]);\n"
  "          // private uint3 p = vload3(jj, sur->p);\n"
  "          //\n"
  "          // y[0] = vload3(p.x, sur->v);\n"
  "          // y[1] = vload3(p.y, sur->v);\n"
  "          // y[2] = vload3(p.z, sur->v);\n"
  "\n"
  "          real sum = r_zero;\n"
  "\n"
  "          for(uint q = 0; q < quad->nq; ++q)\n"
  "          {\n"
  "            real t = quad->qx[q];\n"
  "            real s = quad->qx[q + quad->nq];\n"
  "\n"
  "            xx = (r_one - t) * x[0] + (t - s) * x[1] + s * x[2];\n"
  "\n"
  "            t = quad->qy[q];\n"
  "            s = quad->qy[q + quad->nq];\n"
  "\n"
  "            yy = (r_one - t) * y[0] + (t - s) * y[1] + s * y[2];\n"
  "\n"
  "            sum += quad->w[q] * laplace3d(xx, yy);\n"
  "          }\n"
  "\n"
  "          res += gcii->xt_tmp[j] * gc * gr * sum;\n"
  "        }\n"
  "      }\n"
  "    }\n"
  "  }\n"
  "\n"
  "    // const uint loops = gcii->cidx_size / WRK_GRP_SIZE0;\n"
  "    //\n"
  "    // uint l = 0;\n"
  "    //\n"
  "    // for(l = 0; l < loops; ++l)\n"
  "    // {\n"
  "    //   gcii->xt_tmp[lid0] = gcii->xt[l * WRK_GRP_SIZE0 + lid0];\n"
  "    //\n"
  "    //   if((i + lid0) < gcii->ridx_size)\n"
  "    //   {\n"
  "    //     for(uint j = 0; j < WRK_GRP_SIZE0; ++j)\n"
  "    //     {\n"
  "    //       const uint jj = gcii->cidx[l * WRK_GRP_SIZE0 + j]; // Bank conflict\n"
  "    //       const real gc = sur->g[jj];\n"
  "    //\n"
  "    //       private uint3 p = vload3(jj, sur->p);\n"
  "    //\n"
  "    //       y[0] = vload3(p.x, sur->v);\n"
  "    //       y[1] = vload3(p.y, sur->v);\n"
  "    //       y[2] = vload3(p.z, sur->v);\n"
  "    //\n"
  "    //       real sum = r_zero;\n"
  "    //\n"
  "    //       for(uint q = 0; q < quad->nq; ++q)\n"
  "    //       {\n"
  "    //         real t = quad->qx[q];\n"
  "    //         real s = quad->qx[q + quad->nq];\n"
  "    //\n"
  "    //         xx = (r_one - t) * x[0] + (t - s) * x[1] + s * x[2];\n"
  "    //\n"
  "    //         t = quad->qy[q];\n"
  "    //         s = quad->qy[q + quad->nq];\n"
  "    //\n"
  "    //         yy = (r_one - t) * y[0] + (t - s) * y[1] + s * y[2];\n"
  "    //\n"
  "    //         sum += quad->w[q] * laplace3d(xx, yy);\n"
  "    //       }\n"
  "    //\n"
  "    //       res += gcii->xt[l * WRK_GRP_SIZE0 + j] * gc * gr * sum;\n"
  "    //     }\n"
  "    //   }\n"
  "    // }\n"
  "    //\n"
  "    // if(lid0 < (gcii->cidx_size - l * WRK_GRP_SIZE0))\n"
  "    // {\n"
  "    //   gcii->xt_tmp[lid0] = gcii->xt[l * WRK_GRP_SIZE0 + lid0];\n"
  "    // }\n"
  "    //\n"
  "    // barrier(CLK_LOCAL_MEM_FENCE);\n"
  "    //\n"
  "    // if((i + lid0) < gcii->ridx_size)\n"
  "    // {\n"
  "    //   for(uint j = 0; j < (gcii->cidx_size - l * WRK_GRP_SIZE0); ++j)\n"
  "    //   {\n"
  "    //     const uint jj = gcii->cidx[l * WRK_GRP_SIZE0 + j]; // Bank conflict\n"
  "    //     const real gc = sur->g[jj];\n"
  "    //\n"
  "    //     private uint3 p = vload3(jj, sur->p);\n"
  "    //\n"
  "    //     y[0] = vload3(p.x, sur->v);\n"
  "    //     y[1] = vload3(p.y, sur->v);\n"
  "    //     y[2] = vload3(p.z, sur->v);\n"
  "    //\n"
  "    //     real sum = r_zero;\n"
  "    //\n"
  "    //     for(uint q = 0; q < quad->nq; ++q)\n"
  "    //     {\n"
  "    //       real t = quad->qx[q];\n"
  "    //       real s = quad->qx[q + quad->nq];\n"
  "    //\n"
  "    //       xx = (r_one - t) * x[0] + (t - s) * x[1] + s * x[2];\n"
  "    //\n"
  "    //       t = quad->qy[q];\n"
  "    //       s = quad->qy[q + quad->nq];\n"
  "    //\n"
  "    //       yy = (r_one - t) * y[0] + (t - s) * y[1] + s * y[2];\n"
  "    //\n"
  "    //       sum += quad->w[q] * laplace3d(xx, yy);\n"
  "    //     }\n"
  "    //\n"
  "    //     res += gcii->xt[l * WRK_GRP_SIZE0 + j] * gc * gr * sum;\n"
  "    //   }\n"
  "    // }\n"
  "\n"
  "  if((i + lid0) < gcii->ridx_size)\n"
  "    gcii->yt[i + lid0] += res;\n"
  "}\n"
  "\n"
  "kernel void\n"
  "fastaddeval_h2matrix_avector_1(       const uint dim,\n"
  "                                      const uint n,\n"
  "                               global const real *vs,\n"
  "                               global const uint *p,\n"
  "                               global const real *g,\n"
  "                                      const uint nq,\n"
  "                               constant     real *qx,\n"
  "                               constant     real *qy,\n"
  "                               constant     real *w,\n"
  "                                      const uint num_row_leafs,\n"
  "                               global const uint *num_h2_leafs_per_cluster,\n"
  "                               global const uint *idx_offs,\n"
  "                               global const uint *ridx_sizes,\n"
  "                               global const uint *cidx_sizes,\n"
  "                               global const uint *ridx_offs,\n"
  "                               global const uint *cidx_offs,\n"
  "                               global const uint *ridxs,\n"
  "                               global const uint *cidxs,\n"
  "                                      const real alpha,\n"
  "                               global const uint *xt_offs,\n"
  "                               global const uint *yt_offs,\n"
  "                               global const real *xts,\n"
  "                               global       real *yt)\n"
  "{\n"
  "  const size_t grpid0 = get_group_id(0);\n"
  "\n"
  "  if(grpid0 >= num_row_leafs)\n"
  "    return;\n"
  "  else\n"
  "  {\n"
  "    local real g_tmp[WRK_GRP_SIZE0], xt_tmp[WRK_GRP_SIZE0];\n"
  "    local real3 v_tmp[WRK_GRP_SIZE0][3];\n"
  "\n"
  "    gcidxinfo gcii;\n"
  "    geom      sur;\n"
  "    quad      q;\n"
  "\n"
  "    init_geom(&sur, dim, n, vs, p, g, (local void *) v_tmp, g_tmp);\n"
  "    init_quad(&q, nq, qx, qy, w);\n"
  "    init_row_gcidxinfo(&gcii,\n"
  "                       num_h2_leafs_per_cluster[grpid0],\n"
  "                       idx_offs[grpid0],\n"
  "                       ridx_sizes[grpid0],\n"
  "                       ridx_offs[grpid0],\n"
  "                       ridxs,\n"
  "                       yt_offs[grpid0],\n"
  "                       yt,\n"
  "                       xt_tmp);\n"
  "\n"
  "    private real3 x[3], y[3], xx, yy;\n"
  "\n"
  "    uint i = 0;\n"
  "\n"
  "    for(i = 0; i < gcii.ridx_size; i += WRK_GRP_SIZE0)\n"
  "      mvm_sub_row(&sur,\n"
  "                  &q,\n"
  "                  &gcii,\n"
  "                  i,\n"
  "                  cidx_sizes,\n"
  "                  cidx_offs,\n"
  "                  cidxs,\n"
  "                  xt_offs,\n"
  "                  xts,\n"
  "                  x,\n"
  "                  y,\n"
  "                  xx,\n"
  "                  yy);\n"
  "  }\n"
  "}\n"
  "\n"
  "void\n"
  "clPrefixSum(local real* x, const size_t n)\n"
  "{\n"
  "  const size_t lid0 = get_local_id(0);\n"
  "\n"
  "  for(uint i = 2; i <= n; i <<= 1)\n"
  "  {\n"
  "    if((lid0 % i) == (i - 1))\n"
  "      x[lid0] += x[lid0 - (i >> 1)];\n"
  "\n"
  "    barrier(CLK_LOCAL_MEM_FENCE);\n"
  "  }\n"
  "\n"
  "  for(uint i = (n >> 1); i >= 2; i >>= 1)\n"
  "  {\n"
  "    if((lid0 >= i) && ((lid0 % i) == ((i >> 1) - 1)))\n"
  "      x[lid0] += x[lid0 - (i >> 1)];\n"
  "\n"
  "    barrier(CLK_LOCAL_MEM_FENCE);\n"
  "  }\n"
  "}\n"
  "\n"
  "kernel void\n"
  "fastaddeval_h2matrix_avector_2(       const uint dim,\n"
  "                                      const uint n,\n"
  "                               global const real *vs,\n"
  "                               global const uint *p,\n"
  "                               global const real *g,\n"
  "                                      const uint nq,\n"
  "                               constant     real *qx,\n"
  "                               constant     real *qy,\n"
  "                               constant     real *w,\n"
  "                                      const uint num_row_leafs,\n"
  "                               global const uint *num_h2_leafs_per_cluster,\n"
  "                               global const uint *idx_offs,\n"
  "                               global const uint *ridx_sizes,\n"
  "                               global const uint *cidx_sizes,\n"
  "                               global const uint *ridx_offs,\n"
  "                               global const uint *cidx_offs,\n"
  "                               global const uint *ridxs,\n"
  "                               global const uint *cidxs,\n"
  "                                      const real alpha,\n"
  "                               global const uint *xt_offs,\n"
  "                               global const uint *yt_offs,\n"
  "                               global const real *xts,\n"
  "                               global       real *yt)\n"
  "{\n"
  "  const size_t grpid0 = get_group_id(0);\n"
  "\n"
  "  if(grpid0 >= num_row_leafs)\n"
  "    return;\n"
  "  else\n"
  "  {\n"
  "    const uint   num_h2_leafs = num_h2_leafs_per_cluster[grpid0];\n"
  "    const size_t lid0         = get_local_id(0);\n"
  "\n"
  "    gcidxinfo gcii;\n"
  "    geom      sur;\n"
  "    quad      quad;\n"
  "\n"
  "    local   real  sum[WRK_GRP_SIZE0], yt_tmp[WRK_GRP_SIZE0];\n"
  "    private real3 x[3], y[3], xx, yy;\n"
  "\n"
  "    event_t event = 0;\n"
  "\n"
  "    init_geom(&sur, dim, n, vs, p, g, 0, 0);\n"
  "    init_quad(&quad, nq, qx, qy, w);\n"
  "\n"
  "    /* Same information for all threads in a work group */\n"
  "    init_row_gcidxinfo(&gcii,\n"
  "                       num_h2_leafs,\n"
  "                       idx_offs[grpid0],\n"
  "                       ridx_sizes[grpid0],\n"
  "                       ridx_offs[grpid0],\n"
  "                       ridxs,\n"
  "                       yt_offs[grpid0],\n"
  "                       yt,\n"
  "                       0);\n"
  "\n"
  "    /* Different information for all threads in a work group */\n"
  "    if(lid0 < num_h2_leafs)\n"
  "      set_column_info_gcidxinfo(&gcii,\n"
  "                                cidx_sizes[gcii.idx_off + lid0],\n"
  "                                cidx_offs[gcii.idx_off + lid0],\n"
  "                                cidxs,\n"
  "                                xt_offs[gcii.idx_off + lid0],\n"
  "                                xts);\n"
  "\n"
  "    for(uint k = 0; k < gcii.ridx_size; k += WRK_GRP_SIZE0)\n"
  "    {\n"
  "      const uint loops = (k + WRK_GRP_SIZE0) < gcii.ridx_size\n"
  "                           ? WRK_GRP_SIZE0\n"
  "                           : gcii.ridx_size - k;\n"
  "\n"
  "      if(lid0 < num_h2_leafs)\n"
  "      {\n"
  "        for(uint i = 0; i < loops; ++i)\n"
  "        {\n"
  "          const uint ii = gcii.ridx[k + i];\n"
  "          const real gr = sur.g[ii];\n"
  "\n"
  "          private uint3 p = vload3(ii, sur.p);\n"
  "\n"
  "          x[0] = vload3(p.x, sur.v);\n"
  "          x[1] = vload3(p.y, sur.v);\n"
  "          x[2] = vload3(p.z, sur.v);\n"
  "\n"
  "          real res = r_zero;\n"
  "\n"
  "          for(uint j = 0; j < gcii.cidx_size; ++j)\n"
  "          {\n"
  "            const uint jj = gcii.cidx[j];\n"
  "            const real gc = sur.g[jj];\n"
  "\n"
  "            private uint3 p = vload3(jj, sur.p);\n"
  "\n"
  "            y[0] = vload3(p.x, sur.v);\n"
  "            y[1] = vload3(p.y, sur.v);\n"
  "            y[2] = vload3(p.z, sur.v);\n"
  "\n"
  "            real sum = r_zero;\n"
  "\n"
  "            for(uint q = 0; q < quad.nq; ++q)\n"
  "            {\n"
  "              real t = quad.qx[q];\n"
  "              real s = quad.qx[q + quad.nq];\n"
  "\n"
  "              xx = (r_one - t) * x[0] + (t - s) * x[1] + s * x[2];\n"
  "\n"
  "              t = quad.qy[q];\n"
  "              s = quad.qy[q + quad.nq];\n"
  "\n"
  "              yy = (r_one - t) * y[0] + (t - s) * y[1] + s * y[2];\n"
  "\n"
  "              sum += quad.w[q] * laplace3d(xx, yy);\n"
  "            }\n"
  "\n"
  "            res += gcii.xt[j] * gc * gr * sum;\n"
  "          }\n"
  "\n"
  "          /* Parallel (sum) reduction to get the entry of yt */\n"
  "          sum[lid0] = res;\n"
  "\n"
  "          uint num_leafs = gcii.num_h2_leafs;\n"
  "\n"
  "          while(num_leafs > 1)\n"
  "          {\n"
  "            const uint corrector = num_leafs % 2;\n"
  "\n"
  "            num_leafs >>= 1;\n"
  "\n"
  "            barrier(CLK_LOCAL_MEM_FENCE);\n"
  "\n"
  "            if((lid0 >= corrector) && (lid0 < (num_leafs + corrector)))\n"
  "              sum[lid0] += sum[lid0 + num_leafs];\n"
  "\n"
  "            num_leafs += corrector;\n"
  "          }\n"
  "\n"
  "          barrier(CLK_LOCAL_MEM_FENCE);\n"
  "\n"
  "          if(lid0 == 0)\n"
  "            yt_tmp[i] = alpha * sum[0];\n"
  "          //  gcii.yt[k + i] += alpha * sum[0];\n"
  "        }\n"
  "      }\n"
  "\n"
  "      wait_group_events(1, &event);\n"
  "      event = async_work_group_copy(gcii.yt + k,\n"
  "                                    yt_tmp,\n"
  "                                    loops,\n"
  "                                    0);\n"
  "    }\n"
  "\n"
  "    wait_group_events(1, &event);\n"
  "    // if(lid0 < num_h2_leafs)\n"
  "    // {\n"
  "    //   for(uint i = 0; i < gcii.ridx_size; ++i)\n"
  "    //   {\n"
  "    //     const uint ii = gcii.ridx[i];\n"
  "    //     const real gr = sur.g[ii];\n"
  "    //\n"
  "    //     private uint3 p = vload3(ii, sur.p);\n"
  "    //\n"
  "    //     x[0] = vload3(p.x, sur.v);\n"
  "    //     x[1] = vload3(p.y, sur.v);\n"
  "    //     x[2] = vload3(p.z, sur.v);\n"
  "    //\n"
  "    //     real res = r_zero;\n"
  "    //\n"
  "    //     for(uint j = 0; j < gcii.cidx_size; ++j)\n"
  "    //     {\n"
  "    //       const uint jj = gcii.cidx[j];\n"
  "    //       const real gc = sur.g[jj];\n"
  "    //\n"
  "    //       private uint3 p = vload3(jj, sur.p);\n"
  "    //\n"
  "    //       y[0] = vload3(p.x, sur.v);\n"
  "    //       y[1] = vload3(p.y, sur.v);\n"
  "    //       y[2] = vload3(p.z, sur.v);\n"
  "    //\n"
  "    //       real sum = r_zero;\n"
  "    //\n"
  "    //       for(uint q = 0; q < quad.nq; ++q)\n"
  "    //       {\n"
  "    //         real t = quad.qx[q];\n"
  "    //         real s = quad.qx[q + quad.nq];\n"
  "    //\n"
  "    //         xx = (r_one - t) * x[0] + (t - s) * x[1] + s * x[2];\n"
  "    //\n"
  "    //         t = quad.qy[q];\n"
  "    //         s = quad.qy[q + quad.nq];\n"
  "    //\n"
  "    //         yy = (r_one - t) * y[0] + (t - s) * y[1] + s * y[2];\n"
  "    //\n"
  "    //         sum += quad.w[q] * laplace3d(xx, yy);\n"
  "    //       }\n"
  "    //\n"
  "    //       res += gcii.xt[j] * gc * gr * sum;\n"
  "    //     }\n"
  "    //\n"
  "    //     tmp[lid0] = res;\n"
  "    //\n"
  "    //     uint num_leafs = gcii.num_h2_leafs;\n"
  "    //\n"
  "    //     while(num_leafs > 1)\n"
  "    //     {\n"
  "    //       const uint corrector = num_leafs % 2;\n"
  "    //\n"
  "    //       num_leafs >>= 1;\n"
  "    //\n"
  "    //       barrier(CLK_LOCAL_MEM_FENCE);\n"
  "    //\n"
  "    //       if((lid0 >= corrector) && (lid0 < (num_leafs + corrector)))\n"
  "    //         tmp[lid0] += tmp[lid0 + num_leafs];\n"
  "    //\n"
  "    //       num_leafs += corrector;\n"
  "    //     }\n"
  "    //\n"
  "    //     barrier(CLK_LOCAL_MEM_FENCE);\n"
  "    //\n"
  "    //     if(lid0 == 0)\n"
  "    //       gcii.yt[i] += alpha * tmp[0];\n"
  "    //   }\n"
  "    // }\n"
  "  }\n"
  "}\n"
  "\n"
  "kernel void\n"
  "fastaddeval_h2matrix_avector_3(       const uint dim,\n"
  "                                      const uint n,\n"
  "                               global const real *vs,\n"
  "                               global const uint *p,\n"
  "                               global const real *g,\n"
  "                                      const uint nq,\n"
  "                               constant     real *qx,\n"
  "                               constant     real *qy,\n"
  "                               constant     real *w,\n"
  "                                      const uint num_row_leafs,\n"
  "                               global const uint *num_h2_leafs_per_cluster,\n"
  "                               global const uint *idx_offs,\n"
  "                               global const uint *ridx_sizes,\n"
  "                               global const uint *cidx_sizes,\n"
  "                               global const uint *ridx_offs,\n"
  "                               global const uint *cidx_offs,\n"
  "                               global const uint *ridxs,\n"
  "                               global const uint *cidxs,\n"
  "                                      const real alpha,\n"
  "                               global const uint *xt_offs,\n"
  "                               global const uint *yt_offs,\n"
  "                               global const real *xts,\n"
  "                               global       real *yt)\n"
  "{\n"
  "  const size_t grpid0 = get_group_id(0);\n"
  "\n"
  "  if(grpid0 >= num_row_leafs)\n"
  "    return;\n"
  "  else\n"
  "  {\n"
  "    const uint   num_h2_leafs = num_h2_leafs_per_cluster[grpid0];\n"
  "    const size_t lid0         = get_local_id(0);\n"
  "\n"
  "    gcidxinfo gcii;\n"
  "    geom      sur;\n"
  "    quad      quad;\n"
  "\n"
  "    local   real  sum[WRK_GRP_SIZE0], yt_tmp[WRK_GRP_SIZE0];\n"
  "    private real3 x[3], y[3], xx, yy;\n"
  "\n"
  "    event_t event = 0;\n"
  "\n"
  "    init_geom(&sur, dim, n, vs, p, g, 0, 0);\n"
  "    init_quad(&quad, nq, qx, qy, w);\n"
  "\n"
  "    /* Same information for all threads in a work group */\n"
  "    init_row_gcidxinfo(&gcii,\n"
  "                       num_h2_leafs,\n"
  "                       idx_offs[grpid0],\n"
  "                       ridx_sizes[grpid0],\n"
  "                       ridx_offs[grpid0],\n"
  "                       ridxs,\n"
  "                       yt_offs[grpid0],\n"
  "                       yt,\n"
  "                       0);\n"
  "\n"
  "    /* Different information for all threads in a work group */\n"
  "    if(lid0 < num_h2_leafs)\n"
  "      set_column_info_gcidxinfo(&gcii,\n"
  "                                cidx_sizes[gcii.idx_off + lid0],\n"
  "                                cidx_offs[gcii.idx_off + lid0],\n"
  "                                cidxs,\n"
  "                                xt_offs[gcii.idx_off + lid0],\n"
  "                                xts);\n"
  "\n"
  "    for(uint k = 0; k < gcii.ridx_size; k += WRK_GRP_SIZE0)\n"
  "    {\n"
  "      const uint loops = (k + WRK_GRP_SIZE0) < gcii.ridx_size\n"
  "                           ? WRK_GRP_SIZE0\n"
  "                           : gcii.ridx_size - k;\n"
  "\n"
  "      if(lid0 < num_h2_leafs)\n"
  "      {\n"
  "        for(uint i = 0; i < loops; ++i)\n"
  "        {\n"
  "          const uint ii = gcii.ridx[k + i];\n"
  "          const real gr = sur.g[ii];\n"
  "\n"
  "          private uint3 p = vload3(ii, sur.p);\n"
  "\n"
  "          x[0] = vload3(p.x, sur.v);\n"
  "          x[1] = vload3(p.y, sur.v);\n"
  "          x[2] = vload3(p.z, sur.v);\n"
  "\n"
  "          real res = r_zero;\n"
  "\n"
  "          for(uint j = 0; j < gcii.cidx_size; ++j)\n"
  "          {\n"
  "            if(j != (gcii.cidx_size - 1))\n"
  "            {\n"
  "              const uint jj = gcii.cidx[j + 1];\n"
  "\n"
  "              prefetch(sur.g + jj, 1);\n"
  "\n"
  "              private uint3 p = vload3(jj, sur.p);\n"
  "\n"
  "              prefetch(sur.v + sur.dim * p.x, 1);\n"
  "              prefetch(sur.v + sur.dim * p.y, 1);\n"
  "              prefetch(sur.v + sur.dim * p.z, 1);\n"
  "            }\n"
  "\n"
  "            const uint jj = gcii.cidx[j];\n"
  "            const real gc = sur.g[jj];\n"
  "\n"
  "            private uint3 p = vload3(jj, sur.p);\n"
  "\n"
  "            y[0] = vload3(p.x, sur.v);\n"
  "            y[1] = vload3(p.y, sur.v);\n"
  "            y[2] = vload3(p.z, sur.v);\n"
  "\n"
  "            real sum = r_zero;\n"
  "\n"
  "            for(uint q = 0; q < quad.nq; ++q)\n"
  "            {\n"
  "              real t = quad.qx[q];\n"
  "              real s = quad.qx[q + quad.nq];\n"
  "\n"
  "              xx = (r_one - t) * x[0] + (t - s) * x[1] + s * x[2];\n"
  "\n"
  "              t = quad.qy[q];\n"
  "              s = quad.qy[q + quad.nq];\n"
  "\n"
  "              yy = (r_one - t) * y[0] + (t - s) * y[1] + s * y[2];\n"
  "\n"
  "              sum += quad.w[q] * laplace3d(xx, yy);\n"
  "            }\n"
  "\n"
  "            res += gcii.xt[j] * gc * gr * sum;\n"
  "          }\n"
  "\n"
  "          /* Parallel (sum) reduction to get the entry of yt */\n"
  "          sum[lid0] = res;\n"
  "\n"
  "          uint num_leafs = gcii.num_h2_leafs;\n"
  "\n"
  "          while(num_leafs > 1)\n"
  "          {\n"
  "            const uint corrector = num_leafs % 2;\n"
  "\n"
  "            num_leafs >>= 1;\n"
  "\n"
  "            barrier(CLK_LOCAL_MEM_FENCE);\n"
  "\n"
  "            if((lid0 >= corrector) && (lid0 < (num_leafs + corrector)))\n"
  "              sum[lid0] += sum[lid0 + num_leafs];\n"
  "\n"
  "            num_leafs += corrector;\n"
  "          }\n"
  "\n"
  "          barrier(CLK_LOCAL_MEM_FENCE);\n"
  "\n"
  "          if(lid0 == 0)\n"
  "            yt_tmp[i] = alpha * sum[0];\n"
  "          //  gcii.yt[k + i] += alpha * sum[0];\n"
  "        }\n"
  "      }\n"
  "\n"
  "      wait_group_events(1, &event);\n"
  "      event = async_work_group_copy(gcii.yt + k,\n"
  "                                    yt_tmp,\n"
  "                                    loops,\n"
  "                                    0);\n"
  "    }\n"
  "\n"
  "    wait_group_events(1, &event);\n"
  "  }\n"
  "}\n"
  "\n"
  "#endif\n"
};

/** @} */

#endif // CLGREENCROSS_CL