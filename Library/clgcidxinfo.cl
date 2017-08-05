/* ------------------------------------------------------------
 * This is the file "clgcidxinfo.cl" of this master thesis.
 * All rights reserved, Steffen Boerm 2009
 * ------------------------------------------------------------ */

/**
 * @file      cl/clgcidxinfo.cl
 * @author    Bennet Carstensen
 * @date      2017
 * @copyright All rights reserved, Bennet Carstensen 2017
 */

#ifndef CLGCIDXINFO_CL
#define CLGCIDXINFO_CL

/** @addtogroup gcidxinfo
 *  @{ */

/*  @brief The source code to performe operations on the module
 *         @ref gcidxinfo via OpenCL as a string. */
static const char clgcidxinfo_src[] =
{
  "\n"
  "typedef struct _gcidxinfo gcidxinfo;\n"
  "\n"
  "typedef gcidxinfo *pgcidxinfo;\n"
  "\n"
  "typedef const gcidxinfo *pcgcidxinfo;\n"
  "\n"
  "struct _gcidxinfo\n"
  "{\n"
  "               uint num_h2_leafs;\n"
  "\n"
  "               uint idx_off;\n"
  "\n"
  "               uint ridx_size;\n"
  "\n"
  "               uint ridx_off;\n"
  "\n"
  "  global const uint *ridx;\n"
  "\n"
  "               uint yt_off;\n"
  "\n"
  "  global       real *yt;\n"
  "\n"
  "               uint cidx_size;\n"
  "\n"
  "               uint cidx_off;\n"
  "\n"
  "  global const uint *cidx;\n"
  "\n"
  "  local        uint *cidx_tmp;\n"
  "\n"
  "               uint xt_off;\n"
  "\n"
  "  global const real *xt;\n"
  "\n"
  "  local        real *xt_tmp;\n"
  "};\n"
  "\n"
  "pgcidxinfo\n"
  "new_row_gcidxinfo(       const uint num_h2_leafs,\n"
  "                         const uint idx_off,\n"
  "                         const uint ridx_size,\n"
  "                         const uint ridx_off,\n"
  "                  global const uint *ridxs,\n"
  "                         const uint yt_off,\n"
  "                  global       real *yts,\n"
  "                  local        uint *cidx_tmp,\n"
  "                  local        real *xt_tmp)\n"
  "{\n"
  "  gcidxinfo gcii;\n"
  "\n"
  "  gcii.num_h2_leafs = num_h2_leafs;\n"
  "  gcii.idx_off      = idx_off;\n"
  "  gcii.ridx_size    = ridx_size;\n"
  "  gcii.ridx_off     = ridx_off;\n"
  "  gcii.ridx         = ridxs + ridx_off;\n"
  "  gcii.yt_off       = yt_off;\n"
  "  gcii.yt           = yts + yt_off;\n"
  "\n"
  "  gcii.cidx_size    = 0;\n"
  "  gcii.cidx_off     = 0;\n"
  "  gcii.cidx         = 0;\n"
  "  gcii.cidx_tmp     = cidx_tmp;\n"
  "  gcii.xt_off       = 0;\n"
  "  gcii.xt           = 0;\n"
  "  gcii.xt_tmp       = xt_tmp;\n"
  "\n"
  "  return &gcii;\n"
  "}\n"
  "\n"
  "void\n"
  "set_column_info_gcidxinfo(       pgcidxinfo gcii,\n"
  "                                 const uint cidx_size,\n"
  "                                 const uint cidx_off,\n"
  "                          global const uint *cidxs,\n"
  "                                 const uint xt_off,\n"
  "                          global const real *xts)\n"
  "{\n"
  "  gcii->cidx_size = cidx_size;\n"
  "  gcii->cidx_off  = cidx_off;\n"
  "  gcii->cidx      = cidxs + cidx_off;\n"
  "  gcii->xt_off    = xt_off;\n"
  "  gcii->xt        = xts + xt_off;\n"
  "}\n"
};

/** @} */

#endif // CLGCIDXINFO_CL