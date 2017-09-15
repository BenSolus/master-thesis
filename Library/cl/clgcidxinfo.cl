/* ------------------------------------------------------------
 * This is the file "clgcidxinfo.cl" of this master thesis.
 * All rights reserved, Bennet Carstensen 2017
 * ------------------------------------------------------------ */

/**
 * @file      cl/clgcidxinfo.cl
 * @author    Bennet Carstensen
 * @date      2017
 * @copyright All rights reserved, Bennet Carstensen 2017
 */

#ifndef IDXINFO_CL
#define IDXINFO_CL

#include "clbasic.cl"

/** @brief @ref quad is just an abbreviation for the struct @ref _quad. */
typedef struct _gcidxinfo gcidxinfo;

/** @brief Pointer to a @ref _quad "quad" object. */
typedef gcidxinfo *pgcidxinfo;

/** @brief Pointer to a constant @ref _quad "quad" object. */
typedef const gcidxinfo *pcgcidxinfo;

struct _gcidxinfo
{
               uint num_h2_leafs;

               uint idx_off;

               uint ridx_size;

               uint ridx_off;

  global const uint *ridx;

               uint yt_off;

  global       real *yt;

  local        real *ytl;

               uint cidx_size;

               uint cidx_off;

  global const uint *cidx;

               uint xt_off;

  global const real *xt;
};

void
init_row_gcidxinfo(       pgcidxinfo gcii,
                          const uint num_h2_leafs,
                          const uint idx_off,
                          const uint ridx_size,
                          const uint ridx_off,
                   global const uint *ridxs,
                          const uint yt_off,
                   global       real *yts,
                   local        real *ytl)
{
  gcii->num_h2_leafs = num_h2_leafs;
  gcii->idx_off      = idx_off;
  gcii->ridx_size    = ridx_size;
  gcii->ridx_off     = ridx_off;
  gcii->ridx         = ridxs + ridx_off;
  gcii->yt_off       = yt_off;
  gcii->yt           = yts + yt_off;
  gcii->ytl          = ytl;

  gcii->cidx_size    = 0;
  gcii->cidx_off     = 0;
  gcii->cidx         = 0;
  gcii->xt_off       = 0;
  gcii->xt           = 0;
}

void
set_column_info_gcidxinfo(       pgcidxinfo gcii,
                                 const uint cidx_size,
                                 const uint cidx_off,
                          global const uint *cidxs,
                                 const uint xt_off,
                          global const real *xts)
{
  gcii->cidx_size = cidx_size;
  gcii->cidx_off  = cidx_off;
  gcii->cidx      = cidxs + cidx_off;
  gcii->xt_off    = xt_off;
  gcii->xt        = xts + xt_off;
}

#endif // IDXINFO_CL
