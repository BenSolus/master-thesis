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

kernel void
fastaddeval_h2matrix_avector(
                             global const void *x,
                             global const void *p,
                             global const real *g,
                                    const uint num_row_leafs,
                             global const uint *num_h2_leafs_per_cluster,
                             global const uint *idx_offs,
                             global const uint *buf_offs,
                             global const uint *ridxs_sizes,
                             global const uint *cidxs_sizes,
                             global const uint *ridxs_offs,
                             global const uint *cidxs_offs,
                             global const uint *ridxs,
                             global const uint *cidxs,
                             global const uint *xt_offs,
                             global const uint *yt_offs,
                             global const real *xts,
                             global       real *yts)
{
  const size_t i = get_global_id(0);
  const size_t j = get_global_id(1);

  const size_t size = get_local_size(0);

  if((i >= num_row_leafs) || (j >= num_h2_leafs_per_cluster[i]))
    return;
  else
  {
    const uint idx_off      = idx_offs[i];
    const uint buf_off      = buf_offs[idx_off + j];

    const uint ridx_size    = ridxs_sizes[idx_off + j];
    const uint cidx_size    = cidxs_sizes[idx_off + j];

    global const uint* ridx = ridxs + ridxs_offs[idx_off + j];
    global const uint* cidx = cidxs + cidxs_offs[idx_off + j];

    global const real* xt   = xts + xt_offs[idx_off + j];
    global const real* yt   = yts + yt_offs[idx_off + j];

    // Calculate...
  }
}

#endif
