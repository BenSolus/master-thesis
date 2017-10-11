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

#include "clfastaddevalgca.cl"
#include "clintegralinfo.cl"

#ifdef USE_FLOAT
typedef float8  real8;
typedef float16 real16;
#else
typedef double8  real8;
typedef double16 real16;
#endif

kernel void
fastaddeval_h2matrix_avector_0(         const uint dim,
                                        const uint n,
                               global   const real *vs,
                               global   const uint *p,
                               global   const real *g,
                                        const uint nq,
                               global   const real *xqs,
                               global   const real *yqs,
                               global   const real *wqs,
                               global   const real *bases,
                                        const uint num_row_leafs,
                               global   const uint *rows_this_device,
                               global   const uint *num_h2_leafs_per_cluster,
                               global   const uint *idx_offs,
                               global   const uint *ridx_sizes,
                               global   const uint *cidx_sizes,
                               global   const uint *ridx_offs,
                               global   const uint *cidx_offs,
                               global   const uint *ridxs,
                               global   const uint *cidxs,
                               global   const uint *xt_offs,
                               global   const uint *yt_offs,
                                        const real bem_alpha,
                                        const real alpha,
                               global   const real *xts,
                               global         real *yt)
{
  const size_t grpid0 = get_group_id(0);

  if(grpid0 >= num_row_leafs)
    /* Stop work-groups which won't have any writing cluster assigned to. */
    return;
  else
  {
    const real kernel_factor = r_four_pi; // TODO: Move this to an kernel parameter

    /* Index of the writing cluster for this group. */
    const uint wcidx        = rows_this_device[grpid0];

    /* Get the number of leaf H^2-matrices which the writing cluster is a part
     * of. */
    const uint num_h2_leafs = num_h2_leafs_per_cluster[wcidx];

    const size_t lid0 = get_local_id(0);

    local real xql[2 * QUADRATUR_ORDER];
    local real yql[2 * QUADRATUR_ORDER];
    local real wql[QUADRATUR_ORDER];

    for(uint i = 0; i < nq; i += get_local_size(0))
    {
      if((i + lid0) < nq)
      {
        xql[i + lid0]      = xqs[i + lid0];
        xql[nq + i + lid0] = xqs[nq + i + lid0];
        yql[i + lid0]      = yqs[i + lid0];
        yql[nq + i + lid0] = yqs[nq + i + lid0];
        wql[i + lid0]      = wqs[i + lid0];
      }
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    if(lid0 >= num_h2_leafs)
      /* Stop work-items which won't have any leaf H^2-matrix assigned to. */
      return; //

    /* Information about the indices of basis functions, their number and their
     * position relative to the whole H^2-matrix for each leaf H^2-matrix. */
    gcidxinfo gcii;

    /* Geometry informations needed to perform the quadrature. */
    geom      sur;

    /* The quadrature rule needed to reconstruct the farfield. */
    singquadl sq;

    init_singquadl(&sq, dim, nq, xql, yql, wql, bases);
    init_geom(&sur, dim, n, vs, p, g, 0, 0, 0);

    /* Same H2-matrix farfield information for all threads in a work group. */
    init_row_gcidxinfo(&gcii,
                       num_h2_leafs,
                       idx_offs[wcidx],
                       ridx_sizes[wcidx],
                       ridx_offs[wcidx],
                       ridxs,
                       yt_offs[wcidx],
                       yt,
                       0);

    /* Different farfield H2-matrix information for all threads in a work group */
    set_column_info_gcidxinfo(&gcii,
                              cidx_sizes[gcii.idx_off + lid0],
                              cidx_offs[gcii.idx_off + lid0],
                              cidxs,
                              xt_offs[gcii.idx_off + lid0],
                              xts);

    /* Finally reconstruct the farfield and perform the corresponding MVM. */
    fastaddeval_farfield(&gcii, &sur, &sq, bem_alpha, kernel_factor, alpha);
  }
}

kernel void
fastaddeval_nf_common(         const uint dim,
                               const uint n,
                      global   const real *vs,
                      global   const uint *p,
                      global   const real *g,
                               const uint nq,
                      global   const real *xqs,
                      global   const real *yqs,
                      global   const real *wqs,
                      global   const real *bases,
                               const uint num_nf_writing_clusters,
                      global   const uint *nf_writings_this_device,
                      global   const uint *num_nf_h2_leafs_per_cluster,
                      global   const uint *nf_idx_offs,
                      global   const uint *nf_ridx_sizes,
                      global   const uint *nf_cidx_sizes,
                      global   const uint *nf_ridx_offs,
                      global   const uint *nf_cidx_offs,
                      global   const uint *nf_ridxs,
                      global   const uint *nf_cidxs,
                      global   const uint *nf_xt_offs,
                      global   const uint *nf_yt_offs,
                               const real bem_alpha,
                               const real alpha,
                      global   const real *xts,
                      global         real *yt)
{
  const size_t grpid0 = get_group_id(0);

  if(grpid0 >= num_nf_writing_clusters)
    return;
  else
  {
    const size_t lid0          = get_local_id(0);

    const uint   row_this_group = nf_writings_this_device[grpid0];
    const uint   num_h2_leafs   = num_nf_h2_leafs_per_cluster[row_this_group];
    const real   kernel_factor  = r_four_pi;

    gcidxinfo gcii;
    geom      sur;
    singquadg sq;

    init_singquadg(&sq, dim, nq, xqs, yqs, wqs, bases, 0, 0, 0);
    init_geom(&sur, dim, n, vs, p, g, 0, 0, 0);

    init_row_gcidxinfo(&gcii,
                       num_h2_leafs,
                       nf_idx_offs[row_this_group],
                       nf_ridx_sizes[row_this_group],
                       nf_ridx_offs[row_this_group],
                       nf_ridxs,
                       nf_yt_offs[row_this_group],
                       yt,
                       0);

    if(lid0 >= gcii.num_h2_leafs)
      return;

    /* Different nearfield H2-matrix information for all threads in a work
     * group. */
    set_column_info_gcidxinfo(&gcii,
                              nf_cidx_sizes[gcii.idx_off + lid0],
                              nf_cidx_offs[gcii.idx_off + lid0],
                              nf_cidxs,
                              nf_xt_offs[gcii.idx_off + lid0],
                              xts);

    mvm_on_the_fly_gca(&gcii, &sur, &sq, bem_alpha, kernel_factor, alpha);
  }
}

kernel void
fastaddeval_nf_min_vert(         const uint dim,
                                 const uint n,
                        global   const real *vs,
                        global   const uint *p,
                        global   const real *g,
                                 const uint nq_min_vert,
                        constant       real *xqs_min_vert,
                        constant       real *yqs_min_vert,
                        constant       real *wqs_min_vert,
                        global   const real *bases_min_vert,
                                 const uint num_nf_writing_clusters,
                        global   const uint *nf_writings_this_device,
                        global   const uint *num_nf_h2_leafs_per_cluster,
                        global   const uint *nf_idx_offs,
                        global   const uint *nf_ridx_sizes,
                        global   const uint *nf_cidx_sizes,
                        global   const uint *nf_ridx_offs,
                        global   const uint *nf_cidx_offs,
                        global   const uint *nf_ridxs,
                        global   const uint *nf_cidxs,
                        global   const uint *nf_xt_offs,
                        global   const uint *nf_yt_offs,
                        global   const uint *num_min_vert,
                        global   const uint *idx_off_min_vert,
                        global   const uint *rows_min_vert,
                        global   const uint *cols_min_vert,
                        global   const uint *cidx_min_vert,
                                 const real bem_alpha,
                                 const real alpha,
                        global   const real *xt,
                        global         real *yt)
{
  const size_t grpid0 = get_group_id(0);

  if(grpid0 >= num_nf_writing_clusters)
    return;
  else
  {
    const real   kernel_factor = r_four_pi;

    uint   row_this_group      = nf_writings_this_device[grpid0];

    uint   num_h2_leafs        = num_nf_h2_leafs_per_cluster[row_this_group];

    gcidxinfo  gcii;
    geom       sur;
    intinfo    iinfo;
    singquadc  sq;

    local real ytl[1];

    init_row_gcidxinfo(&gcii,
                       num_h2_leafs,
                       nf_idx_offs[row_this_group],
                       nf_ridx_sizes[row_this_group],
                       nf_ridx_offs[row_this_group],
                       nf_ridxs,
                       nf_yt_offs[row_this_group],
                       yt,
                       ytl);



    init_intinfo(&iinfo,
                 num_min_vert[row_this_group],
                 idx_off_min_vert[row_this_group],
                 rows_min_vert,
                 cols_min_vert,
                 cidx_min_vert);

    init_singquadc_min_vert(&sq,
                            dim,
                            nq_min_vert,
                            xqs_min_vert,
                            yqs_min_vert,
                            wqs_min_vert,
                            bases_min_vert);

    init_geom(&sur, dim, n, vs, p, g, 0, 0, 0);

    eval_integrals(&gcii,
                   &iinfo,
                   &sur,
                   &sq,
                   bem_alpha,
                   kernel_factor,
                   alpha,
                   xt);
  }
}

kernel void
fastaddeval_nf_min_edge(         const uint dim,
                                 const uint n,
                        global   const real *vs,
                        global   const uint *p,
                        global   const real *g,
                                 const uint nq,
                        constant       real *xqs,
                        constant       real *yqs,
                        constant       real *wqs,
                        global   const real *bases,
                                 const uint num_nf_writing_clusters,
                        global   const uint *nf_writings_this_device,
                        global   const uint *num_nf_h2_leafs_per_cluster,
                        global   const uint *nf_idx_offs,
                        global   const uint *nf_ridx_sizes,
                        global   const uint *nf_cidx_sizes,
                        global   const uint *nf_ridx_offs,
                        global   const uint *nf_cidx_offs,
                        global   const uint *nf_ridxs,
                        global   const uint *nf_cidxs,
                        global   const uint *nf_xt_offs,
                        global   const uint *nf_yt_offs,
                        global   const uint *num_min_vert,
                        global   const uint *idx_off_min_vert,
                        global   const uint *rows_min_vert,
                        global   const uint *cols_min_vert,
                        global   const uint *cidx_min_vert,
                                 const real bem_alpha,
                                 const real alpha,
                        global   const real *xt,
                        global         real *yt)
{
  const size_t grpid0 = get_group_id(0);

  if(grpid0 >= num_nf_writing_clusters)
    return;
  else
  {
    const real   kernel_factor = r_four_pi;

    uint   row_this_group      = nf_writings_this_device[grpid0];

    uint   num_h2_leafs        = num_nf_h2_leafs_per_cluster[row_this_group];

    gcidxinfo  gcii;
    geom       sur;
    intinfo    iinfo;
    singquadc  sq;

    local real ytl[1];

    init_row_gcidxinfo(&gcii,
                       num_h2_leafs,
                       nf_idx_offs[row_this_group],
                       nf_ridx_sizes[row_this_group],
                       nf_ridx_offs[row_this_group],
                       nf_ridxs,
                       nf_yt_offs[row_this_group],
                       yt,
                       ytl);



    init_intinfo(&iinfo,
                 num_min_vert[row_this_group],
                 idx_off_min_vert[row_this_group],
                 rows_min_vert,
                 cols_min_vert,
                 cidx_min_vert);

    init_singquadc_min_edge(&sq,
                            dim,
                            nq,
                            xqs,
                            yqs,
                            wqs,
                            bases);

    init_geom(&sur, dim, n, vs, p, g, 0, 0, 0);

    eval_integrals(&gcii,
                   &iinfo,
                   &sur,
                   &sq,
                   bem_alpha,
                   kernel_factor,
                   alpha,
                   xt);
  }
}

kernel void
fastaddeval_nf_uncommon(         const uint dim,
                                 const uint n,
                        global   const real *vs,
                        global   const uint *p,
                        global   const real *g,
                        global   const uint *nq,
                        constant       real *xqs,
                        constant       real *yqs,
                        constant       real *wqs,
                        global   const real *bases,
                                 const uint num_nf_writing_clusters,
                        global   const uint *nf_writings_this_device,
                        global   const uint *num_nf_h2_leafs_per_cluster,
                        global   const uint *nf_idx_offs,
                        global   const uint *nf_ridx_sizes,
                        global   const uint *nf_cidx_sizes,
                        global   const uint *nf_ridx_offs,
                        global   const uint *nf_cidx_offs,
                        global   const uint *nf_ridxs,
                        global   const uint *nf_cidxs,
                        global   const uint *nf_xt_offs,
                        global   const uint *nf_yt_offs,
                        global   const uint *num_min_vert,
                        global   const uint *idx_off_min_vert,
                        global   const uint *rows_min_vert,
                        global   const uint *cols_min_vert,
                        global   const uint *cidx_min_vert,
                                 const real bem_alpha,
                                 const real alpha,
                        global   const real *xt,
                        global         real *yt)
{
  const size_t grpid0 = get_group_id(0);

  if(grpid0 >= num_nf_writing_clusters)
    return;
  else
  {
    const real   kernel_factor = r_four_pi;

    uint   row_this_group      = nf_writings_this_device[grpid0];

    uint   num_h2_leafs        = num_nf_h2_leafs_per_cluster[row_this_group];

    gcidxinfo  gcii;
    geom       sur;
    intinfo    iinfo;
    singquadc  sq;

    local real ytl[1];

    init_row_gcidxinfo(&gcii,
                       num_h2_leafs,
                       nf_idx_offs[row_this_group],
                       nf_ridx_sizes[row_this_group],
                       nf_ridx_offs[row_this_group],
                       nf_ridxs,
                       nf_yt_offs[row_this_group],
                       yt,
                       ytl);

    init_intinfo(&iinfo,
                 num_min_vert[row_this_group],
                 idx_off_min_vert[row_this_group],
                 rows_min_vert,
                 cols_min_vert,
                 cidx_min_vert);

    init_singquadc_uncommon(&sq,
                            dim,
                            nq,
                            xqs,
                            yqs,
                            wqs,
                            bases);

    init_geom(&sur, dim, n, vs, p, g, 0, 0, 0);

    eval_integrals(&gcii,
                   &iinfo,
                   &sur,
                   &sq,
                   bem_alpha,
                   kernel_factor,
                   alpha,
                   xt);
  }
}

#endif // GREENCROSS_CL
