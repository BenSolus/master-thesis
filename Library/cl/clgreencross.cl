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
                                        const uint nq_min_vert,
                               constant       real *xqs_min_vert,
                               constant       real *yqs_min_vert,
                               constant       real *wqs_min_vert,
                               global   const real *bases_min_vert,
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
                               global   const real *xts,
                               global         real *yt)
{
  const size_t grpid0 = get_group_id(0);

  if((grpid0 >= num_row_leafs) && (grpid0 >= num_nf_writing_clusters))
    return;
  else
  {
    const size_t lid0          = get_local_id(0);

    const real   kernel_factor = r_four_pi;

    uint   row_this_group      = rows_this_device[grpid0];

    uint   num_h2_leafs        = num_h2_leafs_per_cluster[row_this_group];

    gcidxinfo gcii;
    geom      sur;
    singquadg sq;

    // local real3 vl[SIZE][3];
    // local uint3 pl[SIZE];
    // local real  gl[SIZE];

    // local real xql[8 * QUADRATUR_ORDER];
    // local real yql[8 * QUADRATUR_ORDER];
    // local real wql[4 * QUADRATUR_ORDER];

    local real ytl[SIZE + 1];

    // event_t events[3];
    //
    // events[0] = async_work_group_copy(xql, xqs, 8 * QUADRATUR_ORDER, 0);
    // events[1] = async_work_group_copy(yql, yqs, 8 * QUADRATUR_ORDER, 0);
    // events[2] = async_work_group_copy(wql, wqs, 4 * QUADRATUR_ORDER, 0);

    init_geom(&sur, dim, n, vs, p, g, 0, 0, 0);
    init_singquadg(&sq, dim, nq, xqs, yqs, wqs, bases, 0, 0, 0);

    if(grpid0 < num_row_leafs)
    {
      /* Same H2-matrix farfield information for all threads in a work group */
      init_row_gcidxinfo(&gcii,
                         num_h2_leafs,
                         idx_offs[row_this_group],
                         ridx_sizes[row_this_group],
                         ridx_offs[row_this_group],
                         ridxs,
                         yt_offs[row_this_group],
                         yt,
                         0);

      /* Different farfield H2-matrix information for all threads in a work group */
      if(lid0 < gcii.num_h2_leafs)
        set_column_info_gcidxinfo(&gcii,
                                  cidx_sizes[gcii.idx_off + lid0],
                                  cidx_offs[gcii.idx_off + lid0],
                                  cidxs,
                                  xt_offs[gcii.idx_off + lid0],
                                  xts);

      fastaddeval_farfield(&gcii, &sur, &sq, bem_alpha, kernel_factor, alpha);
      //mvm_on_the_fly_gca(&gcii, &sur, &sq, bem_alpha, kernel_factor, alpha);
    }

    barrier(CLK_GLOBAL_MEM_FENCE);

    /* Comput the first quadrature points for all nearfield matrices, e.g.,
     * distant polygons have the smallest number of quadrature points (say "n")
     * thus we compute the first n quadrature points for all entries, e.g. those
     * with 0 to 3 identical vertices in the row and column polygon. */

    if(grpid0 >= num_nf_writing_clusters)
      return;

    row_this_group = nf_writings_this_device[grpid0];

    num_h2_leafs   = num_nf_h2_leafs_per_cluster[row_this_group];

    /* Same H2-matrix nearfield information for all threads in a work group. */
    init_row_gcidxinfo(&gcii,
                       num_h2_leafs,
                       nf_idx_offs[row_this_group],
                       nf_ridx_sizes[row_this_group],
                       nf_ridx_offs[row_this_group],
                       nf_ridxs,
                       nf_yt_offs[row_this_group],
                       yt,
                       ytl);

    if(lid0 < gcii.num_h2_leafs)
    {
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

    barrier(CLK_GLOBAL_MEM_FENCE);

    intinfo   iinfo;
    singquadc sq_min_vert;

    init_intinfo(&iinfo,
                 num_min_vert[row_this_group],
                 idx_off_min_vert[row_this_group],
                 rows_min_vert,
                 cols_min_vert,
                 cidx_min_vert);

    init_singquadc_min_vert(&sq_min_vert,
                            dim,
                            nq_min_vert,
                            xqs_min_vert,
                            yqs_min_vert,
                            wqs_min_vert,
                            bases_min_vert);

    // eval_integrals(&gcii,
    //                &iinfo,
    //                &sur,
    //                &sq_min_vert,
    //                bem_alpha,
    //                kernel_factor,
    //                alpha,
    //                xts);
  }
}

kernel void
fastaddeval_nf_common(         const uint dim,
                               const uint n,
                      global   const real *vs,
                      global   const uint *p,
                      global   const real *g)
{

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

    local real ytl[SIZE + 1];

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

    local real ytl[SIZE + 1];

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

#endif // GREENCROSS_CL
