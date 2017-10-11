/* ------------------------------------------------------------
 * This is the file "laplace3d.c" of this master thesis.
 * All rights reserved, Bennet Carstensen 2017
 * ------------------------------------------------------------ */

/**
 * @file      laplace3d.c
 * @author    Bennet Carstensen
 * @date      2017
 * @copyright All rights reserved, Bennet Carstensen 2017
 */

#include "laplace3dgca.h"

#include "laplacebem3d.h"
#include "ocl_system.h"

static inline field
slp_kernel_laplacebem3d(const real * x, const real * y,
                        const real * nx, const real * ny, void *data)
{
  real      dist[3];
  real      norm;

  field     res;

  (void) nx;
  (void) ny;
  (void) data;

  dist[0] = x[0] - y[0];
  dist[1] = x[1] - y[1];
  dist[2] = x[2] - y[2];
  norm = REAL_NORMSQR3(dist[0], dist[1], dist[2]);

  res = REAL_RSQRT(norm);

  return res;
}

static inline void
#ifdef __AVX__
#ifdef USE_FLOAT
slp_kernel_simd_laplace3dgca(const __m256 *x, const __m256 *y,
                             const __m256 *nx, const __m256 * ny, void *data,
                             __m256 *res_re, __m256 *res_im)
{
  (void) nx;
  (void) ny;
  (void) data;
  (void) res_im;

  const __m256 dist[3] = { _mm256_sub_ps(x[0], y[0]),
                           _mm256_sub_ps(x[1], y[1]),
                           _mm256_sub_ps(x[2], y[2])};

  *res_re = _mm256_mul_ps(dist[0], dist[0]);
  *res_re = _mm256_add_ps(_mm256_mul_ps(dist[1], dist[1]), *res_re);
  *res_re = _mm256_add_ps(_mm256_mul_ps(dist[2], dist[2]), *res_re);

  __m256 s = _mm256_mul_ps(*res_re, _mm256_set1_ps(0.5f));
  *res_re  = _mm256_rsqrt_ps(*res_re);

  s = _mm256_mul_ps(s, *res_re);

  const __m256 t = _mm256_mul_ps(*res_re, *res_re);

  *res_re = _mm256_mul_ps(_mm256_set1_ps(1.5f), *res_re);

  *res_re = _mm256_sub_ps(_mm256_mul_ps(s, t), *res_re);
}
#endif // USE_FLOAT
#endif // __AVX__

pgreencross
new_greencross_laplace3d(psurface3d gr, uint res, uint q, uint m, real accur)
{
  pclustergeometry cg;
  pgreencross      gc;

  /* Allocate and init object. */

  gc = (pgreencross) allocmem(sizeof(greencross));

  init_greencross(gc, 3);

  gc->geom = (void *) gr;

  /* Upload geometry data to OpenCL devices. */

  for(uint i = 0; i < ocl_system.num_devices; ++i)
  {
    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           3 * gr->vertices,
                           sizeof(real),
                           gr->x,
                           NULL,
                           &gc->buf_x[i]);

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           3 * gr->triangles,
                           sizeof(uint),
                           gr->t,
                           NULL,
                           &gc->buf_p[i]);

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           gr->triangles,
                           sizeof(real),
                           gr->g,
                           NULL,
                           &gc->buf_g[i]);

  }

  gc->n         = gr->triangles;

  gc->m         = m;

  gc->K         = 4 * m;

  gc->aca_accur = accur;

  gc->bem       = (void*) new_slp_laplace_bem3d(gr,
                                                q,
                                                q + 2,
                                                BASIS_CONSTANT_BEM3D,
                                                BASIS_CONSTANT_BEM3D);

  gc->sq_gca    = build_from_singquad2d(((pbem3d) gc->bem)->sq);

  gc->sq_uncommon = build_uncommon_from_singquad2d(((pbem3d) gc->bem)->sq);

  gc->sq_partial_min_vert =
    build_min_vert_from_singquad2d(((pbem3d) gc->bem)->sq);

  gc->sq_partial_min_edge =
    build_min_edge_from_singquad2d(((pbem3d) gc->bem)->sq);

  gc->kernel_3d = slp_kernel_laplacebem3d;

#ifdef __AVX__
  gc->kernel_simd_3d = slp_kernel_simd_laplace3dgca;
#endif

  /* Set and get cluster(basis) from geometry. */

  cg     = build_bem3d_const_clustergeometry((pbem3d) gc->bem, &gc->idx);

  gc->rc = build_adaptive_cluster(cg, gc->n, gc->idx, res);
  gc->cc = build_adaptive_cluster(cg, gc->n, gc->idx, res);

  gc->rb = build_from_cluster_clusterbasis(gc->rc);
  gc->cb = build_from_cluster_clusterbasis(gc->cc);

  setup_h2matrix_aprx_greenhybrid_bem3d((pbem3d) gc->bem,
                                        gc->rb,
                                        gc->cb,
                                        NULL,
                                        gc->m,
                                        1,
                                        1.0,
                                        accur,
                                        build_bem3d_cube_quadpoints);

  /* Clean up */

  del_clustergeometry(cg);

  return gc;
}
