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

#include "laplace3d.h"

#include "laplacebem3d.h"
#include "ocl_system.h"

/*
 * @brief Substructure used for approximating @ref _hmatrix "h-", @ref
 * _uniformhmatrix "uniformh-" and @ref _h2matrix "h2matrices".
 *
 * This struct contains many needed parameters for the different approximation
 * techniques such as interpolation or green based methods.
 */

pgreencross
new_greencross_laplace3d(psurface3d gr, uint res, uint q, uint m, real accur)
{
  pclustergeometry cg;
  pgreencross      gc;

  cl_int           err = 0;

  /* Allocate and init object. */

  gc = (pgreencross) allocmem(sizeof(greencross));

  init_greencross(gc, 3);

  /* Set and get cluster from geometry. */

  gc->geom = (void *) gr;

  for(uint i = 0; i < ocl_system.num_devices; ++i)
  {
    gc->buf_x[i] = clCreateBuffer(ocl_system.contexts[i],
                                  CL_MEM_READ_ONLY,
                                  ocl_system.max_package_size,
                                  NULL,
                                  &err);

    CL_CHECK(err);

    CL_CHECK(clEnqueueWriteBuffer
               (ocl_system.queues[i * ocl_system.queues_per_device],
                gc->buf_x[i],
                false,
                0,
                3 * gr->vertices * sizeof(real),
                gr->x,
                0,
                NULL,
                NULL));

    gc->buf_p[i] = clCreateBuffer(ocl_system.contexts[i],
                                  CL_MEM_READ_ONLY,
                                  3 * gr->triangles * sizeof(uint),
                                  NULL,
                                  &err);

    CL_CHECK(err);

    CL_CHECK(clEnqueueWriteBuffer
               (ocl_system.queues[i * ocl_system.queues_per_device],
                gc->buf_p[i],
                false,
                0,
                3 * gr->triangles * sizeof(uint),
                gr->t,
                0,
                NULL,
                &gc->event_p[i]));

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           gr->triangles,
                           sizeof(real),
                           gr->g,
                           NULL,
                           &gc->buf_g[i]);
  }

  gc->n  = gr->triangles;

  gc->m  = m;

  gc->K  = 4 * m;

  gc->bem = (void*) new_slp_laplace_bem3d(gr,
                                          q,
                                          q + 2,
                                          BASIS_CONSTANT_BEM3D,
                                          BASIS_CONSTANT_BEM3D);

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

  for(uint i = 0; i < ocl_system.num_devices; ++i)
  {
    gc->buf_qx[i] = clCreateBuffer(ocl_system.contexts[i],
                                   CL_MEM_READ_ONLY,
                                   2 * ((pbem3d) gc->bem)->sq->n_dist * sizeof(real),
                                   NULL,
                                   &err);

    CL_CHECK(err);

    CL_CHECK(clEnqueueWriteBuffer
               (ocl_system.queues[i * ocl_system.queues_per_device],
                gc->buf_qx[i],
                false,
                0,
                2 * ((pbem3d) gc->bem)->sq->n_dist * sizeof(real),
                ((pbem3d) gc->bem)->sq->x_dist,
                0,
                NULL,
                NULL));

    gc->buf_qy[i] = clCreateBuffer(ocl_system.contexts[i],
                                   CL_MEM_READ_ONLY,
                                   2 * ((pbem3d) gc->bem)->sq->n_dist * sizeof(real),
                                   NULL,
                                   &err);

    CL_CHECK(clEnqueueWriteBuffer
               (ocl_system.queues[i * ocl_system.queues_per_device],
                gc->buf_qy[i],
                false,
                0,
                2 * ((pbem3d) gc->bem)->sq->n_dist * sizeof(real),
                ((pbem3d) gc->bem)->sq->y_dist,
                0,
                NULL,
                NULL));

    gc->buf_w[i]  = clCreateBuffer(ocl_system.contexts[i],
                                   CL_MEM_READ_ONLY,
                                   ((pbem3d) gc->bem)->sq->n_dist * sizeof(real),
                                   NULL,
                                   &err);

    CL_CHECK(err);

    CL_CHECK(clEnqueueWriteBuffer
               (ocl_system.queues[i * ocl_system.queues_per_device],
                gc->buf_w[i],
                false,
                0,
                ((pbem3d) gc->bem)->sq->n_dist * sizeof(real),
                ((pbem3d) gc->bem)->sq->w_dist + 9 * ((pbem3d) gc->bem)->sq->n_dist,
                0,
                NULL,
                NULL));
  }

  /* Clean up */

  del_clustergeometry(cg);

  return gc;
}
