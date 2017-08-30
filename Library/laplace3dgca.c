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

  /* Upload quadrature data to OpenCL devices. */

  for(uint i = 0; i < ocl_system.num_devices; ++i)
  {
    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           2 * ((pbem3d) gc->bem)->sq->n_dist,
                           sizeof(real),
                           ((pbem3d) gc->bem)->sq->x_dist,
                           NULL,
                           &gc->buf_qx[i]);

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           2 * ((pbem3d) gc->bem)->sq->n_dist,
                           sizeof(real),
                           ((pbem3d) gc->bem)->sq->y_dist,
                           NULL,
                           &gc->buf_qy[i]);

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           ((pbem3d) gc->bem)->sq->n_dist,
                           sizeof(real),
                           ((pbem3d) gc->bem)->sq->w_dist +
                           9 * ((pbem3d) gc->bem)->sq->n_dist,
                           NULL,
                           &gc->buf_w[i]);
  }

  /* Clean up */

  del_clustergeometry(cg);

  return gc;
}
