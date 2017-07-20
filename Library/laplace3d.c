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

pgreencross
new_greencross_laplace3d(psurface3d gr, uint res, uint q, uint m)
{
  pclustergeometry cg;
  pgreencross      gc;

  /* Allocate and init object. */

  gc = (pgreencross) allocmem(sizeof(greencross));

  init_greencross(gc, 3);

  /* Set and get cluster from geometry. */

  gc->geom = (void *) gr;

  gc->n  = gr->triangles;

  gc->m  = m;

  gc->K  = 4 * m;

  cg     = build_clustergeometry_greencross((const void *) gr, 3, &gc->idx);

  gc->rc = build_adaptive_cluster(cg, gc->n, gc->idx, res);
  gc->cc = build_adaptive_cluster(cg, gc->n, gc->idx, res);

  gc->rb = build_from_cluster_clusterbasis(gc->rc);
  gc->cb = build_from_cluster_clusterbasis(gc->cc);

  gc->bem = new_slp_laplace_bem3d(gr,
                                  q,
                                  q + 2,
                                  BASIS_CONSTANT_BEM3D,
                                  BASIS_CONSTANT_BEM3D);

  setup_h2matrix_aprx_greenhybrid_bem3d((pbem3d) gc->bem,
                                        gc->rb,
                                        gc->cb,
                                        NULL,
                                        gc->m,
                                        1,
                                        0.5,
                                        10e-07,
                                        build_bem3d_cube_quadpoints);

  /* Clean up */

  del_clustergeometry(cg);

  return gc;
}
