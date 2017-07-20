/* ------------------------------------------------------------
 * This is the file "laplace2d.c" of this master thesis.
 * All rights reserved, Bennet Carstensen 2017
 * ------------------------------------------------------------ */

/**
 * @file      laplace2d.c
 * @author    Bennet Carstensen
 * @date      2017
 * @copyright All rights reserved, Bennet Carstensen 2017
 */

#include "laplace2d.h"

#include "laplacebem2d.h"
#include "singquad1d.h"

pgreencross
new_greencross_laplace2d(pcurve2d gr, uint res, uint q, uint m)
{
  pclustergeometry cg;
  pgreencross      gc;

  /* Allocate and init object. */

  gc = (pgreencross) allocmem(sizeof(greencross));

  init_greencross(gc, 2);

  /* Set and get cluster from geometry. */

  gc->geom = (void *) gr;

  gc->n  = gr->edges;

  gc->m  = m;

  gc->K  = 4 * m;

  cg     = build_clustergeometry_greencross((const void *) gr, 2, &gc->idx);

  gc->rc = build_adaptive_cluster(cg, gc->n, gc->idx, res);
  gc->cc = build_adaptive_cluster(cg, gc->n, gc->idx, res);

  gc->rb = build_from_cluster_clusterbasis(gc->rc);
  gc->cb = build_from_cluster_clusterbasis(gc->cc);

  gc->bem = new_slp_laplace_bem2d(gr, q, BASIS_CONSTANT_BEM2D);

  setup_h2matrix_aprx_greenhybrid_bem2d((pbem2d) gc->bem,
                                        gc->rb,
                                        gc->cb,
                                        NULL,
                                        m,
                                        1,
                                        0.5,
                                        10e-07,
                                        build_bem2d_rect_quadpoints);

  /* Clean up */

  del_clustergeometry(cg);

  return gc;
}
