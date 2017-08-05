/* ------------------------------------------------------------
 * This is the file "clgeom.cl" of this master thesis.
 * All rights reserved, Steffen Boerm 2009
 * ------------------------------------------------------------ */

/**
 * @file      cl/clgeom.cl
 * @author    Bennet Carstensen
 * @date      2017
 * @copyright All rights reserved, Bennet Carstensen 2017
 */

#include "clbasic.cl"

/** @brief @ref quad is just an abbreviation for the struct @ref _quad. */
typedef struct _geom geom;

/** @brief Pointer to a @ref _quad "quad" object. */
typedef geom *pgeom;

/** @brief Pointer to a constant @ref _quad "quad" object. */
typedef const geom *pcgeom;

struct _geom
{
          uint dim;

          uint n;

   global const real *v;

   global const uint *p;

   global const real *g;
};

pgeom
new_geom(       const uint dim,
                const uint n,
         global const real *v,
         global const uint *p,
         global const real *g)
{
  geom gr;

  gr.dim = dim;
  gr.n   = n;
  gr.v   = v;
  gr.p   = p;
  gr.g   = g;

  return &gr;
}
