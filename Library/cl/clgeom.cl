/* ------------------------------------------------------------
 * This is the file "clgeom.cl" of this master thesis.
 * All rights reserved, Bennet Carstensen 2017
 * ------------------------------------------------------------ */

/**
 * @file      cl/clgeom.cl
 * @author    Bennet Carstensen
 * @date      2017
 * @copyright All rights reserved, Bennet Carstensen 2017
 */

#include "clbasic.cl"

#ifndef GEOM_CL
#define GEOM_CL

#ifdef USE_FLOAT
typedef float2 real2;
typedef float3 real3;
#else
typedef double2 real2;
typedef double3 real3;
#endif

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

   local        real3 (*vl)[3];
   local        uint3 *pl;
   local        real  *gl;
};

void
init_geom(       pgeom       sur,
                 const uint  dim,
                 const uint  n,
          global const real  *v,
          global const uint  *p,
          global const real  *g,
          local        real3 (*vl)[3],
          local        uint3 *pl,
          local        real  *gl)
{
  sur->dim = dim;
  sur->n   = n;

  sur->v   = v;
  sur->p   = p;
  sur->g   = g;

  sur->vl  = vl;
  sur->pl  = pl;
  sur->gl  = gl;
}

#endif // GEOM_CL
