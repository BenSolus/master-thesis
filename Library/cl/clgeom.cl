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

#ifdef USE_FLOAT
typedef float2 real2;
typedef float3 real3;
#else
typedef double3 real2;
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

   local        void *v_tmp;

   global const uint *p;

   global const real *g;

   local        real *g_tmp;
};

void
init_geom(       pgeom      sur,
                 const uint dim,
                 const uint n,
          global const real *v,
          global const uint *p,
          global const real *g,
          local        void *v_tmp,
          local        real *g_tmp)
{
  sur->dim   = dim;
  sur->n     = n;
  sur->g     = g;
  sur->g_tmp = g_tmp;
  sur->v     = v;
  sur->v_tmp = v_tmp;
  sur->p     = p;
}
