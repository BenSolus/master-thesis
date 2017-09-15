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

#ifndef CLGEOM_CL
#define CLGEOM_CL

/** @addtogroup geom
 *  @{ */

/*  @brief The source code to performe operations on the module
 *         @ref geom via OpenCL as a string. */
static const char clgeom_src[] =
{
  "\n"
  "#ifndef GEOM_CL\n"
  "#define GEOM_CL\n"
  "\n"
  "#ifdef USE_FLOAT\n"
  "typedef float2 real2;\n"
  "typedef float3 real3;\n"
  "#else\n"
  "typedef double2 real2;\n"
  "typedef double3 real3;\n"
  "#endif\n"
  "\n"
  "typedef struct _geom geom;\n"
  "\n"
  "typedef geom *pgeom;\n"
  "\n"
  "typedef const geom *pcgeom;\n"
  "\n"
  "struct _geom\n"
  "{\n"
  "          uint dim;\n"
  "\n"
  "          uint n;\n"
  "\n"
  "   global const real *v;\n"
  "   global const uint *p;\n"
  "   global const real *g;\n"
  "\n"
  "   local        real3 (*vl)[3];\n"
  "   local        uint3 *pl;\n"
  "   local        real  *gl;\n"
  "};\n"
  "\n"
  "void\n"
  "init_geom(       pgeom       sur,\n"
  "                 const uint  dim,\n"
  "                 const uint  n,\n"
  "          global const real  *v,\n"
  "          global const uint  *p,\n"
  "          global const real  *g,\n"
  "          local        real3 (*vl)[3],\n"
  "          local        uint3 *pl,\n"
  "          local        real  *gl)\n"
  "{\n"
  "  sur->dim = dim;\n"
  "  sur->n   = n;\n"
  "\n"
  "  sur->v   = v;\n"
  "  sur->p   = p;\n"
  "  sur->g   = g;\n"
  "\n"
  "  sur->vl  = vl;\n"
  "  sur->pl  = pl;\n"
  "  sur->gl  = gl;\n"
  "}\n"
  "\n"
  "#endif // GEOM_CL\n"
};

/** @} */

#endif // CLGEOM_CL