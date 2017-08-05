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

#ifndef CLGEOM_CL
#define CLGEOM_CL

/** @addtogroup geom
 *  @{ */

/*  @brief The source code to performe operations on the module
 *         @ref geom via OpenCL as a string. */
static const char clgeom_src[] =
{
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
  "\n"
  "   global const uint *p;\n"
  "\n"
  "   global const real *g;\n"
  "};\n"
  "\n"
  "pgeom\n"
  "new_geom(       const uint dim,\n"
  "                const uint n,\n"
  "         global const real *v,\n"
  "         global const uint *p,\n"
  "         global const real *g)\n"
  "{\n"
  "  geom gr;\n"
  "\n"
  "  gr.dim = dim;\n"
  "  gr.n   = n;\n"
  "  gr.v   = v;\n"
  "  gr.p   = p;\n"
  "  gr.g   = g;\n"
  "\n"
  "  return &gr;\n"
  "}\n"
};

/** @} */

#endif // CLGEOM_CL