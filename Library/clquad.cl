/* ------------------------------------------------------------
 * This is the file "clquad.cl" of this master thesis.
 * All rights reserved, Steffen Boerm 2009
 * ------------------------------------------------------------ */

/**
 * @file      cl/clquad.cl
 * @author    Bennet Carstensen
 * @date      2017
 * @copyright All rights reserved, Bennet Carstensen 2017
 */

#ifndef CLQUAD_CL
#define CLQUAD_CL

/** @addtogroup quad
 *  @{ */

/*  @brief The source code to performe operations on the module
 *         @ref quad via OpenCL as a string. */
static const char clquad_src[] =
{
  "\n"
  "typedef struct _quad quad;\n"
  "\n"
  "typedef quad *pquad;\n"
  "\n"
  "typedef const quad *pcquad;\n"
  "\n"
  "struct _quad\n"
  "{\n"
  "           uint nq;\n"
  "\n"
  "  constant real *qx;\n"
  "\n"
  "  constant real *qy;\n"
  "\n"
  "  constant real *w;\n"
  "};\n"
  "\n"
  "void\n"
  "init_quad(         pquad      q,\n"
  "                   const uint nq,\n"
  "          constant       real *qx,\n"
  "          constant       real *qy,\n"
  "          constant       real *w)\n"
  "{\n"
  "  q->nq = nq;\n"
  "  q->qx = qx;\n"
  "  q->qy = qy;\n"
  "  q->w  = w;\n"
  "}\n"
};

/** @} */

#endif // CLQUAD_CL