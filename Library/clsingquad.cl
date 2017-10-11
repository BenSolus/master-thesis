/* ------------------------------------------------------------
 * This is the file "clsingquad.cl" of this master thesis.
 * All rights reserved, Bennet Carstensen 2017
 * ------------------------------------------------------------ */

/**
 * @file      cl/clsingquad.cl
 * @author    Bennet Carstensen
 * @date      2017
 * @copyright All rights reserved, Bennet Carstensen 2017
 */

#ifndef CLSINGQUAD_CL
#define CLSINGQUAD_CL

/** @addtogroup singquad
 *  @{ */

/*  @brief The source code to performe operations on the module
 *         @ref singquad via OpenCL as a string. */
static const char clsingquad_src[] =
{
  "\n"
  "#ifndef SINGQUAD_CL\n"
  "#define SINGQUAD_CL\n"
  "\n"
  "#ifdef USE_FLOAT\n"
  "typedef float4  real4;\n"
  "#else\n"
  "typedef double4  real4;\n"
  "#endif\n"
  "\n"
  "typedef struct _singquadl singquadl;\n"
  "\n"
  "typedef singquadl *psingquadl;\n"
  "\n"
  "typedef const singquadl *pcsingquadl;\n"
  "\n"
  "typedef struct _singquadc singquadc;\n"
  "\n"
  "typedef singquadc *psingquadc;\n"
  "\n"
  "typedef const singquadc *pcsingquadc;\n"
  "\n"
  "typedef struct _singquadg singquadg;\n"
  "\n"
  "typedef singquadg *psingquadg;\n"
  "\n"
  "typedef const singquadg *pcsingquadg;\n"
  "\n"
  "struct _singquadl\n"
  "{\n"
  "  uint nq;\n"
  "\n"
  "  local const real *xqs;\n"
  "  local const real *yqs;\n"
  "  local const real *wqs;\n"
  "\n"
  "  real4 bases;\n"
  "};\n"
  "\n"
  "struct _singquadc\n"
  "{\n"
  "  uint4    nq;\n"
  "\n"
  "  constant real *xqs;\n"
  "  constant real *yqs;\n"
  "  constant real *wqs;\n"
  "\n"
  "  real4 bases;\n"
  "\n"
  "  int   offset;\n"
  "};\n"
  "\n"
  "struct _singquadg\n"
  "{\n"
  "  uint nq;\n"
  "\n"
  "  global const real *xqs;\n"
  "  global const real *yqs;\n"
  "  global const real *wqs;\n"
  "\n"
  "  real4 bases;\n"
  "\n"
  "  int   offset;\n"
  "\n"
  "  local real *xqs_buf;\n"
  "  local real *yqs_buf;\n"
  "  local real *wqs_buf;\n"
  "\n"
  "  event_t events[10];\n"
  "};\n"
  "\n"
  "void\n"
  "init_singquadl(             psingquadl sq,\n"
  "                      const uint       dim,\n"
  "                      const uint       nq,\n"
  "               local  const real       *xqs,\n"
  "               local  const real       *yqs,\n"
  "               local  const real       *wqs,\n"
  "               global const real       *bases)\n"
  "{\n"
  "  sq->nq    = nq;\n"
  "\n"
  "  sq->xqs   = xqs;\n"
  "  sq->yqs   = yqs;\n"
  "  sq->wqs   = wqs;\n"
  "\n"
  "  sq->bases = vload4(0, bases);\n"
  "}\n"
  "\n"
  "void\n"
  "init_singquadc_common(               psingquadc sq,\n"
  "                               const uint       dim,\n"
  "                               const uint       nq,\n"
  "                      constant       real       *xqs,\n"
  "                      constant       real       *yqs,\n"
  "                      constant       real       *wqs,\n"
  "                      global   const real       *bases)\n"
  "{\n"
  "  sq->nq    = nq;\n"
  "\n"
  "  sq->xqs   = xqs;\n"
  "  sq->yqs   = yqs;\n"
  "  sq->wqs   = wqs;\n"
  "\n"
  "  sq->bases = vload4(0, bases);\n"
  "\n"
  "  sq->offset = 0;\n"
  "}\n"
  "\n"
  "void\n"
  "init_singquadc_uncommon(               psingquadc sq,\n"
  "                                 const uint       dim,\n"
  "                        global   const uint       *nq,\n"
  "                        constant       real       *xqs,\n"
  "                        constant       real       *yqs,\n"
  "                        constant       real       *wqs,\n"
  "                        global   const real       *bases)\n"
  "{\n"
  "  sq->nq.x  =  nq[0];\n"
  "  sq->nq.y  =  nq[1];\n"
  "  sq->nq.z  =  nq[2];\n"
  "  sq->nq.w  =  0;\n"
  "\n"
  "  sq->xqs   = xqs;\n"
  "  sq->yqs   = yqs;\n"
  "  sq->wqs   = wqs;\n"
  "\n"
  "  sq->bases = vload4(0, bases);\n"
  "\n"
  "  sq->offset = -1;\n"
  "}\n"
  "\n"
  "void\n"
  "init_singquadc_min_vert(               psingquadc sq,\n"
  "                                 const uint       dim,\n"
  "                                 const uint       nq,\n"
  "                        constant       real       *xqs,\n"
  "                        constant       real       *yqs,\n"
  "                        constant       real       *wqs,\n"
  "                        global   const real       *bases)\n"
  "{\n"
  "  sq->nq    = nq;\n"
  "\n"
  "  sq->xqs   = xqs;\n"
  "  sq->yqs   = yqs;\n"
  "  sq->wqs   = wqs;\n"
  "\n"
  "  sq->bases = vload4(0, bases);\n"
  "\n"
  "  sq->offset = -1;\n"
  "}\n"
  "\n"
  "void\n"
  "init_singquadc_min_edge(               psingquadc sq,\n"
  "                                 const uint       dim,\n"
  "                                 const uint       nq,\n"
  "                        constant       real       *xqs,\n"
  "                        constant       real       *yqs,\n"
  "                        constant       real       *wqs,\n"
  "                        global   const real       *bases)\n"
  "{\n"
  "  sq->nq    = nq;\n"
  "\n"
  "  sq->xqs   = xqs;\n"
  "  sq->yqs   = yqs;\n"
  "  sq->wqs   = wqs;\n"
  "\n"
  "  sq->bases = vload4(0, bases);\n"
  "\n"
  "  sq->offset = -2;\n"
  "}\n"
  "\n"
  "void\n"
  "init_singquadg(             psingquadg sq,\n"
  "                      const uint       dim,\n"
  "                      const uint       nq,\n"
  "               global const real       *xqs,\n"
  "               global const real       *yqs,\n"
  "               global const real       *wqs,\n"
  "               global const real       *bases,\n"
  "               local        real       *xqs_buf,\n"
  "               local        real       *yqs_buf,\n"
  "               local        real       *wqs_buf)\n"
  "{\n"
  "  sq->nq      = nq;\n"
  "\n"
  "  sq->xqs     = xqs;\n"
  "  sq->yqs     = yqs;\n"
  "  sq->wqs     = wqs;\n"
  "\n"
  "  sq->bases   = vload4(0, bases);\n"
  "\n"
  "  sq->xqs_buf = xqs_buf;\n"
  "  sq->yqs_buf = yqs_buf;\n"
  "  sq->wqs_buf = wqs_buf;\n"
  "}\n"
  "\n"
  "void\n"
  "init_singquadg_min_vert(             psingquadg sq,\n"
  "                               const uint       dim,\n"
  "                               const uint       nq,\n"
  "                        global const real       *xqs,\n"
  "                        global const real       *yqs,\n"
  "                        global const real       *wqs,\n"
  "                        global const real       *bases)\n"
  "{\n"
  "  sq->nq    = nq;\n"
  "\n"
  "  sq->xqs   = xqs;\n"
  "  sq->yqs   = yqs;\n"
  "  sq->wqs   = wqs;\n"
  "\n"
  "  sq->bases = vload4(0, bases);\n"
  "\n"
  "  sq->offset = -1;\n"
  "}\n"
  "\n"
  "uint\n"
  "select_quadraturel(              psingquadl sq,\n"
  "                   private const uint3      x,\n"
  "                   private const uint3      y,\n"
  "                   private       uint       *xp,\n"
  "                   private       uint       *yp,\n"
  "                   local   const real       **xq,\n"
  "                   local   const real       **yq,\n"
  "                   local   const real       **wq,\n"
  "                   private       real       *base)\n"
  "{\n"
  "  uint p = (x.x == y.x) + (x.x == y.y) + (x.x == y.z) +\n"
  "           (x.y == y.x) + (x.y == y.y) + (x.y == y.z) +\n"
  "           (x.z == y.x) + (x.z == y.y) + (x.z == y.z);\n"
  "\n"
  "  xp[0] = 0; xp[1] = 1; xp[2] = 2; yp[0] = 0; yp[1] = 1; yp[2] = 2;\n"
  "\n"
  "  *xq   = sq->xqs + (2 * p * sq->nq);\n"
  "  *yq   = sq->yqs + (2 * p * sq->nq);\n"
  "  *wq   = sq->wqs + (p * sq->nq);\n"
  "  *base = ((private real*) &sq->bases)[p];\n"
  "\n"
  "  p = 0;\n"
  "\n"
  "  for(uint i = 0; i < 3; ++i)\n"
  "  {\n"
  "    for(uint j = 0; j < 3; ++j)\n"
  "    {\n"
  "      if(((private const uint*) &x)[i] == ((private const uint*) &y)[j])\n"
  "      {\n"
  "        xp[p] = i;\n"
  "        yp[p] = j;\n"
  "        ++p;\n"
  "        break;\n"
  "      }\n"
  "    }\n"
  "  }\n"
  "\n"
  "  // printf(\"%u (%u %u %u) (%u %u %u)\\n\", p, xp[0], xp[1], xp[2], yp[0], yp[1], yp[2]);\n"
  "\n"
  "  uint q = p;\n"
  "\n"
  "  for(uint i = 0; i < 3; ++i)\n"
  "  {\n"
  "    uint j = 0;\n"
  "\n"
  "    for(j = 0; (j < q) &&\n"
  "               (((private const uint*) &x)[i] !=\n"
  "                ((private const uint*) &x)[xp[j]]); ++j)\n"
  "      ;\n"
  "\n"
  "    if(j == q)\n"
  "      xp[q++] = i;\n"
  "  }\n"
  "\n"
  "  q = p;\n"
  "\n"
  "  for(uint i = 0; i < 3; ++i)\n"
  "  {\n"
  "    uint j = 0;\n"
  "\n"
  "    for(j = 0; (j < q) &&\n"
  "               (((private const uint*) &y)[i] !=\n"
  "                ((private const uint*) &y)[yp[j]]); ++j)\n"
  "      ;\n"
  "\n"
  "    if(j == q)\n"
  "      yp[q++] = i;\n"
  "  }\n"
  "\n"
  "  return p;\n"
  "}\n"
  "\n"
  "uint\n"
  "select_quadraturec(               psingquadc sq,\n"
  "                   private  const uint3      x,\n"
  "                   private  const uint3      y,\n"
  "                   private        uint       *xp,\n"
  "                   private        uint       *yp,\n"
  "                   private        uint       *nq,\n"
  "                   constant       real       **xq,\n"
  "                   constant       real       **yq,\n"
  "                   constant       real       **wq,\n"
  "                   private        real       *base)\n"
  "{\n"
  "  uint p = (x.x == y.x) + (x.x == y.y) + (x.x == y.z) +\n"
  "           (x.y == y.x) + (x.y == y.y) + (x.y == y.z) +\n"
  "           (x.z == y.x) + (x.z == y.y) + (x.z == y.z);\n"
  "\n"
  "  xp[0] = 0; xp[1] = 1; xp[2] = 2; yp[0] = 0; yp[1] = 1; yp[2] = 2;\n"
  "\n"
  "  *nq   = ((private uint*) &sq->nq)[p + sq->offset];\n"
  "  *xq   = sq->xqs + (2 * (p + sq->offset) * *nq);\n"
  "  *yq   = sq->yqs + (2 * (p + sq->offset) * *nq);\n"
  "  *wq   = sq->wqs + (1 * (p + sq->offset) * *nq);\n"
  "  *base = ((private real*) &sq->bases)[p + sq->offset];\n"
  "\n"
  "  p = 0;\n"
  "\n"
  "  for(uint i = 0; i < 3; ++i)\n"
  "  {\n"
  "    for(uint j = 0; j < 3; ++j)\n"
  "    {\n"
  "      if(((private const uint*) &x)[i] == ((private const uint*) &y)[j])\n"
  "      {\n"
  "        xp[p] = i;\n"
  "        yp[p] = j;\n"
  "        ++p;\n"
  "        break;\n"
  "      }\n"
  "    }\n"
  "  }\n"
  "\n"
  "  uint q = p;\n"
  "\n"
  "  for(uint i = 0; i < 3; ++i)\n"
  "  {\n"
  "    uint j = 0;\n"
  "\n"
  "    for(j = 0; (j < q) &&\n"
  "               (((private const uint*) &x)[i] !=\n"
  "                ((private const uint*) &x)[xp[j]]); ++j)\n"
  "      ;\n"
  "\n"
  "    if(j == q)\n"
  "      xp[q++] = i;\n"
  "  }\n"
  "\n"
  "  q = p;\n"
  "\n"
  "  for(uint i = 0; i < 3; ++i)\n"
  "  {\n"
  "    uint j = 0;\n"
  "\n"
  "    for(j = 0; (j < q) &&\n"
  "               (((private const uint*) &y)[i] !=\n"
  "                ((private const uint*) &y)[yp[j]]); ++j)\n"
  "      ;\n"
  "\n"
  "    if(j == q)\n"
  "      yp[q++] = i;\n"
  "  }\n"
  "\n"
  "  return p;\n"
  "}\n"
  "\n"
  "uint\n"
  "select_uncommon_quadraturec(               psingquadc sq,\n"
  "                            private  const uint3      x,\n"
  "                            private  const uint3      y,\n"
  "                            private        uint       *xp,\n"
  "                            private        uint       *yp,\n"
  "                            private        uint       *nq,\n"
  "                            constant       real       **xq,\n"
  "                            constant       real       **yq,\n"
  "                            constant       real       **wq,\n"
  "                            private        real       *base)\n"
  "{\n"
  "  uint p = (x.x == y.x) + (x.x == y.y) + (x.x == y.z) +\n"
  "           (x.y == y.x) + (x.y == y.y) + (x.y == y.z) +\n"
  "           (x.z == y.x) + (x.z == y.y) + (x.z == y.z);\n"
  "\n"
  "  xp[0] = 0; xp[1] = 1; xp[2] = 2; yp[0] = 0; yp[1] = 1; yp[2] = 2;\n"
  "\n"
  "  *nq   = ((private uint*) &sq->nq)[p + sq->offset];\n"
  "  *xq   = sq->xqs;\n"
  "  *yq   = sq->yqs;\n"
  "  *wq   = sq->wqs;\n"
  "  *base = ((private real*) &sq->bases)[p + sq->offset];\n"
  "\n"
  "  for(uint i = 0; i < (p + sq->offset); ++i)\n"
  "  {\n"
  "    *xq += 2 * ((private uint*) &sq->nq)[i];\n"
  "    *yq += 2 * ((private uint*) &sq->nq)[i];\n"
  "    *wq += 1 * ((private uint*) &sq->nq)[i];\n"
  "  }\n"
  "\n"
  "  p = 0;\n"
  "\n"
  "  for(uint i = 0; i < 3; ++i)\n"
  "  {\n"
  "    for(uint j = 0; j < 3; ++j)\n"
  "    {\n"
  "      if(((private const uint*) &x)[i] == ((private const uint*) &y)[j])\n"
  "      {\n"
  "        xp[p] = i;\n"
  "        yp[p] = j;\n"
  "        ++p;\n"
  "        break;\n"
  "      }\n"
  "    }\n"
  "  }\n"
  "\n"
  "  uint q = p;\n"
  "\n"
  "  for(uint i = 0; i < 3; ++i)\n"
  "  {\n"
  "    uint j = 0;\n"
  "\n"
  "    for(j = 0; (j < q) &&\n"
  "               (((private const uint*) &x)[i] !=\n"
  "                ((private const uint*) &x)[xp[j]]); ++j)\n"
  "      ;\n"
  "\n"
  "    if(j == q)\n"
  "      xp[q++] = i;\n"
  "  }\n"
  "\n"
  "  q = p;\n"
  "\n"
  "  for(uint i = 0; i < 3; ++i)\n"
  "  {\n"
  "    uint j = 0;\n"
  "\n"
  "    for(j = 0; (j < q) &&\n"
  "               (((private const uint*) &y)[i] !=\n"
  "                ((private const uint*) &y)[yp[j]]); ++j)\n"
  "      ;\n"
  "\n"
  "    if(j == q)\n"
  "      yp[q++] = i;\n"
  "  }\n"
  "\n"
  "  return p;\n"
  "}\n"
  "\n"
  "uint\n"
  "select_quadratureg(              psingquadg sq,\n"
  "                   private const uint3      x,\n"
  "                   private const uint3      y,\n"
  "                   private       uint       *xp,\n"
  "                   private       uint       *yp,\n"
  "                   global  const real       **xq,\n"
  "                   global  const real       **yq,\n"
  "                   global  const real       **wq,\n"
  "                   private       real       *base)\n"
  "{\n"
  "  uint p = (x.x == y.x) + (x.x == y.y) + (x.x == y.z) +\n"
  "           (x.y == y.x) + (x.y == y.y) + (x.y == y.z) +\n"
  "           (x.z == y.x) + (x.z == y.y) + (x.z == y.z);\n"
  "\n"
  "  xp[0] = 0; xp[1] = 1; xp[2] = 2; yp[0] = 0; yp[1] = 1; yp[2] = 2;\n"
  "\n"
  "  *xq   = sq->xqs + (2 * p * sq->nq);\n"
  "  *yq   = sq->yqs + (2 * p * sq->nq);\n"
  "  *wq   = sq->wqs + (p * sq->nq);\n"
  "  *base = ((private real*) &sq->bases)[p];\n"
  "\n"
  "  p = 0;\n"
  "\n"
  "  for(uint i = 0; i < 3; ++i)\n"
  "  {\n"
  "    for(uint j = 0; j < 3; ++j)\n"
  "    {\n"
  "      if(((private const uint*) &x)[i] == ((private const uint*) &y)[j])\n"
  "      {\n"
  "        xp[p] = i;\n"
  "        yp[p] = j;\n"
  "        ++p;\n"
  "        break;\n"
  "      }\n"
  "    }\n"
  "  }\n"
  "\n"
  "  // printf(\"%u (%u %u %u) (%u %u %u)\\n\", p, xp[0], xp[1], xp[2], yp[0], yp[1], yp[2]);\n"
  "\n"
  "  uint q = p;\n"
  "\n"
  "  for(uint i = 0; i < 3; ++i)\n"
  "  {\n"
  "    uint j = 0;\n"
  "\n"
  "    for(j = 0; (j < q) &&\n"
  "               (((private const uint*) &x)[i] !=\n"
  "                ((private const uint*) &x)[xp[j]]); ++j)\n"
  "      ;\n"
  "\n"
  "    if(j == q)\n"
  "      xp[q++] = i;\n"
  "  }\n"
  "\n"
  "  q = p;\n"
  "\n"
  "  for(uint i = 0; i < 3; ++i)\n"
  "  {\n"
  "    uint j = 0;\n"
  "\n"
  "    for(j = 0; (j < q) &&\n"
  "               (((private const uint*) &y)[i] !=\n"
  "                ((private const uint*) &y)[yp[j]]); ++j)\n"
  "      ;\n"
  "\n"
  "    if(j == q)\n"
  "      yp[q++] = i;\n"
  "  }\n"
  "\n"
  "  return p;\n"
  "}\n"
  "\n"
  "#endif // SINGQUAD_CL\n"
};

/** @} */

#endif // CLSINGQUAD_CL