/* ------------------------------------------------------------
 This is the file "clavector.cl" of this master thesis.
 All rights reserved, Steffen Boerm 2009
 ------------------------------------------------------------ */

 /**
  * @file      cl/clavector.cl
  * @author    Bennet Carstensen
  * @date      2017
  * @copyright All rights reserved, Steffen Boerm 2009
  */

#ifndef CLAVECTOR_CL
#define CLAVECTOR_CL

/** @addtogroup avector
 *  @{ */

/*  @brief The source code to performe operations on the module
 *         @ref avector via OpenCL as a string. */
static const char clavector_src[] =
{
  "#ifndef CLAVECTOR_H\n"
  "#define CLAVECTOR_H\n"
  "\n"
  "\n"
  "typedef struct _avector avector;\n"
  "\n"
  "typedef avector *pavector;\n"
  "\n"
  "typedef const avector *pcavector;\n"
  "\n"
  "\n"
  "struct _avector {\n"
  "  /** @brief Vector coefficients. */\n"
  "  global field *v;\n"
  "\n"
  "  /** @brief Vector dimension. */\n"
  "  size_t dim;\n"
  "};\n"
  "\n"
  "pavector\n"
  "new_avector(global field *v, size_t dim)\n"
  "{\n"
  "  avector a;\n"
  "\n"
  "  a.v   = v;\n"
  "  a.dim = dim;\n"
  "\n"
  "  return &a;\n"
  "}\n"
  "\n"
  "void\n"
  "print_avector(pcavector v)\n"
  "{\n"
  "  const size_t dim = v->dim;\n"
  "\n"
  "  if((get_local_id(0) == 0) && (get_local_id(1) == 0) && (get_local_id(2) == 0))\n"
  "  {\n"
  "    (void) printf(\"Group (%u, %u, %u): avector(%u)\\n\", get_group_id(0),\n"
  "                                                       get_group_id(1),\n"
  "                                                       get_group_id(2),\n"
  "                                                       dim);\n"
  "    if(dim == 0)\n"
  "      return;\n"
  "\n"
  "    (void) printf(\"  (%.5e\", v->v[0]);\n"
  "    for(size_t i = 1; i < dim; i++)\n"
  "      (void) printf(\" %.5e\", v->v[i]);\n"
  "    (void) printf(\")\\n\");\n"
  "  }\n"
  "}\n"
  "\n"
  "\n"
  "#endif // CLAVECTOR_H\n"
};

/** @} */

#endif // CLAVECTOR_CL