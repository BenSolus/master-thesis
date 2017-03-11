/* ------------------------------------------------------------
 * This is the file "objects.c" of this master thesis.
 * All rights reserved, Bennet Carstensen 2017
 * ------------------------------------------------------------ */

#include "objects.h"

void
init_objects(pobjects objs)
{
  assert(objs->x   == NULL);
  assert(objs->mem == NULL);

  const uint n   = objs->n;
  const uint dim = objs->dim;

  real **x  = (real **) allocmem(sizeof(real *) * n);
  real *mem = (real *)  allocreal(n * dim);

  for(uint i = 0; i < n; ++i)
  {
    x[i] = mem;
    mem += dim;
  }

  objs->x   = x;
  objs->mem = mem;
}

void
uninit_objects(pobjects objs)
{
  if(objs->x != NULL)
  {
    freemem(objs->x);
    objs->x = NULL;
  }
  if(objs->mem != NULL)
  {
    freemem(objs->mem);
    objs->mem = NULL;
  }
}

pobjects
new_objects(uint n, uint dim)
{
  pobjects objs = (pobjects) allocmem(sizeof(objects));

  objs->x   = NULL;
  objs->n   = n;
  objs->dim = dim;
  objs->mem = NULL;

  init_objects(objs);

  return objs;
}

void
del_objects(pobjects objs)
{
  if(objs != NULL)
  {
    uninit_objects(objs);

    freemem(objs);
    objs = NULL;
  }
}
