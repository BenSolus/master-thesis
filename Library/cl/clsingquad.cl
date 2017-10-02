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

#include "clbasic.cl"

#ifndef SINGQUAD_CL
#define SINGQUAD_CL

#ifdef USE_FLOAT
typedef float4  real4;
#else
typedef double4  real4;
#endif

/** @brief @ref quadl is just an abbreviation for the struct @ref _quadl. */
typedef struct _singquadl singquadl;

/** @brief Pointer to a @ref _quadc "quadl" object. */
typedef singquadl *psingquadl;

/** @brief Pointer to a constant @ref _quadc "quadl" object. */
typedef const singquadl *pcsingquadl;

/** @brief @ref quadc is just an abbreviation for the struct @ref _quadc. */
typedef struct _singquadc singquadc;

/** @brief Pointer to a @ref _quadc "quadc" object. */
typedef singquadc *psingquadc;

/** @brief Pointer to a constant @ref _quadc "quadc" object. */
typedef const singquadc *pcsingquadc;

/** @brief @ref quadg is just an abbreviation for the struct @ref _quadg. */
typedef struct _singquadg singquadg;

/** @brief Pointer to a @ref _quadg "quadg" object. */
typedef singquadg *psingquadg;

/** @brief Pointer to a constant @ref _quadg "quadg" object. */
typedef const singquadg *pcsingquadg;

struct _singquadl
{
  uint nq;

  local const real *xqs;
  local const real *yqs;
  local const real *wqs;

  real4 bases;
};

struct _singquadc
{
  uint nq;

  constant real *xqs;
  constant real *yqs;
  constant real *wqs;

  real4 bases;

  int   offset;
};

struct _singquadg
{
  uint nq;

  global const real *xqs;
  global const real *yqs;
  global const real *wqs;

  real4 bases;

  int   offset;

  local real *xqs_buf;
  local real *yqs_buf;
  local real *wqs_buf;

  event_t events[10];
};

void
init_singquadl(             psingquadl sq,
                      const uint       dim,
                      const uint       nq,
               local  const real       *xqs,
               local  const real       *yqs,
               local  const real       *wqs,
               global const real       *bases)
{
  sq->nq    = nq;

  sq->xqs   = xqs;
  sq->yqs   = yqs;
  sq->wqs   = wqs;

  sq->bases = vload4(0, bases);
}

void
init_singquadc_any(               psingquadc sq,
                            const uint       dim,
                            const uint       nq,
                   constant       real       *xqs,
                   constant       real       *yqs,
                   constant       real       *wqs,
                   global   const real       *bases)
{
  sq->nq    = nq;

  sq->xqs   = xqs;
  sq->yqs   = yqs;
  sq->wqs   = wqs;

  sq->bases = vload4(0, bases);

  sq->offset = 0;
}

void
init_singquadc_min_vert(               psingquadc sq,
                                 const uint       dim,
                                 const uint       nq,
                        constant       real       *xqs,
                        constant       real       *yqs,
                        constant       real       *wqs,
                        global   const real       *bases)
{
  sq->nq    = nq;

  sq->xqs   = xqs;
  sq->yqs   = yqs;
  sq->wqs   = wqs;

  sq->bases = vload4(0, bases);

  sq->offset = -1;
}

void
init_singquadc_min_edge(               psingquadc sq,
                                 const uint       dim,
                                 const uint       nq,
                        constant       real       *xqs,
                        constant       real       *yqs,
                        constant       real       *wqs,
                        global   const real       *bases)
{
  sq->nq    = nq;

  sq->xqs   = xqs;
  sq->yqs   = yqs;
  sq->wqs   = wqs;

  sq->bases = vload4(0, bases);

  sq->offset = -2;
}

void
init_singquadg(             psingquadg sq,
                      const uint       dim,
                      const uint       nq,
               global const real       *xqs,
               global const real       *yqs,
               global const real       *wqs,
               global const real       *bases,
               local        real       *xqs_buf,
               local        real       *yqs_buf,
               local        real       *wqs_buf)
{
  sq->nq      = nq;

  sq->xqs     = xqs;
  sq->yqs     = yqs;
  sq->wqs     = wqs;

  sq->bases   = vload4(0, bases);

  sq->xqs_buf = xqs_buf;
  sq->yqs_buf = yqs_buf;
  sq->wqs_buf = wqs_buf;
}

void
init_singquadg_min_vert(             psingquadg sq,
                               const uint       dim,
                               const uint       nq,
                        global const real       *xqs,
                        global const real       *yqs,
                        global const real       *wqs,
                        global const real       *bases)
{
  sq->nq    = nq;

  sq->xqs   = xqs;
  sq->yqs   = yqs;
  sq->wqs   = wqs;

  sq->bases = vload4(0, bases);

  sq->offset = -1;
}

uint
select_quadraturel(              psingquadl sq,
                   private const uint3      x,
                   private const uint3      y,
                   private       uint       *xp,
                   private       uint       *yp,
                   local   const real       **xq,
                   local   const real       **yq,
                   local   const real       **wq,
                   private       real       *base)
{
  uint p = (x.x == y.x) + (x.x == y.y) + (x.x == y.z) +
           (x.y == y.x) + (x.y == y.y) + (x.y == y.z) +
           (x.z == y.x) + (x.z == y.y) + (x.z == y.z);

  xp[0] = 0; xp[1] = 1; xp[2] = 2; yp[0] = 0; yp[1] = 1; yp[2] = 2;

  *xq   = sq->xqs + (2 * p * sq->nq);
  *yq   = sq->yqs + (2 * p * sq->nq);
  *wq   = sq->wqs + (p * sq->nq);
  *base = ((private real*) &sq->bases)[p];

  p = 0;

  for(uint i = 0; i < 3; ++i)
  {
    for(uint j = 0; j < 3; ++j)
    {
      if(((private const uint*) &x)[i] == ((private const uint*) &y)[j])
      {
        xp[p] = i;
        yp[p] = j;
        ++p;
        break;
      }
    }
  }

  // printf("%u (%u %u %u) (%u %u %u)\n", p, xp[0], xp[1], xp[2], yp[0], yp[1], yp[2]);

  uint q = p;

  for(uint i = 0; i < 3; ++i)
  {
    uint j = 0;

    for(j = 0; (j < q) &&
               (((private const uint*) &x)[i] !=
                ((private const uint*) &x)[xp[j]]); ++j)
      ;

    if(j == q)
      xp[q++] = i;
  }

  q = p;

  for(uint i = 0; i < 3; ++i)
  {
    uint j = 0;

    for(j = 0; (j < q) &&
               (((private const uint*) &y)[i] !=
                ((private const uint*) &y)[yp[j]]); ++j)
      ;

    if(j == q)
      yp[q++] = i;
  }

  return p;
}

uint
select_quadraturec(               psingquadc sq,
                   private  const uint3      x,
                   private  const uint3      y,
                   private        uint       *xp,
                   private        uint       *yp,
                   constant       real       **xq,
                   constant       real       **yq,
                   constant       real       **wq,
                   private        real       *base)
{
  uint p = (x.x == y.x) + (x.x == y.y) + (x.x == y.z) +
           (x.y == y.x) + (x.y == y.y) + (x.y == y.z) +
           (x.z == y.x) + (x.z == y.y) + (x.z == y.z);

  xp[0] = 0; xp[1] = 1; xp[2] = 2; yp[0] = 0; yp[1] = 1; yp[2] = 2;

  *xq   = sq->xqs + (2 * (p + sq->offset) * sq->nq);
  *yq   = sq->yqs + (2 * (p + sq->offset) * sq->nq);
  *wq   = sq->wqs + ((p + sq->offset) * sq->nq);
  *base = ((private real*) &sq->bases)[p + sq->offset];

  p = 0;

  for(uint i = 0; i < 3; ++i)
  {
    for(uint j = 0; j < 3; ++j)
    {
      if(((private const uint*) &x)[i] == ((private const uint*) &y)[j])
      {
        xp[p] = i;
        yp[p] = j;
        ++p;
        break;
      }
    }
  }

  // printf("%u (%u %u %u) (%u %u %u)\n", p, xp[0], xp[1], xp[2], yp[0], yp[1], yp[2]);

  uint q = p;

  for(uint i = 0; i < 3; ++i)
  {
    uint j = 0;

    for(j = 0; (j < q) &&
               (((private const uint*) &x)[i] !=
                ((private const uint*) &x)[xp[j]]); ++j)
      ;

    if(j == q)
      xp[q++] = i;
  }

  q = p;

  for(uint i = 0; i < 3; ++i)
  {
    uint j = 0;

    for(j = 0; (j < q) &&
               (((private const uint*) &y)[i] !=
                ((private const uint*) &y)[yp[j]]); ++j)
      ;

    if(j == q)
      yp[q++] = i;
  }

  return p;
}

uint
select_quadratureg(              psingquadg sq,
                   private const uint3      x,
                   private const uint3      y,
                   private       uint       *xp,
                   private       uint       *yp,
                   global  const real       **xq,
                   global  const real       **yq,
                   global  const real       **wq,
                   private       real       *base)
{
  uint p = (x.x == y.x) + (x.x == y.y) + (x.x == y.z) +
           (x.y == y.x) + (x.y == y.y) + (x.y == y.z) +
           (x.z == y.x) + (x.z == y.y) + (x.z == y.z);

  xp[0] = 0; xp[1] = 1; xp[2] = 2; yp[0] = 0; yp[1] = 1; yp[2] = 2;

  *xq   = sq->xqs + (2 * p * sq->nq);
  *yq   = sq->yqs + (2 * p * sq->nq);
  *wq   = sq->wqs + (p * sq->nq);
  *base = ((private real*) &sq->bases)[p];

  p = 0;

  for(uint i = 0; i < 3; ++i)
  {
    for(uint j = 0; j < 3; ++j)
    {
      if(((private const uint*) &x)[i] == ((private const uint*) &y)[j])
      {
        xp[p] = i;
        yp[p] = j;
        ++p;
        break;
      }
    }
  }

  // printf("%u (%u %u %u) (%u %u %u)\n", p, xp[0], xp[1], xp[2], yp[0], yp[1], yp[2]);

  uint q = p;

  for(uint i = 0; i < 3; ++i)
  {
    uint j = 0;

    for(j = 0; (j < q) &&
               (((private const uint*) &x)[i] !=
                ((private const uint*) &x)[xp[j]]); ++j)
      ;

    if(j == q)
      xp[q++] = i;
  }

  q = p;

  for(uint i = 0; i < 3; ++i)
  {
    uint j = 0;

    for(j = 0; (j < q) &&
               (((private const uint*) &y)[i] !=
                ((private const uint*) &y)[yp[j]]); ++j)
      ;

    if(j == q)
      yp[q++] = i;
  }

  return p;
}

#endif // SINGQUAD_CL
