/* ------------------------------------------------------------
 * This is the file "clsingquad.cl" of this master thesis.
 * All rights reserved, Steffen Boerm 2009
 * ------------------------------------------------------------ */

/**
 * @file      cl/clsingquad.cl
 * @author    Bennet Carstensen
 * @date      2017
 * @copyright All rights reserved, Bennet Carstensen 2017
 */

#include "clbasic.cl"

#ifdef USE_FLOAT
typedef float4  real4;
#else
typedef double4  real4;
#endif

// /** @brief @ref quadc is just an abbreviation for the struct @ref _quadc. */
// typedef struct _quadc quadc;
//
// /** @brief Pointer to a @ref _quadc "quadc" object. */
// typedef quadc *pquadc;
//
// /** @brief Pointer to a constant @ref _quadc "quadc" object. */
// typedef const quadc *pcquadc;

/** @brief @ref quadg is just an abbreviation for the struct @ref _quadg. */
typedef struct _singquadg singquadg;

/** @brief Pointer to a @ref _quadg "quadg" object. */
typedef singquadg *psingquadg;

/** @brief Pointer to a constant @ref _quadg "quadg" object. */
typedef const singquadg *pcsingquadg;

// struct _quadc
// {
//            uint nq;
//
//   constant real *qx;
//
//   constant real *qy;
//
//   constant real *w;
// };

struct _singquadg
{
  uint nq;

  global const real *xqs;
  global const real *yqs;
  global const real *wqs;

  real4 bases;
};

// void
// init_quadc(              pquadc q,
//                    const uint   nq,
//            constant      real   *qx,
//            constant      real   *qy,
//            constant      real   *w)
// {
//   q->nq = nq;
//   q->qx = qx;
//   q->qy = qy;
//   q->w  = w;
// }

void
init_singquadg(             psingquadg sq,
                      const uint       dim,
                      const uint       nq,
               global const real       *xqs,
               global const real       *yqs,
               global const real       *wqs,
               global const real       *bases)
{
  sq->nq        = nq;

  sq->xqs       = xqs;
  sq->yqs       = yqs;
  sq->wqs       = wqs;

  sq->bases     = vload4(0, bases);
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
