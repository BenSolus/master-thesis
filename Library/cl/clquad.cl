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

#include "clbasic.cl"

/** @brief @ref quad is just an abbreviation for the struct @ref _quad. */
typedef struct _quad quad;

/** @brief Pointer to a @ref _quad "quad" object. */
typedef quad *pquad;

/** @brief Pointer to a constant @ref _quad "quad" object. */
typedef const quad *pcquad;

struct _quad
{
           uint nq;

  constant real *qx;

  constant real *qy;

  constant real *w;
};

void
init_quad(         pquad      q,
                   const uint nq,
          constant       real *qx,
          constant       real *qy,
          constant       real *w)
{
  q->nq = nq;
  q->qx = qx;
  q->qy = qy;
  q->w  = w;
}
