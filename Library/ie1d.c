/* ------------------------------------------------------------
 * This is a shortened version of the
 *  "ie1d.c" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2016
 * ------------------------------------------------------------ */

 /** @file ie1d.c
 *  @author Steffen B&ouml;rm
 */

#include <stdio.h>

#include "ie1d.h"
#include "basic.h"
#include "cluster.h"
#include "gaussquad.h"

/* ------------------------------------------------------------
 * Constructor and destructor
 * ------------------------------------------------------------ */

pie1d
new_ie1d(uint n, uint q)
{
  pie1d     ie;

  /* Allocate storage */
  ie = (pie1d) allocmem(sizeof(ie1d));

  /* Initialize entries */
  ie->n = n;
  ie->mmax = 0;
  ie->xi = 0;
  ie->m = 0;

  /* Initialize callback functions */
  ie->nearfield = 0;
  ie->fillrk = 0;
  ie->fillV = 0;
  ie->fillE = 0;
  ie->fillS = 0;

  /* Initialize Gauss quadrature */
  ie->q = q;
  if(q > 0){
    ie->xq = allocreal(q);
    ie->wq = allocreal(q);
    assemble_gauss(q, ie->xq, ie->wq);
  }

  return ie;
}

void
del_ie1d(pie1d ie)
{
  uint      m;

  if(ie->q > 0){
  /* Release quadrature points */
  freemem(ie->wq);
  freemem(ie->xq);
  }
  /* Release interpolation points */
  if (ie->xi) {
    for (m = 1; m <= ie->mmax; m++)
      freemem(ie->xi[m]);
    freemem(ie->xi);

    /* Make sure the released storage cannot be used by mistake */
    ie->xi = 0;
  }

  /* Release object */
  freemem(ie);
}

/* ------------------------------------------------------------
 * Construct a cluster tree
 * ------------------------------------------------------------ */

static pcluster
buildcluster(uint size, uint off, uint * idx, real h, uint leafsize)
{
  pcluster  t;
  uint      size0, size1;

  t = 0;
  if (size <= leafsize)
    t = new_cluster(size, idx + off, 0, 1);
  else {
    t = new_cluster(size, idx + off, 2, 1);

    size0 = size / 2;
    t->son[0] = buildcluster(size0, off, idx, h, leafsize);

    size1 = size - size0;
    t->son[1] = buildcluster(size1, off + size0, idx, h, leafsize);
  }

  t->bmin[0] = off * h;
  t->bmax[0] = (off + size) * h;

  update_cluster(t);

  return t;
}

pcluster
build_ie1d_cluster(pcie1d ie, uint leafsize)
{
  uint      n = ie->n;
  uint     *idx;
  uint      i;

  idx = (uint *) allocmem(sizeof(uint) * n);
  for (i = 0; i < n; i++)
    idx[i] = i;

  return buildcluster(n, 0, idx, 1.0 / n, leafsize);
}


/* ------------------------------------------------------------
 * Nearfield integration
 * ------------------------------------------------------------ */

static    real
antiderivative(real x, real y)
{
  if (x == y)
    return 0.0;
  else
    return 0.25 * REAL_SQR(x - y) * (3.0 - 2.0 * REAL_LOG(REAL_ABS(x - y)));
}

void
nearfield_ie1d(const uint * ridx, const uint * cidx, pcie1d ie, pamatrix N)
{
  uint      rows = N->rows;
  uint      cols = N->cols;
  real      h = 1.0 / ie->n;
  uint      i, j, ii, jj;
  real      x0, x1, y0, y1, val;

  for (j = 0; j < cols; j++) {
    jj = (cidx ? cidx[j] : j);
    y0 = jj * h;
    y1 = (jj + 1) * h;

    for (i = 0; i < rows; i++) {
      ii = (ridx ? ridx[i] : i);
      x0 = ii * h;
      x1 = (ii + 1) * h;

      val =
	antiderivative(x1, y0) + antiderivative(x0,
						y1) - 2.0 * antiderivative(x0,
									   y0);

      setentry_amatrix(N, i, j, val);
    }
  }
}


/* ------------------------------------------------------------
 *  Determine approximation orders
 * ------------------------------------------------------------ */

static    real
  maxorder(real m0, real m1, pccluster t)
  {
    real      mt;
    uint      i;

    mt = m0;

    if (t->sons > 0) {
      mt = 0.0;
      for (i = 0; i < t->sons; i++)
        mt = REAL_MAX(mt, maxorder(m0, m1, t->son[i]));

      mt += m1;
    }

    return mt;
  }

static void
  setorder(real mt, real m0, real m1, uint tname, pccluster t, uint * m)
  {
    uint      tname1;
    uint      i;

    m[tname] = (uint) mt;

    tname1 = tname + 1;
    for (i = 0; i < t->sons; i++) {
      setorder(mt - m1, m0, m1, tname1, t->son[i], m);

      tname1 += t->son[i]->desc;
    }
    assert(tname1 == tname + t->desc);
  }


void
  prepare_orders_ie1d(pie1d ie, pccluster root, real m0, real m1)
  {
    real      mroot;
    uint      mmax;
    uint      i, m;

    assert(ie->xi == 0);

    mroot = maxorder(m0, m1, root);

    ie->mmax = mmax = (uint) mroot;
    ie->xi = (preal *) allocmem(sizeof(preal) * (mmax + 1));

    ie->xi[0] = 0;
    for (m = 1; m <= mmax; m++) {
      ie->xi[m] = allocreal(m);
      for (i = 0; i < m; i++)
        ie->xi[m][i] = REAL_COS(M_PI * (m - i - 0.5) / m);
    }

    ie->m = (uint *) allocmem(sizeof(uint) * root->desc);
    setorder(mroot, m0, m1, 0, root, ie->m);
  }

/* ------------------------------------------------------------
 * Approximation by Taylor expansion
 * ------------------------------------------------------------ */


void
  setup_aprx_taylor_ie1d(pie1d ie, pccluster root, real m0, real m1)
  {
    ie->nearfield = nearfield_ie1d;
    ie->fillrk = fillrk_taylor_ie1d;

    prepare_orders_ie1d(ie, root, m0, m1);
  }


void
fillrk_taylor_ie1d(pccluster rc, uint rname,
		   pccluster cc, uint cname, pcie1d ie, prkmatrix r)
  {
    /* Fill in your version */

    /* Which is the smaller cluster? */
    const bool rcLeqcc = rc->size <= cc->size;
    /* 1 / n */
    const real h       = 1.0 / ie->n;
    /* Central point of the smaller cluster */
    const real xt      = rcLeqcc ? (rc->bmax[0] + rc->bmin[0]) * 0.5
                                 : (cc->bmax[0] + cc->bmin[0]) * 0.5;

    setrank_rkmatrix(r, ie->m[rcLeqcc ? rname : cname]);

    /* Rank of r */
    const uint k       = r->k;

    /* Index set used for the integrals of the monomials */
    const uint *monomIdx = rcLeqcc ? rc->idx : cc->idx;
    /* Index set used for the integrals of the derivations */
    const uint *deriIdx  = rcLeqcc ? cc->idx : rc->idx;

    /* Matrix used to save the integrals of the monomials */
    pamatrix monomMatrix = rcLeqcc ? &r->A   : &r->B;
    /* Matrix used to save the integrals of the derivations */
    pamatrix deriMatrix  = rcLeqcc ? &r->B   : &r->A;

    /* Number of rows for the current matrix */
    uint n;

    /*
     * First calculate all integrals of the monomials and save them in the
     * appropriate matrix
     */
    n = monomMatrix->rows;

    for(uint i = 0; i < n; ++i)
    {
      const real ii     = (real) (monomIdx ? monomIdx[i] : i);
      const real imxt   = ((ii * h) - xt);
      const real ip1mxt = (((ii + 1.0) * h) - xt);

      real ddend, dsor, minend, subend;

      dsor   = 1.0;
      minend = 1.0;
      subend = 1.0;

      for(uint v = 0; v < k; ++v)
      {
        /* divisor = (v + 1)! */
        dsor *= (real) (v + 1);

        /* minuend = (((i + 1) / n) - xt)^{v} */
        minend *= ip1mxt;
        /* subtrahend = ((i / n) - xt)^{v} */
        subend *= imxt;

        /* dividend = (((i + 1) / n) - xt)^{v} - ((i / n) - xt)^{v} */
        ddend   = minend - subend;

        /*
         *                      (((i + 1) / n) - xt)^{v} - ((i / n) - xt)^{v}
         * monomMatrix_{i, v} = ---------------------------------------------
         *                                            (v + 1)!
         */
        setentry_amatrix(monomMatrix, i, v, ddend / dsor);
      }
    }

    /*
     * Second calculate all integrals of the derivations and save them in the
     * appropriate matrix
     */
    n = deriMatrix->rows;

    for(uint j = 0; j < n; ++j)
    {
      const real jj = (real) (deriIdx ? deriIdx[j] : j);

      const real xtmj   = xt - (jj * h);
      const real xtmjp1 = xt - ((jj + 1.0) * h);

      /*
       * deriMatrix_{i, 1} = (xt - ((i + 1) / n)) * ln|xt - ((i + 1) / n)| -
       *                     (xt - (i / n))       * ln|xt - (i / n)|
       */
      setentry_amatrix(deriMatrix,
                       j,
                       0,
                       xtmjp1 * (REAL_LOG(REAL_ABS(xtmjp1)) - 1.0) -
                       xtmj   * (REAL_LOG(REAL_ABS(xtmj))   - 1.0));

      /* deriMatrix_{i, 1} = ln|xt - ((i + 1) / n)| - ln|xt - (i / n)| */
      setentry_amatrix(deriMatrix,
                       j,
                       1,
                       REAL_LOG(REAL_ABS(xtmjp1)) -
                       REAL_LOG(REAL_ABS(xtmj)));

      real ddend, dsor, facVm2, minend, subend;

      facVm2 = 1.0;
      minend = 1.0;
      subend = 1.0;

      for(uint v = 2; v < k; ++v)
      {
        /* minuend = (xt - (i / n))^{v - 1} */
        minend *= xtmj;
        /* minuend = (xt - ((i + 1) / n))^{v - 1} */
        subend *= xtmjp1;

        /* dividend = (-1)^v * (v - 2)! * ((xt - (i / n))^{v - 1} - (xt - ((i + 1) / n))^{v - 1}) */
        ddend = facVm2 * (minend - subend);
        /* divisor = (xt - ((i + 1) / n))^{v - 1} * (xt - (i / n))^{v - 1} */
        dsor  = minend * subend;

        /*
         *                                         (xt - (i / n))^{v - 1} - (xt - ((i + 1) / n))^{v - 1}
         * deriMatrix_{j, v} = (-1)^v * (v - 2)! * -----------------------------------------------------
         *                                         (xt - ((i + 1) / n))^{v - 1} * (xt - (i / n))^{v - 1}
         */
        setentry_amatrix(deriMatrix, j, v, ddend / dsor);

        /* facVm2 = (-1.0)^v * (v - 2)! (for the next iteration) */
        facVm2 *= -1.0 * (real) (v - 1.0);
      }
    }
  }



/* ------------------------------------------------------------
 * Fill a hierarchical matrix
 * ------------------------------------------------------------ */

static void
fill_hmatrix(phmatrix G, uint rname, uint cname, pcie1d ie)
{
  uint      rsons, csons;
  uint      rname1, cname1;
  uint      i, j;

  if (G->son) {
    /* Handle submatrices */
    rsons = G->rsons;
    csons = G->csons;

    /* Compute name of first column son */
    cname1 = (G->son[0]->cc == G->cc ? cname : cname + 1);

    for (j = 0; j < csons; j++) {
      /* Compute name of first row son */
      rname1 = (G->son[0]->rc == G->rc ? rname : rname + 1);

      for (i = 0; i < rsons; i++) {
	/* Fill submatrices */
	fill_hmatrix(G->son[i + j * rsons], rname1, cname1, ie);

	/* Switch to next row son */
	rname1 += G->son[i]->rc->desc;
      }
      assert(rname1 == rname + G->rc->desc);

      /* Switch to next column son */
      cname1 += G->son[j * rsons]->cc->desc;
    }
    assert(cname1 == cname + G->cc->desc);
  }
  else if (G->f) {
    /* Fill inadmissible submatrix */
    ie->nearfield(G->rc->idx, G->cc->idx, ie, G->f);
  }
  else if (G->r) {
    /* Fill admissible submatrix */
    ie->fillrk(G->rc, rname, G->cc, cname, ie, G->r);
  }
}

void
fill_hmatrix_ie1d(pcie1d ie, phmatrix G)
{
  fill_hmatrix(G, 0, 0, ie);
}


static    field
  eval_kernel(real x, real y)
  {
    return -REAL_LOG(REAL_ABS(x - y));
  }

static    real
  eval_lagrange(uint nu, uint m, pcreal xi, real x)
  {
    real      val;
    uint      i;

    val = 1.0;

    for (i = 0; i < nu; i++)
      val *= (x - xi[i]) / (xi[nu] - xi[i]);

    for (i = nu + 1; i < m; i++)
      val *= (x - xi[i]) / (xi[nu] - xi[i]);

    return val;
  }

void
  fillrk_interpolation_ie1d(pccluster rc, uint rname,
                            pccluster cc, uint cname, pcie1d ie, prkmatrix r)
  {
    /* rc < cc ?*/
    const bool  rcLeqcc  = rc->size <= cc->size;

    /*
     * Set rank of the rkmatrix to the approximation order of the smaller
     * cluster
     */
    setrank_rkmatrix(r, ie->m[rcLeqcc ? rname : cname]);

    const uint  k        = r->k;        /* approximation order */
    const uint  q        = ie->q;       /* quadrature order */
    const real  h        = 1.0 / ie->n; /* step size */

    /* Index set used for the integrals of the Lagrange polynomials */
    const uint* lagrangeIdx = rcLeqcc ? rc->idx : cc->idx;
    /* Index set used for the integrals of the kernel function */
    const uint* kernelIdx   = rcLeqcc ? cc->idx : rc->idx;

    pcreal xi = ie->xi[k]; /* interpolation points for [-1, 1] */
    pcreal xq = ie->xq;    /* quadrature points for [-1, 1] */
    pcreal w  = ie->wq;    /* quadrature weights */

    /* Matrix used to save the integrals of the Lagrange polynomials */
    pamatrix    lagrangeMat = rcLeqcc ? &r->A   : &r->B;
    /* Matrix used to save the integrals of the kernel function */
    pamatrix    kernelMat   = rcLeqcc ? &r->B   : &r->A;

    /* Transformed interpolation points */
    preal phixi = (preal) malloc(sizeof(real) * k);
    /* transformed quadrature points */
    preal phixq = (preal) malloc(sizeof(real) * q);

    /* Number of rows of the current matrix */
    uint  n;

    /*
     * First calculate all integrals of the Lagrange polynomials and save them
     * in the appropriate matrix
     */
    n = lagrangeMat->rows;

    /* Factors for point-mapping */
    real bma  = ((lagrangeIdx ? lagrangeIdx[n - 1]   : n - 1) + 1
                 - (lagrangeIdx ? lagrangeIdx[0]     : 0))         * h * 0.5;
    real bpa  = ((lagrangeIdx ? lagrangeIdx[0]       : 0)
                 + (lagrangeIdx ? lagrangeIdx[n - 1] : n - 1) + 1) * h * 0.5;

    /* Transform the interpolation points */
    for(uint i = 0; i < k; ++i)
      phixi[i] = bpa + bma * xi[i];

    bma = h * 0.5;

    for(uint i = 0; i < n; ++i)
    {
      const real ii    = lagrangeIdx ? lagrangeIdx[i] : i;
      bpa   = (ii + 0.5) * h;

      /* Transform the quadrature points */
      for(uint j = 0; j < q; ++j)
        phixq[j] = bpa + bma * xq[j];

      for(uint v = 0; v < k; ++v)
      {
        real y = 0.0;

        for(uint j = 0; j < q; ++j)
          y += w[j] * eval_lagrange(v, k, phixi, phixq[j]);

        /* lagrangeMat_{i, v} = (1 / (2 * n)) * sum_{j = 0}^q w_j L_v(xq_j) */
        setentry_amatrix(lagrangeMat, i, v, bma * y);
      }
    }

    /*
     * Second calculate all integrals of the kernel function and save them in the
     * appropriate matrix
     */

    n = kernelMat->rows;

    for(uint j = 0; j < n; ++j)
    {
      const real jj  = (real) (kernelIdx ? kernelIdx[j] : j);
      bpa = (jj + 0.5) * h;

      /* Transform the quadrature points */
      for(uint i = 0; i < q; ++i)
        phixq[i] = bpa + bma * xq[i];

      for(uint v = 0; v < k; ++v)
      {
        real y = 0.0;

        for(uint i = 0; i < q; ++i)
          y += w[i] * eval_kernel(phixi[v], phixq[i]);

        /* kernelMat_{j, v} = (1 / (2 * n)) * sum_{i = 0}^q w_i g(xi_v, xq_i) */
        setentry_amatrix(kernelMat, j, v, bma * y);
      }
    }
  }

 /* ----------------------------------------------
    NEW - STUFF
    ---------------------------------------------- */

void
fillV_interpolation_ie1d(pccluster t, uint tname, pcie1d ie, pamatrix V)
{

  /* Fill cluster basis leaf case */

  const uint k = ie->m[tname]; /* approximation order */
  const uint q = ie->q; /* quadrature order */
  const uint n = t->size; /* matrix row number */
  const real h = 1.0 / ie->n; /* step size */

  /* Index set used for the integrals of the Lagrange polynomials */
  const uint* lagrangeIdx = t->idx;

  pcreal xi = ie->xi[k]; /* interpolation points for [-1, 1] */
  pcreal xq = ie->xq;    /* quadrature points for [-1, 1] */
  pcreal w = ie->wq;     /* quadrature weights */

  /* Transformed interpolation points */
  preal phixi = (preal) malloc(sizeof(real) * k);
  /* Transformed quadrature points */
  preal phixq = (preal) malloc(sizeof(real) * q);

  /* Resize matrix V */
  resize_amatrix(V, n, k);

  /* Factors for point-mapping */
  real bma = ((lagrangeIdx ? lagrangeIdx[n - 1] : n - 1) + 1
             - (lagrangeIdx ? lagrangeIdx[0] : 0)) * h * 0.5;
  real bpa = ((lagrangeIdx ? lagrangeIdx[0] : 0)
             + (lagrangeIdx ? lagrangeIdx[n - 1] : n - 1) + 1) * h * 0.5;

  /* Transform the interpolation points */
  for(uint i = 0; i < k; ++i)
    phixi[i] = bpa + bma * xi[i];

  bma = h * 0.5;

  for(uint i = 0; i < n; ++i) {
    const real ii = lagrangeIdx ? lagrangeIdx[i] : i;
    bpa = (ii + 0.5) * h;

    /* Transform the quadrature points */
    for(uint j = 0; j < q; ++j)
      phixq[j] = bpa + bma * xq[j];

    for(uint v = 0; v < k; ++v) {
      real y = 0.0;
      for(uint j = 0; j < q; ++j)
        y += w[j] * eval_lagrange(v, k, phixi, phixq[j]);

      setentry_amatrix(V, i, v, bma * y);
    }
  }
}

void
fillE_interpolation_ie1d(pccluster s, uint sname,
			 pccluster f, uint fname, pcie1d ie, pamatrix E)
{

  /* Fill cluster basis transfer case, 's' is the son, 'f' the father */

  const uint sk = ie->m[sname];
  const uint fk = ie->m[fname]; /* approximation orders */
  const uint sn = s->size;
  const uint fn = f->size; /* cluster sizes */
  const real h = 1.0 / ie->n; /* step size */

  /* Index set used for the integrals of the Lagrange polynomials */
  const uint* sidx = s->idx;
  const uint* fidx = f->idx;

  pcreal sxi = ie->xi[sk];
  pcreal fxi = ie->xi[fk]; /* interpolation points for [-1, 1] */

  /* Transformed interpolation points */
  preal sphixi = (preal) malloc(sizeof(real) * sk);
  preal fphixi = (preal) malloc(sizeof(real) * fk);

  /* Resize matrix E */
  resize_amatrix(E, sk, fk);

  /* Factors for point-mapping */
  real sbma = ((sidx ? sidx[sn - 1] : sn - 1) + 1
             - (sidx ? sidx[0] : 0)) * h * 0.5;
  real sbpa = ((sidx ? sidx[0] : 0)
             + (sidx ? sidx[sn - 1] : sn - 1) + 1) * h * 0.5;
  real fbma = ((fidx ? fidx[fn - 1] : fn - 1) + 1
             - (fidx ? fidx[0] : 0)) * h * 0.5;
  real fbpa = ((fidx ? fidx[0] : 0)
             + (fidx ? fidx[fn - 1] : fn - 1) + 1) * h * 0.5;

  /* Transform the interpolation points */
  for(uint i = 0; i < sk; ++i)
    sphixi[i] = sbpa + sbma * sxi[i];
  for(uint j = 0; j < fk; ++j)
    fphixi[j] = fbpa + fbma * fxi[j];

  for(uint i = 0; i < sk; ++i) {
    for(uint j = 0; j < fk; ++j)
      setentry_amatrix(E, i, j, eval_lagrange(j, fk, fphixi, sphixi[i]));
  }
}

void
fillS_interpolation_ie1d(pccluster rc, uint rname,
			 pccluster cc, uint cname, pcie1d ie, pamatrix S)
{

  /* Fill coupling matrix corresponding to rc[rname] \times cc[cname] */

  const uint* ridx = rc->idx;
  const uint* cidx = cc->idx; /* index sets */
  const uint rk = ie->m[rname];
  const uint ck = ie->m[cname]; /* approximation orders */
  const uint rn = rc->size;
  const uint cn = cc->size; /* cluster sizes */
  const real h = 1.0 / ie->n; /* step size */
  pcreal rxi = ie->xi[rk];
  pcreal cxi = ie->xi[ck]; /* interpolation points for [-1, 1] */

  /* Transformed interpolation points */
  preal rphixi = (preal) malloc(sizeof(real) * rk);
  preal cphixi = (preal) malloc(sizeof(real) * ck);

  /* Resize matrix S */
  resize_amatrix(S, rk, ck);

  /* Factors for point-mapping */
  real rbma  = ((ridx ? ridx[rn - 1] : rn - 1) + 1
               - (ridx ? ridx[0] : 0)) * h * 0.5;
  real rbpa  = ((ridx ? ridx[0] : 0)
               + (ridx ? ridx[rn - 1] : rn - 1) + 1) * h * 0.5;
  real cbma  = ((cidx ? cidx[cn - 1] : cn - 1) + 1
               - (cidx ? cidx[0] : 0)) * h * 0.5;
  real cbpa  = ((cidx ? cidx[0] : 0)
               + (cidx ? cidx[cn - 1] : cn - 1) + 1) * h * 0.5;

  /* Transform the interpolation points */
  for(uint i = 0; i < rk; ++i)
    rphixi[i] = rbpa + rbma * rxi[i];
  for(uint j = 0; j < ck; ++j)
    cphixi[j] = cbpa + cbma * cxi[j];

  for(uint i = 0; i < rk; ++i) {
    for(uint j = 0; j < ck; ++j)
      setentry_amatrix(S, i, j, eval_kernel(rphixi[i], cphixi[j]));
  }
}

void
  setup_aprx_interpolation_ie1d(pie1d ie, pccluster root, real m0, real m1)
  {
    ie->nearfield = nearfield_ie1d;
    ie->fillrk = fillrk_interpolation_ie1d;

  ie->fillV = fillV_interpolation_ie1d;
  ie->fillE = fillE_interpolation_ie1d;
  ie->fillS = fillS_interpolation_ie1d;

    prepare_orders_ie1d(ie, root, m0, m1);
  }

/*---------------------------------------------------
   H2-Matrix
  --------------------------------------------------- */



/* ------------------------------------------------------------
 * Fill a cluster basis
 * ------------------------------------------------------------ */

static void
fill_clusterbasis(pclusterbasis cb, uint tname, pcie1d ie)
{
  uint      i, tname1;

  if (cb->son) {
    /* Compute name of first son */
    tname1 = tname + 1;
    for (i = 0; i < cb->sons; i++) {
      /* Prepare son's cluster basis */
      fill_clusterbasis(cb->son[i], tname1, ie);

      /* Switch to next son */
      tname1 += cb->son[i]->t->desc;
    }
    assert(tname1 == tname + cb->t->desc);

    /* Set the rank */
    setrank_clusterbasis(cb, ie->m[tname]);

    /* Compute name of first son */
    tname1 = tname + 1;
    for (i = 0; i < cb->sons; i++) {
      /* Fill transfer matrix */
      ie->fillE(cb->son[i]->t, tname1, cb->t, tname,
		ie, getE_clusterbasis(cb->son[i]));

      /* Switch to next son */
      tname1 += cb->son[i]->t->desc;
    }
    assert(tname1 == tname + cb->t->desc);
  }
  else {
    /* Set the rank */
    setrank_clusterbasis(cb, ie->m[tname]);

    /* Fill leaf matrix */
    ie->fillV(cb->t, tname, ie, getV_clusterbasis(cb));
  }
}

void
fill_clusterbasis_ie1d(pcie1d ie, pclusterbasis cb)
{
  fill_clusterbasis(cb, 0, ie);
}

/* ------------------------------------------------------------
 * Fill an H2-matrix
 * ------------------------------------------------------------ */

static void
fill_h2matrix(ph2matrix G, uint rname, uint cname, pcie1d ie)
{
  uint      rsons, csons;
  uint      rname1, cname1;
  uint      i, j;

  if (G->son) {
    /* Handle submatrices */
    rsons = G->rsons;
    csons = G->csons;

    /* Compute name of first column son */
    cname1 = (G->son[0]->cb == G->cb ? cname : cname + 1);

    for (j = 0; j < csons; j++) {
      /* Compute name of first row son */
      rname1 = (G->son[0]->rb == G->rb ? rname : rname + 1);

      for (i = 0; i < rsons; i++) {
	/* Fill submatrices */
	fill_h2matrix(G->son[i + j * rsons], rname1, cname1, ie);

	/* Switch to next row son */
	rname1 += G->son[i]->rb->t->desc;
      }
      assert(rname1 == rname + G->rb->t->desc);

      /* Switch to next column son */
      cname1 += G->son[j * rsons]->cb->t->desc;
    }
    assert(cname1 == cname + G->cb->t->desc);
  }
  else if (G->f) {
    /* Fill inadmissible submatrix */
    ie->nearfield(G->rb->t->idx, G->cb->t->idx, ie, G->f);
  }
  else if (G->u) {
    /* Fill admissible submatrix */
    ie->fillS(G->rb->t, rname, G->cb->t, cname, ie, getS_uniform(G->u));
  }
}

void
fill_h2matrix_ie1d(pcie1d ie, ph2matrix G)
{
  fill_h2matrix(G, 0, 0, ie);
}
