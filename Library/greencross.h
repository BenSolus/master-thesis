/* ------------------------------------------------------------
 * This is the file "greencross.h" of this master thesis.
 * All rights reserved, Bennet Carstensen 2017
 * ------------------------------------------------------------ */

/**
 * @file      greencross.h
 * @author    Bennet Carstensen
 * @date      2017
 * @copyright All rights reserved, Bennet Carstensen 2017
 */

#ifndef GREENCROSS_H
#define GREENCROSS_H

#include "h2matrix.h"

#ifdef __APPLE__
    #include "OpenCL/opencl.h"
#else
    #include "CL/cl.h"
#endif

/** @defgroup greencross greencross
 *  @brief Algorithms and functions to perform the Green cross approximation
 *         method.
 *
 * @{*/

/** @brief @ref greencross is just an abbreviation for the struct @ref
                _greencross. */
typedef struct _greencross greencross;

/** @brief Pointer to a @ref _greencross "greencross" object. */
typedef greencross *pgreencross;

/** @brief Pointer to a constant @ref _greencross "greencross" object. */
typedef const greencross *pcgreencross;

/** @brief @ref gcopencl is just an abbreviation for the struct @ref
                _gcopencl. */
typedef struct _gcopencl gcopencl;

/** @brief Pointer to a @ref gcopencl object. */
typedef gcopencl *pgcopencl;

/** @brief Pointer to a constant @ref gcopencl object. */
typedef const gcopencl *pcgcopencl;

/** @brief @ref gcopencl is just an abbreviation for the struct @ref
                _gcopencl. */
typedef struct _oclwork oclwork;

/** @brief Pointer to a @ref gcopencl object. */
typedef oclwork *poclwork;

/** @brief Pointer to a constant @ref gcopencl object. */
typedef const oclwork *pcoclwork;

static const uint greencross_min_dim = 2;
static const uint greencross_max_dim = 3;

/** @brief Main container object for performing the Green cross approximation
 *         method. */
struct _greencross
{
  /** @brief Dimension of the problem. */
  uint dim;

  /** @brief Geometry of the problem. */
  void *geom;

  cl_mem *buf_x;

  cl_mem *buf_p;

  cl_event *event_p;

  cl_mem *buf_g;

  /** @brief Number of basis functions */
  uint n;

  /** @brief Index set */
  uint *idx;

  /** @brief Approximation order of Green's formula. */
  uint m;

  uint K;

  real accur;

  /** @brief Row cluster */
  pcluster rc;

  /** @brief Column cluster */
  pcluster cc;

  pclusterbasis rb;

  pclusterbasis cb;

  void *bem;

  cl_mem *buf_qx;

  cl_mem *buf_qy;

  cl_mem *buf_w;

  pgcopencl gcocl;

  poclwork  oclwrk;
};

/** @brief Initilize components related to the dimension of the problem and
 *         zero initializes all other components of the given @ref _greencross
 *         "greencross" object @p gc.
 *
 *  @note Should be called at the beginning of the initialization of each
 *        greencross object.
 *  @note Should always be matched with a call of @ref uninit_greencross.
 *
 *  @param[in, out] gc  The greencross object which needs to be initialized.
 *  @param[in]      dim The dimension of the problem. */
HEADER_PREFIX void
init_greencross(pgreencross gc, uint dim);

/** @brief Uninitilize a @ref _greencross "greencross" object.
 *
 * Invalidates member pointers, freeing corresponding storage if appropriate,
 * and prepares the object for deletion.
 *
 *  @note  Should always be matched with a call of @ref init_greencross.
 *
 *  @param[out] gc  The greencross object which needs to be initialized. */
HEADER_PREFIX void
uninit_greencross(pgreencross gc);

/** @brief Free memory and set pointer to NULL of the corresponding
 *         @ref _greencross greencross object @p gc.
 *
 *  @param gc The greencross object which needs to be deleted. */
HEADER_PREFIX void
del_greencross(pgreencross gc);

/** @brief Creats a <a href="http://www.h2lib.org/doc/d8/da5/struct__clustergeometry.html">
 *         clustergeometry</a> object from a triangular mesh.
 *
 *  The geometry is taken from @p data interpreted as different objects
 *  depending on the dimension of the problem @p dim.
 *
 *  @attention For @p dim = 2, @p data is expected to be an object of type
 *             <a href="http://www.h2lib.org/doc/d8/df9/struct__curve2d.html">
 *             curve2d</a>.
 *
 *  @param[in]  data Pointer to an object containing the informations to
 *                   construct the cluster geometry.
 *  @param[in]  dim  Dimension of the cluster geometry.
 *  @param[out] idx  Index set which will be allocated by this method to index
 *                   each triangle in data. The entries of idx will be
 *                   0, 1, ..., N-1, where N is equal to the number of
 *                   triangles.
 *  @return          A valid cluster geometry object which can be used to
 *                   construct a cluster tree along with the index set @p idx.*/
HEADER_PREFIX pclustergeometry
build_clustergeometry_greencross(const void *data, const uint dim, uint **idx);


/** @brief Constructs the matrix resulting from the Galerkin discretization of a
 *         variational formulation described in @p gc.
 *
 * The arrays @p ridx and @p cidx serve as the indices for the set of base
 * functions described through the geometry of the problem
 * @ref _greencross.geom "greencross::geom". Further choices like quadrature for
 * solving the integral are taken from @p gc.
 *
 * @param[in]  gc    Object containing the problem describtion like base
 *                   functions and Quadratur choices.
 * @param[in]  rsize Number of base functions used from @ref _greencross.geom
 *                   "greencross::geom" for the row of the Galerkin
 *                   discretization matrix.
 * @param[in]  ridx  Indices of the base functions used from @ref
 *                   _greencross.geom "greencross::geom" for the row of the
 *                   Galerkin discretization matrix.
 * @param[in]  csize Number of base functions used from @ref _greencross.geom
 *                   "greencross::geom" for the column of the Galerkin
 *                   discretization matrix.
 * @param[in]  cidx  Indices of the base functions used from @ref
 *                   _greencross.geom "greencross::geom" for the column of the
 *                   Galerkin discretization matrix.
 * @param[out] G     Matrix which will contain the Galerkin discretization. */
HEADER_PREFIX void
assemble_amatrix_greencross(pcgreencross gc,
                            const uint   *ridx,
                            const uint   *cidx,
                            bool         ntrans,
                            pamatrix     G);

HEADER_PREFIX pamatrix
build_amatrix_greencross(pcgreencross gc,
                         const uint rows,
                         const uint *ridx,
                         const uint cols,
                         const uint *cidx,
                         bool       ntrans);

HEADER_PREFIX void
assemble_green_rkmatrix_greencross(pcgreencross gc,
                                   pccluster    row,
                                   pccluster    col,
                                   prkmatrix    AB);

HEADER_PREFIX prkmatrix
build_green_rkmatrix_greencross(pcgreencross gc, pccluster row, pccluster col);

HEADER_PREFIX void
assemble_green_hmatrix_greencross(pgreencross gc, pcblock b, phmatrix H);

HEADER_PREFIX phmatrix
build_green_hmatrix_greencross(pgreencross gc, void *eta);

uint
assemble_green_cross_leaf_clusterbasis_greencross(pcgreencross  gc,
                                                  const bool    row,
                                                  pclusterbasis c,
                                                  uint          **pidx);

/** @brief Constructs a Green cross approximation of a problem given by a
 *         @ref _greencross "greencross", represented as
 *         an <a href="http://www.h2lib.org/doc/d7/ddd/group__h2matrix.html">
 *         @f$\mathcal{H}^{2}@f$-Matrix</a> \cite Börm2016GCA.
 *
 * Currently, if we find an admissible leaf by the clusters @p t and @p s given
 * by the @ref _greencross "greencross" object @p gc, we calculate the first
 * factor of Green's representation formula as a degenerated approximation
 * @f$A_{t} \in \mathbb{\hat{t} \times K}@f$ with @f$\hat{t}@f$ and K given
 * by @p t and @p gc respectively. We approximate this matrix further with the
 * adaptive cross approximation using full pivoting, resulting in a low-rank
 * approximation @f$A_{t} \approx CD^{*}@f$ and the row permutation matrix
 * @f$P@f$. Finaly, with @f$V_{t} := C(PC)^{-1}@f$, @f$S_{ts} :=
 * PG|_{\hat{t} \times \hat{s}}@f$ and @f$V_{s} := I_{\hat{s}}f$, with
 * @f$I_{\hat{s}} \in \mathbb{R}^{\hat{s} \times \hat{s}}@f$ the identity
 * matrix and @f$\hat{s}@f$ given by @p s, we have an approximation
 * @f$G|_{\hat{t} \times \hat{s}} \approx V_{t} S_{ts} V_{s}^{*}.
 *
 * @param[in] gc The @ref _greencross "greencross" object contains the problem
 *                description and several data to perform the approximation.
 * @result    H2 <a href="http://www.h2lib.org/doc/d7/ddd/group__h2matrix.html">
 *               @f$\mathcal{H}^{2}@f$-Matrix</a> object which will be filled
 *               according to the Green cross approximation method. */
HEADER_PREFIX ph2matrix
build_green_cross_h2matrix_greencross(pgreencross gc, void *eta);

HEADER_PREFIX void
addeval_h2matrix_avector_greencross(pgreencross gc,
                                    field       alpha,
                                    pch2matrix  H2,
                                    pcavector   x,
                                    pavector    y);

HEADER_PREFIX void
addevaltrans_h2matrix_avector_greencross(pgreencross gc,
                                         field       alpha,
                                         pch2matrix  H2,
                                         pcavector   x,
			                                   pavector    y);

HEADER_PREFIX void
mvm_h2matrix_avector_greencross(pgreencross gc,
                                field       alpha,
                                bool        h2trans,
                                pch2matrix  H2,
                                pcavector   x,
		                            pavector    y);

/** @} */

#endif // GREENCROSS_H
