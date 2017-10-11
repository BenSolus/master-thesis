/* ------------------------------------------------------------
 * This is the file "singquadgca.c" of this master thesis.
 * All rights reserved, Bennet Carstensen 2017
 * ------------------------------------------------------------ */

/**
 * @file      singquadgca.c
 * @author    Bennet Carstensen
 * @date      2017
 * @copyright All rights reserved, Bennet Carstensen 2017
 */

#include "singquadgca.h"

#include "basic.h"
#include "ocl_system.h"

#include <string.h>

void
init_singquadgca(psingquadgca sq_gca)
{
  if(sq_gca != NULL)
  {
    sq_gca->offset    = 0;

    memset(sq_gca->nq, 0, sizeof(uint) * 3);

    sq_gca->xqs       = NULL;
    sq_gca->yqs       = NULL;
    sq_gca->wqs       = NULL;
    sq_gca->bases     = NULL;
    sq_gca->buf_nq    = NULL;
    sq_gca->buf_xqs   = NULL;
    sq_gca->buf_yqs   = NULL;
    sq_gca->buf_wqs   = NULL;
    sq_gca->buf_bases = NULL;
  }
}

void
uninit_singquadgca(psingquadgca sq_gca)
{
  if(sq_gca != NULL)
  {
    if(sq_gca->xqs != NULL)
      freemem(sq_gca->xqs);

    if(sq_gca->yqs != NULL)
      freemem(sq_gca->yqs);

    if(sq_gca->wqs != NULL)
      freemem(sq_gca->wqs);

    if(sq_gca->buf_nq != NULL)
    {
      for(uint i = 0; i < ocl_system.num_devices; ++i)
        CL_CHECK(clReleaseMemObject(sq_gca->buf_nq[i]));

      freemem(sq_gca->buf_nq);
    }

    if(sq_gca->buf_xqs != NULL)
    {
      for(uint i = 0; i < ocl_system.num_devices; ++i)
        CL_CHECK(clReleaseMemObject(sq_gca->buf_xqs[i]));

      freemem(sq_gca->buf_xqs);
    }

    if(sq_gca->buf_yqs != NULL)
    {
      for(uint i = 0; i < ocl_system.num_devices; ++i)
      CL_CHECK(clReleaseMemObject(sq_gca->buf_yqs[i]));

      freemem(sq_gca->buf_yqs);
    }

    if(sq_gca->buf_wqs != NULL)
    {
      for(uint i = 0; i < ocl_system.num_devices; ++i)
      CL_CHECK(clReleaseMemObject(sq_gca->buf_wqs[i]));

      freemem(sq_gca->buf_wqs);
    }

    if(sq_gca->buf_bases != NULL)
    {
      for(uint i = 0; i < ocl_system.num_devices; ++i)
      CL_CHECK(clReleaseMemObject(sq_gca->buf_bases[i]));

      freemem(sq_gca->buf_bases);
    }
  }

  init_singquadgca(sq_gca);
}

void
del_singquadgca(psingquadgca sq_gca)
{
  uninit_singquadgca(sq_gca);

  if(sq_gca != NULL)
    freemem(sq_gca);
}

psingquadgca
build_from_singquad2d(psingquad2d sq)
{
  psingquadgca sq_gca = (psingquadgca) allocmem(sizeof(singquadgca));

  init_singquadgca(sq_gca);

  const uint nq = sq_gca->nq[0] = sq->n_dist;

  sq_gca->xqs   = allocreal(8 * nq);
  sq_gca->yqs   = allocreal(8 * nq);
  sq_gca->wqs   = allocreal(4 * nq);
  sq_gca->bases = allocreal(4);

  uint vnq = ROUNDUP(sq->n_dist, VREAL);

  real *xq = sq->x_dist;
  real *yq = sq->y_dist;
  real *wq = sq->w_dist + (9 * vnq);

  memcpy(sq_gca->xqs, xq, nq * sizeof(real));
  memcpy(sq_gca->xqs + nq, xq + vnq, nq * sizeof(real));
  memcpy(sq_gca->yqs, yq, nq * sizeof(real));
  memcpy(sq_gca->yqs + nq, yq + vnq, nq * sizeof(real));
  memcpy(sq_gca->wqs, wq, nq * sizeof(real));

  sq_gca->bases[0]  = sq->base_dist;

  vnq = ROUNDUP(sq->n_vert, VREAL);

  xq = sq->x_vert;
  yq = sq->y_vert;
  wq = sq->w_vert + (9 * vnq) + (sq->n_vert - nq);

  memcpy(sq_gca->xqs + (2 * nq), xq + (sq->n_vert - nq), nq * sizeof(real));
  memcpy(sq_gca->xqs + (3 * nq), xq + vnq + (sq->n_vert - nq), nq * sizeof(real));
  memcpy(sq_gca->yqs + (2 * nq), yq + (sq->n_vert - nq), nq * sizeof(real));
  memcpy(sq_gca->yqs + (3 * nq), yq + vnq + (sq->n_vert - nq), nq * sizeof(real));
  memcpy(sq_gca->wqs + nq, wq, nq * sizeof(real));

  sq_gca->bases[1]  = sq->base_vert;

  vnq = ROUNDUP(sq->n_edge, VREAL);

  xq = sq->x_edge;
  yq = sq->y_edge;
  wq = sq->w_edge + (9 * vnq) + (sq->n_edge - nq);

  memcpy(sq_gca->xqs + (4 * nq), xq + (sq->n_edge - nq), nq * sizeof(real));
  memcpy(sq_gca->xqs + (5 * nq), xq + vnq + (sq->n_edge - nq), nq * sizeof(real));
  memcpy(sq_gca->yqs + (4 * nq), yq + (sq->n_edge - nq), nq * sizeof(real));
  memcpy(sq_gca->yqs + (5 * nq), yq + vnq + (sq->n_edge - nq), nq * sizeof(real));
  memcpy(sq_gca->wqs + (2 * nq), wq, nq * sizeof(real));

  sq_gca->bases[2]  = sq->base_edge;

  vnq = ROUNDUP(sq->n_id, VREAL);

  xq = sq->x_id;
  yq = sq->y_id;
  wq = sq->w_id + (9 * vnq) + (sq->n_id - nq);

  memcpy(sq_gca->xqs + (6 * nq), xq + (sq->n_id - nq), nq * sizeof(real));
  memcpy(sq_gca->xqs + (7 * nq), xq + vnq + (sq->n_id - nq), nq * sizeof(real));
  memcpy(sq_gca->yqs + (6 * nq), yq + (sq->n_id - nq), nq * sizeof(real));
  memcpy(sq_gca->yqs + (7 * nq), yq + vnq + (sq->n_id - nq), nq * sizeof(real));
  memcpy(sq_gca->wqs + (3 * nq), wq, nq * sizeof(real));

  sq_gca->bases[3]  = sq->base_id;

  sq_gca->offset    = 0;

  sq_gca->buf_xqs   = (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));
  sq_gca->buf_yqs   = (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));
  sq_gca->buf_wqs   = (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));
  sq_gca->buf_bases = (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));

  for(uint i = 0; i < ocl_system.num_devices; ++i)
  {
    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           8 * nq,
                           sizeof(real),
                           sq_gca->xqs,
                           NULL,
                           &sq_gca->buf_xqs[i]);

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           8 * nq,
                           sizeof(real),
                           sq_gca->yqs,
                           NULL,
                           &sq_gca->buf_yqs[i]);

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           4 * nq,
                           sizeof(real),
                           sq_gca->wqs,
                           NULL,
                           &sq_gca->buf_wqs[i]);

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           4,
                           sizeof(real),
                           sq_gca->bases,
                           NULL,
                           &sq_gca->buf_bases[i]);
  }

  return sq_gca;
}

psingquadgca
build_min_vert_from_singquad2d(psingquad2d sq)
{
  psingquadgca sq_gca = (psingquadgca) allocmem(sizeof(singquadgca));

  init_singquadgca(sq_gca);

  const uint nq = sq_gca->nq[0] = sq->n_vert - sq->n_dist;

  sq_gca->xqs   = allocreal(6 * nq);
  sq_gca->yqs   = allocreal(6 * nq);
  sq_gca->wqs   = allocreal(3 * nq);
  sq_gca->bases = allocreal(3);

  uint vnq = ROUNDUP(sq->n_vert, VREAL);

  real *xq = sq->x_vert;
  real *yq = sq->y_vert;
  real *wq = sq->w_vert + (9 * vnq);

  memcpy(sq_gca->xqs,      xq,       nq * sizeof(real));
  memcpy(sq_gca->xqs + nq, xq + vnq, nq * sizeof(real));
  memcpy(sq_gca->yqs,      yq,       nq * sizeof(real));
  memcpy(sq_gca->yqs + nq, yq + vnq, nq * sizeof(real));
  memcpy(sq_gca->wqs,      wq,       nq * sizeof(real));

  sq_gca->bases[0]  = sq->base_vert;

  vnq = ROUNDUP(sq->n_edge, VREAL);

  xq = sq->x_edge;
  yq = sq->y_edge;
  wq = sq->w_edge + (9 * vnq) + (sq->n_edge - sq->n_vert);

  uint offset_src = sq->n_edge - sq->n_vert;

  memcpy(sq_gca->xqs + (2 * nq), xq + offset_src,       nq * sizeof(real));
  memcpy(sq_gca->xqs + (3 * nq), xq + vnq + offset_src, nq * sizeof(real));
  memcpy(sq_gca->yqs + (2 * nq), yq + offset_src,       nq * sizeof(real));
  memcpy(sq_gca->yqs + (3 * nq), yq + vnq + offset_src, nq * sizeof(real));
  memcpy(sq_gca->wqs + (1 * nq), wq,                    nq * sizeof(real));

  sq_gca->bases[1]  = sq->base_edge;

  vnq = ROUNDUP(sq->n_id, VREAL);

  xq = sq->x_id;
  yq = sq->y_id;
  wq = sq->w_id + (9 * vnq) + (sq->n_id - sq->n_vert);

  offset_src = sq->n_id - sq->n_vert;

  memcpy(sq_gca->xqs + (4 * nq), xq + offset_src,       nq * sizeof(real));
  memcpy(sq_gca->xqs + (5 * nq), xq + vnq + offset_src, nq * sizeof(real));
  memcpy(sq_gca->yqs + (4 * nq), yq + offset_src,       nq * sizeof(real));
  memcpy(sq_gca->yqs + (5 * nq), yq + vnq + offset_src, nq * sizeof(real));
  memcpy(sq_gca->wqs + (2 * nq), wq,                    nq * sizeof(real));

  sq_gca->bases[2]  = sq->base_id;

  sq_gca->offset    = -1;

  sq_gca->buf_xqs   = (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));
  sq_gca->buf_yqs   = (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));
  sq_gca->buf_wqs   = (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));
  sq_gca->buf_bases = (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));

  for(uint i = 0; i < ocl_system.num_devices; ++i)
  {
    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           6 * nq,
                           sizeof(real),
                           sq_gca->xqs,
                           NULL,
                           &sq_gca->buf_xqs[i]);

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           6 * nq,
                           sizeof(real),
                           sq_gca->yqs,
                           NULL,
                           &sq_gca->buf_yqs[i]);

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           3 * nq,
                           sizeof(real),
                           sq_gca->wqs,
                           NULL,
                           &sq_gca->buf_wqs[i]);

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           3,
                           sizeof(real),
                           sq_gca->bases,
                           NULL,
                           &sq_gca->buf_bases[i]);
  }

  return sq_gca;
}

psingquadgca
build_min_edge_from_singquad2d(psingquad2d sq)
{
  psingquadgca sq_gca = (psingquadgca) allocmem(sizeof(singquadgca));

  init_singquadgca(sq_gca);

  const uint nq = sq_gca->nq[0] = sq->n_edge - sq->n_vert;

  sq_gca->xqs   = allocreal(4 * nq);
  sq_gca->yqs   = allocreal(4 * nq);
  sq_gca->wqs   = allocreal(2 * nq);
  sq_gca->bases = allocreal(2);

  uint vnq = ROUNDUP(sq->n_edge, VREAL);

  real *xq = sq->x_edge;
  real *yq = sq->y_edge;
  real *wq = sq->w_edge + (9 * vnq);

  uint offset_src = 0;

  memcpy(sq_gca->xqs + (0 * nq), xq + offset_src,       nq * sizeof(real));
  memcpy(sq_gca->xqs + (1 * nq), xq + vnq + offset_src, nq * sizeof(real));
  memcpy(sq_gca->yqs + (0 * nq), yq + offset_src,       nq * sizeof(real));
  memcpy(sq_gca->yqs + (1 * nq), yq + vnq + offset_src, nq * sizeof(real));
  memcpy(sq_gca->wqs + (0 * nq), wq,                    nq * sizeof(real));

  sq_gca->bases[0]  = sq->base_edge;

  vnq = ROUNDUP(sq->n_id, VREAL);

  xq = sq->x_id;
  yq = sq->y_id;
  wq = sq->w_id + (9 * vnq) + (sq->n_id - sq->n_edge);

  offset_src = sq->n_id - sq->n_edge;

  memcpy(sq_gca->xqs + (2 * nq), xq + offset_src,       nq * sizeof(real));
  memcpy(sq_gca->xqs + (3 * nq), xq + vnq + offset_src, nq * sizeof(real));
  memcpy(sq_gca->yqs + (2 * nq), yq + offset_src,       nq * sizeof(real));
  memcpy(sq_gca->yqs + (3 * nq), yq + vnq + offset_src, nq * sizeof(real));
  memcpy(sq_gca->wqs + (1 * nq), wq,                    nq * sizeof(real));

  sq_gca->bases[1]  = sq->base_id;

  sq_gca->offset    = -2;

  sq_gca->buf_xqs   = (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));
  sq_gca->buf_yqs   = (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));
  sq_gca->buf_wqs   = (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));
  sq_gca->buf_bases = (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));

  for(uint i = 0; i < ocl_system.num_devices; ++i)
  {
    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           4 * nq,
                           sizeof(real),
                           sq_gca->xqs,
                           NULL,
                           &sq_gca->buf_xqs[i]);

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           4 * nq,
                           sizeof(real),
                           sq_gca->yqs,
                           NULL,
                           &sq_gca->buf_yqs[i]);

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           2 * nq,
                           sizeof(real),
                           sq_gca->wqs,
                           NULL,
                           &sq_gca->buf_wqs[i]);

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           2,
                           sizeof(real),
                           sq_gca->bases,
                           NULL,
                           &sq_gca->buf_bases[i]);
  }

  return sq_gca;
}

psingquadgca
build_uncommon_from_singquad2d(psingquad2d sq)
{
  psingquadgca sq_gca = (psingquadgca) allocmem(sizeof(singquadgca));

  init_singquadgca(sq_gca);

  sq_gca->nq[0] = sq->n_vert - sq->n_dist;
  sq_gca->nq[1] = sq->n_edge - sq->n_dist;
  sq_gca->nq[2] = sq->n_id   - sq->n_dist;

  sq_gca->xqs   = allocreal(2 * sq_gca->nq[0] + 2 * sq_gca->nq[1] + 2 * sq_gca->nq[2]);
  sq_gca->yqs   = allocreal(2 * sq_gca->nq[0] + 2 * sq_gca->nq[1] + 2 * sq_gca->nq[2]);
  sq_gca->wqs   = allocreal(1 * sq_gca->nq[0] + 1 * sq_gca->nq[1] + 1 * sq_gca->nq[2]);
  sq_gca->bases = allocreal(3);

  uint vnq = ROUNDUP(sq->n_vert, VREAL);

  real *xq = sq->x_vert;
  real *yq = sq->y_vert;
  real *wq = sq->w_vert + (9 * vnq);

  memcpy(sq_gca->xqs + (0 * sq_gca->nq[0]), xq,       sq_gca->nq[0] * sizeof(real));
  memcpy(sq_gca->xqs + (1 * sq_gca->nq[0]), xq + vnq, sq_gca->nq[0] * sizeof(real));
  memcpy(sq_gca->yqs + (0 * sq_gca->nq[0]), yq,       sq_gca->nq[0] * sizeof(real));
  memcpy(sq_gca->yqs + (1 * sq_gca->nq[0]), yq + vnq, sq_gca->nq[0] * sizeof(real));
  memcpy(sq_gca->wqs + (0 * sq_gca->nq[0]), wq,       sq_gca->nq[0] * sizeof(real));

  sq_gca->bases[0]  = sq->base_vert;

  vnq = ROUNDUP(sq->n_edge, VREAL);

  xq = sq->x_edge;
  yq = sq->y_edge;
  wq = sq->w_edge + (9 * vnq);

  memcpy(sq_gca->xqs + (2 * sq_gca->nq[0]),
         xq,
         sq_gca->nq[1] * sizeof(real));
  memcpy(sq_gca->xqs + (2 * sq_gca->nq[0] + sq_gca->nq[1]),
         xq + vnq,
         sq_gca->nq[1] * sizeof(real));
  memcpy(sq_gca->yqs + (2 * sq_gca->nq[0]),
         yq,
         sq_gca->nq[1] * sizeof(real));
  memcpy(sq_gca->yqs + (2 * sq_gca->nq[0] + sq_gca->nq[1]),
         yq + vnq,
         sq_gca->nq[1] * sizeof(real));
  memcpy(sq_gca->wqs + (1 * sq_gca->nq[0]),
         wq,
         sq_gca->nq[1] * sizeof(real));

  sq_gca->bases[1]  = sq->base_edge;

  vnq = ROUNDUP(sq->n_id, VREAL);

  xq = sq->x_id;
  yq = sq->y_id;
  wq = sq->w_id + (9 * vnq);

  memcpy(sq_gca->xqs + (2 * sq_gca->nq[0] + 2 * sq_gca->nq[1]),
         xq,
         sq_gca->nq[2] * sizeof(real));
  memcpy(sq_gca->xqs + (2 * sq_gca->nq[0] + 2 * sq_gca->nq[1] + sq_gca->nq[2]),
         xq + vnq,
         sq_gca->nq[2] * sizeof(real));
  memcpy(sq_gca->yqs + (2 * sq_gca->nq[0] + 2 * sq_gca->nq[1]),
         yq,
         sq_gca->nq[2] * sizeof(real));
  memcpy(sq_gca->yqs + (2 * sq_gca->nq[0] + 2 * sq_gca->nq[1] + sq_gca->nq[2]),
         yq + vnq,
         sq_gca->nq[2] * sizeof(real));
  memcpy(sq_gca->wqs + (sq_gca->nq[0] + sq_gca->nq[1]),
         wq,
         sq_gca->nq[2] * sizeof(real));

  sq_gca->bases[2]  = sq->base_id;

  sq_gca->offset    = -1;

  sq_gca->buf_nq    = (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));
  sq_gca->buf_xqs   = (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));
  sq_gca->buf_yqs   = (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));
  sq_gca->buf_wqs   = (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));
  sq_gca->buf_bases = (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));

  for(uint i = 0; i < ocl_system.num_devices; ++i)
  {
    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           3,
                           sizeof(real),
                           sq_gca->nq,
                           NULL,
                           &sq_gca->buf_nq[i]);

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           2 * sq_gca->nq[0] + 2 * sq_gca->nq[1] + 2 * sq_gca->nq[2],
                           sizeof(real),
                           sq_gca->xqs,
                           NULL,
                           &sq_gca->buf_xqs[i]);

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           2 * sq_gca->nq[0] + 2 * sq_gca->nq[1] + 2 * sq_gca->nq[2],
                           sizeof(real),
                           sq_gca->yqs,
                           NULL,
                           &sq_gca->buf_yqs[i]);

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           sq_gca->nq[0] + sq_gca->nq[1] + sq_gca->nq[2],
                           sizeof(real),
                           sq_gca->wqs,
                           NULL,
                           &sq_gca->buf_wqs[i]);

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           3,
                           sizeof(real),
                           sq_gca->bases,
                           NULL,
                           &sq_gca->buf_bases[i]);
  }

  return sq_gca;
}

void
select_quadrature_singquadgca(psingquadgca sq_gca,
                              const uint   *tv,
                              const uint   *sv,
                              uint         *tp,
                              uint         *sp,
                              real         **x,
                              real         **y,
                              real         **w,
                              real         *base)
{
  uint p = (uint) (tv[0] == sv[0]) + (tv[0] == sv[1]) + (tv[0] == sv[2]) +
                  (tv[1] == sv[0]) + (tv[1] == sv[1]) + (tv[1] == sv[2]) +
                  (tv[2] == sv[0]) + (tv[2] == sv[1]) + (tv[2] == sv[2]);

  if(p > 3)
  {
    fprintf(stderr, "error: Unknown quadrature situation!\n");
    abort();
  }

  tp[0] = 0, tp[1] = 1, tp[2] = 2;
  sp[0] = 0, sp[1] = 1, sp[2] = 2;

  *x    = sq_gca->xqs + (2 * (p + sq_gca->offset) * sq_gca->nq[0]);
  *y    = sq_gca->yqs + (2 * (p + sq_gca->offset) * sq_gca->nq[0]);
  *w    = sq_gca->wqs + ((p + sq_gca->offset) * sq_gca->nq[0]);
  *base = sq_gca->bases[p + sq_gca->offset];

  p = 0;

  for(uint i = 0; i < 3; ++i)
  {
    for(uint j = 0; j < 3; ++j)
    {
      if(tv[i] == sv[j])
      {
        tp[p] = i;
        sp[p] = j;
        p++;
        break;
      }
    }
  }

  uint q = p;

  for(uint i = 0; i < 3; i++)
  {
    uint j = 0;

    for(j = 0; j < q && tv[i] != tv[tp[j]]; j++)
      ;
    if(j == q)
      tp[q++] = i;
  }

  assert(q == 3);

  q = p;

  for(uint i = 0; i < 3; i++)
  {
    uint j = 0;

    for(j = 0; j < q && sv[i] != sv[sp[j]]; j++)
      ;
    if(j == q)
      sp[q++] = i;
  }

  assert(q == 3);

  assert(p <= 3);
}

void
select_quadrature_uncommon_singquadgca(psingquadgca sq_gca,
                                       const uint   *tv,
                                       const uint   *sv,
                                       uint         *tp,
                                       uint         *sp,
                                       uint         *nq,
                                       real         **x,
                                       real         **y,
                                       real         **w,
                                       real         *base)
{
  uint p = (uint) (tv[0] == sv[0]) + (tv[0] == sv[1]) + (tv[0] == sv[2]) +
           (tv[1] == sv[0]) + (tv[1] == sv[1]) + (tv[1] == sv[2]) +
           (tv[2] == sv[0]) + (tv[2] == sv[1]) + (tv[2] == sv[2]);

  if(p > 3)
  {
    fprintf(stderr, "error: Unknown quadrature situation!\n");
    abort();
  }

  tp[0] = 0, tp[1] = 1, tp[2] = 2;
  sp[0] = 0, sp[1] = 1, sp[2] = 2;

  *nq   = sq_gca->nq[p + sq_gca->offset];
  *x    = sq_gca->xqs;
  *y    = sq_gca->yqs;
  *w    = sq_gca->wqs;

  for(uint i = 0; i < (p + sq_gca->offset); ++i)
  {
    *x += 2 * sq_gca->nq[i];
    *y += 2 * sq_gca->nq[i];
    *w += 1 * sq_gca->nq[i];
  }

  *base = sq_gca->bases[p + sq_gca->offset];

  p = 0;

  for(uint i = 0; i < 3; ++i)
  {
    for(uint j = 0; j < 3; ++j)
    {
      if(tv[i] == sv[j])
      {
        tp[p] = i;
        sp[p] = j;
        p++;
        break;
      }
    }
  }

  uint q = p;

  for(uint i = 0; i < 3; i++)
  {
    uint j = 0;

    for(j = 0; j < q && tv[i] != tv[tp[j]]; j++)
      ;
    if(j == q)
      tp[q++] = i;
  }

  assert(q == 3);

  q = p;

  for(uint i = 0; i < 3; i++)
  {
    uint j = 0;

    for(j = 0; j < q && sv[i] != sv[sp[j]]; j++)
      ;
    if(j == q)
      sp[q++] = i;
  }

  assert(q == 3);

  assert(p <= 3);
}
