/* ------------------------------------------------------------
 * This is the file "gcopencl.c" of this master thesis.
 * All rights reserved, Bennet Carstensen 2017
 * ------------------------------------------------------------ */

/**
 * @file      gcopencl.c
 * @author    Bennet Carstensen
 * @date      2017
 * @copyright All rights reserved, Bennet Carstensen 2017
 */

#include "gcopencl.h"

#include "ocl_system.h"

void
init_gcopencl(pgcopencl gcocl)
{
  assert(gcocl != NULL);

  cl_int res;

  gcocl->num_row_leafs            = 0;
  gcocl->max_num_h2_leafs_per_row = 0;
  gcocl->size_ridx                = 0;
  gcocl->size_cidx                = 0;
  gcocl->roff                     = 0;
  gcocl->coff                     = 0;

  gcocl->names_row_leafs          = NULL;
  gcocl->num_h2_leafs_per_row     = NULL;
  gcocl->buf_num_h2_leafs_per_row =
    (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));
  gcocl->workload_per_row         = NULL;
  gcocl->idx_off                  = NULL;
  gcocl->buf_idx_off              =
    (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));
  gcocl->col_names_per_row        = NULL;
  gcocl->h2_leafs_per_row         = NULL;
  gcocl->ridx_sizes               = NULL;
  gcocl->buf_ridx_sizes           =
    (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));
  gcocl->ridx_off                 = NULL;
  gcocl->buf_ridx_off             =
    (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));
  gcocl->host_ridx                = NULL;
  gcocl->buf_ridx                 =
    (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));
  gcocl->cidx_sizes               = NULL;
  gcocl->buf_cidx_sizes           =
    (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));
  gcocl->cidx_off                 = NULL;
  gcocl->host_cidx_off            = NULL;
  gcocl->buf_cidx_off             =
    (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));
  gcocl->host_cidx                = NULL;
  gcocl->buf_cidx                 =
    (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));
  gcocl->xtoffs                   = NULL;
  gcocl->host_xtoffs              = NULL;
  gcocl->buf_xtoffs               =
    (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));
  gcocl->buf_xt                   =
    (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));
  gcocl->xt_events                =
    (cl_event *) calloc(ocl_system.num_devices, sizeof(cl_event));
  gcocl->ytoffs                   = NULL;
  gcocl->buf_ytoffs               =
    (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));
  gcocl->buf_yt                   =
    (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));

  gcocl->yt_events                =
    (cl_event *) calloc(ocl_system.num_devices, sizeof(cl_event));

  for(uint i = 0; i < ocl_system.num_devices; ++i)
  {
    gcocl->buf_ridx[i] = clCreateBuffer(ocl_system.contexts[i],
                                        CL_MEM_READ_ONLY,
                                        ocl_system.max_package_size,
                                        NULL,
                                        &res);

    CL_CHECK(res);

    gcocl->buf_cidx[i] = clCreateBuffer(ocl_system.contexts[i],
                                        CL_MEM_READ_ONLY,
                                        ocl_system.max_package_size,
                                        NULL,
                                        &res);

    CL_CHECK(res);

    gcocl->buf_xt[i] = clCreateBuffer(ocl_system.contexts[i],
                                      CL_MEM_READ_ONLY,
                                      2 * ocl_system.max_package_size,
                                      NULL,
                                      &res);

    CL_CHECK(res);

    gcocl->buf_yt[i] = clCreateBuffer(ocl_system.contexts[i],
                                      CL_MEM_READ_WRITE,
                                      2 * ocl_system.max_package_size,
                                      NULL,
                                      &res);

    CL_CHECK(res);

//    CL_CHECK(clEnqueueFillBuffer(ocl_system.queues[i],
//                                 gcocl->buf_yt[i],
//                                 &r_zero,
//                                 sizeof(real),
//                                 0,
//                                 ocl_system.max_package_size,
//                                 0,
//                                 NULL,
//                                 NULL));
  }
}

void
uninit_gcopencl(pgcopencl gcocl)
{
  if(gcocl->names_row_leafs != NULL)
  {
    freemem(gcocl->names_row_leafs);
    gcocl->names_row_leafs = NULL;
  }

  if(gcocl->num_h2_leafs_per_row != NULL)
  {
    freemem(gcocl->num_h2_leafs_per_row);
    gcocl->num_h2_leafs_per_row = NULL;
  }

  if(gcocl->buf_num_h2_leafs_per_row != NULL)
  {
    for(uint i = 0; i < ocl_system.num_devices; ++i)
      if(gcocl->buf_num_h2_leafs_per_row[i] != NULL)
      {
        CL_CHECK(clReleaseMemObject(gcocl->buf_num_h2_leafs_per_row[i]));
        gcocl->buf_num_h2_leafs_per_row[i] = NULL;
      }

    freemem(gcocl->buf_num_h2_leafs_per_row);
    gcocl->buf_num_h2_leafs_per_row = NULL;
  }

  if(gcocl->workload_per_row != NULL)
  {
    freemem(gcocl->workload_per_row);
    gcocl->workload_per_row = NULL;
  }

  if(gcocl->idx_off != NULL)
  {
    freemem(gcocl->idx_off);
    gcocl->idx_off = NULL;
  }

  if(gcocl->buf_idx_off != NULL)
  {
    for(uint i = 0; i < ocl_system.num_devices; ++i)
      if(gcocl->buf_idx_off[i] != NULL)
      {
        CL_CHECK(clReleaseMemObject(gcocl->buf_idx_off[i]));
        gcocl->buf_idx_off[i] = NULL;
      }

    freemem(gcocl->buf_idx_off);
    gcocl->buf_idx_off = NULL;
  }

  if(gcocl->h2_leafs_per_row != NULL)
  {
    for(uint i = 0; i < gcocl->num_row_leafs; ++i)
      if(gcocl->h2_leafs_per_row[i] != NULL)
      {
        freemem(gcocl->h2_leafs_per_row[i]);
        gcocl->h2_leafs_per_row[i] = NULL;
      }

    freemem(gcocl->h2_leafs_per_row);
    gcocl->h2_leafs_per_row = NULL;
  }

  if(gcocl->col_names_per_row != NULL)
  {
    for(uint i = 0; i < gcocl->num_row_leafs; ++i)
      if(gcocl->col_names_per_row[i] != NULL)
      {
        freemem(gcocl->col_names_per_row[i]);
        gcocl->col_names_per_row[i] = NULL;
      }

    freemem(gcocl->col_names_per_row);
    gcocl->col_names_per_row = NULL;
  }

  if(gcocl->buf_ridx_sizes != NULL)
  {
    freemem(gcocl->ridx_sizes);
    gcocl->buf_ridx_sizes = NULL;
  }

  if(gcocl->buf_ridx_sizes != NULL)
  {
    for(uint i = 0; i < ocl_system.num_devices; ++i)
      if(gcocl->buf_ridx_sizes[i] != NULL)
      {
        CL_CHECK(clReleaseMemObject(gcocl->buf_ridx_sizes[i]));
        gcocl->buf_ridx_sizes[i] = NULL;
      }

    freemem(gcocl->buf_ridx_sizes);
    gcocl->buf_ridx_sizes = NULL;
  }

  if(gcocl->ridx_off != NULL)
  {
    freemem(gcocl->ridx_off);
    gcocl->ridx_off = NULL;
  }

  if(gcocl->buf_ridx_off != NULL)
  {
    for(uint i = 0; i < ocl_system.num_devices; ++i)
      if(gcocl->buf_ridx_off[i] != NULL)
      {
        CL_CHECK(clReleaseMemObject(gcocl->buf_ridx_off[i]));
        gcocl->buf_ridx_off[i] = NULL;
      }

    freemem(gcocl->buf_ridx_off);
    gcocl->buf_ridx_off = NULL;
  }

  if(gcocl->host_ridx != NULL)
  {
    freemem(gcocl->host_ridx);
    gcocl->host_ridx = NULL;
  }

  if(gcocl->buf_ridx != NULL)
  {
    for(uint i = 0; i < ocl_system.num_devices; ++i)
      if(gcocl->buf_ridx[i] != NULL)
      {
        CL_CHECK(clReleaseMemObject(gcocl->buf_ridx[i]));
        gcocl->buf_ridx[i] = NULL;
      }

    freemem(gcocl->buf_ridx);
    gcocl->buf_ridx = NULL;
  }

  if(gcocl->cidx_sizes != NULL)
  {
    freemem(gcocl->cidx_sizes);
    gcocl->cidx_sizes = NULL;
  }

  if(gcocl->buf_cidx_sizes != NULL)
  {
    for(uint i = 0; i < ocl_system.num_devices; ++i)
      if(gcocl->buf_cidx_sizes[i] != NULL)
      {
        CL_CHECK(clReleaseMemObject(gcocl->buf_cidx_sizes[i]));
        gcocl->buf_cidx_sizes[i] = NULL;
      }

    freemem(gcocl->buf_cidx_sizes);
    gcocl->buf_cidx_sizes = NULL;
  }

  if(gcocl->cidx_off != NULL)
  {
    for(uint i = 0; i < gcocl->num_row_leafs; ++i)
      if(gcocl->cidx_off[i] != NULL)
      {
        freemem(gcocl->cidx_off[i]);
        gcocl->cidx_off[i] = NULL;
      }

    freemem(gcocl->cidx_off);
    gcocl->cidx_off = NULL;
  }

  if(gcocl->host_cidx_off != NULL)
  {
    freemem(gcocl->host_cidx_off);
    gcocl->host_cidx_off = NULL;
  }

  if(gcocl->buf_cidx_off != NULL)
  {
    for(uint i = 0; i < ocl_system.num_devices; ++i)
      if(gcocl->buf_cidx_off[i] != NULL)
      {
        CL_CHECK(clReleaseMemObject(gcocl->buf_cidx_off[i]));
        gcocl->buf_cidx_off[i] = NULL;
      }

    freemem(gcocl->buf_cidx_off);
    gcocl->buf_cidx_off = NULL;
  }

  if(gcocl->host_cidx != NULL)
  {
    freemem(gcocl->host_cidx);
    gcocl->host_cidx = NULL;
  }

  if(gcocl->buf_cidx != NULL)
  {
    for(uint i = 0; i < ocl_system.num_devices; ++i)
      if(gcocl->buf_cidx[i] != NULL)
      {
        CL_CHECK(clReleaseMemObject(gcocl->buf_cidx[i]));
        gcocl->buf_cidx[i] = NULL;
      }

    freemem(gcocl->buf_cidx);
    gcocl->buf_cidx = NULL;
  }

  if(gcocl->xtoffs != NULL)
  {
    for(uint i = 0; i < gcocl->num_row_leafs; ++i)
      if(gcocl->xtoffs[i] != NULL)
      {
        freemem(gcocl->xtoffs[i]);
        gcocl->xtoffs[i] = NULL;
      }

    freemem(gcocl->xtoffs);
    gcocl->xtoffs = NULL;
  }

  if(gcocl->host_xtoffs != NULL)
  {
    freemem(gcocl->host_xtoffs);
    gcocl->host_xtoffs = NULL;
  }

  if(gcocl->buf_xtoffs != NULL)
  {
    for(uint i = 0; i < ocl_system.num_devices; ++i)
      if(gcocl->buf_xtoffs[i] != NULL)
      {
        CL_CHECK(clReleaseMemObject(gcocl->buf_xtoffs[i]));
        gcocl->buf_xtoffs[i] = NULL;
      }

    freemem(gcocl->buf_xtoffs);
    gcocl->buf_xtoffs = NULL;
  }

  if(gcocl->buf_xt != NULL)
  {
    for(uint i = 0; i < ocl_system.num_devices; ++i)
      CL_CHECK(clReleaseMemObject(gcocl->buf_xt[i]));

    freemem(gcocl->buf_xt);
    gcocl->buf_xt = NULL;
  }

  if(gcocl->xt_events != NULL)
  {
    for(uint i = 0; i < ocl_system.num_devices; ++i)
      CL_CHECK(clReleaseEvent(gcocl->xt_events[i]));

    freemem(gcocl->xt_events);
    gcocl->xt_events = NULL;
  }

  if(gcocl->ytoffs != NULL)
  {
    freemem(gcocl->ytoffs);
    gcocl->ytoffs = NULL;
  }

  if(gcocl->buf_ytoffs != NULL)
  {
    for(uint i = 0; i < ocl_system.num_devices; ++i)
      if(gcocl->buf_ytoffs[i] != NULL)
      {
        CL_CHECK(clReleaseMemObject(gcocl->buf_ytoffs[i]));
        gcocl->buf_ytoffs[i] = NULL;
      }

    freemem(gcocl->buf_ytoffs);
    gcocl->buf_ytoffs = NULL;
  }

  if(gcocl->buf_yt != NULL)
  {
    for(uint i = 0; i < ocl_system.num_devices; ++i)
      if(gcocl->buf_yt[i] != NULL)
      {
        CL_CHECK(clReleaseMemObject(gcocl->buf_yt[i]));
        gcocl->buf_yt[i] = NULL;
      }

    freemem(gcocl->buf_yt);
    gcocl->buf_yt = NULL;
  }

  if(gcocl->yt_events != NULL)
  {
    for(uint i = 0; i < ocl_system.num_devices; ++i)
      CL_CHECK(clReleaseEvent(gcocl->yt_events[i]));

    freemem(gcocl->yt_events);
    gcocl->yt_events = NULL;
  }
}

void
del_gcopencl(pgcopencl gcocl)
{
  if(gcocl != NULL)
  {
    uninit_gcopencl(gcocl);

    freemem(gcocl);
  }
}

void
del_gcopencls(const uint num_gcocls, pgcopencl *gcocls)
{
  if(gcocls != NULL)
  {
    for(uint i = 0; i < num_gcocls; ++i)
    {
      del_gcopencl(gcocls[i]);
      gcocls[i] = NULL;
    }

    freemem(gcocls);
  }
}

size_t
getsize_gcopencl(pcgcopencl gcocl)
{
  size_t sz = 0;

  sz = sizeof(gcocl);

  /* names_row_leafs, num_h2_leafs_per_row, workload_per_row, idx_off,
   * ridx_sizes, ridx_off, ytoffs */
  sz += gcocl->num_row_leafs * 7 * sizeof(uint);

  /* host_ridx, host_cidx */
  sz += (gcocl->size_ridx + gcocl->size_cidx) * sizeof(uint);

  /* host_cidx_off, host_xtoffs */
  sz += (gcocl->idx_off[gcocl->num_row_leafs - 1] +
         gcocl->num_h2_leafs_per_row[gcocl->num_row_leafs - 1]) *
         2 * sizeof(uint);

  /* buf_num_h2_leafs_per_row, buf_idx_off, buf_ridx_size, buf_ridx_off,
   * buf_ridx, buf_cidx_size, buf_cidx_off, buf_cidx, buf_xtoffs, buf_xt,
   * buf_ytoffs, buf_yt */
  sz += 12 * ocl_system.num_devices * sizeof(cl_mem);

  /* yt_events */
  sz += ocl_system.num_devices * sizeof(cl_event);

  for(uint i = 0; i < gcocl->num_row_leafs; ++i)
  {
    /* h2_leafs_per_row */
    sz += gcocl->num_h2_leafs_per_row[i] * sizeof(ph2matrix);

    /* col_names_per_row, cidx_off, xtoffs */
    sz += 3 * gcocl->num_h2_leafs_per_row[i] * sizeof(uint);
  }

  return sz;
}

static bool
are_equal_clusters(pccluster a, pccluster b)
{
  if(a->size != b->size)
    return false;

  for(uint i = 0; i < a->size; ++i)
    if(a->idx[i] != b->idx[i])
      return false;

  return true;
}

static void
iterate_recursively_h2matrix(ph2matrix  H2,
                             const uint xtoff,
                             const uint ytoff,
                             const uint level,
                             uint       *num_gcocls,
                             pgcopencl  **gcocls)
{
  if(level >= *num_gcocls)
  {
    *num_gcocls += 1;
    *gcocls = realloc(*gcocls, *num_gcocls * sizeof(gcopencl));

    if(*gcocls == NULL)
    {
      fprintf(stderr, "error: failed to realloc array of gcopencls!");
      exit(1);
    }

    pgcopencl gcocl = *gcocls[*num_gcocls - 1];

    init_gcopencl(gcocl);

    gcocl->num_row_leafs           = 1;
    gcocl->num_h2_leafs_per_row    = calloc(1, sizeof(uint));
    gcocl->num_h2_leafs_per_row[0] = 1;
    gcocl->h2_leafs_per_row        = calloc(1, sizeof(ph2matrix*));
    gcocl->h2_leafs_per_row[0]     = calloc(1, sizeof(ph2matrix));
    gcocl->h2_leafs_per_row[0][0]  = H2;
    gcocl->xtoffs                  = calloc(1, sizeof(uint*));
    gcocl->xtoffs[0]               = calloc(1, sizeof(uint));
    gcocl->xtoffs[0][0]            = xtoff;
    gcocl->ytoffs[0]               = ytoff;
  }
  else
  {
    pgcopencl gcocl = *gcocls[level];

    bool row_already_stored = false;

    for(uint i = 0; i < gcocl->num_row_leafs; ++i)
      if(are_equal_clusters(H2->rb->t, gcocl->h2_leafs_per_row[level][i]->rb->t))
      {
        row_already_stored = true;

        const uint j = gcocl->num_h2_leafs_per_row[i];

        gcocl->num_h2_leafs_per_row[i] += 1;

        gcocl->h2_leafs_per_row[i]      =
          realloc(gcocl->h2_leafs_per_row[i],
                  gcocl->num_h2_leafs_per_row[i] *sizeof(ph2matrix));

        if(gcocl->h2_leafs_per_row[i] == NULL)
        {
          fprintf(stderr, "error: failed to realloc array of H2-matrices!");
          exit(1);
        }

        gcocl->h2_leafs_per_row[i][j]   = H2;

        gcocl->xtoffs[i]                =
          realloc(gcocl->xtoffs[i],
                  gcocl->num_h2_leafs_per_row[i] * sizeof(uint));

        if(gcocl->xtoffs[i] == NULL)
        {
          fprintf(stderr, "error: failed to realloc array of offsets!");
          exit(1);
        }

        gcocl->xtoffs[i][j]             = xtoff;

      }

    if(!row_already_stored)
    {

    }
  }

  if(H2->son)
  {
    pcclusterbasis rb = H2->rb;
    pcclusterbasis cb = H2->cb;

    const uint rsons  = H2->rsons;
    const uint csons  = H2->csons;

    uint xtoff1        = cb->k;

    for(uint j = 0; j < csons; ++j)
    {
      assert(csons == 1 || cb->sons > 0);

      uint ytoff1 = rb->k;

      for(uint i = 0; i < rsons; ++i)
      {
        assert(rsons == 1 || rb->sons > 0);

        iterate_recursively_h2matrix
          (H2->son[i + j * rsons],
           xtoff + xtoff1 - (cb->sons > 0 ? 0 : cb->k),
           ytoff + ytoff1 - (rb->sons > 0 ? 0 : rb->k),
           level + 1,
           num_gcocls,
           gcocls);

        ytoff1 += (rb->sons > 0 ? rb->son[i]->ktree : rb->t->size);
      }

      assert(ytoff1 == rb->ktree);

      xtoff1 += (cb->sons > 0 ? cb->son[j]->ktree : cb->t->size);
    }

    assert(xtoff1 == cb->ktree);
  }
}

//void
//get_ocl_tree_informations_gcocl(ph2matrix H2,
//                                uint      *num_levels,
//                                pgcopencl **gcocls)
//{
//  del_gcopencls(*num_levels, *gcocls);
//
//  *num_levels = 0;
//  *gcocls     = (pgcopencl *) calloc(0, sizeof(gcopencl));
//
//  iterate_recursively_h2matrix(H2, 0, 0, 0, num_levels, gcocls);
//
//}