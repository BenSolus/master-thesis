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

#include "clgreencross.cl"
#include "clgeom.cl"
#include "clquad.cl"
#include "clgcidxinfo.cl"

const char *src_code_strs[] = { clgeom_src, clquad_src, clgcidxinfo_src, clgreencross_src };
const char *kernel_names[1]  = { "fastaddeval_h2matrix_avector_3" };

void
init_gcopencl(pgcopencl gcocl)
{
  assert(gcocl != NULL);

  cl_int res;

  gcocl->num_row_leafs = 0;
  gcocl->size_ridx     = 0;
  gcocl->size_cidx     = 0;
  gcocl->roff          = 0;
  gcocl->coff          = 0;

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
  gcocl->ytoffs                   = NULL;
  gcocl->buf_ytoffs               =
    (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));
  gcocl->buf_yt                   =
    (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));
  gcocl->kernels                  = (cl_kernel *) calloc(1, sizeof(cl_kernel));
  gcocl->row_clusters             = NULL;
  gcocl->cur_h2matrix_leafs       = NULL;

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
                                      ocl_system.max_package_size,
                                      NULL,
                                      &res);

    CL_CHECK(res);

    gcocl->buf_yt[i] = clCreateBuffer(ocl_system.contexts[i],
                                      CL_MEM_READ_WRITE,
                                      ocl_system.max_package_size,
                                      NULL,
                                      &res);

    CL_CHECK(res);

    CL_CHECK(clEnqueueFillBuffer(ocl_system.queues[i],
                                 gcocl->buf_yt[i],
                                 &r_zero,
                                 sizeof(real),
                                 0,
                                 ocl_system.max_package_size,
                                 0,
                                 NULL,
                                 NULL));
  }

  char add_flags[50];

  int n = sprintf(add_flags,
                  "-cl-nv-verbose -DWRK_GRP_SIZE0=%u",
                  64);

  assert(n > 0);

  setup_kernels_fix(4,
                    src_code_strs,
                    add_flags,
                    1,
                    kernel_names,
                    &gcocl->kernels);
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

  if(gcocl->kernels != NULL)
  {
    freemem(gcocl->kernels);
    gcocl->kernels = NULL;
  }

  if(gcocl->row_clusters != NULL)
  {
    freemem(gcocl->row_clusters);
    gcocl->row_clusters = NULL;
  }

  if(gcocl->cur_h2matrix_leafs != NULL)
  {
    freemem(gcocl->cur_h2matrix_leafs);
    gcocl->cur_h2matrix_leafs = NULL;
  }

  if(gcocl->kernels != NULL)
  {
    for(uint k = 0; k < ocl_system.num_devices; ++k)
      for(uint i = 0; i < 1; ++i)
        for(uint j = 0; j < ocl_system.queues_per_device; ++j)
          if(gcocl->kernels[j + k * ocl_system.queues_per_device +
                            i * ocl_system.num_devices *
                            ocl_system.queues_per_device] != NULL)
          {
            CL_CHECK(clReleaseKernel
                       (gcocl->kernels[j + k * ocl_system.queues_per_device +
                                       i * ocl_system.num_devices *
                                       ocl_system.queues_per_device]));
          }

    freemem(gcocl->kernels);
    gcocl->kernels = NULL;
  }
}
