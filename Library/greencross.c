/* ------------------------------------------------------------
 * This is the file "greencross.c" of this master thesis.
 * All rights reserved, Bennet Carstensen 2017
 * ------------------------------------------------------------ */

#include "greencross.h"

#include "ocl_system.h"

#include "clamatrix.cl"
#include "greencross.cl"

static greencross_ocl ocl_gc =
{
  .num_kernels = 0
};

void
init_greencross_ocl()
{
  const char *kernel_names[] =
  {
    "eval_leaf_forward"
  };

  ocl_gc.num_kernels = sizeof(kernel_names) / sizeof(kernel_names[0]);

  cl_device_id *devices;
  cl_int       res;
  cl_uint      ndevices;
  cl_uint      nqueues;
  cl_uint      nthreads;

  const char *src_strs[2] = { clamatrix_src, greencross_src };

  get_opencl_devices(&devices, &ndevices);

  set_opencl_devices(devices, ndevices, 1);

  ocl_gc.kernels = (cl_kernel *) allocmem(sizeof(cl_kernel) *
                                          ocl_gc.num_kernels);

  setup_kernels_fix(2, src_strs, 1, kernel_names, &ocl_gc.kernels);

  nqueues  = ocl_system.queues_per_device;
  nthreads = ndevices * nqueues;

  ocl_gc.mem_coefs = (cl_mem *) allocmem(nthreads * sizeof(cl_mem));

  for(uint i = 0; i < ndevices; ++i)
    for(uint j = 0; j < nthreads; ++j)
    {
      ocl_gc.mem_coefs[j + i * nqueues] = clCreateBuffer
      (
        ocl_system.contexts[i],
        CL_MEM_READ_WRITE,
        ocl_system.max_package_size,
        NULL,
        &res
      );

      CL_CHECK(res);
    }
}

void
uninit_greencross_ocl()
{
  cl_uint ndevices = ocl_system.num_devices;
  cl_uint nqueues  = ocl_system.queues_per_device;

  for(uint i = 0; i < ndevices; ++i)
    for(uint j = 0; j < nqueues; ++j)
    {
      clReleaseMemObject(ocl_gc.mem_coefs[j + i * nqueues]);
    }

  freemem(ocl_gc.mem_coefs);
  ocl_gc.mem_coefs = NULL;
}

void
forward_greencross(pclusterbasis cb, pcavector x, pavector xt)
{
  if(cb->sons)
    for(uint i = 0; i < cb->sons; ++i)
      forward_greencross(cb->son[i], x, xt);
  else
  {
    cl_command_queue queue;
    cl_kernel        kernel;

    cl_int res;

    size_t global_off[]  = {  0,  0 };
    size_t local_size[]  = { 16, 16 };
    size_t global_size[] = { 16, 16 };

    queue = ocl_system.queues[0];
    kernel = ocl_gc.kernels[0];

    res = clEnqueueNDRangeKernel(queue,
                                 kernel,
                                 2,
                                 global_off,
                                 global_size,
                                 local_size,
                                 0,
                                 NULL,
                                 NULL);

    CL_CHECK(res);
  }
}

// INLINE_PREFIX void
// get_leaf_infos_clusterbasis(pclusterbasis cb, uint* n, uint* sizeAll)
// {
//   if(cb->sons)
//     for(uint i = 0; i < cb->sons; ++i)
//       get_leaf_infos_clusterbasis(cb->son[i], n, sizeAll);
//   else
//   {
//     ++*n;
//     *sizeAll += cb->V.rows + cb->V.cols;
//   }
// }
//
// INLINE_PREFIX void
// extract_leaf_content_clusterbasis(pclusterbasis cb,
//                                   real* Vs,
//                                   uint* rows,
//                                   uint* cols,
//                                   uint* i,
//                                   uint* off)
// {
//   if(cb->sons)
//     for(uint j = 0; j < cb->sons; ++j)
//       extract_leaf_content_clusterbasis(cb->son[j], Vs, rows, cols, i, off);
//   else
//   {
//     const uint nrows = cb->V.rows;
//     const uint ncols = cb->V.cols;
//
//     memcpy(Vs + *off, cb->V.a, nrows * ncols);
//
//     off += nrows * ncols;
//
//     rows[*i] = nrows;
//     cols[*i] = ncols;
//
//     ++*i;
//   }
// }
//
// void
// eval_leafs_forward(pclusterbasis cb, pcavector x, pavector xt)
// {
//
//   const char *kernel_names[] = { "eval_leaf_forward" };
//
//   cl_mem           bufCoefsV;
//   cl_command_queue queue;
//   cl_device_id     *devices;
//   cl_event         h2d[2], calc;
//   cl_kernel        *kernels;
//   cl_int           res;
//   cl_uint          ndevices;
//
//   get_opencl_devices(&devices, &ndevices);
//
//   set_opencl_devices(devices, ndevices, 1);
//
//   kernels = (cl_kernel *) allocmem(sizeof(cl_kernel));
//
//   setup_kernels(greencross_src, 1, kernel_names, &kernels);
//
//   queue = ocl_system.queues[0];
//
//   bufCoefsV = clCreateBuffer(ocl_system.contexts[0],
//                              CL_MEM_READ_ONLY,
//                              (size_t) sizeof(real) * cb->V.rows * cb->V.cols,
//                              cb->V.a,
//                              &res);
//
//   CL_CHECK(res);
//
//
//   size_t global_offset[] = {  0,  0 };
//   size_t local_size[]    = { 32, 32 };
//   size_t global_size[]   = { 32, 32 };
//
//   res = clEnqueueNDRangeKernel(ocl_system.queues[1],
//                                kernels[1],
//                                2,
//                                global_offset,
//                                global_size,
//                                local_size,
//                                1,
//                                h2d,
//                                &calc);
//
//   CL_CHECK(res);
// }
