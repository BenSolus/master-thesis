/* ------------------------------------------------------------
 * This is the file "ocl_system.c" of this master thesis.
 * All rights reserved, Bennet Carstensen 2015
 * ------------------------------------------------------------ */

#include "clsettings.h"
#include "ocl_system.h"

#include <string.h>

const char *default_flags = "-Werror -cl-mad-enable -cl-fast-relaxed-math"
#ifdef USE_COMPLEX
                            " -DUSE_COMPLEX"
#endif
#ifdef USE_FLOAT
                            " -DUSE_FLOAT"
#endif
;

void
setup_kernels_fix(const uint n,
                  const char **src_strs,
                  const char *add_flags,
                  const uint num_kernels,
                  const char **kernel_names,
                  cl_kernel **kernels)
{
  cl_uint   num_devices = ocl_system.num_devices;
  cl_uint   num_queues = ocl_system.queues_per_device;
  const char **srcs;
  char      *flags;
  size_t    size_src_str[n + 1];
  cl_int    res;
  cl_uint    i, j, k;
  cl_program program;

  size_src_str[0] = strlen(clsettings_src);
  for(i = 0;  i < n; ++i)
    size_src_str[i + 1] = strlen(src_strs[i]);
  srcs = (const char **) allocmem((n + 1) * sizeof(char *));
  srcs[0] = clsettings_src;
  for(i = 0; i < n; ++i)
    srcs[i + 1] = src_strs[i];

  /****************************************************
   * Create OpenCL program object.
   ****************************************************/

  *kernels =
    (cl_kernel *) allocmem(num_queues * num_kernels * num_devices *
                           sizeof(cl_kernel));

  flags = (char *) malloc(strlen(default_flags) + strlen(add_flags) + 2);
  strcpy(flags, default_flags);
  strcat(flags, " ");
  strcat(flags, add_flags);

  for (k = 0; k < num_devices; ++k) {

    program = clCreateProgramWithSource(ocl_system.contexts[k], n + 1, srcs,
          size_src_str, &res);
    CL_CHECK(res)

    /****************************************************
     * Compile OpenCL program object.
     ****************************************************/
      res = clBuildProgram(program,
                           1,
                           ocl_system.devices + k,
                           flags,
                           NULL,
                           NULL);

    if (res != CL_SUCCESS) {
      cl_device_id dev_id = ocl_system.devices[k];
      size_t    len;
      char      buffer[204800];
      cl_build_status bldstatus;
      printf("\nError %d: Failed to build program executable [ %s ]\n", res,
       get_error_string(res));
      res = clGetProgramBuildInfo(program, dev_id, CL_PROGRAM_BUILD_STATUS,
          sizeof(bldstatus), (void *) &bldstatus,
          &len);
      if (res != CL_SUCCESS) {
        printf("Build Status error %d: %s\n", res, get_error_string(res));
        exit(1);
      }
      if (bldstatus == CL_BUILD_SUCCESS)
        printf("Build Status: CL_BUILD_SUCCESS\n");
      if (bldstatus == CL_BUILD_NONE)
        printf("Build Status: CL_BUILD_NONE\n");
      if (bldstatus == CL_BUILD_ERROR)
        printf("Build Status: CL_BUILD_ERROR\n");
      if (bldstatus == CL_BUILD_IN_PROGRESS)
        printf("Build Status: CL_BUILD_IN_PROGRESS\n");
      res = clGetProgramBuildInfo(program, dev_id, CL_PROGRAM_BUILD_OPTIONS,
          sizeof(buffer), buffer, &len);
      if (res != CL_SUCCESS) {
        printf("Build Options error %d: %s\n", res, get_error_string(res));
        exit(1);
      }
      printf("Build Options: %s\n", buffer);
      res = clGetProgramBuildInfo(program, dev_id, CL_PROGRAM_BUILD_LOG,
          sizeof(buffer), buffer, &len);
      if (res != CL_SUCCESS) {
        printf("Build Log error %d: %s\n", res, get_error_string(res));
        exit(1);
      }
      printf("Build Log:\n%s\n", buffer);
      abort();
    }

    cl_device_id dev_id = ocl_system.devices[k];
    size_t       len;
    char         buffer[204800];

    if(clGetProgramBuildInfo(program, dev_id, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len) != CL_SUCCESS)
    {
      printf("Build Log error %d: %s\n", res, get_error_string(res));
      exit(1);
    }

    printf("Build Log:\n%s\n", buffer);

    /****************************************************
     * Build kernels
     ****************************************************/

    for (i = 0; i < num_kernels; ++i) {
      for (j = 0; j < num_queues; ++j) {
        (*kernels)[j + k * num_queues + i * num_devices * num_queues] =
          clCreateKernel(program, kernel_names[i], &res);
        CL_CHECK(res);
      }
    }

    clReleaseProgram(program);
  }

  freemem(flags);

  freemem(srcs);
}

void
create_and_fill_buffer(cl_context       context,
                       cl_mem_flags     flags,
                       cl_command_queue queue,
                       size_t           num,
                       size_t           size,
                       void             *src,
                       cl_event         *event,
                       cl_mem           *buffer)
{
  cl_int res;

  *buffer = clCreateBuffer(context, flags, num * size, NULL, &res);

  CL_CHECK(res);

  CL_CHECK(clEnqueueWriteBuffer(queue,
                                *buffer,
                                false,
                                0,
                                num * size, src,
                                0,
                                NULL,
                                event));
}

cl_kernel*
delete_kernels(const uint num_kernels, cl_kernel** kernels)
{
  if(*kernels != NULL)
  {
    for(uint k = 0; k < ocl_system.num_devices; ++k)
    {
      for (uint i = 0; i < num_kernels; ++i)
      {
        for (uint j = 0; j < ocl_system.queues_per_device; ++j)
        {
          if(*kernels[j + k * ocl_system.queues_per_device +
                      i * ocl_system.num_devices *
                      ocl_system.queues_per_device] != NULL)
          {
            CL_CHECK(clReleaseKernel
                       (*kernels[j + k * ocl_system.queues_per_device +
                                 i * ocl_system.num_devices *
                                 ocl_system.queues_per_device]));
          }
        }
      }
    }

    freemem(*kernels);
    *kernels = NULL;
  }

  return *kernels;
}