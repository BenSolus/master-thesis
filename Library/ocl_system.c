/* ------------------------------------------------------------
 * This is the file "ocl_system.c" of this master thesis.
 * All rights reserved, Bennet Carstensen 2015
 * ------------------------------------------------------------ */

#include "clsettings.h"
#include "ocl_system.h"

#include <string.h>

void
setup_kernels_fix(const uint n,
                  const char **src_strs,
                  const uint num_kernels,
                  const char **kernel_names,
                  cl_kernel **kernels)
{
  cl_uint   num_devices = ocl_system.num_devices;
  cl_uint   num_queues = ocl_system.queues_per_device;
  const char **srcs;
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

  for (k = 0; k < num_devices; ++k) {

    program = clCreateProgramWithSource(ocl_system.contexts[k], n + 1, srcs,
          size_src_str, &res);
    CL_CHECK(res)

    /****************************************************
     * Compile OpenCL program object.
     ****************************************************/
      res = clBuildProgram(program, 1, ocl_system.devices + k,
         "-Werror -cl-mad-enable -cl-fast-relaxed-math"
#ifdef USE_COMPLEX
         " -DUSE_COMPLEX"
#endif
#ifdef USE_FLOAT
         " -DUSE_FLOAT"
#endif
         , NULL, NULL);

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

  freemem(srcs);
}
