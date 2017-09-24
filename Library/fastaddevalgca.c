/* ------------------------------------------------------------
 * This is the file "fastaddevalgca.c" of this master thesis.
 * All rights reserved, Bennet Carstensen 2017
 * ------------------------------------------------------------ */

/**
 * @file      fastaddevalgca.c
 * @author    Bennet Carstensen
 * @date      2017
 * @copyright All rights reserved, Bennet Carstensen 2017
 */

#include "fastaddevalgca.h"

#include "ocl_system.h"

void
init_fastaddevalgca(pfastaddevalgca feval)
{
  if(feval != NULL)
  {
    feval->num_wrk_pkgs = 0;

    feval->buf_xt             = NULL;
    feval->events_xt          = NULL;
#ifndef USE_OPENMP
    feval->buf_yt             = NULL;
#else
    feval->buf_yt_ff          = NULL;
    feval->buf_yt_nf_common   = NULL;
    feval->buf_yt_nf_min_vert = NULL;
    feval->buf_yt_nf_min_edge = NULL;
#endif
    feval->events_yt          = NULL;
  }
}

void
uninit_fastaddevalgca(pfastaddevalgca feval)
{
  if(feval != NULL)
  {
    if(feval->buf_xt != NULL)
    {
      for(uint i = 0; i < ocl_system.num_devices; ++i)
        CL_CHECK(clReleaseMemObject(feval->buf_xt[i]));

      freemem(feval->buf_xt);
    }

    if(feval->events_xt != NULL)
    {
      for(uint i = 0; i < ocl_system.num_devices; ++i)
        CL_CHECK(clReleaseEvent(feval->events_xt[i]));

      freemem(feval->events_xt);
    }

#ifndef USE_OPENMP
    if(feval->buf_yt != NULL)
    {
      for(uint i = 0; i < ocl_system.num_devices; ++i)
        CL_CHECK(clReleaseMemObject(feval->buf_yt[i]));

      freemem(feval->buf_yt);
    }
#else
    if(feval->buf_yt_ff != NULL)
    {
      for(uint i = 0; i < ocl_system.num_devices; ++i)
      CL_CHECK(clReleaseMemObject(feval->buf_yt_ff[i]));

      freemem(feval->buf_yt_ff);
    }

    if(feval->buf_yt_nf_common != NULL)
    {
      for(uint i = 0; i < ocl_system.num_devices; ++i)
      CL_CHECK(clReleaseMemObject(feval->buf_yt_nf_common[i]));

      freemem(feval->buf_yt_nf_common);
    }

    if(feval->buf_yt_nf_min_vert != NULL)
    {
      for(uint i = 0; i < ocl_system.num_devices; ++i)
      CL_CHECK(clReleaseMemObject(feval->buf_yt_nf_min_vert[i]));

      freemem(feval->buf_yt_nf_min_vert);
    }

    if(feval->buf_yt_nf_min_edge != NULL)
    {
      for(uint i = 0; i < ocl_system.num_devices; ++i)
      CL_CHECK(clReleaseMemObject(feval->buf_yt_nf_min_edge[i]));

      freemem(feval->buf_yt_nf_min_edge);
    }
#endif

    if(feval->events_yt != NULL)
    {
      for(uint i = 0; i < ocl_system.num_devices; ++i)
        CL_CHECK(clReleaseEvent(feval->events_yt[i]));

      freemem(feval->events_yt);
    }

    init_fastaddevalgca(feval);
  }
}

pfastaddevalgca
new_fastaddevalgca(pcclusterbasis rb,
                   pcclusterbasis cb,
                   const uint     num_wrk_pkgs)
{
  pfastaddevalgca feval = (pfastaddevalgca) allocmem(sizeof(fastaddevalgca));

  init_fastaddevalgca(feval);

  feval->num_wrk_pkgs       = num_wrk_pkgs;
  feval->buf_xt             =
    (cl_mem *)   calloc(ocl_system.num_devices, sizeof(cl_mem));
  feval->events_xt          =
    (cl_event *) calloc(feval->num_wrk_pkgs,    sizeof(cl_event));
#ifndef USE_OPENMP
  feval->buf_yt             =
    (cl_mem *)   calloc(ocl_system.num_devices, sizeof(cl_mem));
  feval->events_yt          =
    (cl_event *) calloc(feval->num_wrk_pkgs,    sizeof(cl_event));
#else
  feval->buf_yt_ff          =
    (cl_mem *)   calloc(ocl_system.num_devices, sizeof(cl_mem));
  feval->buf_yt_nf_common   =
    (cl_mem *)   calloc(ocl_system.num_devices, sizeof(cl_mem));
  feval->buf_yt_nf_min_vert =
    (cl_mem *)   calloc(ocl_system.num_devices, sizeof(cl_mem));
  feval->buf_yt_nf_min_edge =
    (cl_mem *)   calloc(ocl_system.num_devices, sizeof(cl_mem));
  feval->events_yt          =
    (cl_event *) calloc(4 * feval->num_wrk_pkgs,    sizeof(cl_event));
#endif

  for(uint i = 0; i < ocl_system.num_devices; ++i)
  {
    cl_int res;

    feval->buf_xt[i] = clCreateBuffer(ocl_system.contexts[i],
                                      CL_MEM_READ_ONLY,
                                      cb->ktree * sizeof(real),
                                      NULL,
                                      &res);

    CL_CHECK(res);

#ifndef USE_OPENMP
    feval->buf_yt[i] = clCreateBuffer(ocl_system.contexts[i],
                                      CL_MEM_READ_ONLY,
                                      rb->ktree * sizeof(real),
                                      NULL,
                                      &res);

    CL_CHECK(res);
#else
    feval->buf_yt_ff[i] = clCreateBuffer(ocl_system.contexts[i],
                                      CL_MEM_READ_ONLY,
                                      rb->ktree * sizeof(real),
                                      NULL,
                                      &res);

    CL_CHECK(res);

    feval->buf_yt_nf_common[i] = clCreateBuffer(ocl_system.contexts[i],
                                      CL_MEM_READ_ONLY,
                                      rb->ktree * sizeof(real),
                                      NULL,
                                      &res);

    CL_CHECK(res);

    feval->buf_yt_nf_min_vert[i] = clCreateBuffer(ocl_system.contexts[i],
                                                  CL_MEM_READ_ONLY,
                                                  rb->ktree * sizeof(real),
                                                  NULL,
                                                 &res);

    CL_CHECK(res);

    feval->buf_yt_nf_min_edge[i] = clCreateBuffer(ocl_system.contexts[i],
                                                  CL_MEM_READ_ONLY,
                                                  rb->ktree * sizeof(real),
                                                  NULL,
                                                  &res);

    CL_CHECK(res);
#endif
  }

  return feval;
}

void
del_fastaddevalgca(pfastaddevalgca feval)
{
  uninit_fastaddevalgca(feval);

  if(feval != NULL)
    freemem(feval);
}