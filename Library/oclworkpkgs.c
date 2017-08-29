/* ------------------------------------------------------------
 * This is the file "greencross.c" of this master thesis.
 * All rights reserved, Bennet Carstensen 2017
 * ------------------------------------------------------------ */

#include "oclworkpkgs.h"

#include "ocl_system.h"

void
init_oclwork(poclworkpgs oclwrk)
{
  oclwrk->num_wrk_pkgs     = 0;

  oclwrk->wrk_per_pkg      = NULL;

  oclwrk->num_rows_per_pkg = NULL;

  oclwrk->rows_per_pkg     = NULL;

  oclwrk->buf_rows_this_device =
    (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));
}

void
uninit_oclwork(poclworkpgs oclwrk)
{
  if(oclwrk->wrk_per_pkg != NULL)
  {
    freemem(oclwrk->wrk_per_pkg);
    oclwrk->wrk_per_pkg = NULL;
  }

  if(oclwrk->num_rows_per_pkg)
  {
    freemem(oclwrk->num_rows_per_pkg);
    oclwrk->num_rows_per_pkg = NULL;
  }

  if(oclwrk->rows_per_pkg != NULL)
  {
    for(uint i = 0; i < oclwrk->num_wrk_pkgs; ++i)
      if(oclwrk->rows_per_pkg[i] != NULL)
      {
        freemem(oclwrk->rows_per_pkg[i]);
        oclwrk->rows_per_pkg[i] = NULL;
      }

    freemem(oclwrk->rows_per_pkg);
    oclwrk->rows_per_pkg = NULL;
  }
}
