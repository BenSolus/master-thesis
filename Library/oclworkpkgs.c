/* ------------------------------------------------------------
 * This is the file "greencross.c" of this master thesis.
 * All rights reserved, Bennet Carstensen 2017
 * ------------------------------------------------------------ */

#include "oclworkpkgs.h"

#include "ocl_system.h"

void
init_oclwork(poclworkpgs oclwrk)
{
  if(oclwrk != NULL)
  {
    oclwrk->num_wrk_pkgs         = 0;

    oclwrk->wrk_per_pkg          = NULL;
    oclwrk->first_idx_of_pkgs    = NULL;
    oclwrk->last_idx_of_pkg      = NULL;
    oclwrk->num_rows_per_pkg     = NULL;
    oclwrk->rows_per_pkg         = NULL;
    oclwrk->buf_rows_this_device = NULL;
  }
}

void
uninit_oclwork(poclworkpgs oclwrk)
{
  if(oclwrk != NULL)
  {
    if(oclwrk->wrk_per_pkg != NULL)
      freemem(oclwrk->wrk_per_pkg);

    if(oclwrk->first_idx_of_pkgs != NULL)
      freemem(oclwrk->first_idx_of_pkgs);

    if(oclwrk->last_idx_of_pkg != NULL)
      freemem(oclwrk->first_idx_of_pkgs);

    if(oclwrk->num_rows_per_pkg != NULL)
      freemem(oclwrk->num_rows_per_pkg);

    if(oclwrk->rows_per_pkg != NULL)
    {
      for(uint i = 0; i < oclwrk->num_wrk_pkgs; ++i)
        if(oclwrk->rows_per_pkg[i] != NULL)
        {
          freemem(oclwrk->rows_per_pkg[i]);
          oclwrk->rows_per_pkg[i] = NULL;
        }

      freemem(oclwrk->rows_per_pkg);
    }

    if(oclwrk->buf_rows_this_device != NULL)
    {
      for(uint i = 0; i < ocl_system.num_devices; ++i)
        CL_CHECK(clReleaseMemObject(oclwrk->buf_rows_this_device[i]));

      freemem(oclwrk->buf_rows_this_device);
    }

    init_oclwork(oclwrk);
  }
}

poclworkpgs
new_equidistant_distributed_oclwork(pcgcopencl gcocl,
                                    ph2matrix  H2,
                                    const uint num_packages)
{
//  if(num_packages > ocl_system.num_devices)
//  {
//    fprintf(stderr, "error: can't create more packages than devices!\n");
//    abort();
//  }

  poclworkpgs oclwrk = (poclworkpgs) allocmem(sizeof(oclworkpgs));;

  init_oclwork(oclwrk);

  oclwrk->num_wrk_pkgs      = num_packages;
  oclwrk->wrk_per_pkg       =
    (uint *) calloc(oclwrk->num_wrk_pkgs, sizeof(uint));
  oclwrk->num_rows_per_pkg  =
    (uint *) calloc(oclwrk->num_wrk_pkgs, sizeof(uint));
  oclwrk->first_idx_of_pkgs =
    (uint *) calloc(oclwrk->num_wrk_pkgs, sizeof(uint));
  oclwrk->last_idx_of_pkg   =
    (uint *) calloc(oclwrk->num_wrk_pkgs, sizeof(uint));
  oclwrk->num_rows_per_pkg  =
    (uint *) calloc(oclwrk->num_wrk_pkgs, sizeof(uint));
  oclwrk->rows_per_pkg      =
    (uint **) calloc(oclwrk->num_wrk_pkgs, sizeof(uint*));

  for(uint i = 0; i < oclwrk->num_wrk_pkgs; ++i)
    oclwrk->rows_per_pkg[i] = (uint *) calloc(0, sizeof(uint));

  oclwrk->buf_rows_this_device =
    (cl_mem *) calloc(oclwrk->num_wrk_pkgs, sizeof(cl_mem));

  const uint ktree = H2->rb->ktree; // TODO: Needs to be changed for transposed matrices

  for(uint i = 0; i < oclwrk->num_wrk_pkgs; ++i)
  {
    oclwrk->first_idx_of_pkgs[i] = ktree * i / oclwrk->num_wrk_pkgs;
    oclwrk->last_idx_of_pkg[i]   = ktree * (i + 1) / oclwrk->num_wrk_pkgs;
  }

  for(uint i = 0; i < gcocl->num_row_leafs; ++i)
  {
    for(uint j = 0; j < oclwrk->num_wrk_pkgs; ++j)
    {
      const bool is_responsible =
        (oclwrk->first_idx_of_pkgs[j] <= gcocl->ytoffs[i]) &&
        (oclwrk->last_idx_of_pkg[j]   >
         (gcocl->ytoffs[i] + gcocl->ridx_sizes[i] - 1));

      if(is_responsible)
      {
        oclwrk->wrk_per_pkg[j]      += gcocl->workload_per_row[i];

        oclwrk->num_rows_per_pkg[j] += 1;

        oclwrk->rows_per_pkg[j] =
          realloc(oclwrk->rows_per_pkg[j],
                  oclwrk->num_rows_per_pkg[j] * sizeof(uint));

        if(oclwrk->rows_per_pkg[j] == NULL)
        {
          fprintf(stderr, "error: failed updating work package!");
          abort();
        }

        oclwrk->rows_per_pkg[j][oclwrk->num_rows_per_pkg[j] - 1] = i;
      }
      else
      {
        const bool crosses_left_border =
          (gcocl->ytoffs[i] < oclwrk->first_idx_of_pkgs[j]) &&
           ((gcocl->ytoffs[i] + gcocl->ridx_sizes[i] - 1) >=
            oclwrk->first_idx_of_pkgs[j]);

        const bool crosses_right_bordrer =
          (gcocl->ytoffs[i] < oclwrk->last_idx_of_pkg[j]) &&
           ((gcocl->ytoffs[i] + gcocl->ridx_sizes[i] - 1) >=
            oclwrk->last_idx_of_pkg[j]);

        if(crosses_left_border || crosses_right_bordrer)
        {
          uint num_including_idxs = gcocl->ytoffs[i] + gcocl->ridx_sizes[i];

          if(crosses_left_border)
            num_including_idxs -= oclwrk->first_idx_of_pkgs[j];
          else
            num_including_idxs -= oclwrk->last_idx_of_pkg[j];

          if(((int) gcocl->ridx_sizes[i] - (int) num_including_idxs) > 0)
          {
            int k = -1;

            for(uint l = 0; l < oclwrk->num_wrk_pkgs; ++l)
            {
              if(l != j)
              {
                if(crosses_left_border)
                {
                  const bool same_border = oclwrk->last_idx_of_pkg[l] ==
                                           oclwrk->first_idx_of_pkgs[j];

                  if(same_border)
                  {
                    k = l;

                    oclwrk->first_idx_of_pkgs[j] =
                    oclwrk->last_idx_of_pkg[l]   = gcocl->ytoffs[i];
                  }
                }
                else
                {
                  const bool same_border = oclwrk->last_idx_of_pkg[j] ==
                                           oclwrk->first_idx_of_pkgs[l];

                  if(same_border)
                  {
                    k = l;

                    oclwrk->first_idx_of_pkgs[l] =
                    oclwrk->last_idx_of_pkg[j]   = gcocl->ytoffs[i] +
                                                   gcocl->ridx_sizes[i] + 1;
                  }
                }
              }
            }

            if(k < 0)
            {
              fprintf(stderr, "error: Haven't found an adjent work package!");
              abort();
            }

            oclwrk->wrk_per_pkg[j]      += gcocl->workload_per_row[i];

            oclwrk->num_rows_per_pkg[j] += 1;

            oclwrk->rows_per_pkg[j] =
              realloc(oclwrk->rows_per_pkg[j],
                      oclwrk->num_rows_per_pkg[j] * sizeof(uint));

            if(oclwrk->rows_per_pkg[j] == NULL)
            {
              fprintf(stderr, "error: failed updating work package!");
              abort();
            }

            oclwrk->rows_per_pkg[j][oclwrk->num_rows_per_pkg[j] - 1] = i;
          }
          else
          {
            int k = -1;

            for(uint l = 0; l < oclwrk->num_wrk_pkgs; ++l)
            {
              if(l != j)
              {
                if(crosses_left_border)
                {
                  const bool same_border = oclwrk->last_idx_of_pkg[l] ==
                                           oclwrk->first_idx_of_pkgs[j];

                  if(same_border)
                  {
                    k = l;

                    oclwrk->first_idx_of_pkgs[j] =
                    oclwrk->last_idx_of_pkg[l]   = gcocl->ytoffs[i] +
                                                   gcocl->ridx_sizes[i] + 1;
                  }
                }
                else
                {
                  const bool same_border = oclwrk->last_idx_of_pkg[j] ==
                                           oclwrk->first_idx_of_pkgs[l];

                  if(same_border)
                  {
                    k = l;

                    oclwrk->first_idx_of_pkgs[l] =
                    oclwrk->last_idx_of_pkg[j]   = gcocl->ytoffs[i];
                  }
                }
              }
            }

            if(k < 0)
            {
              fprintf(stderr, "error: Haven't found an adjent work package!");
              abort();
            }

            oclwrk->wrk_per_pkg[k]      += gcocl->workload_per_row[i];

            oclwrk->num_rows_per_pkg[k] += 1;

            oclwrk->rows_per_pkg[k] =
              realloc(oclwrk->rows_per_pkg[k],
                      oclwrk->num_rows_per_pkg[k] * sizeof(uint));

            if(oclwrk->rows_per_pkg[k] == NULL)
            {
              fprintf(stderr, "error: failed updating work package!");
              abort();
            }

            oclwrk->rows_per_pkg[k][oclwrk->num_rows_per_pkg[k] - 1] = i;
          }
        }
      }
    }
  }

//  for(uint i = 0; i < oclwrk->num_wrk_pkgs; ++i)
//  {
//    printf(
//      "Package %u with %u writing clusters: First index: %u, Last index: %u\n",
//      i, oclwrk->num_rows_per_pkg[i], oclwrk->first_idx_of_pkgs[i],
//      oclwrk->last_idx_of_pkg[i]);
//
//    for(uint j = 0; j < oclwrk->num_rows_per_pkg[i]; ++j)
//      printf("  Writing cluster %u: First index: %u, Last index: %u\n", oclwrk->rows_per_pkg[i][j], gcocl->ytoffs[oclwrk->rows_per_pkg[i][j]], gcocl->ytoffs[oclwrk->rows_per_pkg[i][j]] + gcocl->ridx_sizes[oclwrk->rows_per_pkg[i][j]]);
//  }
  /* Write the indices of the clusters lists a device is responsible for to the
   * corresponding device. */
  for(uint i = 0; i < oclwrk->num_wrk_pkgs; ++i)
  {
    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           oclwrk->num_rows_per_pkg[i],
                           sizeof(uint),
                           oclwrk->rows_per_pkg[i],
                           NULL,
                           &oclwrk->buf_rows_this_device[i]);
  }

  return oclwrk;
}

void
del_oclwrk(poclworkpgs oclwrk)
{
  uninit_oclwork(oclwrk);

  if(oclwrk != NULL)
    freemem(oclwrk);
}