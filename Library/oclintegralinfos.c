/* ------------------------------------------------------------
 * This is the file "oclintegralinfos.c" of this master thesis.
 * All rights reserved, Bennet Carstensen 2017
 * ------------------------------------------------------------ */

/**
 * @file      oclintegralinfos.c
 * @author    Bennet Carstensen
 * @date      2017
 * @copyright All rights reserved, Bennet Carstensen 2017
 */

#include "oclintegralinfos.h"

#include "ocl_system.h"
#include "singquad2d.h"

#include <memory.h>

void
init_integralinfos(pintegralinfos iinfos)
{
  if(iinfos != NULL)
  {
    iinfos->num_integral_grps = 0;
    iinfos->num_integrals     = NULL;
    iinfos->buf_num_integrals = NULL;
    iinfos->idx_off           = NULL;
    iinfos->buf_idx_off       = NULL;
    iinfos->rows              = NULL;
    iinfos->host_rows         = NULL;
    iinfos->buf_rows          = NULL;
    iinfos->cols              = NULL;
    iinfos->host_cols         = NULL;
    iinfos->cidx              = NULL;
    iinfos->host_cidx         = NULL;
    iinfos->buf_cidx          = NULL;
  }
}

void
uninit_integralinfos(pintegralinfos iinfos)
{
  if(iinfos != NULL)
  {
    if(iinfos->num_integrals != NULL)
      freemem(iinfos->num_integrals);

    if(iinfos->buf_num_integrals != NULL)
    {
      for(uint i = 0; i < ocl_system.num_devices; ++i)
        CL_CHECK(clReleaseMemObject(iinfos->buf_num_integrals[i]));

      freemem(iinfos->buf_num_integrals);
    }

    if(iinfos->idx_off != NULL)
      freemem(iinfos->idx_off);

    if(iinfos->buf_idx_off != NULL)
    {
      for(uint i = 0; i < ocl_system.num_devices; ++i)
        CL_CHECK(clReleaseMemObject(iinfos->buf_idx_off[i]));

      freemem(iinfos->buf_idx_off);
    }

    if(iinfos->rows != NULL)
    {
      for(uint i = 0; i < iinfos->num_integral_grps; ++i)
        if(iinfos->rows[i] != NULL)
          freemem(iinfos->rows[i]);

      freemem(iinfos->rows);
    }

    if(iinfos->host_rows != NULL)
      freemem(iinfos->host_rows);

    if(iinfos->buf_rows != NULL)
    {
      for(uint i = 0; i < ocl_system.num_devices; ++i)
      CL_CHECK(clReleaseMemObject(iinfos->buf_rows[i]));

      freemem(iinfos->buf_rows);
    }

    if(iinfos->cols != NULL)
    {
      for(uint i = 0; i < iinfos->num_integral_grps; ++i)
        if(iinfos->cols[i] != NULL)
          freemem(iinfos->cols[i]);

      freemem(iinfos->cols);
    }

    if(iinfos->host_cols != NULL)
      freemem(iinfos->host_cols);

    if(iinfos->buf_cols != NULL)
    {
      for(uint i = 0; i < ocl_system.num_devices; ++i)
      CL_CHECK(clReleaseMemObject(iinfos->buf_cols[i]));

      freemem(iinfos->buf_cols);
    }

    if(iinfos->cidx != NULL)
    {
      for(uint i = 0; i < iinfos->num_integral_grps; ++i)
        if(iinfos->cidx[i] != NULL)
          freemem(iinfos->cidx[i]);

      freemem(iinfos->cidx);
    }

    if(iinfos->host_cidx != NULL)
      freemem(iinfos->host_cidx);

    if(iinfos->buf_cidx != NULL)
    {
      for(uint i = 0; i < ocl_system.num_devices; ++i)
      CL_CHECK(clReleaseMemObject(iinfos->buf_cidx[i]));

      freemem(iinfos->buf_cidx);
    }
  }
}

pintegralinfos
new_integralinfos()
{
  pintegralinfos iinfos = (pintegralinfos) allocmem(sizeof(integralinfos));

  iinfos->buf_num_integrals =
    (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));

  iinfos->buf_idx_off =
    (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));

  iinfos->buf_rows =
    (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));

  iinfos->buf_cols =
    (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));

  iinfos->buf_cidx =
    (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));

  return iinfos;
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
push_back_uint(uint **dest, const uint size, uint entry)
{
  *dest = realloc(*dest, (size + 1) * sizeof(uint*));

  if(*dest == NULL)
  {
    fprintf(stderr, "error: failed to push back a uint!\n");
    abort();
  }

  (*dest)[size] = entry;
}

static void
iterate_recursively_h2matrix(pcgcopencl     idxinfos,
                             pch2matrix     H2,
                             uint           xtoff,
                             const uint     dim,
                             void           *poly_idx,
                             const uint     min_num_id_verts,
                             pintegralinfos iinfos)
{
  if(H2->f)
  {
    xtoff += H2->cb->k;

    /* Find the writing cluster of the current H^2-matrix inside the index info
     * object.*/
    int  i = -1;

    for(uint k = 0; k < idxinfos->num_row_leafs; ++k)
      if(are_equal_clusters(H2->rb->t, idxinfos->h2_leafs_per_row[k][0]->rb->t))
      {
        i = k;
        break;
      }

    /* We face an error if we haven't found any matching writing cluster. */
    if(i < 0)
    {
      fprintf(stderr, "error: Unknown error in "
        "iterate_recursively_h2matrix!\n");

      abort();
    }

    uint j = 0;

    /* Find the reading cluster of the current H^2-matrix inside the index info
     * object.*/
    for(; j < idxinfos->num_h2_leafs_per_row[i]; ++j)
      if(are_equal_clusters(H2->cb->t, idxinfos->h2_leafs_per_row[i][j]->cb->t))
        break;

    /* We face an error if we haven't found any matching reading cluster. */
    if(j >= idxinfos->num_h2_leafs_per_row[i])
    {
      fprintf(stderr, "error: Unknown error in "
        "iterate_recursively_h2matrix!\n");

      abort();
    }

    /* Get the necessary informations of all entries of this matrix where the
     * row and column polygon share at least one common vertex. */
    if(H2->f)
    {
      const uint idx_off = idxinfos->idx_off[i];
      const uint *ridx   = idxinfos->host_ridx + idxinfos->ridx_off[i];
      const uint *cidx   = idxinfos->host_cidx + idxinfos->cidx_off[i][j];

      uint (*polys)[dim] = (uint(*)[dim]) poly_idx;

      for(uint k = 0; k < idxinfos->cidx_sizes[idx_off + j]; ++k)
      {
        const uint kk = cidx[k];

        for(uint l = 0; l < idxinfos->ridx_sizes[i]; ++l)
        {
          if(dim == 3)
          {
            if(fast_select_quadrature((uint(*)[3]) polys, kk, ridx[l]) > min_num_id_verts)
            {
              const uint old_size = iinfos->num_integrals[i];

              iinfos->num_integrals[i] += 1;

              push_back_uint(&iinfos->rows[i], old_size, l);

              push_back_uint(&iinfos->cols[i], old_size, xtoff + k);

              push_back_uint(&iinfos->cidx[i], old_size, kk);
            }
          }
          else
            fprintf(stderr, "2D not implemented yet in "
              "iterate_recursively_h2matrix");
        }
      }
    }
  }
  else if(H2->son)
  {
    /* Recursivly call the children. */

    pcclusterbasis rb = H2->rb;
    pcclusterbasis cb = H2->cb;

    const uint rsons  = H2->rsons;
    const uint csons  = H2->csons;

    uint xtoff1 = cb->k;

    for(uint j = 0; j < csons; ++j)
    {
      assert(csons == 1 || cb->sons > 0);

      uint ytoff1 = rb->k;

      for(uint i = 0; i < rsons; ++i)
      {
        assert(rsons == 1 || rb->sons > 0);

        iterate_recursively_h2matrix(idxinfos,
                                     H2->son[i + j * rsons],
                                     xtoff + xtoff1 - (cb->sons > 0 ? 0 : cb->k),
                                     dim,
                                     poly_idx,
                                     min_num_id_verts,
                                     iinfos);

        ytoff1 += (rb->sons > 0 ? rb->son[i]->ktree : rb->t->size);
      }

      assert(ytoff1 == rb->ktree);

      xtoff1 += (cb->sons > 0 ? cb->son[j]->ktree : cb->t->size);
    }

    assert(xtoff1 == cb->ktree);
  }
}

pintegralinfos
build_from_idxinfos_integralinfos(pcgcopencl idxinfos,
                                  pch2matrix H2,
                                  const uint dim,
                                  void       *poly_idxs,
                                  const uint min_num_id_verts)
{
  pintegralinfos iinfos = new_integralinfos();

  iinfos->num_integral_grps = idxinfos->num_row_leafs;

  iinfos->num_integrals = (uint*)   calloc(iinfos->num_integral_grps, sizeof(uint));
  iinfos->idx_off       = (uint*)   calloc(iinfos->num_integral_grps, sizeof(uint));
  iinfos->rows          = (uint **) calloc(iinfos->num_integral_grps, sizeof(uint*));
  iinfos->cols          = (uint **) calloc(iinfos->num_integral_grps, sizeof(uint*));
  iinfos->cidx          = (uint **) calloc(iinfos->num_integral_grps, sizeof(uint*));

  for(uint i = 0; i < iinfos->num_integral_grps; ++i)
  {
    iinfos->rows[i] = (uint *) calloc(0, sizeof(uint));
    iinfos->cols[i] = (uint *) calloc(0, sizeof(uint));
    iinfos->cidx[i] = (uint *) calloc(0, sizeof(uint));
  }

  iterate_recursively_h2matrix(idxinfos, H2, 0, dim, poly_idxs, min_num_id_verts, iinfos);

  iinfos->idx_off = (uint*) calloc(iinfos->num_integral_grps, sizeof(uint));

  for(uint i = 0; i < iinfos->num_integral_grps; ++i)
  {
    /* Get index offsets for informations about integrals where row and column
     * polygon share at least min_num_id_verts common vertices. */
    if(i > 0)
      iinfos->idx_off[i] = iinfos->idx_off[i - 1] + iinfos->num_integrals[i - 1];
  }

  /* Allocate memory for nondimensionalize versions of the arrays. */

  iinfos->host_rows =
    (uint *) calloc(iinfos->idx_off[iinfos->num_integral_grps - 1] +
                    iinfos->num_integrals[iinfos->num_integral_grps - 1],
                    sizeof(uint));

  iinfos->host_cols =
    (uint *) calloc(iinfos->idx_off[iinfos->num_integral_grps - 1] +
                    iinfos->num_integrals[iinfos->num_integral_grps - 1],
                    sizeof(uint));

  iinfos->host_cidx =
    (uint *) calloc(iinfos->idx_off[iinfos->num_integral_grps - 1] +
                    iinfos->num_integrals[iinfos->num_integral_grps - 1],
                    sizeof(uint));

  /* Write in the nondimensionalize versions of the arrays. */
  for(uint i = 0; i < iinfos->num_integral_grps; ++i)
  {
    memcpy(iinfos->host_rows + iinfos->idx_off[i],
           iinfos->rows[i],
           iinfos->num_integrals[i] * sizeof(uint));

    memcpy(iinfos->host_cols + iinfos->idx_off[i],
           iinfos->cols[i],
           iinfos->num_integrals[i] * sizeof(uint));

    memcpy(iinfos->host_cidx + iinfos->idx_off[i],
           iinfos->cidx[i],
           iinfos->num_integrals[i] * sizeof(uint));
  }

  /* Write everything to the devices. */
  for(uint i = 0; i < ocl_system.num_devices; ++i)
  {
    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           iinfos->num_integral_grps,
                           sizeof(uint),
                           iinfos->num_integrals,
                           NULL,
                           &iinfos->buf_num_integrals[i]);

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           iinfos->num_integral_grps,
                           sizeof(uint),
                           iinfos->idx_off,
                           NULL,
                           &iinfos->buf_idx_off[i]);

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           iinfos->idx_off[iinfos->num_integral_grps - 1] +
                           iinfos->num_integrals[iinfos->num_integral_grps - 1],
                           sizeof(uint),
                           iinfos->host_rows,
                           NULL,
                           &iinfos->buf_rows[i]);

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           iinfos->idx_off[iinfos->num_integral_grps - 1] +
                           iinfos->num_integrals[iinfos->num_integral_grps - 1],
                           sizeof(uint),
                           iinfos->host_cols,
                           NULL,
                           &iinfos->buf_cols[i]);

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           iinfos->idx_off[iinfos->num_integral_grps - 1] +
                           iinfos->num_integrals[iinfos->num_integral_grps - 1],
                           sizeof(uint),
                           iinfos->host_cidx,
                           NULL,
                           &iinfos->buf_cidx[i]);
  }

  return iinfos;
}
