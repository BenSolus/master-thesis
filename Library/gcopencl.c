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
#include "singquad2d.h"

#include <string.h>

void
init_gcopencl(pgcopencl gcocl)
{
  if(gcocl != NULL)
  {
    gcocl->is_farfield              = false;

    gcocl->num_row_leafs            = 0;
    gcocl->max_num_h2_leafs_per_row = 0;
    gcocl->size_ridx                = 0;
    gcocl->size_cidx                = 0;
    gcocl->roff                     = 0;
    gcocl->coff                     = 0;

    gcocl->names_row_leafs          = NULL;
    gcocl->num_h2_leafs_per_row     = NULL;
    gcocl->buf_num_h2_leafs_per_row = NULL;
    gcocl->workload_per_row         = NULL;
    gcocl->idx_off                  = NULL;
    gcocl->buf_idx_off              = NULL;
    gcocl->col_names_per_row        = NULL;
    gcocl->h2_leafs_per_row         = NULL;
    gcocl->ridx_sizes               = NULL;
    gcocl->buf_ridx_sizes           = NULL;
    gcocl->ridx_off                 = NULL;
    gcocl->buf_ridx_off             = NULL;
    gcocl->host_ridx                = NULL;
    gcocl->buf_ridx                 = NULL;
    gcocl->cidx_sizes               = NULL;
    gcocl->buf_cidx_sizes           = NULL;
    gcocl->cidx_off                 = NULL;
    gcocl->host_cidx_off            = NULL;
    gcocl->buf_cidx_off             = NULL;
    gcocl->host_cidx                = NULL;
    gcocl->buf_cidx                 = NULL;
    gcocl->num_any_id               = NULL;
    gcocl->buf_num_any_id           = NULL;
    gcocl->idx_off_any_id           = NULL;
    gcocl->buf_idx_off_any_id       = NULL;
    gcocl->rows_any_id              = NULL;
    gcocl->host_rows_any_id         = NULL;
    gcocl->buf_rows_any_id          = NULL;
    gcocl->cols_any_id              = NULL;
    gcocl->host_cols_any_id         = NULL;
    gcocl->buf_cols_any_id          = NULL;
    gcocl->cidx_any_id              = NULL;
    gcocl->host_cidx_any_id         = NULL;
    gcocl->buf_cidx_any_id          = NULL;
    gcocl->xtoffs                   = NULL;
    gcocl->host_xtoffs              = NULL;
    gcocl->buf_xtoffs               = NULL;
    gcocl->ytoffs                   = NULL;
    gcocl->buf_ytoffs               = NULL;
  }
}

void
uninit_gcopencl(pgcopencl gcocl)
{
  if(gcocl != NULL)
  {
    /* Free all memory/release all objects. */
    if (gcocl->names_row_leafs != NULL)
      freemem(gcocl->names_row_leafs);

    if (gcocl->num_h2_leafs_per_row != NULL)
      freemem(gcocl->num_h2_leafs_per_row);

    if (gcocl->buf_num_h2_leafs_per_row != NULL) {
      for (uint i = 0; i < ocl_system.num_devices; ++i) CL_CHECK(
        clReleaseMemObject(gcocl->buf_num_h2_leafs_per_row[i]));

      freemem(gcocl->buf_num_h2_leafs_per_row);
    }

    if (gcocl->workload_per_row != NULL)
      freemem(gcocl->workload_per_row);

    if (gcocl->idx_off != NULL)
      freemem(gcocl->idx_off);

    if (gcocl->buf_idx_off != NULL) {
      for (uint i = 0; i < ocl_system.num_devices; ++i) CL_CHECK(
        clReleaseMemObject(gcocl->buf_idx_off[i]));

      freemem(gcocl->buf_idx_off);
    }

    if (gcocl->h2_leafs_per_row != NULL) {
      for (uint i = 0; i < gcocl->num_row_leafs; ++i)
        if (gcocl->h2_leafs_per_row[i] != NULL) {
          freemem(gcocl->h2_leafs_per_row[i]);
          gcocl->h2_leafs_per_row[i] = NULL;
        }

      freemem(gcocl->h2_leafs_per_row);
    }

    if (gcocl->col_names_per_row != NULL) {
      for (uint i = 0; i < gcocl->num_row_leafs; ++i)
        if (gcocl->col_names_per_row[i] != NULL) {
          freemem(gcocl->col_names_per_row[i]);
          gcocl->col_names_per_row[i] = NULL;
        }

      freemem(gcocl->col_names_per_row);
    }

    if (gcocl->buf_ridx_sizes != NULL)
      freemem(gcocl->ridx_sizes);

    if (gcocl->buf_ridx_sizes != NULL) {
      for (uint i = 0; i < ocl_system.num_devices; ++i)
        if (gcocl->buf_ridx_sizes[i] != NULL) {
          CL_CHECK(clReleaseMemObject(gcocl->buf_ridx_sizes[i]));
          gcocl->buf_ridx_sizes[i] = NULL;
        }

      freemem(gcocl->buf_ridx_sizes);
    }

    if (gcocl->ridx_off != NULL)
      freemem(gcocl->ridx_off);

    if (gcocl->buf_ridx_off != NULL) {
      for (uint i = 0; i < ocl_system.num_devices; ++i)
        if (gcocl->buf_ridx_off[i] != NULL) {
          CL_CHECK(clReleaseMemObject(gcocl->buf_ridx_off[i]));
          gcocl->buf_ridx_off[i] = NULL;
        }

      freemem(gcocl->buf_ridx_off);
    }

    if (gcocl->host_ridx != NULL)
      freemem(gcocl->host_ridx);

    if (gcocl->buf_ridx != NULL) {
      for (uint i = 0; i < ocl_system.num_devices; ++i) CL_CHECK(
        clReleaseMemObject(gcocl->buf_ridx[i]));

      freemem(gcocl->buf_ridx);
    }

    if (gcocl->cidx_sizes != NULL)
      freemem(gcocl->cidx_sizes);

    if (gcocl->buf_cidx_sizes != NULL) {
      for (uint i = 0; i < ocl_system.num_devices; ++i) CL_CHECK(
        clReleaseMemObject(gcocl->buf_cidx_sizes[i]));

      freemem(gcocl->buf_cidx_sizes);
    }

    if (gcocl->cidx_off != NULL) {
      for (uint i = 0; i < gcocl->num_row_leafs; ++i)
        if (gcocl->cidx_off[i] != NULL) {
          freemem(gcocl->cidx_off[i]);
          gcocl->cidx_off[i] = NULL;
        }

      freemem(gcocl->cidx_off);
    }

    if (gcocl->host_cidx_off != NULL)
      freemem(gcocl->host_cidx_off);

    if (gcocl->buf_cidx_off != NULL) {
      for (uint i = 0; i < ocl_system.num_devices; ++i) CL_CHECK(
        clReleaseMemObject(gcocl->buf_cidx_off[i]));

      freemem(gcocl->buf_cidx_off);
    }

    if (gcocl->host_cidx != NULL)
      freemem(gcocl->host_cidx);

    if(gcocl->buf_cidx != NULL)
    {
      for(uint i = 0; i < ocl_system.num_devices; ++i) CL_CHECK(
        clReleaseMemObject(gcocl->buf_cidx[i]));

      freemem(gcocl->buf_cidx);
    }

    if(gcocl->num_any_id != NULL)
      freemem(gcocl->num_any_id);

    if(gcocl->buf_num_any_id != NULL)
    {
      for(uint i = 0; i < ocl_system.num_devices; ++i) CL_CHECK(
        clReleaseMemObject(gcocl->buf_num_any_id[i]));

      freemem(gcocl->buf_num_any_id);
    }

    if(gcocl->idx_off_any_id != NULL)
      freemem(gcocl->idx_off_any_id);

    if(gcocl->buf_idx_off_any_id != NULL)
    {
      for(uint i = 0; i < ocl_system.num_devices; ++i) CL_CHECK(
        clReleaseMemObject(gcocl->buf_idx_off_any_id[i]));

      freemem(gcocl->buf_idx_off_any_id);
    }

    if(gcocl->rows_any_id != NULL)
    {
      for(uint i = 0; i < gcocl->num_row_leafs; ++i)
        if(gcocl->rows_any_id[i] != NULL)
        {
          freemem(gcocl->rows_any_id[i]);
          gcocl->rows_any_id[i] = NULL;
        }

      freemem(gcocl->rows_any_id);
    }

    if(gcocl->host_rows_any_id != NULL)
      freemem(gcocl->host_rows_any_id);

    if(gcocl->buf_rows_any_id != NULL)
    {
      for(uint i = 0; i < ocl_system.num_devices; ++i) CL_CHECK(
        clReleaseMemObject(gcocl->buf_rows_any_id[i]));

      freemem(gcocl->buf_rows_any_id);
    }

    if(gcocl->cols_any_id != NULL)
    {
      for(uint i = 0; i < gcocl->num_row_leafs; ++i)
        if(gcocl->cols_any_id[i] != NULL)
        {
          freemem(gcocl->cols_any_id[i]);
          gcocl->cols_any_id[i] = NULL;
        }

      freemem(gcocl->cols_any_id);
    }

    if(gcocl->host_cols_any_id != NULL)
      freemem(gcocl->host_cols_any_id);

    if(gcocl->buf_cols_any_id != NULL)
    {
      for(uint i = 0; i < ocl_system.num_devices; ++i) CL_CHECK(
        clReleaseMemObject(gcocl->buf_cols_any_id[i]));

      freemem(gcocl->buf_cols_any_id);
    }

    if(gcocl->cidx_any_id != NULL)
    {
      for(uint i = 0; i < gcocl->num_row_leafs; ++i)
        if(gcocl->cidx_any_id[i] != NULL)
        {
          freemem(gcocl->cidx_any_id[i]);
          gcocl->cidx_any_id[i] = NULL;
        }

      freemem(gcocl->cidx_any_id);
    }

    if(gcocl->host_cidx_any_id != NULL)
      freemem(gcocl->host_cidx_any_id);

    if(gcocl->buf_cidx_any_id != NULL)
    {
      for(uint i = 0; i < ocl_system.num_devices; ++i) CL_CHECK(
        clReleaseMemObject(gcocl->buf_cidx_any_id[i]));

      freemem(gcocl->buf_cidx_any_id);
    }

    if(gcocl->xtoffs != NULL)
    {
      for(uint i = 0; i < gcocl->num_row_leafs; ++i)
        if (gcocl->xtoffs[i] != NULL)
        {
          freemem(gcocl->xtoffs[i]);
          gcocl->xtoffs[i] = NULL;
        }

      freemem(gcocl->xtoffs);
    }

    if (gcocl->host_xtoffs != NULL)
      freemem(gcocl->host_xtoffs);

    if (gcocl->buf_xtoffs != NULL) {
      for (uint i = 0; i < ocl_system.num_devices; ++i) CL_CHECK(
        clReleaseMemObject(gcocl->buf_xtoffs[i]));

      freemem(gcocl->buf_xtoffs);
    }

    if (gcocl->ytoffs != NULL)
      freemem(gcocl->ytoffs);

    if (gcocl->buf_ytoffs != NULL) {
      for (uint i = 0; i < ocl_system.num_devices; ++i) CL_CHECK(
        clReleaseMemObject(gcocl->buf_ytoffs[i]));

      freemem(gcocl->buf_ytoffs);
    }

    /* Set all variables to 0/NULL. */
    init_gcopencl(gcocl);
  }
}

pgcopencl
new_gcopencl(const information_t info)
{
  pgcopencl gcocl = (pgcopencl) allocmem(sizeof(gcopencl));

  init_gcopencl(gcocl);

  gcocl->buf_num_h2_leafs_per_row =
    (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));

  gcocl->buf_idx_off =
    (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));

  gcocl->buf_ridx_sizes =
    (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));

  gcocl->buf_ridx_off =
    (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));

  gcocl->buf_ridx =
    (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));

  gcocl->buf_cidx_sizes =
    (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));

  gcocl->buf_cidx_off =
    (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));

  gcocl->buf_cidx =
    (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));

  if(info == NEARFIELD)
  {
    gcocl->buf_num_any_id =
      (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));

    gcocl->buf_idx_off_any_id =
      (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));

    gcocl->buf_rows_any_id =
      (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));

    gcocl->buf_cols_any_id =
      (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));

    gcocl->buf_cidx_any_id =
      (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));
  }

  gcocl->buf_xtoffs =
    (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));

  gcocl->buf_ytoffs =
    (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));

  return gcocl;
}

pgcopencl
new_nearfield_gcopencl(ph2matrix H2, const uint dim, void *poly_idx)
{
  pgcopencl gcocl    = new_gcopencl(NEARFIELD);

  get_ocl_informations(H2, NEARFIELD, dim, poly_idx, gcocl);

  return gcocl;
}

void
del_gcopencl(pgcopencl gcocl)
{
  uninit_gcopencl(gcocl);

  if(gcocl != NULL)
    freemem(gcocl);
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

static void
push_back_h2matrix(ph2matrix **dest, const uint size, ph2matrix H2)
{
  *dest = realloc(*dest, (size + 1) * sizeof(ph2matrix*));

  if(*dest == NULL)
  {
    fprintf(stderr, "error: failed to push back a H^2-matrix!\n");
    abort();
  }

  (*dest)[size] = H2;
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
push_back_memcpy_uint(uint       **dest,
                      const uint size_dest,
                      uint       *src,
                      const uint size_src)
{
  *dest = realloc(*dest, (size_dest + size_src) * sizeof(uint));

  if(*dest == NULL)
  {
    fprintf(stderr, "error: failed to push back a uint array!\n");
    abort();
  }

  memcpy(*dest + size_dest, src, size_src * sizeof(uint));
}

void
get_ocl_nearfild_informations(ph2matrix H2,
                              uint      mname,
                              uint      rname,
                              uint      cname,
                              uint      pardepth,
                              void      *data)
{
  (void) mname;
  (void) pardepth;

  pgcopencl gcocl = (pgcopencl) data;

  if(H2->f)
  {
    bool row_already_stored = false;

    for(uint i = 0; i < gcocl->num_row_leafs; ++i)
      if(gcocl->names_row_leafs[i] == rname)
      {
        row_already_stored = true;

        bool column_already_used = false;

        const uint old_size = gcocl->num_h2_leafs_per_row[i];

        gcocl->num_h2_leafs_per_row[i] += 1;

        push_back_h2matrix(&gcocl->h2_leafs_per_row[i], old_size, H2);

        push_back_uint(&gcocl->col_names_per_row[i], old_size, cname);

        for(uint j = 0; j < gcocl->num_row_leafs; ++j)
          for(uint k = 0; k < gcocl->num_h2_leafs_per_row[j]; ++k)
            if(i != j || ((gcocl->num_h2_leafs_per_row[i] - 1) != k))
              if(gcocl->col_names_per_row[j][k] == cname)
              {
                column_already_used = true;

                push_back_uint(&gcocl->cidx_off[i], old_size, gcocl->cidx_off[j][k]);

                break;
              }

        if(!column_already_used)
        {
          push_back_uint(&gcocl->cidx_off[i], old_size, gcocl->coff);

          push_back_memcpy_uint(&gcocl->host_cidx,
                                gcocl->coff,
                                H2->cb->t->idx,
                                H2->cb->t->size);

//          gcocl->host_cidx = realloc
//            (gcocl->host_cidx, (gcocl->coff + H2->cb->t->size) * sizeof(uint));
//
//          assert(gcocl->host_cidx != NULL);
//
//          memcpy(gcocl->host_cidx + gcocl->coff,
//                 H2->cb->t->idx,
//                 H2->cb->t->size * sizeof(uint));

          gcocl->coff += H2->cb->t->size;
        }

        break;
      }

    if(!row_already_stored)
    {
      bool column_already_used = false;

      const uint old_size = gcocl->num_row_leafs;

      gcocl->num_row_leafs += 1;

      push_back_uint(&gcocl->names_row_leafs, old_size, rname);

      push_back_uint(&gcocl->num_h2_leafs_per_row, old_size, 1);

      /*************************** h2_leafs_per_row ***************************/

      gcocl->h2_leafs_per_row = realloc(gcocl->h2_leafs_per_row,
                                        gcocl->num_row_leafs * sizeof(ph2matrix*));

      assert(gcocl->h2_leafs_per_row != NULL);

      gcocl->h2_leafs_per_row[gcocl->num_row_leafs - 1] =
        calloc(1, sizeof(ph2matrix*));

      gcocl->h2_leafs_per_row[gcocl->num_row_leafs - 1][0] = H2;


      /*************************** col_names_per_row **************************/

      gcocl->col_names_per_row = realloc(gcocl->col_names_per_row,
                                         gcocl->num_row_leafs * sizeof(uint*));

      assert(gcocl->col_names_per_row != NULL);

      gcocl->col_names_per_row[gcocl->num_row_leafs - 1] =
        calloc(1, sizeof(uint*));

      gcocl->col_names_per_row[gcocl->num_row_leafs - 1][0] = cname;

      push_back_uint(&gcocl->ridx_off, old_size, gcocl->roff);

      push_back_memcpy_uint(&gcocl->host_ridx,
                            gcocl->roff,
                            H2->rb->t->idx,
                            H2->rb->t->size);

//      gcocl->host_ridx = realloc
//        (gcocl->host_ridx, (gcocl->roff + H2->rb->t->size) * sizeof(uint));
//
//      assert(gcocl->host_ridx != NULL);
//
//      memcpy(gcocl->host_ridx + gcocl->roff,
//             H2->rb->t->idx,
//             H2->rb->t->size * sizeof(uint));

      gcocl->roff += H2->rb->t->size;

      for(int i = 0; i < (int) gcocl->num_row_leafs - 2; ++i)
        for(uint j = 0; j < gcocl->num_h2_leafs_per_row[i]; ++j)
          if(gcocl->col_names_per_row[i][j] == cname)
          {
            column_already_used = true;

            /**************************** cidx_off ****************************/

            gcocl->cidx_off = realloc(gcocl->cidx_off,
                                      gcocl->num_row_leafs * sizeof(uint*));

            assert(gcocl->cidx_off != NULL);

            gcocl->cidx_off[gcocl->num_row_leafs - 1] =
              calloc(1, sizeof(uint*));

            gcocl->cidx_off[gcocl->num_row_leafs - 1][0] =
              gcocl->cidx_off[i][j];

            break;
          }

      if(!column_already_used)
      {
        /**************************** cidx_off ****************************/

        gcocl->cidx_off = realloc(gcocl->cidx_off,
                                  gcocl->num_row_leafs * sizeof(uint*));

        assert(gcocl->cidx_off != NULL);

        gcocl->cidx_off[gcocl->num_row_leafs - 1] =
          calloc(1, sizeof(uint*));

        gcocl->cidx_off[gcocl->num_row_leafs - 1][0] = gcocl->coff;

        push_back_memcpy_uint(&gcocl->host_cidx,
                              gcocl->coff,
                              H2->cb->t->idx,
                              H2->cb->t->size);

        gcocl->coff += H2->cb->t->size;
      }
    }
  }
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
                             uint       xtoff,
                             uint       ytoff,
                             const uint dim,
                             void       *poly_idx,
                             pgcopencl  gcocl)
{
  if((H2->u && gcocl->is_farfield) || (H2->f && !gcocl->is_farfield))
  {
    /* Store offsets if we have reached a leaf. */

    /* Nearfield matrices have an additional offset. */
    if(H2->f && !gcocl->is_farfield)
    {
      xtoff += H2->cb->k;
      ytoff += H2->rb->k;
    }

    /* Find the writing cluster of the current H^2-matrix inside the gcopencl
     * object.*/
    int  i = -1;

    for(uint k = 0; k < gcocl->num_row_leafs; ++k)
      if(are_equal_clusters(H2->rb->t, gcocl->h2_leafs_per_row[k][0]->rb->t))
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

    /* Write the offsets in the corresponding entries. */
    gcocl->ytoffs[i] = ytoff;

    uint j = 0;

    for(; j < gcocl->num_h2_leafs_per_row[i]; ++j)
      if(are_equal_clusters(H2->cb->t, gcocl->h2_leafs_per_row[i][j]->cb->t))
      {
        gcocl->xtoffs[i][j] = xtoff;
        break;
      }

    /* We face an error if we haven't found any matching reading cluster. */
    if(j >= gcocl->num_h2_leafs_per_row[i])
    {
      fprintf(stderr, "error: Unknown error in "
        "iterate_recursively_h2matrix!\n");

      abort();
    }

    /* Get the necessary informations of all entries of this matrix where the
     * row and column polygon share at least one common vertex. */
    if(H2->f)
    {
      const uint idx_off = gcocl->idx_off[i];
      const uint *ridx   = gcocl->host_ridx + gcocl->ridx_off[i];
      const uint *cidx   = gcocl->host_cidx + gcocl->cidx_off[i][j];

      uint (*polys)[dim] = (uint(*)[dim]) poly_idx;

      for(uint k = 0; k < gcocl->cidx_sizes[idx_off + j]; ++k)
      {
        const uint kk = cidx[k];

        for(uint l = 0; l < gcocl->ridx_sizes[i]; ++l)
        {
          if(dim == 3)
          {
            if(fast_select_quadrature((uint(*)[3]) polys, kk, ridx[l]) > 0)
            {
              const uint old_size = gcocl->num_any_id[i];

              gcocl->num_any_id[i] += 1;

              push_back_uint(&gcocl->rows_any_id[i], old_size, l);

              push_back_uint(&gcocl->cols_any_id[i], old_size, xtoff + k);

              push_back_uint(&gcocl->cidx_any_id[i], old_size, kk);
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

        iterate_recursively_h2matrix(H2->son[i + j * rsons],
                                     xtoff + xtoff1 - (cb->sons > 0 ? 0 : cb->k),
                                     ytoff + ytoff1 - (rb->sons > 0 ? 0 : rb->k),
                                     dim,
                                     poly_idx,
                                     gcocl);

        ytoff1 += (rb->sons > 0 ? rb->son[i]->ktree : rb->t->size);
      }

      assert(ytoff1 == rb->ktree);

      xtoff1 += (cb->sons > 0 ? cb->son[j]->ktree : cb->t->size);
    }

    assert(xtoff1 == cb->ktree);
  }
}

void
get_ocl_informations(ph2matrix           H2,
                     const information_t info,
                     const uint          dim,
                     void                *poly_idx,
                     pgcopencl           gcocl)
{
  gcocl->num_row_leafs            = 0;
  gcocl->max_num_h2_leafs_per_row = 0;
  gcocl->size_ridx                = 0;
  gcocl->size_cidx                = 0;
  gcocl->roff                     = 0;
  gcocl->coff                     = 0;

  gcocl->names_row_leafs          = (uint *) calloc(0, sizeof(uint));
  gcocl->num_h2_leafs_per_row     = (uint *) calloc(0, sizeof(uint));
  gcocl->h2_leafs_per_row         = (ph2matrix **) calloc(0, sizeof(ph2matrix*));
  gcocl->col_names_per_row        = (uint **) calloc(0, sizeof(uint*));
  gcocl->ridx_off                 = (uint *) calloc(0, sizeof(uint));
  gcocl->host_ridx                = (uint *) calloc(0, sizeof(uint));
  gcocl->cidx_off                 = (uint **) calloc(0, sizeof(uint*));
  gcocl->host_cidx                = (uint *) calloc(0, sizeof(uint));

  switch(info)
  {
    case NEARFIELD: gcocl->is_farfield = false;
                    iterate_h2matrix(H2,
                                     0,
                                     0,
                                     0,
                                     0,
                                     get_ocl_nearfild_informations,
                                     NULL,
                                     (void *) gcocl);
                    break;
    default:        fprintf(stderr, "error: given case for "
                                    "get_ocl_informations not implemented yet");
                    abort();
  }

  gcocl->size_ridx   = gcocl->roff;
  gcocl->size_cidx   = gcocl->coff;
  gcocl->workload_per_row = (uint *) calloc(gcocl->num_row_leafs, sizeof(uint));
  gcocl->idx_off     = (uint *) calloc(gcocl->num_row_leafs, sizeof(uint));
  gcocl->ridx_sizes  = (uint *) calloc(gcocl->num_row_leafs, sizeof(uint));
  gcocl->num_any_id  = (uint *) calloc(gcocl->num_row_leafs, sizeof(uint));
  gcocl->rows_any_id = (uint **) calloc(gcocl->num_row_leafs, sizeof(uint*));
  gcocl->cols_any_id = (uint **) calloc(gcocl->num_row_leafs, sizeof(uint*));
  gcocl->cidx_any_id = (uint **) calloc(gcocl->num_row_leafs, sizeof(uint*));

  for(uint i = 0; i < gcocl->num_row_leafs; ++i)
  {
    gcocl->rows_any_id[i] = (uint*) calloc(0, sizeof(uint));
    gcocl->cols_any_id[i] = (uint*) calloc(0, sizeof(uint));
    gcocl->cidx_any_id[i] = (uint*) calloc(0, sizeof(uint));

    /* Get maximum number of H2-matrix leafs which share the same writing
     * cluster. */
    if(gcocl->num_h2_leafs_per_row[i] > gcocl->max_num_h2_leafs_per_row)
      gcocl->max_num_h2_leafs_per_row = gcocl->num_h2_leafs_per_row[i];

    /* Get workload of the current writing cluster. */
    for(uint j = 0; j < gcocl->num_h2_leafs_per_row[i]; ++j)
      gcocl->workload_per_row[i] +=
        gcocl->is_farfield
        ? gcocl->h2_leafs_per_row[i][j]->u->S.rows *
          gcocl->h2_leafs_per_row[i][j]->u->S.cols
        : gcocl->h2_leafs_per_row[i][j]->f->rows *
          gcocl->h2_leafs_per_row[i][j]->f->cols;

    /* Get the index offset for the current writing cluster needed by the
     * devices. */
    gcocl->idx_off[i] =
      (i == 0 ? 0 : gcocl->idx_off[i - 1] + gcocl->num_h2_leafs_per_row[i - 1]);

    /* Sum up the number of indices needed by the devices. */
    gcocl->ridx_sizes[i] = gcocl->h2_leafs_per_row[i][0]->rb->t->size;
    // TODO: ^Differs between farfield and nearfield
  }

  gcocl->cidx_sizes =
    (uint *) calloc(gcocl->idx_off[gcocl->num_row_leafs - 1] +
                    gcocl->num_h2_leafs_per_row[gcocl->num_row_leafs - 1],
                    sizeof(uint));

  gcocl->host_cidx_off =
    (uint *) calloc(gcocl->idx_off[gcocl->num_row_leafs - 1] +
                    gcocl->num_h2_leafs_per_row[gcocl->num_row_leafs - 1],
                    sizeof(uint));

  /* Get number of indices for all H2-matrix leafs, ready to be written to
   * OpenCL devices. */
  for (uint i = 0; i < gcocl->num_row_leafs; ++i)
  {
    for(uint j = 0; j < gcocl->num_h2_leafs_per_row[i]; ++j)
      gcocl->cidx_sizes[gcocl->idx_off[i] + j] =
        gcocl->h2_leafs_per_row[i][j]->cb->t->size;
        // TODO: ^Differs between farfield and nearfield

    memcpy(gcocl->host_cidx_off + gcocl->idx_off[i],
           gcocl->cidx_off[i],
           gcocl->num_h2_leafs_per_row[i] * sizeof(uint));
  }

  /* Write H2-matrix informations to the devices. */
  for(uint i = 0; i < ocl_system.num_devices; ++i)
  {
    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           gcocl->num_row_leafs,
                           sizeof(uint),
                           gcocl->num_h2_leafs_per_row,
                           NULL,
                           &gcocl->buf_num_h2_leafs_per_row[i]);

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           gcocl->num_row_leafs,
                           sizeof(uint),
                           gcocl->idx_off,
                           NULL,
                           &gcocl->buf_idx_off[i]);

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           gcocl->num_row_leafs,
                           sizeof(uint),
                           gcocl->ridx_sizes,
                           NULL,
                           &gcocl->buf_ridx_sizes[i]);

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           gcocl->num_row_leafs,
                           sizeof(uint),
                           gcocl->ridx_off,
                           NULL,
                           &gcocl->buf_ridx_off[i]);

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           gcocl->size_ridx,
                           sizeof(uint),
                           gcocl->host_ridx,
                           NULL,
                           &gcocl->buf_ridx[i]);

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           gcocl->idx_off[gcocl->num_row_leafs - 1] +
                           gcocl->num_h2_leafs_per_row[gcocl->num_row_leafs - 1],
                           sizeof(uint),
                           gcocl->cidx_sizes,
                           NULL,
                           &gcocl->buf_cidx_sizes[i]);

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           gcocl->idx_off[gcocl->num_row_leafs - 1] +
                           gcocl->num_h2_leafs_per_row[gcocl->num_row_leafs - 1],
                           sizeof(uint),
                           gcocl->host_cidx_off,
                           NULL,
                           &gcocl->buf_cidx_off[i]);

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           gcocl->size_cidx,
                           sizeof(uint),
                           gcocl->host_cidx,
                           NULL,
                           &gcocl->buf_cidx[i]);
  }

  gcocl->xtoffs             = (uint **) calloc(gcocl->num_row_leafs,
                                               sizeof(uint*));

  gcocl->ytoffs             = (uint *) calloc(gcocl->num_row_leafs,
                                              sizeof(uint));

  for(uint i = 0; i < gcocl->num_row_leafs; ++i)
    gcocl->xtoffs[i]       =
      (uint *) calloc(gcocl->num_h2_leafs_per_row[i], sizeof(uint));

  iterate_recursively_h2matrix(H2, 0, 0, dim, poly_idx, gcocl);

  gcocl->idx_off_any_id = (uint*) calloc(gcocl->num_row_leafs, sizeof(uint));

  gcocl->host_xtoffs    =
    (uint *) calloc(gcocl->idx_off[gcocl->num_row_leafs - 1] +
                    gcocl->num_h2_leafs_per_row[gcocl->num_row_leafs - 1],
                    sizeof(uint));


  for(uint i = 0; i < gcocl->num_row_leafs; ++i)
  {
    /* Get index offsets for informations about integrals where row and column
     * polygon share at least one common vertex. */
    if(i > 0)
      gcocl->idx_off_any_id[i] = gcocl->idx_off_any_id[i - 1] +
                                 gcocl->num_any_id[i - 1];

    /* Make the offsets for xt ready to be written to OpenCL devices. */
    memcpy(gcocl->host_xtoffs + gcocl->idx_off[i],
           gcocl->xtoffs[i],
           gcocl->num_h2_leafs_per_row[i] * sizeof(uint));
  }

  gcocl->host_rows_any_id =
    (uint *) calloc(gcocl->idx_off_any_id[gcocl->num_row_leafs - 1] +
                    gcocl->num_any_id[gcocl->num_row_leafs - 1],
                    sizeof(uint));

  gcocl->host_cols_any_id =
    (uint *) calloc(gcocl->idx_off_any_id[gcocl->num_row_leafs - 1] +
                    gcocl->num_any_id[gcocl->num_row_leafs - 1],
                    sizeof(uint));

  gcocl->host_cidx_any_id =
    (uint *) calloc(gcocl->idx_off_any_id[gcocl->num_row_leafs - 1] +
                    gcocl->num_any_id[gcocl->num_row_leafs - 1],
                    sizeof(uint));

  for(uint i = 0; i < gcocl->num_row_leafs; ++i)
  {
    memcpy(gcocl->host_rows_any_id + gcocl->idx_off_any_id[i],
           gcocl->rows_any_id[i],
           gcocl->num_any_id[i] * sizeof(uint));

    memcpy(gcocl->host_cols_any_id + gcocl->idx_off_any_id[i],
           gcocl->cols_any_id[i],
           gcocl->num_any_id[i] * sizeof(uint));

    memcpy(gcocl->host_cidx_any_id + gcocl->idx_off_any_id[i],
           gcocl->cidx_any_id[i],
           gcocl->num_any_id[i] * sizeof(uint));
  }

  /* Write offsets to the devices. */
  for(uint i = 0; i < ocl_system.num_devices; ++i)
  {
    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           gcocl->num_row_leafs,
                           sizeof(uint),
                           gcocl->num_any_id,
                           NULL,
                           &gcocl->buf_num_any_id[i]);

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           gcocl->num_row_leafs,
                           sizeof(uint),
                           gcocl->idx_off_any_id,
                           NULL,
                           &gcocl->buf_idx_off_any_id[i]);

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           gcocl->idx_off_any_id[gcocl->num_row_leafs - 1] +
                           gcocl->num_any_id[gcocl->num_row_leafs - 1],
                           sizeof(uint),
                           gcocl->host_rows_any_id,
                           NULL,
                           &gcocl->buf_rows_any_id[i]);

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           gcocl->idx_off_any_id[gcocl->num_row_leafs - 1] +
                           gcocl->num_any_id[gcocl->num_row_leafs - 1],
                           sizeof(uint),
                           gcocl->host_cols_any_id,
                           NULL,
                           &gcocl->buf_cols_any_id[i]);

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           gcocl->idx_off_any_id[gcocl->num_row_leafs - 1] +
                           gcocl->num_any_id[gcocl->num_row_leafs - 1],
                           sizeof(uint),
                           gcocl->host_cidx_any_id,
                           NULL,
                           &gcocl->buf_cidx_any_id[i]);

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           gcocl->idx_off[gcocl->num_row_leafs - 1] +
                           gcocl->num_h2_leafs_per_row[gcocl->num_row_leafs - 1],
                           sizeof(uint),
                           gcocl->host_xtoffs,
                           NULL,
                           &gcocl->buf_xtoffs[i]);

    create_and_fill_buffer(ocl_system.contexts[i],
                           CL_MEM_READ_ONLY,
                           ocl_system.queues[i * ocl_system.queues_per_device],
                           gcocl->num_row_leafs,
                           sizeof(uint),
                           gcocl->ytoffs,
                           NULL,
                           &gcocl->buf_ytoffs[i]);
  }
}
