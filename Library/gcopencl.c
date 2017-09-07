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

    if (gcocl->xtoffs != NULL) {
      for (uint i = 0; i < gcocl->num_row_leafs; ++i)
        if (gcocl->xtoffs[i] != NULL) {
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
new_gcopencl()
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

  gcocl->buf_xtoffs =
    (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));

  gcocl->buf_ytoffs =
    (cl_mem *) calloc(ocl_system.num_devices, sizeof(cl_mem));

  return gcocl;
}

pgcopencl
new_nearfield_gcopencl(ph2matrix H2)
{
  pgcopencl gcocl    = new_gcopencl();

  get_ocl_informations(H2, NEARFIELD, gcocl);

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

        // Update H2-matrix leaf counter for this cluster
        gcocl->num_h2_leafs_per_row[i] += 1;

        /************************* h2_leafs_per_row ***************************/

        gcocl->h2_leafs_per_row[i] =
          realloc(gcocl->h2_leafs_per_row[i],
                  gcocl->num_h2_leafs_per_row[i] * sizeof(ph2matrix*));

        gcocl->h2_leafs_per_row[i][gcocl->num_h2_leafs_per_row[i] - 1] = H2;

        /************************* col_names_per_row **************************/

        gcocl->col_names_per_row[i] =
          realloc(gcocl->col_names_per_row[i],
                  gcocl->num_h2_leafs_per_row[i] * sizeof(ph2matrix*));

        gcocl->col_names_per_row[i][gcocl->num_h2_leafs_per_row[i] - 1] = cname;

        for(uint j = 0; j < gcocl->num_row_leafs; ++j)
          for(uint k = 0; k < gcocl->num_h2_leafs_per_row[j]; ++k)
            if(i != j || ((gcocl->num_h2_leafs_per_row[i] - 1) != k))
              if(gcocl->col_names_per_row[j][k] == cname)
              {
                column_already_used = true;

                /********************** buf_cidx_off ************************/

                gcocl->cidx_off[i] =
                  realloc(gcocl->cidx_off[i],
                          gcocl->num_h2_leafs_per_row[i] * sizeof(ph2matrix*));

                gcocl->cidx_off[i][gcocl->num_h2_leafs_per_row[i] - 1] =
                  gcocl->cidx_off[j][k];

                break;
              }

        if(!column_already_used)
        {
          /************************** buf_cidx_off ****************************/

          gcocl->cidx_off[i] =
            realloc(gcocl->cidx_off[i],
                    gcocl->num_h2_leafs_per_row[i] * sizeof(ph2matrix*));

          gcocl->cidx_off[i][gcocl->num_h2_leafs_per_row[i] - 1] = gcocl->coff;

          /**************************** host_cidx *****************************/

          gcocl->host_cidx = realloc
            (gcocl->host_cidx, (gcocl->coff + H2->cb->t->size) * sizeof(uint));

          assert(gcocl->host_cidx != NULL);

          memcpy(gcocl->host_cidx + gcocl->coff,
                 H2->cb->t->idx,
                 H2->cb->t->size * sizeof(uint));

          gcocl->coff += H2->cb->t->size;
        }

        break;
      }

    if(!row_already_stored)
    {
      bool column_already_used = false;

      // Update number of leafs
      gcocl->num_row_leafs += 1;

      /**************************** names_row_leafs ***************************/

      gcocl->names_row_leafs = realloc(gcocl->names_row_leafs,
                                       gcocl->num_row_leafs * sizeof(uint));

      assert(gcocl->names_row_leafs != NULL);

      gcocl->names_row_leafs[gcocl->num_row_leafs - 1] = rname;

      /************************* num_h2_leafs_per_row *************************/

      gcocl->num_h2_leafs_per_row = realloc(gcocl->num_h2_leafs_per_row,
                                            gcocl->num_row_leafs * sizeof(uint));

      assert(gcocl->num_h2_leafs_per_row != NULL);

      gcocl->num_h2_leafs_per_row[gcocl->num_row_leafs - 1] = 1;

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

      /******************************* ridx_off *******************************/

      gcocl->ridx_off = realloc(gcocl->ridx_off,
                                gcocl->num_row_leafs * sizeof(uint));

      assert(gcocl->ridx_off != NULL);

      gcocl->ridx_off[gcocl->num_row_leafs - 1] = gcocl->roff;

      /****************************** host_ridx *******************************/

      gcocl->host_ridx = realloc
        (gcocl->host_ridx, (gcocl->roff + H2->rb->t->size) * sizeof(uint));

      assert(gcocl->host_ridx != NULL);

      memcpy(gcocl->host_ridx + gcocl->roff,
             H2->rb->t->idx,
             H2->rb->t->size * sizeof(uint));

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

        /**************************** host_cidx ****************************/

        gcocl->host_cidx = realloc
          (gcocl->host_cidx, (gcocl->coff + H2->cb->t->size) * sizeof(uint));

        assert(gcocl->host_cidx != NULL);

        memcpy(gcocl->host_cidx + gcocl->coff,
               H2->cb->t->idx,
               H2->cb->t->size * sizeof(uint));

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
iterate_recursively_h2matrix(ph2matrix H2,
                             uint      xtoff,
                             uint      ytoff,
                             pgcopencl gcocl)
{
  if((H2->u && gcocl->is_farfield) || (H2->f && !gcocl->is_farfield))
  {
    /* Write offsets if we have reached a leaf. */

    /* Nearfield matrices have an additional offset. */
    if(H2->f && !gcocl->is_farfield)
    {
      xtoff += H2->cb->k;
      ytoff += H2->rb->k;
    }

    /* Find the writing cluster of the current H^2-matrix inside the gcopencl
     * object.*/
    int  i = -1;
    uint j;

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

    for(j = 0; j < gcocl->num_h2_leafs_per_row[i]; ++j)
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
get_ocl_informations(ph2matrix H2, const information_t info, pgcopencl gcocl)
{
  gcocl->num_row_leafs            = 0;
  gcocl->max_num_h2_leafs_per_row = 0;
  gcocl->size_ridx                = 0;
  gcocl->size_cidx                = 0;
  gcocl->roff                     = 0;
  gcocl->coff                     = 0;

  gcocl->names_row_leafs          = (uint *) calloc(0, sizeof(uint));
  gcocl->num_h2_leafs_per_row     = (uint *) calloc(0, sizeof(uint));
  gcocl->h2_leafs_per_row         = (ph2matrix **) calloc(0, sizeof(ph2matrix *));
  gcocl->col_names_per_row        = (uint **) calloc(0, sizeof(uint *));
  gcocl->ridx_off                 = (uint *) calloc(0, sizeof(uint));
  gcocl->host_ridx                = (uint *) calloc(0, sizeof(uint));
  gcocl->cidx_off                 = (uint **) calloc(0, sizeof(uint *));
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

  gcocl->size_ridx  = gcocl->roff;
  gcocl->size_cidx  = gcocl->coff;
  gcocl->workload_per_row = (uint *) calloc(gcocl->num_row_leafs, sizeof(uint));
  gcocl->idx_off    = (uint *) calloc(gcocl->num_row_leafs, sizeof(uint));
  gcocl->ridx_sizes = (uint *) calloc(gcocl->num_row_leafs, sizeof(uint));

  for(uint i = 0; i < gcocl->num_row_leafs; ++i)
  {
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

  /* Get offsets of individual writing clusters relative to the whole
   * H^2-matrix. */

  gcocl->xtoffs             = (uint **) calloc(gcocl->num_row_leafs,
                                               sizeof(uint*));

  gcocl->ytoffs             = (uint *) calloc(gcocl->num_row_leafs,
                                              sizeof(uint));

  for(uint i = 0; i < gcocl->num_row_leafs; ++i)
    gcocl->xtoffs[i]       =
      (uint *) calloc(gcocl->num_h2_leafs_per_row[i], sizeof(uint));

  iterate_recursively_h2matrix(H2, 0, 0, gcocl);

  gcocl->host_xtoffs   =
    (uint *) calloc(gcocl->idx_off[gcocl->num_row_leafs - 1] +
                    gcocl->num_h2_leafs_per_row[gcocl->num_row_leafs - 1],
                    sizeof(uint));

  /* Make the offsets ready to be written to OpenCL devices. */
  for(uint i = 0; i < gcocl->num_row_leafs; ++i)
  {
    memcpy(gcocl->host_xtoffs + gcocl->idx_off[i],
           gcocl->xtoffs[i],
           gcocl->num_h2_leafs_per_row[i] * sizeof(uint));
  }

  /* Write offsets to the devices. */
  for(uint i = 0; i < ocl_system.num_devices; ++i)
  {
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
