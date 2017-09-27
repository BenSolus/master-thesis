
#include "laplace3dgca.h"

#include "bem3d.h"
#include "macrosurface3d.h"
#include "opencl.h"

#include <ctime>
#include <iostream>
#include <limits>

static const real eps = std::numeric_limits<real>::epsilon();

static real eta   = 1.0;       // Parameter for the accuracy of hierarchical
                               // clustering
static uint n     = 2048;     // Problem size
static uint m     = 8;         // Approximation order
static uint q     = 2;         // Quadratur order
static uint res   = 16;     // Cluster resolution
static real aca   = 1e-5;     // ACA resolution
static uint tests = 1;

size_t
getcouplingsize_h2matrix(pch2matrix H2)
{
  size_t sz = 0;

  if(H2->u)
    sz += getsize_amatrix(&H2->u->S);

  for(uint j = 0; j < H2->csons; ++j)
    for(uint i = 0; i < H2->rsons; ++i)
      sz += getfarsize_h2matrix(H2->son[i + j * H2->rsons]);

  return sz;
}

int
main(int argc, char *argv[])
{
  pgreencross gc;
  ph2matrix   H2;
  pavector    x;
  pavector    y_ref;
  pmacrosurface3d mg;
  psurface3d      gr;

  init_h2lib(&argc, &argv);

  std::cout << "\nMachine epsilon: " << eps << "\n";

  std::cout << "\nProblem size: " << n << "\n";

  std::cout << "\nCluster resolution: " << res << "\n";

  std::cout << "\nNumber of tests for benchmarks performed: "
            << tests
            << "\n\n";

  cl_device_id *devices;
  cl_uint      ndevices;

  get_opencl_devices(&devices, &ndevices);

  set_opencl_devices(devices, ndevices, 1);

  for(uint i = 0; i < ocl_system.num_devices; ++i)
  {
    cl_uint  num_compute_units, preferred_double_vector, preferred_float_vector;
    cl_ulong local_size;
    char     device_name[255];

    CL_CHECK(clGetDeviceInfo(ocl_system.devices[i],
                             CL_DEVICE_NAME,
                             255,
                             &device_name,
                             NULL));

    CL_CHECK(clGetDeviceInfo(ocl_system.devices[i],
                             CL_DEVICE_LOCAL_MEM_SIZE,
                             sizeof(cl_ulong),
                             &local_size,
                             NULL));

    CL_CHECK(clGetDeviceInfo(ocl_system.devices[i],
                             CL_DEVICE_MAX_COMPUTE_UNITS,
                             sizeof(cl_uint),
                             &num_compute_units,
                             NULL));

    CL_CHECK(clGetDeviceInfo(ocl_system.devices[i],
                             CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT,
                             sizeof(cl_uint),
                             &preferred_float_vector,
                             NULL));

    CL_CHECK(clGetDeviceInfo(ocl_system.devices[i],
                             CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE,
                             sizeof(cl_uint),
                             &preferred_double_vector,
                             NULL));

    printf("\n%s: CL_DEVICE_LOCAL_MEM_SIZE: %lu\n", device_name, local_size);
    printf("%s: CL_DEVICE_MAX_COMPUTE_UNITS: %u\n",
           device_name,
           num_compute_units);
    printf("%s: CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT: %u\n",
           device_name,
           preferred_float_vector);
    printf("%s: CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE: %u\n",
           device_name,
           preferred_double_vector);
  }

  std::clock_t    time(0);

  mg      = new_sphere_macrosurface3d();

  gr      = build_from_macrosurface3d_surface3d(mg, REAL_SQRT((real) n * 0.125));

  gc      = new_greencross_laplace3d(gr, res, q, m, aca);

  H2      = build_green_cross_h2matrix_greencross(gc, (void *) &eta);

//  psingquad2d sq2d = ((pbem3d) gc->bem)->sq;

//  const uint vnq = ROUNDUP(sq2d->n_id, VREAL);

//  psingquadgca sq = gc->sq_partial_min_vert;

//  for(uint i = 0; i < sq->nq; ++i)
//    printf("%3.u: xq1: %.5e, xq2: %.5e, yq1: %.5e, yq2: %.5e, wq: %.5e\n", i, (sq2d->x_id + (sq2d->n_id - sq2d->n_vert))[i], ((sq2d->x_id + (sq2d->n_id - sq2d->n_vert)) + vnq)[i], sq2d->y_id[(sq2d->n_id - sq2d->n_vert) + i], sq2d->y_id[(sq2d->n_id - sq2d->n_vert) + vnq + i], (sq2d->w_id + vnq * 9)[(sq2d->n_id - sq2d->n_vert) + i]);
//
//  std::cout << "\n";
//
//  for(uint i = 0; i < sq->nq; ++i)
//    printf("%3.u: xq1: %.5e, xq2: %.5e, yq1: %.5e, yq2: %.5e, wq: %.5e\n", i, sq->xqs[4 * sq->nq + i], sq->xqs[5 * sq->nq + i], sq->yqs[4 * sq->nq + i], sq->yqs[5 * sq->nq + i], sq->wqs[2 * sq->nq + i]);

  printf("\nGeometry:\n"
         "  %u polygons\n"
         "  %u vertices\n"
         "  %f GB\n",
         gr->triangles,
         gr->vertices,
         ((3 * gr->vertices * sizeof(real)) +
         (3 * gr->triangles * sizeof(uint)) +
         (gr->triangles * sizeof(real))) / (1024.0 * 1024.0 * 1024.0));

  printf("\nQuadrature orders:\n"
         "  distant:   %u\n"
         "  vertex:    %u\n"
         "  edge:      %u\n"
         "  identical: %u\n",
         ((pbem3d) gc->bem)->sq->n_dist,
         ((pbem3d) gc->bem)->sq->n_vert,
         ((pbem3d) gc->bem)->sq->n_edge,
         ((pbem3d) gc->bem)->sq->n_id);

  printf("\nH2-Matrix: %u blocks\n"
         "             %u row clusters\n"
         "             %u column clusters\n"
         "             size: %.5f GB\n"
         "             coupling: %.5f GB\n"
         "             neafield: %.5f GB\n",
         H2->desc,
         H2->rb->t->desc,
         H2->cb->t->desc,
         gettotalsize_h2matrix(H2) / (1024.0 * 1024.0 * 1024.0),
         getfarsize_h2matrix(H2) / (1024.0 * 1024.0 * 1024.0),
         getnearsize_h2matrix(H2) / (1024.0 * 1024.0 * 1024.0));

  for(uint i(0); i < tests; ++i)
  {
    pblock    broot = build_strict_block(gc->rc,
                                         gc->cc,
                                         &eta,
                                         admissible_max_cluster);

    ph2matrix H2_feval = build_from_block_h2matrix(broot,
                                                   gc->rb,
                                                   gc->cb);

    std::clock_t begin(std::clock());

    assemble_bem3d_h2matrix_row_clusterbasis((pbem3d) gc->bem, H2->rb);
    assemble_bem3d_h2matrix_col_clusterbasis((pbem3d) gc->bem, H2->cb);
    assemble_bem3d_h2matrix((pbem3d) gc->bem, H2);

    std::clock_t end(std::clock());

    time += end - begin;

    del_h2matrix(H2_feval);
    del_block(broot);
  }

  std::cout << "\nTime needed to construct H2-matrix:\n  "
            << 1000.0 * time / (tests * CLOCKS_PER_SEC)
            << " ms\n\n";

  time = 0;

  for(uint i(0); i < tests; ++i)
  {
    pblock    broot = build_strict_block(gc->rc,
                                         gc->cc,
                                         &eta,
                                         admissible_max_cluster);

    ph2matrix H2_feval = build_from_block_h2matrix(broot,
                                                   gc->rb,
                                                   gc->cb);

    std::clock_t begin(std::clock());

    assemble_bem3d_h2matrix((pbem3d) gc->bem, H2);

    std::clock_t end(std::clock());

    time += end - begin;

    del_h2matrix(H2_feval);
    del_block(broot);
  }

  std::cout << "\nTime needed to construct nearfield and farfield coupling "
            << "matrices:\n  "
            << 1000.0 * time / (tests * CLOCKS_PER_SEC)
            << " ms\n\n";

  /************************** Reference H2Lib H2-MVM **************************/

  time    = 0;

  x       = new_avector(n);
  random_avector(x);

  y_ref   = new_avector(n);

  for(uint i(0); i < tests; ++i)
  {
    clear_avector(y_ref);

    std::clock_t begin(std::clock());

    mvm_h2matrix_avector(1.0, false, H2, x, y_ref);

    std::clock_t end(std::clock());

    time += end - begin;
  }

  std::cout << "Reference H2-MVM:\n  "
            << 1000.0 * time / (tests * CLOCKS_PER_SEC)
            << " ms\n\n";

  /******************************** GCA H2-MVM ********************************/

  pavector y = new_avector(n);

  time = 0;

  real time_ff          = r_zero;
  real time_nf_common   = r_zero;
  real time_min_vert = r_zero;
  real time_nf_min_edge = r_zero;

  for(uint i(0); i < tests; ++i)
  {
    clear_avector(y);

    std::clock_t begin(std::clock());

    mvm_h2matrix_avector_greencross(gc, 1.0, false, H2, x, y);

    std::clock_t end(std::clock());

    time += end - begin;

    cl_ulong start = 0;
    cl_ulong fin   = 0;

    clGetEventProfilingInfo(gc->kernel_events[0],
                            CL_PROFILING_COMMAND_START,
                            sizeof(cl_ulong),
                            &start,
                            NULL);
    clGetEventProfilingInfo(gc->kernel_events[0],
                            CL_PROFILING_COMMAND_END,
                            sizeof(cl_ulong),
                            &fin,
                            NULL);

    time_ff += (real) (fin - start) * (real) 1e-06;

    start = fin = 0;

    clGetEventProfilingInfo(gc->kernel_events[1],
                            CL_PROFILING_COMMAND_START,
                            sizeof(cl_ulong),
                            &start,
                            NULL);
    clGetEventProfilingInfo(gc->kernel_events[1],
                            CL_PROFILING_COMMAND_END,
                            sizeof(cl_ulong),
                            &fin,
                            NULL);

    time_nf_common += (real) (fin - start) * (real) 1e-06;

    start = fin = 0;

    clGetEventProfilingInfo(gc->kernel_events[2],
                            CL_PROFILING_COMMAND_START,
                            sizeof(cl_ulong),
                            &start,
                            NULL);
    clGetEventProfilingInfo(gc->kernel_events[2],
                            CL_PROFILING_COMMAND_END,
                            sizeof(cl_ulong),
                            &fin,
                            NULL);

    time_min_vert += (real) (fin - start) * (real) 1e-06;

    start = fin = 0;

    clGetEventProfilingInfo(gc->kernel_events[3],
                            CL_PROFILING_COMMAND_START,
                            sizeof(cl_ulong),
                            &start,
                            NULL);
    clGetEventProfilingInfo(gc->kernel_events[3],
                            CL_PROFILING_COMMAND_END,
                            sizeof(cl_ulong),
                            &fin,
                            NULL);

    time_nf_min_edge += (real) (fin - start) * (real) 1e-06;
  }

  add_avector(r_minusone, y_ref, y);

  std::cout << "GCA H2-MVM:\n  "
            << 1000.0 * time / (tests * CLOCKS_PER_SEC)
            << " ms, rel. error: " << std::scientific
            << norm2_avector(y) / norm2_avector(y_ref) << std::fixed
            << "\n\n";

  std::cout << "GPU farfield fastaddeval:\n  "
            << time_ff / (real) tests
            << " ms\n\n";

  std::cout << "GPU partial nearfield fastaddeval of all integrals:\n  "
            << time_nf_common / (real) tests
            << " ms\n\n";

  std::cout << "GPU partial nearfield fastaddeval of integrals with at minimum of an identical vertex:\n  "
            << time_min_vert / (real) tests
            << " ms\n\n";

  std::cout << "GPU partial nearfield fastaddeval of integrals with at minimum of an identical edge:\n  "
            << time_nf_min_edge / (real) tests
            << " ms\n\n";

  time = 0;

  pavector xt = new_avector(H2->cb->ktree);
  pavector yt = new_zero_avector(H2->rb->ktree);

  forward_clusterbasis_avector(H2->cb, x, xt);
  for(uint i(0); i < tests; ++i)
  {
    clear_avector(yt);

    std::clock_t begin(std::clock());

    fastaddeval_farfield_cpu_h2matrix_avectors_gca(gc, r_one, xt, yt);

    std::clock_t end(std::clock());

    time += end - begin;
  }

  std::cout << "CPU farfield fastaddeval:\n  "
            << 1000.0 * time / (tests * CLOCKS_PER_SEC)
            << " ms\n\n";

  time = 0;

  for(uint i(0); i < tests; ++i)
  {
    clear_avector(yt);

    std::clock_t begin(std::clock());

    fastaddeval_nearfield_partial_h2matrix_avector_gca(gc, r_one, xt, yt);

    std::clock_t end(std::clock());

    time += end - begin;
  }

  std::cout << "CPU partial nearfield fastaddeval of all integrals:\n  "
            << 1000.0 * time / (tests * CLOCKS_PER_SEC)
            << " ms\n\n";

  time = 0;

  for(uint i(0); i < tests; ++i)
  {
    clear_avector(yt);

    std::clock_t begin(std::clock());

    fastaddeval_nearfield_partial_min_id_vert_gca(gc, r_one, xt, yt);

    std::clock_t end(std::clock());

    time += end - begin;
  }

  std::cout << "CPU partial nearfield fastaddeval of integrals with at minimum of an identical vertex:\n  "
            << 1000.0 * time / (tests * CLOCKS_PER_SEC)
            << " ms\n\n";

  time = 0;

  for(uint i(0); i < tests; ++i)
  {
    clear_avector(yt);

    std::clock_t begin(std::clock());

    fastaddeval_nearfield_partial_min_id_edge_gca(gc, r_one, xt, yt);

    std::clock_t end(std::clock());

    time += end - begin;
  }

  std::cout << "CPU partial nearfield fastaddeval of integrals with at minimum of an identical edge:\n  "
            << 1000.0 * time / (tests * CLOCKS_PER_SEC)
            << " ms\n\n";

  time = 0;

  for(uint i(0); i < tests; ++i)
  {
    random_avector(y);

    std::clock_t begin(std::clock());

    add_avector(f_one, y_ref, y);

    std::clock_t end(std::clock());

    time += end - begin;
  }

  add_avector(r_minusone, y_ref, y);

  std::cout << "add_avector:\n  "
            << 1000.0 * time / (tests * CLOCKS_PER_SEC)
            << " ms\n\n";

  time = 0;

  for(uint i(0); i < tests; ++i)
  {
    std::clock_t begin(std::clock());

    CL_CHECK(clEnqueueReadBuffer(ocl_system.queues[0],
                                 gc->feval->buf_yt[0],
                                 CL_TRUE,
                                 sizeof(real),
                                 y->dim * sizeof(real),
                                 y->v,
                                 0,
                                 NULL,
                                 NULL));

    std::clock_t end(std::clock());

    time += end - begin;
  }

  add_avector(r_minusone, y_ref, y);

  std::cout << "clEnqueueReadBuffer:\n  "
            << 1000.0 * time / (tests * CLOCKS_PER_SEC)
            << " ms\n\n";

  del_avector(yt);
  del_avector(xt);
  del_avector(y);
  del_avector(y_ref);
  del_avector(x);
  del_h2matrix(H2);
  del_greencross(gc);
  del_macrosurface3d(mg);

  for(uint i = 0; i < ocl_system.num_devices; ++i)
    CL_CHECK(clReleaseContext(ocl_system.contexts[i]));

  printf("  %u vectors, %u matrices and %u clusterbasis\n",
         getactives_avector(),
         getactives_amatrix(),
         getactives_clusterbasis());

  uninit_h2lib();
}
