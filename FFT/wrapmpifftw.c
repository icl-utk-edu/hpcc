
#include <stdio.h>
#include <stdlib.h>

#include "hpccfft.h"
#include "wrapmpifftw.h"

hpcc_fftw_mpi_plan
hpcc_fftw_mpi_create_plan(MPI_Comm comm, int n, fftw_direction dir, int flags) {
  hpcc_fftw_mpi_plan p;
  int zero = 0;
  fftw_complex *a = NULL, *b = NULL;
  int rank, size;

  MPI_Comm_size( comm, &size );
  MPI_Comm_rank( comm, &rank );

  p = fftw_malloc( sizeof *p );

  p->n0 = 0;
  p->n = n;
  p->comm = comm;

  p->c = fftw_malloc( ((NDA3+NP)*NBLK+NP) * sizeof(fftw_complex) );
  p->d = fftw_malloc( (NDA3+NP) * sizeof(fftw_complex) );
  p->wx = fftw_malloc( (NDA3/2+NP) * sizeof(fftw_complex) );
  p->wy = fftw_malloc( (NDA3/2+NP) * sizeof(fftw_complex) );
  p->wz = fftw_malloc( (NDA3/2+NP) * sizeof(fftw_complex) );
  p->ww = fftw_malloc( (NDA3*NDA3) * sizeof(fftw_complex) );
  p->www = fftw_malloc( (n/size) * sizeof(fftw_complex) );

  p->dir = dir;
  p->flags = flags;

  pzfft1d_( a, b, p->ww, p->www, &n, &rank, &size, &zero, p, &p->n0 );

  return p;
}

void
hpcc_fftw_mpi_destroy_plan(hpcc_fftw_mpi_plan p) {
  if (!p) return;

  fftw_free( p->www );
  fftw_free( p->ww );
  fftw_free( p->wz );
  fftw_free( p->wy );
  fftw_free( p->wx );
  fftw_free( p->d );
  fftw_free( p->c );
  fftw_free( p );
}

void
hpcc_fftw_mpi(hpcc_fftw_mpi_plan p, int n_fields, fftw_complex *local_data, fftw_complex *work){
  int one = 1, two = 2;
  int rank, size;
  int n;

  MPI_Comm_size( p->comm, &size );
  MPI_Comm_rank( p->comm, &rank );

  n = p->n;

  if (FFTW_FORWARD == p->dir)
    pzfft1d_( local_data, work, p->ww, p->www, &n, &rank, &size, &one, p, &p->n0 );
  else
    pzfft1d_( local_data, work, p->ww, p->www, &n, &rank, &size, &two, p, &p->n0 );
}

void
hpcc_fftw_mpi_local_sizes(hpcc_fftw_mpi_plan p, int *local_n, int *local_start,
  int *local_n_after_transform, int *local_start_after_transform, int *total_local_size) {
  int rank, size;
  int n;
  MPI_Comm_size( p->comm, &size );
  MPI_Comm_rank( p->comm, &rank );
  n = p->n;
  if (local_n) *local_n = n / size;
  if (local_start) *local_start = n / size * rank;
  if (local_n_after_transform) *local_n_after_transform = n / size;
  if (local_start_after_transform) *local_start_after_transform = n / size * rank;
  if (total_local_size) *total_local_size = n / size;
}
