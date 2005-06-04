/* -*- mode: C; tab-width: 2; indent-tabs-mode: nil; fill-column: 79; coding: iso-latin-1-unix -*- */
/* mpifft.c
 */

#include <hpcc.h>

#include "hpccfft.h"
#include "wrapmpifftw.h"

double *HPCC_fft_timings_forward, *HPCC_fft_timings_backward;

static void
MPIFFT0(HPCC_Params *params, int doIO, FILE *outFile, MPI_Comm comm, double *UGflops, s64Int_t *Un,
        double *UmaxErr, int *Ufailure) {
  int commRank, commSize;
  int failure;
  s64Int_t i, n;
  s64Int_t locn, loc0, alocn, aloc0, tls;
  double maxErr, tmp1, tmp2, tmp3, t0, t1, t2, t3, Gflops;
  double deps;
  fftw_complex *inout, *work;
  fftw_mpi_plan p;
  hpcc_fftw_mpi_plan ip;

  failure = 1;
  Gflops = -1.0;
  deps = HPL_dlamch( HPL_MACH_EPS );

  MPI_Comm_size( comm, &commSize );
  MPI_Comm_rank( comm, &commRank );

  /* there are two vectors of size 'n'/'commSize': inout, work, and internal work: 2*'n'/'commSize' */
  n = HPCC_LocalVectorSize( params, 4, sizeof(fftw_complex), 1 );

  /* make sure that the size of the vector size is greater than number of processes squared */
  if (n < commSize)
    n = 2 * commSize;

  n *= commSize; /* global vector size */

  t1 = -MPI_Wtime();
  p = fftw_mpi_create_plan( comm, n, FFTW_FORWARD, FFTW_MEASURE );
  t1 += MPI_Wtime();

  fftw_mpi_local_sizes( p, &locn, &loc0, &alocn, &aloc0, &tls );

  inout = (fftw_complex *)fftw_malloc( tls * (sizeof *inout) );
  work  = (fftw_complex *)fftw_malloc( tls * (sizeof *work) );

  if (! inout || ! work) goto comp_end;

  t0 = -MPI_Wtime();
  HPCC_bcnrand( 2 * tls, 53 * commRank * 2 * tls, inout );
  t0 += MPI_Wtime();

  t2 = -MPI_Wtime();
  fftw_mpi( p, 1, inout, work );
  t2 += MPI_Wtime();

  fftw_mpi_destroy_plan( p );

  ip = HPCC_fftw_mpi_create_plan( comm, n, FFTW_BACKWARD, FFTW_ESTIMATE );
  t3 = -MPI_Wtime();
  HPCC_fftw_mpi( ip, 1, inout, work );
  t3 += MPI_Wtime();
  HPCC_fftw_mpi_destroy_plan( ip );

  HPCC_bcnrand( 2 * tls, 53 * commRank * 2 * tls, work ); /* regenerate data */

  maxErr = 0.0;
  for (i = 0; i < tls; i++) {
    tmp1 = c_re( inout[i] ) - c_re( work[i] );
    tmp2 = c_im( inout[i] ) - c_im( work[i] );
    tmp3 = sqrt( tmp1*tmp1 + tmp2*tmp2 );
    maxErr = maxErr >= tmp3 ? maxErr : tmp3;
  }
  if (maxErr / log(n) / deps < params->test.thrsh) failure = 0;

  if (doIO) {
    fprintf( outFile, "Number of nodes: %d\n", commSize );
    fprintf( outFile, "Vector size: %20.0f\n", tmp1 = (double)n );
    fprintf( outFile, "Generation time: %9.3f\n", t0 );
    fprintf( outFile, "Tuning: %9.3f\n", t1 );
    fprintf( outFile, "Computing: %9.3f\n", t2 );
    fprintf( outFile, "Inverse FFT: %9.3f\n", t3 );
    fprintf( outFile, "max(|x-x0|): %9.3e\n", maxErr );
  }

  if (t2 > 0.0) Gflops = 1e-9 * (5.0 * n * log(n) / log(2.0)) / t2;

  comp_end:

  if (work) fftw_free( work );
  if (inout) fftw_free( inout );

  *UGflops = Gflops;
  *Un = n;
  *UmaxErr = maxErr;
  *Ufailure = failure;
}

int
HPCC_MPIFFT(HPCC_Params *params) {
  int commRank, commSize;
  int i, procPow2, isComputing, doIO, failure;
  s64Int_t n;
  double Gflops = -1.0, maxErr = -1.0;
  MPI_Comm comm;
  FILE *outFile;

  MPI_Comm_size( MPI_COMM_WORLD, &commSize );
  MPI_Comm_rank( MPI_COMM_WORLD, &commRank );

  doIO = commRank == 0 ? 1 : 0;

  if (doIO) {
    outFile = fopen( params->outFname, "a" );
    if (! outFile) outFile = stderr;
  }

  /* Find power of two that is smaller or equal to number of processes */
  for (i = 0, procPow2 = 1;
       procPow2 <= commSize && i < (int)(sizeof(int) * 8 - 1);
       i++, procPow2 <<= 1)
    ;
  i--;
  procPow2 = 1 << i;

  isComputing = commRank < procPow2 ? 1 : 0;

  HPCC_fft_timings_forward = params->MPIFFTtimingsForward;
  HPCC_fft_timings_backward = params->MPIFFTtimingsBackward;

  if (commSize == procPow2)
    comm = MPI_COMM_WORLD;
  else
    MPI_Comm_split( MPI_COMM_WORLD, isComputing ? 0 : MPI_UNDEFINED, commRank, &comm );

  if (isComputing)
    MPIFFT0( params, doIO, outFile, comm, &Gflops, &n, &maxErr, &failure );

  if (commSize != procPow2 && isComputing && comm != MPI_COMM_NULL)
    MPI_Comm_free( &comm );

  params->MPIFFT_N = n;
  params->MPIFFT_maxErr = maxErr;

  MPI_Bcast( &Gflops, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );

  params->MPIFFTGflops = Gflops;

  params->FFTEnblk = FFTE_NBLK;
  params->FFTEnp = FFTE_NP;
  params->FFTEl2size = FFTE_L2SIZE;

  if (doIO) if (outFile != stderr) fclose( outFile );

  return 0;
}
