/* -*- mode: C; tab-width: 2; indent-tabs-mode: nil; fill-column: 79; coding: iso-latin-1-unix -*- */
/* mpifft.c
 */

#include <hpcc.h>

#include "hpccfft.h"
#include "wrapmpifftw.h"

static int
LocalVectorSize(long maxCount) {
  long ln;
  int i, n, maxIntBits;

  /* this is the maximum power of 2 that that can be held in a signed integer (for a 4-byte
     integer, 2**31-1 is the maximum integer, so the maximum power of 2 is 30) */
  maxIntBits = sizeof(int) * 8 - 2;

  /* Find the largest size of a vector. The size of the vector has to be:  a power 2, fit in
     an integer, and small enough to fit in maxCount. */

  for (i = 0, ln = 1; ln <= maxCount && i <= maxIntBits; i++, ln <<= 1)
    ;
  ln >>= 1; /* there was one shift too much */

  if (i <= maxIntBits) n = ln;
  else n = 1 << maxIntBits;

  return n;
}

static void
MPIFFT0(long HPLMaxProcMem, double HPLthshr, int doIO, FILE *outFile, MPI_Comm comm,
        double *UGflops, int *Un, int *Ufailure) {
  int commRank, commSize;
  int i, n, mult, failure = 1, maxVal = 1 << (sizeof(int) * 8 - 2);
  int locn, loc0, alocn, aloc0, tls;
  double maxErr, tmp1, tmp2, tmp3, t0, t1, t2, t3, Gflops = -1.0;
  double deps = HPL_dlamch( HPL_MACH_EPS );
  fftw_complex *inout, *work;
  fftw_mpi_plan p;
  hpcc_fftw_mpi_plan ip;

  MPI_Comm_size( comm, &commSize );
  MPI_Comm_rank( comm, &commRank );

  /* there are three vectors of size 'n'/'commSize': inout, work, and internal work */
  n = LocalVectorSize( HPLMaxProcMem / 3 / sizeof(fftw_complex) );

  /* 'n' should be multiplied by the number of nodes but that may overflow an integer */
  mult = maxVal / n;
  mult = mult <= commSize ? mult : commSize;

  n *= mult;


  t1 = MPI_Wtime();
  p = fftw_mpi_create_plan( comm, n, FFTW_FORWARD, FFTW_MEASURE );
  t1 = MPI_Wtime() - t1;

  fftw_mpi_local_sizes( p, &locn, &loc0, &alocn, &aloc0, &tls );

  inout = fftw_malloc( tls * (sizeof *inout) );
  work  = fftw_malloc( tls * (sizeof *work) );

  if (! inout || ! work) goto comp_end;

  t0 = MPI_Wtime();
  bcnrand( 2 * tls, 0, inout );
  t0 = MPI_Wtime() - t0;

  t2 = MPI_Wtime();
  fftw_mpi( p, 1, inout, work );
  t2 = MPI_Wtime() - t2;

  fftw_mpi_destroy_plan( p );

  ip = hpcc_fftw_mpi_create_plan( comm, n, FFTW_BACKWARD, FFTW_ESTIMATE );
  t3 = MPI_Wtime();
  hpcc_fftw_mpi( ip, 1, inout, work );
  t3 = MPI_Wtime() - t3;
  hpcc_fftw_mpi_destroy_plan( ip );

  bcnrand( 2 * tls, 0, work ); /* regenerate data */

  maxErr = 0.0;
  for (i = 0; i < tls; i++) {
    tmp1 = c_re( inout[i] ) - c_re( work[i] );
    tmp2 = c_im( inout[i] ) - c_im( work[i] );
    tmp3 = sqrt( tmp1*tmp1 + tmp2*tmp2 );
    maxErr = maxErr >= tmp3 ? maxErr : tmp3;
  }
  if (maxErr / log(n) / deps < HPLthshr) failure = 0;

  if (doIO) {
    fprintf( outFile, "Vector size: %d\n", n );
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
  *Ufailure = failure;
}

int
MPIFFT(HPCC_Params *params) {
  int commRank, commSize;
  int i, n, procPow2, isComputing, doIO, failure;
  double Gflops = -1.0;
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

  if (commSize == procPow2)
    comm = MPI_COMM_WORLD;
  else
    MPI_Comm_split( MPI_COMM_WORLD, isComputing ? 0 : MPI_UNDEFINED, commRank, &comm );

  if (isComputing)
    MPIFFT0( params->HPLMaxProcMem, params->test.thrsh, doIO, outFile, comm, &Gflops, &n,
	     &failure );

  params->MPIFFTGflops = Gflops;

  if (doIO) if (outFile != stderr) fclose( outFile );

  MPI_Barrier( MPI_COMM_WORLD );

  return 0;
}
