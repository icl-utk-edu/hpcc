/* -*- mode: C; tab-width: 2; indent-tabs-mode: nil; fill-column: 79; coding: iso-latin-1-unix -*- */
/* tstfft.c
 */

#include <hpcc.h>

#include "hpccfft.h"

static int
TestFFT1(long HPLMaxProcMem, double HPLthshr, int doIO, FILE *outFile,
         double *UGflops, int *Un, int *Ufailure) {
  fftw_complex *in, *out;
  fftw_plan p;
  hpcc_fftw_plan ip;
  double Gflops = -1.0;
  double maxErr, tmp1, tmp2, tmp3, t0, t1, t2, t3;
  int i, n, maxIntBits, failure = 1;
  long ln, maxCount = HPLMaxProcMem/2/sizeof(fftw_complex);
  double deps = HPL_dlamch( HPL_MACH_EPS );

  /* this is the maximum power of 2 that that can be held in a signed integer (for a 4-byte
     integer, 2**31-1 is the maximum integer, so the maximum power of 2 is 30) */
  maxIntBits = sizeof(int) * 8 - 2;

  /* Find the largest size of two (same length) complex vectors. The size of each vector has
     to be:  a power 2 and fit in an integer. Combined size of the vectors has to be small
     enough to fit in main memory (taken from HPL). */

  for (i = 0, ln = 1; ln <= maxCount; i++, ln <<= 1)
    ;
  ln >>= 1; /* there was one shift too much */

  if (i <= maxIntBits) n = ln;
  else n = 1 << maxIntBits;

  /* need to use fftw_malloc() so that the returned pointers will be aligned properly for SSE
     instructions on Intel/AMD systems */
  in  = (fftw_complex *)fftw_malloc( (sizeof *in)  * n );
  out = (fftw_complex *)fftw_malloc( (sizeof *out) * n );

  if (! in || ! out) goto comp_end;

  t0 = -MPI_Wtime();
  HPCC_bcnrand( 2*n, 0, in );
  t0 += MPI_Wtime();

  t1 = -MPI_Wtime();
  p = fftw_create_plan( n, FFTW_FORWARD, FFTW_MEASURE );
  t1 += MPI_Wtime();

  t2 = -MPI_Wtime();
  fftw_one( p, in, out );
  t2 += MPI_Wtime();

  fftw_destroy_plan(p);

  ip = HPCC_fftw_create_plan( n, FFTW_BACKWARD, FFTW_ESTIMATE );
  t3 = -MPI_Wtime();
  HPCC_fftw_one( ip, out, in );
  t3 += MPI_Wtime();
  HPCC_fftw_destroy_plan( ip );

  HPCC_bcnrand( 2*n, 0, out ); /* regenerate data */
  maxErr = 0.0;
  for (i = 0; i < n; i++) {
    tmp1 = c_re( in[i] ) - c_re( out[i] );
    tmp2 = c_im( in[i] ) - c_im( out[i] );
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

  if (out) fftw_free( out );
  if (in)  fftw_free( in );

  *UGflops = Gflops;
  *Un = n;
  *Ufailure = failure;

  return 0;
}

int
HPCC_TestFFT(HPCC_Params *params, int doIO, double *UGflops, int *Un, int *Ufailure) {
  int rv, n, failure = 1;
  double Gflops;
  FILE *outFile;

  if (doIO) {
    outFile = fopen( params->outFname, "a" );
    if (! outFile) {
      outFile = stderr;
      fprintf( outFile, "Cannot open output file.\n" );
      return 1;
    }
  }

  n = 0;
  Gflops = -1.0;
  rv = TestFFT1( params->HPLMaxProcMem, params->test.thrsh, doIO, outFile,
                 &Gflops, &n, &failure );

  if (doIO) {
    fflush( outFile );
    fclose( outFile );
  }

  if (UGflops) *UGflops = Gflops;
  if (Un) *Un = n;
  if (Ufailure) *Ufailure = failure;

  return rv;
}
