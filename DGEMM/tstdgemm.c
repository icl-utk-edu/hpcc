/* -*- mode: C; tab-width: 2; indent-tabs-mode: nil; fill-column: 79; coding: iso-latin-1-unix -*- */
/* tstdgemm.c
 */

#include <hpcc.h>

/* Generates random matrix with entries between 0.0 and 1.0 */
static void
dmatgen(int m, int n, double *a, int lda) {
  int i, j;
  double *a0 = a, rcp = 1.0 / RAND_MAX;
  time_t tt;

  srand( time( &tt ) );

  for (j = 0; j < n; j++) {
    for (i = 0; i < m; i++)
      a0[i] = rcp * rand();

    a0 += lda;
  }
}

/* Copies a matrix */
static void
dmatcpy(int m, int n, double *dst, int ldd, double *src, int lds) {
  int i, j;
  double *dst0 = dst, *src0 = src;

  for (j = 0; j < n; j++) {
    for (i = 0; i < m; i++)
      dst0[i] = src0[i];

    dst0 += ldd;
    src0 += lds;
  }
}

/* Test output of DGEMM by performing 'n_' * 'k_' operations (and some cache trashing).
   In each column of 'c' a random entry is computed (using entry from 'd') and compared to the
   entry in 'c'. Sum of absolute values of differences is used to caluclate a scaled residual.
 */
static void
tst_dgemm(int m_, int n_, int k_, double alpha, double *a, int lda, double *b, int ldb,
    double beta, double *c, int ldc, double *d, int ldd, double *Usres) {
  int i, j, k;
  double *a0, *b0 = b, *c0 = c, *d0 = d, tmp, sum = 0.0, sumc, sumd, norm;
  double deps = HPL_dlamch( HPL_MACH_EPS );

  /* Calculate Frobenius norm of 'c' and 'd' */
  sumc = sumd = 0.0;
  for (j = 0; j < n_; j++) {
    for (i = 0; i < m_; i++) {
      sumc += c0[i] * c0[i];
      sumd += d0[i] * d0[i];
    }
    c0 += ldc;
    d0 += ldd;
  }

  norm = sumc > sumd ? sumc : sumd;
  norm = sqrt(norm);

  c0 = c;
  d0 = d;
  sum = 0.0;
  for (j = 0; j < n_; j++) { /* for each column of 'c' */
    i = rand() % m_;

    a0 = a;
    tmp = 0.0;
    for (k = 0; k < k_; k++) {
      tmp += a0[i] * b0[k];

      a0 += lda;
    }
    tmp = beta * d0[i] + alpha * tmp;
    sum += fabs(tmp - c0[i]); /* this should be zero (or around epsilon) */

    b0 += ldb;
    c0 += ldc;
    d0 += ldd;
  }

  if (Usres) *Usres = sum / norm / ((double)(m_ + n_ + k_) / 3.0) / deps;
}

int
HPCC_TestDGEMM(HPCC_Params *params, int doIO, double *UGflops, int *Un, int *Ufailure) {
  int n, failure = 1;
  double *a, *b, *c0, *c1, alpha, beta, sres;
  double Gflops = 0.0, dn, t0, t1;
  long l_n;
  FILE *outFile;

  if (doIO) {
    outFile = fopen( params->outFname, "a" );
    if (! outFile) {
      outFile = stderr;
      fprintf( outFile, "Cannot open output file.\n" );
      return 1;
    }
  }

  n = (int)sqrt( params->HPLMaxProcMem / sizeof(double) / 4 );
  if (n < 0) n = -n; /* if 'n' has overflown an integer */
  l_n = n;

  a = XMALLOC( double, l_n * l_n );
  b = XMALLOC( double, l_n * l_n );
  c0 = XMALLOC( double, l_n * l_n );
  c1 = XMALLOC( double, l_n * l_n );

  if (! a || ! b || ! c0 || ! c1) {
    goto comp_end;
  }

  dmatgen( n, n, a, n );
  dmatgen( n, n, b, n );
  dmatgen( n, n, c0, n );
  dmatcpy( n, n, c1, n, c0, n );

  alpha = a[n / 2];
  beta  = b[n / 2];

  t0 = MPI_Wtime();
  HPL_dgemm( HplColumnMajor, HplNoTrans, HplNoTrans, n, n, n, alpha, a, n, b, n, beta, c0, n );
  t1 = MPI_Wtime();

  t1 -= t0;
  dn = (double)n;
  if (t1 != 0.0 && t1 != -0.0)
    Gflops = 2.0e-9 * dn * dn * dn / t1;
  else
    Gflops = 0.0;

  tst_dgemm( n, n, n, alpha, a, n, b, n, beta, c0, n, c1, n, &sres );

  if (doIO) fprintf( outFile, "Scaled residual: %g\n", sres );

  if (sres < params->test.thrsh)
    failure = 0;

  comp_end:

  if (c1) free( c1 );
  if (c0) free( c0 );
  if (b) free( b );
  if (a) free( a );

  if (doIO) {
    fflush( outFile );
    fclose( outFile );
  }

  if (UGflops) *UGflops = Gflops;
  if (Un) *Un = n;
  if (Ufailure) *Ufailure = failure;

  return 0;
}
