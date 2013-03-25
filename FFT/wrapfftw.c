
#include <stdio.h>
#include <stdlib.h>

#include "hpccfft.h"

#ifdef _OPENMP
#include <omp.h>
#endif

hpcc_fftw_plan
HPCC_fftw_create_plan(int n, fftw_direction dir, int flags) {
  hpcc_fftw_plan p;
  fftw_complex *a = NULL, *b = NULL;
  size_t w1_size, w2_size, ww1_size, ww2_size, ww3_size, ww4_size;

  p = (hpcc_fftw_plan)fftw_malloc( sizeof *p );
  if (! p) return p;

  w1_size = Mmax( FFTE_NDA2/2 + FFTE_NP, (int)(1.100 * sqrt( n )) );
  w2_size = Mmax( FFTE_NDA2/2 + FFTE_NP, (int)(0.375 * sqrt( n )) );
  ww1_size = Mmax( FFTE_NDA2 + FFTE_NDA4*FFTE_NP + FFTE_NP, (int)(1.0 * sqrt( n )) );
  ww2_size = Mmax( FFTE_NDA2 + FFTE_NDA4*FFTE_NP + FFTE_NP, (int)(3.9 * sqrt( n )) );
  ww3_size = Mmax( FFTE_NDA2 + FFTE_NDA4*FFTE_NP + FFTE_NP, (int)(5.4773 * sqrt( n )) );
  ww4_size = Mmax( FFTE_NDA2 + (1 << 13), (int)(1.0/256.0 * n) );

  p->w1 = (fftw_complex *)fftw_malloc( w1_size * (sizeof *p->w1) );
  p->w2 = (fftw_complex *)fftw_malloc( w2_size * (sizeof *p->w2) );
  p->ww1 = (fftw_complex *)fftw_malloc( ww1_size * (sizeof *p->ww1) );
  p->ww2 = (fftw_complex *)fftw_malloc( ww2_size * (sizeof *p->ww1) );
  p->ww3 = (fftw_complex *)fftw_malloc( ww3_size * (sizeof *p->ww1) );
  p->ww4 = (fftw_complex *)fftw_malloc( ww4_size * (sizeof *p->ww1) );

  p->c_size = Mmax( (FFTE_NDA2+FFTE_NP) * FFTE_NBLK + FFTE_NP, (int)(16.75 * sqrt( n )) );
  p->d_size = Mmax( FFTE_NDA2+FFTE_NP, (int)(1.0 * sqrt( n )) );
#ifdef _OPENMP
#pragma omp parallel
  {
#pragma omp single
    {
      int i;
      i = omp_get_num_threads();
      p->c = (fftw_complex *)fftw_malloc( p->c_size * (sizeof *p->c) * i );
      p->d = (fftw_complex *)fftw_malloc( p->d_size * (sizeof *p->d) * i );
    }
  }
#else
  p->c = (fftw_complex *)fftw_malloc( p->c_size * (sizeof *p->c) );
  p->d = (fftw_complex *)fftw_malloc( p->d_size * (sizeof *p->d) );
#endif

  if (! p->w1 || ! p->w2 || ! p->ww1 || ! p->ww2 || ! p->ww3 || ! p->ww4 || ! p->c || ! p->d) {
    if (p->d) fftw_free( p->d );
    if (p->c) fftw_free( p->c );
    if (p->ww4) fftw_free( p->ww4 );
    if (p->ww3) fftw_free( p->ww3 );
    if (p->ww2) fftw_free( p->ww2 );
    if (p->ww1) fftw_free( p->ww1 );
    if (p->w2) fftw_free( p->w2 );
    if (p->w1) fftw_free( p->w1 );
    fftw_free( p );
    return NULL;
  }

  HPCC_zfft1d( n, a, b, 0, p );

  p->n = n;
  p->dir = dir;
  p->flags = flags;

  return p;
}

void
HPCC_fftw_destroy_plan(hpcc_fftw_plan p) {
  if (! p) return;
  fftw_free( p->d );
  fftw_free( p->c );
  fftw_free( p->ww4 );
  fftw_free( p->ww3 );
  fftw_free( p->ww2 );
  fftw_free( p->ww1 );
  fftw_free( p->w2 );
  fftw_free( p->w1 );
  fftw_free( p );
}

/* Without additional storage of size p->n there is no way to preserve FFTW 2
   semantics (the `in' vector is not modified). But it doesn't matter for the
   calling code: it doesn't rely on this semantics. The change in semantics
   occured while going from FFTE 3.3 to FFTE 4.0. */
void
HPCC_fftw_one(hpcc_fftw_plan p, fftw_complex *in, fftw_complex *out) {
  int i, n;

  if (FFTW_FORWARD == p->dir)
    HPCC_zfft1d( p->n, in, out, -1, p );
  else
    HPCC_zfft1d( p->n, in, out, +1, p );

  n = p->n;
  /* Copy the transform to `out' vector. */
  for (i = 0; i < n; ++i) {
    c_assgn( out[i], in[i] );
  }
}
