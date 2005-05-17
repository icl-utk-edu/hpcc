
#include <stdio.h>
#include <stdlib.h>

#include "hpccfft.h"

hpcc_fftw_plan
HPCC_fftw_create_plan(int n, fftw_direction dir, int flags) {
  hpcc_fftw_plan p;
  fftw_complex *a = NULL, *b = NULL;

  p = fftw_malloc( sizeof *p );

  p->w1 = malloc( (FFTE_NDA2/2 + FFTE_NP) * (sizeof *p->w1) );
  p->w2 = malloc( (FFTE_NDA2/2 + FFTE_NP) * (sizeof *p->w2) );
  p->ww = malloc( ((FFTE_NDA2+FFTE_NP) * 4 + FFTE_NP) * (sizeof *p->ww) );
  p->c = malloc( ((FFTE_NDA2+FFTE_NP) * (FFTE_NBLK + 1) + FFTE_NP) * (sizeof *p->c) );

  HPCC_zfft1d( n, a, b, 0, p );

  p->n = n;
  p->dir = dir;
  p->flags = flags;

  return p;
}

void
HPCC_fftw_destroy_plan(hpcc_fftw_plan p) {
  if (! p) return;
  fftw_free( p->c );
  fftw_free( p->ww );
  fftw_free( p->w2 );
  fftw_free( p->w1 );
  fftw_free( p );
}

void
HPCC_fftw_one(hpcc_fftw_plan p, fftw_complex *in, fftw_complex *out) {
  if (FFTW_FORWARD == p->dir)
    HPCC_zfft1d( p->n, in, out, -1, p );
  else
    HPCC_zfft1d( p->n, in, out, +1, p );
}
