
#include <stdio.h>
#include <stdlib.h>

#include "hpccfft.h"

hpcc_fftw_plan
hpcc_fftw_create_plan(int n, fftw_direction dir, int flags) {
  hpcc_fftw_plan p;
  int zero = 0;
  fftw_complex *a = NULL, *b = NULL;

  p = fftw_malloc( sizeof *p );

  p->n0 = 0;

  p->c = fftw_malloc( ((NDA2+NP) * NBLK + NP) * sizeof(fftw_complex) );
  p->d = fftw_malloc( (NDA2 + NP) * sizeof(fftw_complex) );
  p->w1 = fftw_malloc( (NDA2/2+NP) * sizeof(fftw_complex) );
  p->w2 = fftw_malloc( (NDA2/2+NP) * sizeof(fftw_complex) );
  p->ww1 = fftw_malloc( (NDA2+NDA4*NP+NP) * sizeof(fftw_complex) );
  p->ww2 = fftw_malloc( (NDA2+NDA4*NP+NP) * sizeof(fftw_complex) );
  p->ww3 = fftw_malloc( (NDA2+NDA4*NP+NP) * sizeof(fftw_complex) );
  p->ww4 = fftw_malloc( (NDA2+NDA4*NP+NP) * sizeof(fftw_complex) );

  hpcc_zfft1d_( a, b, &n, &zero, p, &p->n0 );

  p->dir = dir;
  p->flags = flags;

  return p;
}

void hpcc_fftw_destroy_plan(hpcc_fftw_plan p) {
  if (! p) return;
  fftw_free( p->ww4 );
  fftw_free( p->ww3 );
  fftw_free( p->ww2 );
  fftw_free( p->ww1 );
  fftw_free( p->w2 );
  fftw_free( p->w1 );
  fftw_free( p->d );
  fftw_free( p->c );
  fftw_free( p );
}

void
hpcc_fftw_one(hpcc_fftw_plan p, fftw_complex *in, fftw_complex *out) {
  int one = 1, two = 2;

  if (FFTW_FORWARD == p->dir)
    hpcc_zfft1d_( in, out, &p->n0, &one, p, &p->n0 );
  else
    hpcc_zfft1d_( in, out, &p->n0, &two, p, &p->n0 );
}
