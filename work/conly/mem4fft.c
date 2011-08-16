/* -*- mode: C; tab-width: 2; indent-tabs-mode: nil; fill-column: 79; coding: iso-latin-1-unix -*- */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "hpccfft.h"

double G_val;
char G_anon[10] = "ANONYMOUS";

int
loop_start(int upper_bound) {
  if (upper_bound > 1)
    return upper_bound-1;
  return upper_bound;
}

fftw_complex_ptr
copy_ctor(fftw_complex_ptr self, int idx, fftw_complex_ptr other) {
  self->largest_idx = other->largest_idx;
  self->last_idx = other->last_idx;
  self->offset = other->offset + idx;
  self->name = other->name;
  self->size = other->size;

  return other;
}

fftw_complex_ptr
PTR1D(fftw_complex_ptr v, int i, fftw_complex_ptr tmp) {
  return copy_ctor( tmp, i, v );
}

fftw_complex_ptr
PTR2D(fftw_complex_ptr v, int i, int j, int ld, fftw_complex_ptr tmp) {
  return copy_ctor( tmp, i + j * ld, v );
}

void
c_mul3v(fftw_complex v1, fftw_complex v2, fftw_complex v3) {
}

void
fftw_complex_show_stat(fftw_complex_ptr p) {
  printf( "%s[%d] %d\n", p->name, p->size, *(p->largest_idx) );
}

fftw_complex_ptr
fftw_complex_new(char *name, int size) {
  fftw_complex_ptr self;

  self = malloc( sizeof *self );

  self->name = name;
  self->size = size;

  self->largest_idx = &(self->largest_idx_storage);
  self->last_idx = &(self->last_idx_storage);
  *(self->largest_idx) = 0;
  self->offset = 0;

  return self;
}

void
fftw_complex_delete(fftw_complex_ptr self) {
  free( self );
}

static void
plan_init(hpcc_fftw_plan p, int n) {
  p->n = n;
  p->w1 = fftw_complex_new("w1", n);
  p->w2 = fftw_complex_new("w2", n);
  p->ww = fftw_complex_new("ww", n);
  p->ww2 = fftw_complex_new("ww2", n);
  p->ww3 = fftw_complex_new("ww3", n);
  p->ww4 = fftw_complex_new("ww4", n);
  p->c = fftw_complex_new("c", n);
  p->d = fftw_complex_new("d", n);
}

static void
plan_finalize(hpcc_fftw_plan p) {
  fftw_complex_delete( p->w1 );
  fftw_complex_delete( p->w2 );
  fftw_complex_delete( p->ww );
  fftw_complex_delete( p->ww2 );
  fftw_complex_delete( p->ww3 );
  fftw_complex_delete( p->ww4 );
  fftw_complex_delete( p->c );
  fftw_complex_delete( p->d );
}

void
plan_show_stats(hpcc_fftw_plan p) {
  fftw_complex_show_stat( p->w1 );
  fftw_complex_show_stat( p->w2 );
  fftw_complex_show_stat( p->ww );
  fftw_complex_show_stat( p->ww2 );
  fftw_complex_show_stat( p->ww3 );
  fftw_complex_show_stat( p->ww4 );
  fftw_complex_show_stat( p->c );
  fftw_complex_show_stat( p->d );
}

static void
mem4fft(int n) {
  struct hpcc_fftw_plan_s ps;
  fftw_complex_ptr in, out;

  in  = fftw_complex_new("INPUT",  n);
  out = fftw_complex_new("OUTPUT", n);

  plan_init( &ps, n );

  HPCC_zfft1d( n, in, out, 1, &ps );

  fftw_complex_show_stat(in);
  fftw_complex_show_stat(out);

  plan_show_stats( &ps );

  plan_finalize( &ps );

  fftw_complex_delete( out );
  fftw_complex_delete( in );
}

static void
enumerate_all() {
  int i2, i3, i5, n, nlast=1, factors[3];

  for (i5 = 0; i5 <= 13; ++i5)
    for (i3 = 0; i3 <= 19; ++i3)
      for (i2 = 0; i2 <= 30; ++i2) {
        if (i2 == 0 && i3 == 0 && i5 == 0) continue;

        n = HPCC_ipow( 2, i2 ) * HPCC_ipow( 3, i3 ) * HPCC_ipow( 5, i5 );

        /* in case of overflow the loop is aborted with "break":
         * bigger values of "i2" will cause overflow as well */

        /* check for overflow */

        /* overflow created a negative number */
        if (n <= 0) break;

        /* overflow created a smaller number  in subsequent iteration*/
        if (i2 > 0 && n < nlast) break;

        /* overflow created a number not divisible but neither 2, 3, nor 5 */
        if (n % 2 && n % 3 && n % 5) break;

        HPCC_factor235( n, factors );
        /* overflow created a number that has factors other than 2, 3, and 5 */
        if (factors[0] != i2 || factors[1] != i3 || factors[2] != i5) break;

        printf( "F: %d = 2^%d * 3^%d * 5^%d\n", n, factors[0], factors[1], factors[2] );

        mem4fft( n );

        nlast = n;
      }
}

int
main(int argc, char *argv[]) {
  int i, n = 1, cp, nf[3] = {2, 3, 5};

  /* first argument starts with A (as in All)? */
  if (argc == 2 && argv[1][0] == 'A')
    enumerate_all();

  for (i = 0; i < 3; ++i) {
    if (i >= argc-1)
      break;
    if (sscanf( argv[i+1], "%d", &cp ) > 0 && cp > 0) {
      n *= HPCC_ipow( nf[i], cp );
    }
  }

  if (n < 2) n = 2;

  mem4fft( n );

  return 0;
}
