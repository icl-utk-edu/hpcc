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

double *
arraccss(fftw_complex_ptr this, int idx) {
  *(this->last_idx) = this->offset + idx;
  if (*(this->largest_idx) < idx)
    *(this->largest_idx) = this->offset + idx;

  return &G_val;
}

double *
ARR1D(fftw_complex_ptr a, int i) {
  return arraccss( a, i );
}

double *
ARR2D(fftw_complex_ptr a, int i, int j, int ld) {
  return arraccss( a, i + j * ld );
}

double *
ARR3D(fftw_complex_ptr a, int i, int j, int k, int ld1, int ld2) {
  return arraccss( a, i + j * ld1 + k * ld1 * ld2 );
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
c_assgn(char *name, ...) {
  if (name) return;
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

static void mem4fft(int n) {
  struct hpcc_fftw_plan_s ps;
  fftw_complex_ptr in, out;

  in  = fftw_complex_new("INPUT",  n);
  out = fftw_complex_new("OUTPUT", n);

  plan_init( &ps, n );

  HPCC_zfft1d( n, in, out, 1, &ps );

  fftw_complex_show_stat(in);
  fftw_complex_show_stat(out);

  plan_show_stats( &ps );
}

int
main(int argc, char *argv[]) {
  int i, n = 0, cp, nf[3] = {2, 3, 5};

  for (i = 0; i < 3; ++i) {
    if (i >= argc-1)
      break;
    if (sscanf( argv[i+1], "%d", &cp ) > 0 && cp > 0) {
      n += HPCC_ipow( nf[i], cp );
    }
  }

  if (n < 2) n = 2;

  mem4fft( n );

  return 0;
}
