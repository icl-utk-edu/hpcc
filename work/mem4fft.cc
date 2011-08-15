/* -*- mode: C; tab-width: 2; indent-tabs-mode: nil; fill-column: 79; coding: iso-latin-1-unix -*- */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "mem4fft.h"

impl_complex_t G_complex_val;
char G_anon[10] = "ANONYMOUS";

impl_complex_t&
fftw_complex::index_access(int idx) {
  *(this->last_idx) = this->offset + idx;
  if (*(this->largest_idx) < idx)
    *(this->largest_idx) = this->offset + idx;

  return G_complex_val;
}

impl_complex_t& fftw_complex::operator [](long unsigned int idx) {
  return this->index_access(idx);
}

impl_complex_t& fftw_complex::arr2d(int idx) {
  return this->index_access(idx);
}

impl_complex_t& fftw_complex::arr3d(int idx) {
  return this->index_access(idx);
}

fftw_complex *fftw_complex::ptr2d(int idx) {
  this->index_access(idx);
  return this;
}

impl_complex_t& fftw_complex::sqbracket(int idx) {
  return this->index_access(idx);
}

void fftw_complex::ctor(fftw_complex* self, const fftw_complex& other, int offset) {
  self->largest_idx = other.largest_idx;
  self->offset = other.offset + offset;
  self->last_idx = other.last_idx;
  self->name = other.name;
  self->size = other.size;
}

fftw_complex::fftw_complex(const fftw_complex& other) {
  ctor( this, other, other.offset );
}

fftw_complex::fftw_complex(fftw_complex& other, int offset) {
  ctor( this, other, offset );
}

fftw_complex::fftw_complex(fftw_complex* other, int offset) {
  ctor( this, *other, offset );
}

fftw_complex::fftw_complex(const char *name, int size) {
  this->largest_idx = &(this->largest_idx_storage);
  this->last_idx = &(this->last_idx_storage);
  *(this->largest_idx) = 0;
  this->name = name;
  this->size = size;
  this->offset = 0;
}

fftw_complex::fftw_complex(void) {
  this->largest_idx = &(this->largest_idx_storage);
  this->last_idx = &(this->last_idx_storage);
  this->name = G_anon;
  this->size = size;
  this->offset = 0;
}

void c_assgn(fftw_complex& d, fftw_complex&s) { }

void c_assgn(impl_complex_t& d, fftw_complex&s) { }

void c_assgn(impl_complex_t& d, impl_complex_t&s) { }

void c_assgn(fftw_complex& d, impl_complex_t&s) { }

void c_mul3v(fftw_complex& v1, fftw_complex& v2, fftw_complex& v3) {}

void
fftw_complex::show_stat(void) {
  printf( "%s[%d] %d\n", this->name, this->size, *(this->largest_idx) );
}

int
plan_init(hpcc_fftw_plan p, int n) {
  p->n = n;
  p->w1 = new fftw_complex("w1", n);
  p->w2 = new fftw_complex("w2", n);
  p->ww = new fftw_complex("ww", n);
  p->ww2 = new fftw_complex("ww2", n);
  p->ww3 = new fftw_complex("ww3", n);
  p->ww4 = new fftw_complex("ww4", n);
  p->c = new fftw_complex("c", n);
  p->d = new fftw_complex("d", n);

  c_re( ARR2D( p->w1, 0, 0, 0 ) ) = 0.0;

  return 0;
}

void
plan_show_stats(hpcc_fftw_plan p) {
  p->w1->show_stat();
  p->w2->show_stat();
  p->ww->show_stat();
  p->ww2->show_stat();
  p->ww3->show_stat();
  p->ww4->show_stat();
  p->c->show_stat();
  p->d->show_stat();
}

static void mem4fft(int n) {
  struct hpcc_fftw_plan_s ps;
  fftw_complex *in, *out;

  in = new fftw_complex("INPUT", n);
  out = new fftw_complex("INPUT", n);

  plan_init( &ps, n );

  HPCC_zfft1d( n, in, out, 0, &ps );

  in->show_stat();
  out->show_stat();

  plan_show_stats( &ps );
}

int
main(int argc, char *argv[]) {
  int i, j, n = 0, cp, nf[3] = {2, 3, 5};

  for (i = 0; i < 3; ++i) {
    if (i < argc-1 && sscanf( argv[i+1], "%d", &cp ) > 0 && cp > 0) {
      for (j = 1; cp ; --cp)
        j *= nf[i];
      n += j;
    }
  }

  mem4fft( n );

  return 0;
}
