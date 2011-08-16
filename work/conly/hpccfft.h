
#include <math.h>

typedef double floating_point;
typedef long s64Int_t;

typedef struct fftw_complex_s {
  double v;
} fftw_complex;

typedef struct fftw_complex_ptr_s {
  char *name;
  int size, offset;
  int *largest_idx, *last_idx;
  int largest_idx_storage, last_idx_storage;
} *fftw_complex_ptr;

typedef struct hpcc_fftw_plan_s {
 int n;
 fftw_complex_ptr w1, w2, ww, ww2, ww3, ww4, c, d;
} *hpcc_fftw_plan;

extern double G_val;

static double *
arraccss(fftw_complex_ptr this, int idx) {
  *(this->last_idx) = this->offset + idx;
  if (*(this->largest_idx) < idx)
    *(this->largest_idx) = this->offset + idx;

  return &G_val;
}

static double *
ARR1D(fftw_complex_ptr a, int i) {
  return arraccss( a, i );
}

static double *
ARR2D(fftw_complex_ptr a, int i, int j, int ld) {
  return arraccss( a, i + j * ld );
}

static double *
ARR3D(fftw_complex_ptr a, int i, int j, int k, int ld1, int ld2) {
  return arraccss( a, i + j * ld1 + k * ld1 * ld2 );
}

static void
c_assgn_aa(double *d, double *s) {
  *d = *s;
}

static void
c_assgn_av(double *d, fftw_complex s) {
  *d = s.v;
}

static void
c_assgn_va(fftw_complex s, double *d) {
  *d = s.v;
}

extern fftw_complex_ptr PTR1D(fftw_complex_ptr v, int i, fftw_complex_ptr tmp);
extern fftw_complex_ptr PTR2D(fftw_complex_ptr v, int i, int j, int ld, fftw_complex_ptr tmp);

extern void c_mul3v(fftw_complex v1, fftw_complex v2, fftw_complex v3);

extern int HPCC_factor235_8(s64Int_t n, int *ip);
extern int HPCC_zfft1d(int n, fftw_complex_ptr a, fftw_complex_ptr b, int iopt, hpcc_fftw_plan p);
extern int HPCC_fft235(fftw_complex_ptr a, fftw_complex_ptr b, fftw_complex_ptr w, int n, const int *ip);
extern int HPCC_factor235(int n, int *ip);
extern int HPCC_settbl(fftw_complex_ptr w, int n);
extern int HPCC_ipow(int x, int y);

extern int loop_start(int upper_bound);

#define c_re(x) *(x)
#define c_im(x) *(x)

#define V2MIN(r, v) r = (v) < r ? (v) : r

#define FFTE_NP 8
#define FFTE_NBLK 16
#define FFTE_L2SIZE 1048576
#define FFTE_NDA2 65536
