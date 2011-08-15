
#include <math.h>

typedef double floating_point;
typedef long s64Int_t;

typedef struct fftw_complex_s {
  int v;
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

extern double *ARR1D(fftw_complex_ptr, int i);
extern double *ARR2D(fftw_complex_ptr, int i, int j, int ld);
extern double *ARR3D(fftw_complex_ptr, int i, int j, int k, int ld1, int ld2);

extern fftw_complex_ptr PTR1D(fftw_complex_ptr v, int i, fftw_complex_ptr tmp);
extern fftw_complex_ptr PTR2D(fftw_complex_ptr v, int i, int j, int ld, fftw_complex_ptr tmp);

extern void c_assgn(char *format, ...);
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
