
#include <math.h>

#define L2SIZE 524288
#define NP 8
#define NBLK 16
#define NDA2 32768
#define NDA3 1024
#define NDA4 256

typedef int integer;
typedef float real;
typedef double doublereal;
typedef struct {doublereal r, i;} doublecomplex;

#define REAL(x) (x)[0]
#define IMAG(x) (x)[1]
typedef double HPCC_Complex[2];

#ifdef LONG_IS_64BITS
typedef unsigned long u64Int_t;
typedef long s64Int_t;
#else
typedef unsigned long long u64Int_t;
typedef long long s64Int_t;
#endif

#include "wrapfftw.h"

/*
extern void Cfactor(int n, int *ip);
extern void Cfft235a(fftw_complex *a, fftw_complex *b, fftw_complex *w, int n, int *ip);
extern void Cfft235b(fftw_complex *a, fftw_complex *b, fftw_complex *w, int n, int *ip);
extern void Csettbl(fftw_complex *w, int n);
extern void Cfft2(fftw_complex *a, fftw_complex *b, int m);
extern void Cfft3a(fftw_complex *a, fftw_complex *b, fftw_complex *w, int l);
extern void Cfft3b(fftw_complex *a,fftw_complex *b,fftw_complex *w,int m,int l);
extern void Cfft4a(fftw_complex *a,fftw_complex *b,fftw_complex *w, int l);
extern void Cfft4b(fftw_complex *a,fftw_complex *b,fftw_complex *w, int m, int l);
extern void Cfft5a(fftw_complex *a,fftw_complex *b,fftw_complex *w, int l);
extern void Cfft5b(fftw_complex *a,fftw_complex *b,fftw_complex *w, int m, int l);
extern void Cfft8a(fftw_complex *a,fftw_complex *b,fftw_complex *w, int l);
extern void Cfft8b(fftw_complex *a,fftw_complex *b,fftw_complex *w, int m, int l);
*/
extern int hpcc_zfft1d_(void *a, void *b, integer *n, integer *iopt, hpcc_fftw_plan p, integer *n0);

extern int bcnrand(u64Int_t n, u64Int_t a, void *x);
