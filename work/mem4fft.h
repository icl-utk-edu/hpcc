
#define FFTE_NP 8
#define FFTE_NBLK 16
#define FFTE_L2SIZE 1048576
#define FFTE_NDA2 65536

#define ARR2D(a, i, j, lda) (a)->arr2d((i)+(j)*(lda))
#define PTR2D(a, i, j, lda) (a)->ptr2d((i)+(j)*(lda))
#define ARR3D(a, i, j, k, lda1, lda2) (a)->arr3d((i)+(lda1)*((j)+(k)*(lda2)))
#define PTR3D(a, i, j, k, lda1, lda2) (a+(i)+(lda1)*((j)+(k)*(lda2)))
#define ARR4D(a, i, j, k, l, lda1, lda2, lda3) a[(i)+(lda1)*((j)+(lda2)*((k)+(lda3)*(l)))]
#define V3MIN(r, e, v) r = (e); V2MIN(r, v)
#define V2MIN(r, v) r = (v) < r ? (v) : r
#define EMAX(d, v, e) d=(e); d=d>(v)?d:(v)

#define c_re(c)  ((c).real)
#define c_im(c)  ((c).imag)

typedef long unsigned int s64Int_t;

typedef struct impl_complex_s {
  double real, imag;
} impl_complex_t;

class fftw_complex {
  impl_complex_t& index_access(int idx);

public:
  int size, offset;
  int largest_idx_storage, last_idx_storage;
  /* making these pointers creates an easy to spot (out of bounds SEGFAULT) error detection */
  int *largest_idx, *last_idx;
  const char *name;

  void ctor(fftw_complex* self, const fftw_complex& other, int offset);
  fftw_complex(const fftw_complex& other);
  fftw_complex(fftw_complex& other, int offset);
  fftw_complex(fftw_complex* other, int offset);
  fftw_complex(const char *name, int n);
  fftw_complex(void);

  impl_complex_t& operator [](long unsigned int idx);

  impl_complex_t& arr2d(int idx);
  impl_complex_t& arr3d(int idx);
  fftw_complex *ptr2d(int idx);
  impl_complex_t& sqbracket(int idx);

  void show_stat();
};

void c_assgn(fftw_complex& d, fftw_complex& s);
void c_assgn(impl_complex_t& d, impl_complex_t& s);
void c_assgn(fftw_complex& d, impl_complex_t& s);
void c_assgn(impl_complex_t& d, fftw_complex& s);
void c_mul3v(fftw_complex& v1, fftw_complex& v2, fftw_complex& v3);

typedef struct hpcc_fftw_plan_s {
  int n;
  fftw_complex *w1, *w2, *ww, *ww2, *ww3, *ww4, *c, *d;
} *hpcc_fftw_plan;

extern int HPCC_fft235(fftw_complex *a, fftw_complex *b, fftw_complex *w, int n, const int *ip);

extern int HPCC_settbl(fftw_complex *w, int n);
extern int HPCC_factor235(int n, int *ip);
extern int HPCC_zfft1d(int n, fftw_complex *a, fftw_complex *b, int iopt, hpcc_fftw_plan p);
