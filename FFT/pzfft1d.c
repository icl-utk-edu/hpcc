/* pzfft1d.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "hpccfft.h"
#include "wrapmpifftw.h"
#define min(a,b) ((a)<=(b)?(a):(b))
#define max(a,b) ((a)>=(b)?(a):(b))

/* Table of constant values */

static integer c__2 = 2;
static integer c__3 = 3;
static integer c__5 = 5;

static int
pow_ii(int *x_, int *p_) {
  int i, r, x = *x_, p = *p_;

  if (1 == x || 0 == x) return x;
  if (0 == p) return 1;
  if (-1 == x) return (p & 1) ? -1 : 1;
  if (p < 0) return 0;
  r = 1;
  for (i = 0; i < p; i++) r *= x;
  return r;
}
static void
d_cnjg(doublecomplex *r, doublecomplex *v) {
  r->r = v->r; r->i = -v->i;
}

/*     FFTE: A FAST FOURIER TRANSFORM PACKAGE */

/*     (C) COPYRIGHT SOFTWARE, 2000-2004, ALL RIGHTS RESERVED */
/*                BY */
/*         DAISUKE TAKAHASHI */
/*         INSTITURE OF INFORMATION SCIENCES AND ELECTRONICS */
/*         UNIVERSITY OF TSUKUBA */
/*         1-1-1 TENNODAI, TSUKUBA, IBARAKI 305-8573, JAPAN */
/*         E-MAIL: daisuke@is.tsukuba.ac.jp */


/*     1-DIMENSIONAL COMPLEX FFT ROUTINE (MPI VERSION) */

/*     FORTRAN77 SOURCE PROGRAM */

/*     CALL PZFFT1D(A,B,WW,WWW,N,ME,NPU,IOPT) */

/*     A(N/NPU) IS COMPLEX INPUT/OUTPUT VECTOR (COMPLEX*16) */
/* !HPF$ DISTRIBUTE A(CYCLIC) */
/*     B(N/NPU) IS WORK VECTOR (COMPLEX*16) */
/*     WW(NDA3*NDA3) IS SINE/COSINE TABLE (COMPLEX*16) */
/*     WWW(N/NPU) IS SINE/COSINE TABLE (COMPLEX*16) */
/*     N IS THE LENGTH OF THE TRANSFORMS (INTEGER*8) */
/*       ----------------------------------- */
/*         N = (2**IP) * (3**IQ) * (5**IR) */
/*       ----------------------------------- */
/*     ME IS PROCESSOR NUMBER (INTEGER*4) */
/*     NPU IS THE NUMBER OF PROCESSORS (INTEGER*4) */
/*     IOPT = 0 FOR INITIALIZE FFT COEFFICIENTS */
/*          = 1 FOR FORWARD TRANSFORM WHERE INPUT AND OUTPUT ARE IN */
/*              BLOCK DISTRIBUTION */
/*          = 2 FOR INVERSE TRANSFORM WHERE INPUT AND OUTPUT ARE IN */
/*              BLOCK DISTRIBUTION */
/*          = 3 FOR FORWARD TRANSFORM WHERE INPUT AND OUTPUT ARE IN */
/*              CYCLIC DISTRIBUTION */
/*          = 4 FOR INVERSE TRANSFORM WHERE INPUT AND OUTPUT ARE IN */
/*              CYCLIC DISTRIBUTION */
/*          = 5 FOR FORWARD TRANSFORM WHERE INPUT IS IN */
/*              BLOCK DISTRIBUTION AND OUTPUT IS IN CYCLIC DISTRIBUTION */
/*          = 6 FOR INVERSE TRANSFORM WHERE INPUT IS IN */
/*              BLOCK DISTRIBUTION AND OUTPUT IS IN CYCLIC DISTRIBUTION */
/*          = 7 FOR FORWARD TRANSFORM WHERE INPUT IS IN */
/*              CYCLIC DISTRIBUTION AND OUTPUT IS IN BLOCK DISTRIBUTION */
/*          = 8 FOR INVERSE TRANSFORM WHERE INPUT IS IN */
/*              CYCLIC DISTRIBUTION AND OUTPUT IS IN BLOCK DISTRIBUTION */

/*     WRITTEN BY DAISUKE TAKAHASHI */

/* Subroutine */ int pzfft1d_(void *a_, void *b_, 
	void *ww_, void *www_, s64Int_t *n, integer *me, 
	integer *npu, integer *iopt, hpcc_fftw_mpi_plan p, s64Int_t *n0)
{
    /* System generated locals */
    integer i__1, i__2;
    doublecomplex z__1, z__2;

    /* Local variables */
    doublecomplex *c__, *d__;
    integer i__;
    doublereal dn;
    s64Int_t nn;
    integer ip[3], nx, ny, nz;
    doublecomplex *wx, *wy, *wz;
    integer lnn[3], lnx[3], lny[3], lnz[3], lnpu[3];
    doublecomplex *a, *b, *ww, *www;
    extern /* Subroutine */ int pzb2c_(doublecomplex *, doublecomplex *, 
	    doublecomplex *, s64Int_t *, integer *, hpcc_fftw_mpi_plan), pzc2b_(doublecomplex *, 
	    doublecomplex *, doublecomplex *, s64Int_t *, integer *, hpcc_fftw_mpi_plan),
            factor_(integer *, integer *), factor8_(s64Int_t *, integer *),
      settbl_(doublecomplex *, integer *), 
	    settbl3_(void *, integer *, integer *), pzfft1d0_(
	    doublecomplex *, doublecomplex *, doublecomplex *, doublecomplex *
	    , doublecomplex *, doublecomplex *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, doublecomplex *
	    , doublecomplex *, doublecomplex *, integer *, integer *, integer 
	    *, integer *, hpcc_fftw_mpi_plan), settblp_(void *, integer *, integer *, 
	    integer *, integer *, integer *);

    a = (doublecomplex *)a_;
    b = (doublecomplex *)b_;
    ww = (doublecomplex *)ww_;
    www = (doublecomplex *)www_;
    c__ = (doublecomplex *)p->c;
    d__ = (doublecomplex *)p->d;
    wx = (doublecomplex *)p->wx;
    wy = (doublecomplex *)p->wy;
    wz = (doublecomplex *)p->wz;

/*     FFTE: A FAST FOURIER TRANSFORM PACKAGE */

/*     (C) COPYRIGHT SOFTWARE, 2000-2004, ALL RIGHTS RESERVED */
/*                BY */
/*         DAISUKE TAKAHASHI */
/*         INSTITURE OF INFORMATION SCIENCES AND ELECTRONICS */
/*         UNIVERSITY OF TSUKUBA */
/*         1-1-1 TENNODAI, TSUKUBA, IBARAKI 305-8573, JAPAN */
/*         E-MAIL: daisuke@is.tsukuba.ac.jp */


/*     HEADER FILE FOR PARAMETERS */

/*     FORTRAN77 SOURCE PROGRAM */

/*     WRITTEN BY DAISUKE TAKAHASHI */

/* The maximum supported 2-D transform length is 32768. */
/* The maximum supported 3-D transform length is 1024. */
/* The parameter NBLK is a blocking parameter. */
/*      PARAMETER (NBLK=8)  (for PentiumIII and Athlon) */
/*      PARAMETER (NBLK=16) (for Pentium4, Athlon XP, Opteron, Itanium */
/*                           and Itanium2) */
/* The parameter NP is a padding parameter to avoid cache conflicts in */
/* the FFT routines. */
/*      PARAMETER (NP=2) (for PentiumIII) */
/*      PARAMETER (NP=4) (for Athlon, Athlon XP, Opteron and Itanium) */
/*      PARAMETER (NP=8) (for Pentium4 and Itanium2) */
/* Size of L2 cache */
    /* Parameter adjustments */
    --www;
    --ww;
    --b;
    --a;

    /* Function Body */

    p->timings[0] = MPI_Wtime();

    nn = *n / *npu;
    factor_(npu, lnpu);
    factor8_(&nn, lnn);
    for (i__ = 1; i__ <= 3; ++i__) {
	ip[i__ - 1] = lnn[i__ - 1] + lnpu[i__ - 1];
/* Computing MAX */
	i__1 = lnpu[i__ - 1], i__2 = (ip[i__ - 1] + 1) / 3;
	lnz[i__ - 1] = max(i__1,i__2);
/* Computing MAX */
	i__1 = lnpu[i__ - 1], i__2 = (ip[i__ - 1] - lnz[i__ - 1] + 1) / 2;
	lnx[i__ - 1] = max(i__1,i__2);
	lny[i__ - 1] = ip[i__ - 1] - lnx[i__ - 1] - lnz[i__ - 1];
/* L10: */
    }
    nx = pow_ii(&c__2, lnx) * pow_ii(&c__3, &lnx[1]) * pow_ii(&c__5, &lnx[2]);
    ny = pow_ii(&c__2, lny) * pow_ii(&c__3, &lny[1]) * pow_ii(&c__5, &lny[2]);
    nz = pow_ii(&c__2, lnz) * pow_ii(&c__3, &lnz[1]) * pow_ii(&c__5, &lnz[2]);

    if (*n != *n0) {
	settbl_(wx, &nx);
	settbl_(wy, &ny);
	settbl_(wz, &nz);
	settbl3_(&ww[1], &ny, &nz);
	settblp_(&www[1], &nx, &ny, &nz, me, npu);
	*n0 = *n;
    }
    if (*iopt == 0) {
	return 0;
    }

    if (*iopt % 2 == 0) {
	i__1 = nn;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = i__;
	    d_cnjg(&z__1, &a[i__]);
	    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
/* L20: */
	}
    }

    p->timings[1] = MPI_Wtime();
    if (*iopt == 1 || *iopt == 2 || *iopt == 5 || *iopt == 6) {
	pzb2c_(&a[1], &a[1], &b[1], &nn, npu, p );
    }
    p->timings[2] = MPI_Wtime();

/* $OMP PARALLEL PRIVATE(C,D) */
    pzfft1d0_(&a[1], &a[1], &a[1], &b[1], &b[1], c__, c__, d__, wx, wy, wz, &
	    ww[1], &www[1], &nx, &ny, &nz, npu, p );
/* $OMP END PARALLEL */

    p->timings[5] = MPI_Wtime();
    if (*iopt == 1 || *iopt == 2 || *iopt == 7 || *iopt == 8) {
	pzc2b_(&a[1], &a[1], &b[1], &nn, npu, p );
    }
    p->timings[6] = MPI_Wtime();

    if (*iopt % 2 == 0) {
	dn = 1. / (doublereal) (*n);
	i__1 = nn;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = i__;
	    d_cnjg(&z__2, &a[i__]);
	    z__1.r = dn * z__2.r, z__1.i = dn * z__2.i;
	    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
/* L30: */
	}
    }
    p->timings[7] = MPI_Wtime();
    return 0;
} /* pzfft1d_ */

int pzfft1d0_(doublecomplex *a, doublecomplex *ax, 
	doublecomplex *az, doublecomplex *b, doublecomplex *bx, doublecomplex 
	*cy, doublecomplex *cz, doublecomplex *d__, doublecomplex *wx, 
	doublecomplex *wy, doublecomplex *wz, doublecomplex *ww, 
	doublecomplex *www, integer *nx, integer *ny, integer *nz, integer *
	npu, hpcc_fftw_mpi_plan p)
{
    /* System generated locals */
    integer a_dim1, a_dim2, a_offset, ax_dim1, ax_dim2, ax_dim3, ax_offset, 
	    az_dim1, az_dim2, az_offset, b_dim1, b_dim2, b_offset, bx_dim1, 
	    bx_dim2, bx_dim3, bx_offset, cy_dim1, cy_offset, cz_dim1, 
	    cz_offset, ww_dim1, ww_offset, www_dim1, www_dim2, www_offset, 
	    i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8;
    doublecomplex z__1;

    /* Local variables */
    integer i__, j, k, l, ii, jj, kk, lnx[3], lny[3], lnz[3], nnx, nnz;
    extern /* Subroutine */ int fft235a_(doublecomplex *, doublecomplex *, 
	    doublecomplex *, integer *, integer *), factor_(integer *, 
	    integer *), pztrans_(doublecomplex *, doublecomplex *, s64Int_t *, 
	    integer *, hpcc_fftw_mpi_plan);
    s64Int_t nn;


/*     FFTE: A FAST FOURIER TRANSFORM PACKAGE */

/*     (C) COPYRIGHT SOFTWARE, 2000-2004, ALL RIGHTS RESERVED */
/*                BY */
/*         DAISUKE TAKAHASHI */
/*         INSTITURE OF INFORMATION SCIENCES AND ELECTRONICS */
/*         UNIVERSITY OF TSUKUBA */
/*         1-1-1 TENNODAI, TSUKUBA, IBARAKI 305-8573, JAPAN */
/*         E-MAIL: daisuke@is.tsukuba.ac.jp */


/*     HEADER FILE FOR PARAMETERS */

/*     FORTRAN77 SOURCE PROGRAM */

/*     WRITTEN BY DAISUKE TAKAHASHI */

/* The maximum supported 2-D transform length is 32768. */
/* The maximum supported 3-D transform length is 1024. */
/* The parameter NBLK is a blocking parameter. */
/*      PARAMETER (NBLK=8)  (for PentiumIII and Athlon) */
/*      PARAMETER (NBLK=16) (for Pentium4, Athlon XP, Opteron, Itanium */
/*                           and Itanium2) */
/* The parameter NP is a padding parameter to avoid cache conflicts in */
/* the FFT routines. */
/*      PARAMETER (NP=2) (for PentiumIII) */
/*      PARAMETER (NP=4) (for Athlon, Athlon XP, Opteron and Itanium) */
/*      PARAMETER (NP=8) (for Pentium4 and Itanium2) */
/* Size of L2 cache */

    /* Parameter adjustments */
    --d__;
    --wx;
    --wy;
    --wz;
    www_dim1 = *ny;
    www_dim2 = *nx;
    www_offset = 1 + www_dim1 * (1 + www_dim2);
    www -= www_offset;
    cy_dim1 = *ny + 8;
    cy_offset = 1 + cy_dim1;
    cy -= cy_offset;
    b_dim1 = *nx;
    b_dim2 = *ny;
    b_offset = 1 + b_dim1 * (1 + b_dim2);
    b -= b_offset;
    ww_dim1 = *nz;
    ww_offset = 1 + ww_dim1;
    ww -= ww_offset;
    cz_dim1 = *nz + 8;
    cz_offset = 1 + cz_dim1;
    cz -= cz_offset;
    bx_dim1 = *nx / *npu;
    bx_dim2 = *ny;
    bx_dim3 = *nz / *npu;
    bx_offset = 1 + bx_dim1 * (1 + bx_dim2 * (1 + bx_dim3));
    bx -= bx_offset;
    az_dim1 = *nz / *npu;
    az_dim2 = *ny;
    az_offset = 1 + az_dim1 * (1 + az_dim2);
    az -= az_offset;
    ax_dim1 = *nx / *npu;
    ax_dim2 = *ny;
    ax_dim3 = *nz / *npu;
    ax_offset = 1 + ax_dim1 * (1 + ax_dim2 * (1 + ax_dim3));
    ax -= ax_offset;
    a_dim1 = *nx / *npu;
    a_dim2 = *ny;
    a_offset = 1 + a_dim1 * (1 + a_dim2);
    a -= a_offset;

    /* Function Body */
    factor_(nx, lnx);
    factor_(ny, lny);
    factor_(nz, lnz);

    nnx = *nx / *npu;
    nnz = *nz / *npu;
    nn = *nx; nn *= *ny; nn *= *nz; nn /= *npu;

/* $OMP DO */
    i__1 = *ny;
    for (j = 1; j <= i__1; ++j) {
	i__2 = nnx;
	for (ii = 1; ii <= i__2; ii += 16) {
	    i__3 = *nz;
	    for (kk = 1; kk <= i__3; kk += 16) {
/* Computing MIN */
		i__5 = ii + 15;
		i__4 = min(i__5,nnx);
		for (i__ = ii; i__ <= i__4; ++i__) {
/* Computing MIN */
		    i__6 = kk + 15;
		    i__5 = min(i__6,*nz);
		    for (k = kk; k <= i__5; ++k) {
			i__6 = k + (i__ - ii + 1) * cz_dim1;
			i__7 = i__ + (j + k * a_dim2) * a_dim1;
			cz[i__6].r = a[i__7].r, cz[i__6].i = a[i__7].i;
/* L10: */
		    }
/* L20: */
		}
/* L30: */
	    }
/* Computing MIN */
	    i__4 = ii + 15;
	    i__3 = min(i__4,nnx);
	    for (i__ = ii; i__ <= i__3; ++i__) {
		fft235a_(&cz[(i__ - ii + 1) * cz_dim1 + 1], &d__[1], &wz[1], 
			nz, lnz);
/* L40: */
	    }
	    i__3 = *npu;
	    for (l = 1; l <= i__3; ++l) {
		i__4 = nnz;
		for (k = 1; k <= i__4; ++k) {
/* Computing MIN */
		    i__6 = ii + 15;
		    i__5 = min(i__6,nnx);
		    for (i__ = ii; i__ <= i__5; ++i__) {
			i__6 = i__ + (j + (k + l * bx_dim3) * bx_dim2) * 
				bx_dim1;
			i__7 = l + (k - 1) * *npu + (i__ - ii + 1) * cz_dim1;
			i__8 = l + (k - 1) * *npu + j * ww_dim1;
			z__1.r = cz[i__7].r * ww[i__8].r - cz[i__7].i * ww[
				i__8].i, z__1.i = cz[i__7].r * ww[i__8].i + 
				cz[i__7].i * ww[i__8].r;
			bx[i__6].r = z__1.r, bx[i__6].i = z__1.i;
/* L50: */
		    }
/* L60: */
		}
/* L70: */
	    }
/* L80: */
	}
/* L90: */
    }
    p->timings[3] = MPI_Wtime();
/* $OMP SINGLE */
    pztrans_(&bx[bx_offset], &ax[ax_offset], &nn, npu, p );
/* $OMP END SINGLE */
    p->timings[4] = MPI_Wtime();
/* $OMP DO */
    i__1 = nnz;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *npu;
	for (l = 1; l <= i__2; ++l) {
	    i__3 = nnx;
	    for (ii = 1; ii <= i__3; ii += 16) {
		i__4 = *ny;
		for (jj = 1; jj <= i__4; jj += 16) {
/* Computing MIN */
		    i__6 = ii + 15;
		    i__5 = min(i__6,nnx);
		    for (i__ = ii; i__ <= i__5; ++i__) {
/* Computing MIN */
			i__7 = jj + 15;
			i__6 = min(i__7,*ny);
			for (j = jj; j <= i__6; ++j) {
			    i__7 = j + (i__ - ii + 1) * cy_dim1;
			    i__8 = i__ + (j + (k + l * ax_dim3) * ax_dim2) * 
				    ax_dim1;
			    cy[i__7].r = ax[i__8].r, cy[i__7].i = ax[i__8].i;
/* L100: */
			}
/* L110: */
		    }
/* L120: */
		}
/* Computing MIN */
		i__5 = ii + 15;
		i__4 = min(i__5,nnx);
		for (i__ = ii; i__ <= i__4; ++i__) {
		    fft235a_(&cy[(i__ - ii + 1) * cy_dim1 + 1], &d__[1], &wy[
			    1], ny, lny);
/* L130: */
		}
		i__4 = *ny;
		for (j = 1; j <= i__4; ++j) {
/* Computing MIN */
		    i__6 = ii + 15;
		    i__5 = min(i__6,nnx);
		    for (i__ = ii; i__ <= i__5; ++i__) {
			i__6 = l + (i__ - 1) * *npu + (j + k * b_dim2) * 
				b_dim1;
			i__7 = j + (i__ - ii + 1) * cy_dim1;
			i__8 = j + (l + (i__ - 1) * *npu + k * www_dim2) * 
				www_dim1;
			z__1.r = cy[i__7].r * www[i__8].r - cy[i__7].i * www[
				i__8].i, z__1.i = cy[i__7].r * www[i__8].i + 
				cy[i__7].i * www[i__8].r;
			b[i__6].r = z__1.r, b[i__6].i = z__1.i;
/* L140: */
		    }
/* L150: */
		}
/* L160: */
	    }
/* L170: */
	}
	i__2 = *ny;
	for (j = 1; j <= i__2; ++j) {
	    fft235a_(&b[(j + k * b_dim2) * b_dim1 + 1], &d__[1], &wx[1], nx, 
		    lnx);
/* L180: */
	}
/* L190: */
    }
/* $OMP DO */
    i__1 = *nx;
    for (ii = 1; ii <= i__1; ii += 16) {
	i__2 = *ny;
	for (jj = 1; jj <= i__2; jj += 16) {
	    i__3 = nnz;
	    for (kk = 1; kk <= i__3; kk += 16) {
/* Computing MIN */
		i__5 = ii + 15;
		i__4 = min(i__5,*nx);
		for (i__ = ii; i__ <= i__4; ++i__) {
/* Computing MIN */
		    i__6 = jj + 15;
		    i__5 = min(i__6,*ny);
		    for (j = jj; j <= i__5; ++j) {
/* Computing MIN */
			i__7 = kk + 15;
			i__6 = min(i__7,nnz);
			for (k = kk; k <= i__6; ++k) {
			    i__7 = k + (j + i__ * az_dim2) * az_dim1;
			    i__8 = i__ + (j + k * b_dim2) * b_dim1;
			    az[i__7].r = b[i__8].r, az[i__7].i = b[i__8].i;
/* L200: */
			}
/* L210: */
		    }
/* L220: */
		}
/* L230: */
	    }
/* L240: */
	}
/* L250: */
    }
    return 0;
} /* pzfft1d0_ */

/* Subroutine */ int settbl3_(void *w_, integer *ny, integer *nz)
{
    /* System generated locals */
    integer w_dim2, w_offset, i__1, i__2;

    /* Local variables */
    integer j, k;
    doublereal px, pi2;
    doublereal *w = (doublereal *)w_;


    /* Parameter adjustments */
    w_dim2 = *nz;
    w_offset = 1 + 2 * (1 + w_dim2);
    w -= w_offset;

    /* Function Body */
    pi2 = atan(1.) * 8.;
    px = -pi2 / ((doublereal) *ny * (doublereal)*nz);
/* $OMP PARALLEL DO */
    i__1 = *ny;
    for (j = 1; j <= i__1; ++j) {
/* DIR$ VECTOR ALIGNED */
	i__2 = *nz;
	for (k = 1; k <= i__2; ++k) {
	    w[(k + j * w_dim2 << 1) + 1] = cos(px * (doublereal) (k - 1) * (
		    doublereal) (j - 1));
	    w[(k + j * w_dim2 << 1) + 2] = sin(px * (doublereal) (k - 1) * (
		    doublereal) (j - 1));
/* L10: */
	}
/* L20: */
    }
/* $OMP END PARALLEL DO */
    return 0;
} /* settbl3_ */

/* Subroutine */ int settblp_(void *w_, integer *nx, integer *ny, 
	integer *nz, integer *me, integer *npu)
{
    /* System generated locals */
    integer w_dim2, w_dim3, w_offset, i__1, i__2, i__3;

    /* Local variables */
    integer i__, j, k;
    doublereal px, pi2;
    doublereal *w = (doublereal *)w_;


    /* Parameter adjustments */
    w_dim2 = *ny;
    w_dim3 = *nx;
    w_offset = 1 + 2 * (1 + w_dim2 * (1 + w_dim3));
    w -= w_offset;

    /* Function Body */
    pi2 = atan(1.) * 8.;
    px = -pi2 / ((doublereal) *nx * (doublereal) *ny * (doublereal) *nz);
/* $OMP PARALLEL DO */
    i__1 = *nz / *npu;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *nx;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* DIR$ VECTOR ALIGNED */
	    i__3 = *ny;
	    for (j = 1; j <= i__3; ++j) {
		w[(j + (i__ + k * w_dim3) * w_dim2 << 1) + 1] = cos(px * (
			doublereal) (i__ - 1) * ((doublereal) (*me) + (
			doublereal) (k - 1) * (doublereal) (*npu) + (
			doublereal) (j - 1) * (doublereal) (*nz)));
		w[(j + (i__ + k * w_dim3) * w_dim2 << 1) + 2] = sin(px * (
			doublereal) (i__ - 1) * ((doublereal) (*me) + (
			doublereal) (k - 1) * (doublereal) (*npu) + (
			doublereal) (j - 1) * (doublereal) (*nz)));
/* L10: */
	    }
/* L20: */
	}
/* L30: */
    }
/* $OMP END PARALLEL DO */
    return 0;
} /* settblp_ */

