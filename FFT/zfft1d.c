/* zfft1d.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "hpccfft.h"
#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))

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

/*     This has been slightly modified by DHB. */
/*     David H Bailey    13 May 2004 */

/*     FFTE: A FAST FOURIER TRANSFORM PACKAGE */

/*     (C) COPYRIGHT SOFTWARE, 2000-2004, ALL RIGHTS RESERVED */
/*                BY */
/*         DAISUKE TAKAHASHI */
/*         INSTITURE OF INFORMATION SCIENCES AND ELECTRONICS */
/*         UNIVERSITY OF TSUKUBA */
/*         1-1-1 TENNODAI, TSUKUBA, IBARAKI 305-8573, JAPAN */
/*         E-MAIL: daisuke@is.tsukuba.ac.jp */


/*     1-DIMENSIONAL COMPLEX FFT ROUTINE */

/*     FORTRAN77 SOURCE PROGRAM */

/*     CALL ZFFT1D(A,B,N,IOPT) */

/*     A(N) IS COMPLEX INPUT VECTOR (COMPLEX*16) */
/*     B(N) IS COMPLEX OUTPUT VECTOR (COMPLEX*16) */
/*     N IS THE LENGTH OF THE TRANSFORMS (INTEGER*4) */
/*       ----------------------------------- */
/*         N = (2**IP) * (3**IQ) * (5**IR) */
/*       ----------------------------------- */
/*     IOPT = 1 FOR FORWARD TRANSFORM (INTEGER*4) */
/*          = 2 FOR INVERSE TRANSFORM */
/*     DHB: */
/*     IOPT = 0 for initialize w arrays only. */

/*     WRITTEN BY DAISUKE TAKAHASHI */

/* Subroutine */ int zfft1d_(void *a_, void *b_, integer *n, 
	integer *iopt, hpcc_fftw_plan p, integer *n0)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1, z__2;

    /* Local variables */
    doublecomplex *c__, *d__;
    integer i__, n1, n2, m1, m2;
    doublecomplex *w1, *w2;
    doublereal dn;
    integer ip[3], ip1[3], ip2[3];
    doublecomplex *a, *b;
    doublecomplex *ww1, *ww2, *ww3, *ww4;
    extern /* Subroutine */ int fft235b_(doublecomplex *, doublecomplex *, 
	    doublecomplex *, integer *, integer *), factor_(integer *, 
	    integer *), zfft1d0_(doublecomplex *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, doublecomplex *
	    , doublecomplex *, doublecomplex *, doublecomplex *, 
	    doublecomplex *, integer *, integer *, integer *, integer *, 
	    integer *, integer *), settbl_(doublecomplex *, integer *), 
	    settbls_(void *, void *, void *, 
	    void *, integer *, integer *, integer *, integer *);

    a = (doublecomplex *)a_;
    b = (doublecomplex *)b_;
    c__ = (doublecomplex *)p->c;
    d__ = (doublecomplex *)p->d;
    w1 = (doublecomplex *)p->w1;
    w2 = (doublecomplex *)p->w2;
    ww1 = (doublecomplex *)p->ww1;
    ww2 = (doublecomplex *)p->ww2;
    ww3 = (doublecomplex *)p->ww3;
    ww4 = (doublecomplex *)p->ww4;

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
    --b;
    --a;

    /* Function Body */

    factor_(n, ip);

    if (*n <= L2SIZE/16/3 && *n <= NDA2) {
	if (*n != *n0) {
	    settbl_(w1, n);
	    *n0 = *n;
	}
	if (*iopt == 2) {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__2 = i__;
		d_cnjg(&z__1, &a[i__]);
		a[i__2].r = z__1.r, a[i__2].i = z__1.i;
/* L10: */
	    }
	}
	if (*iopt == 0) {
	    return 0;
	}
	fft235b_(&a[1], &b[1], w1, n, ip);
	if (*iopt == 2) {
	    dn = 1. / (doublereal) (*n);
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__2 = i__;
		d_cnjg(&z__2, &b[i__]);
		z__1.r = dn * z__2.r, z__1.i = dn * z__2.i;
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
/* L20: */
	    }
	}
    } else {
	if (*iopt == 2) {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__2 = i__;
		d_cnjg(&z__1, &a[i__]);
		a[i__2].r = z__1.r, a[i__2].i = z__1.i;
/* L30: */
	    }
	}
	for (i__ = 1; i__ <= 3; ++i__) {
	    ip1[i__ - 1] = (ip[i__ - 1] + 1) / 2;
	    ip2[i__ - 1] = ip[i__ - 1] - ip1[i__ - 1];
/* L40: */
	}
	n1 = pow_ii(&c__2, ip1) * pow_ii(&c__3, &ip1[1]) * pow_ii(&c__5, &ip1[
		2]);
	n2 = pow_ii(&c__2, ip2) * pow_ii(&c__3, &ip2[1]) * pow_ii(&c__5, &ip2[
		2]);
	if (pow_ii(&c__2, ip1) < 16 || pow_ii(&c__2, ip2) < 16) {
/* Computing MIN */
	    i__3 = ip1[0] / 2;
	    i__4 = ip1[1] / 2;
	    i__5 = ip1[2] / 2;
	    i__1 = n1, i__2 = pow_ii(&c__2, &i__3) * pow_ii(&c__3, &i__4) * 
		    pow_ii(&c__5, &i__5);
	    m1 = min(i__1,i__2);
/* Computing MIN */
	    i__3 = ip2[0] / 2;
	    i__4 = ip2[1] / 2;
	    i__5 = ip2[2] / 2;
	    i__1 = n2, i__2 = pow_ii(&c__2, &i__3) * pow_ii(&c__3, &i__4) * 
		    pow_ii(&c__5, &i__5);
	    m2 = min(i__1,i__2);
	} else {
/* Computing MIN */
/* Computing MAX */
	    i__5 = ip1[0] / 2;
	    i__3 = 16, i__4 = pow_ii(&c__2, &i__5);
	    i__1 = n1, i__2 = max(i__3,i__4);
	    m1 = min(i__1,i__2);
/* Computing MIN */
/* Computing MAX */
	    i__5 = ip2[0] / 2;
	    i__3 = 16, i__4 = pow_ii(&c__2, &i__5);
	    i__1 = n2, i__2 = max(i__3,i__4);
	    m2 = min(i__1,i__2);
	}

	if (*n != *n0) {
	    settbl_(w1, &n1);
	    settbl_(w2, &n2);
	    settbls_(ww1, ww2, ww3, ww4, &n1, &n2, &m1, &m2);
	    *n0 = *n;
	}
	if (*iopt == 0) {
	    return 0;
	}
/* $OMP PARALLEL PRIVATE(C,D) */
	zfft1d0_(&a[1], &b[1], c__, d__, w1, w2, ww1, ww2, ww3, ww4, &n1, &n2,
		 &m1, &m2, ip1, ip2);
/* $OMP END PARALLEL */
	if (*iopt == 2) {
	    dn = 1. / (doublereal) (*n);
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__2 = i__;
		d_cnjg(&z__2, &b[i__]);
		z__1.r = dn * z__2.r, z__1.i = dn * z__2.i;
		b[i__2].r = z__1.r, b[i__2].i = z__1.i;
/* L50: */
	    }
	}
    }
    return 0;
} /* zfft1d_ */

/* Subroutine */ int zfft1d0_(doublecomplex *a, doublecomplex *b, 
	doublecomplex *c__, doublecomplex *d__, doublecomplex *w1, 
	doublecomplex *w2, doublecomplex *ww1, doublecomplex *ww2, 
	doublecomplex *ww3, doublecomplex *ww4, integer *n1, integer *n2, 
	integer *m1, integer *m2, integer *ip1, integer *ip2)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, ww1_dim1, 
	    ww1_offset, ww2_dim1, ww2_offset, ww3_dim1, ww3_offset, ww4_dim1, 
	    ww4_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8, i__9, 
	    i__10;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Local variables */
    integer i__, j, ii, ij, jj, ik, ir, is, ij0;
    doublecomplex temp;
    extern /* Subroutine */ int fft235a_(doublecomplex *, doublecomplex *, 
	    doublecomplex *, integer *, integer *);


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

/* $OMP DO PRIVATE(IJ,IJ0,IR,J,TEMP) */
    /* Parameter adjustments */
    --d__;
    --w1;
    --w2;
    a_dim1 = *n1;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    c_dim1 = *n2 + 8;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    b_dim1 = *n2;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    ww4_dim1 = *n2 / *m2 + 8;
    ww4_offset = 1 + ww4_dim1;
    ww4 -= ww4_offset;
    ww3_dim1 = *m2 + 8;
    ww3_offset = 1 + ww3_dim1;
    ww3 -= ww3_offset;
    ww2_dim1 = *n2 / *m2 + 8;
    ww2_offset = 1 + ww2_dim1;
    ww2 -= ww2_offset;
    ww1_dim1 = *m2 + 8;
    ww1_offset = 1 + ww1_dim1;
    ww1 -= ww1_offset;
    --ip1;
    --ip2;

    /* Function Body */
    i__1 = *n1;
    for (ii = 1; ii <= i__1; ii += 16) {
	i__2 = *n2;
	for (jj = 1; jj <= i__2; jj += 16) {
/* Computing MIN */
	    i__4 = ii + 15;
	    i__3 = min(i__4,*n1);
	    for (i__ = ii; i__ <= i__3; ++i__) {
/* Computing MIN */
		i__5 = jj + 15;
		i__4 = min(i__5,*n2);
		for (j = jj; j <= i__4; ++j) {
		    i__5 = j + (i__ - ii + 1) * c_dim1;
		    i__6 = i__ + j * a_dim1;
		    c__[i__5].r = a[i__6].r, c__[i__5].i = a[i__6].i;
/* L10: */
		}
/* L20: */
	    }
/* L30: */
	}
/* Computing MIN */
	i__3 = ii + 15;
	i__2 = min(i__3,*n1);
	for (i__ = ii; i__ <= i__2; ++i__) {
	    fft235a_(&c__[(i__ - ii + 1) * c_dim1 + 1], &d__[1], &w2[1], n2, &
		    ip2[1]);
/* L40: */
	}
	if (pow_ii(&c__2, &ip1[1]) < 16 || pow_ii(&c__2, &ip2[1]) < 16) {
	    i__2 = *n2 / *m2;
	    for (is = 1; is <= i__2; ++is) {
		i__3 = *m2;
		for (ik = 1; ik <= i__3; ++ik) {
		    j = ik + (is - 1) * *m2;
/* Computing MIN */
		    i__5 = ii + 15;
		    i__4 = min(i__5,*n1);
		    for (i__ = ii; i__ <= i__4; ++i__) {
			ir = (i__ - 1) / *m1 + 1;
			ij = (i__ - 1) % *m1 + 1;
			i__5 = i__ + j * a_dim1;
			i__6 = j + (i__ - ii + 1) * c_dim1;
			i__7 = ik + ij * ww1_dim1;
			i__8 = is + ij * ww2_dim1;
			z__4.r = ww1[i__7].r * ww2[i__8].r - ww1[i__7].i * 
				ww2[i__8].i, z__4.i = ww1[i__7].r * ww2[i__8]
				.i + ww1[i__7].i * ww2[i__8].r;
			i__9 = ik + ir * ww3_dim1;
			z__3.r = z__4.r * ww3[i__9].r - z__4.i * ww3[i__9].i, 
				z__3.i = z__4.r * ww3[i__9].i + z__4.i * ww3[
				i__9].r;
			i__10 = is + ir * ww4_dim1;
			z__2.r = z__3.r * ww4[i__10].r - z__3.i * ww4[i__10]
				.i, z__2.i = z__3.r * ww4[i__10].i + z__3.i * 
				ww4[i__10].r;
			z__1.r = c__[i__6].r * z__2.r - c__[i__6].i * z__2.i, 
				z__1.i = c__[i__6].r * z__2.i + c__[i__6].i * 
				z__2.r;
			a[i__5].r = z__1.r, a[i__5].i = z__1.i;
/* L50: */
		    }
/* L60: */
		}
/* L70: */
	    }
	} else {
	    ir = (ii - 1) / *m1 + 1;
	    ij0 = (ii - 1) % *m1 + 1;
	    i__2 = *n2 / *m2;
	    for (is = 1; is <= i__2; ++is) {
		i__3 = *m2;
		for (ik = 1; ik <= i__3; ++ik) {
		    i__4 = ik + ir * ww3_dim1;
		    i__5 = is + ir * ww4_dim1;
		    z__1.r = ww3[i__4].r * ww4[i__5].r - ww3[i__4].i * ww4[
			    i__5].i, z__1.i = ww3[i__4].r * ww4[i__5].i + ww3[
			    i__4].i * ww4[i__5].r;
		    temp.r = z__1.r, temp.i = z__1.i;
		    j = ik + (is - 1) * *m2;
		    ij = ij0;
/* Computing MIN */
		    i__5 = ii + 15;
		    i__4 = min(i__5,*n1);
		    for (i__ = ii; i__ <= i__4; ++i__) {
			i__5 = i__ + j * a_dim1;
			i__6 = j + (i__ - ii + 1) * c_dim1;
			i__7 = ik + ij * ww1_dim1;
			i__8 = is + ij * ww2_dim1;
			z__3.r = ww1[i__7].r * ww2[i__8].r - ww1[i__7].i * 
				ww2[i__8].i, z__3.i = ww1[i__7].r * ww2[i__8]
				.i + ww1[i__7].i * ww2[i__8].r;
			z__2.r = z__3.r * temp.r - z__3.i * temp.i, z__2.i = 
				z__3.r * temp.i + z__3.i * temp.r;
			z__1.r = c__[i__6].r * z__2.r - c__[i__6].i * z__2.i, 
				z__1.i = c__[i__6].r * z__2.i + c__[i__6].i * 
				z__2.r;
			a[i__5].r = z__1.r, a[i__5].i = z__1.i;
			++ij;
/* L80: */
		    }
/* L90: */
		}
/* L100: */
	    }
	}
/* L110: */
    }
/* $OMP DO */
    i__1 = *n2;
    for (jj = 1; jj <= i__1; jj += 16) {
/* Computing MIN */
	i__3 = jj + 15;
	i__2 = min(i__3,*n2);
	for (j = jj; j <= i__2; ++j) {
	    fft235a_(&a[j * a_dim1 + 1], &c__[c_offset], &w1[1], n1, &ip1[1]);
/* L120: */
	}
	i__2 = *n1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MIN */
	    i__4 = jj + 15;
	    i__3 = min(i__4,*n2);
	    for (j = jj; j <= i__3; ++j) {
		i__4 = j + i__ * b_dim1;
		i__5 = i__ + j * a_dim1;
		b[i__4].r = a[i__5].r, b[i__4].i = a[i__5].i;
/* L130: */
	    }
/* L140: */
	}
/* L150: */
    }
    return 0;
} /* zfft1d0_ */

/* Subroutine */ int settbls_(void *w1_, void *w2_, void *w3_, 
	void *w4_, integer *n1, integer *n2, integer *m1, integer *m2)
{
    /* System generated locals */
    integer w1_dim2, w1_offset, w2_dim2, w2_offset, w3_dim2, w3_offset, 
	    w4_dim2, w4_offset, i__1, i__2;

    /* Local variables */
    integer j, k, ir, is;
    doublereal px, pi2;
    doublereal *w1 = (doublereal *)w1_, *w2 = (doublereal *)w2_, *w3 = (doublereal *)w3_, *w4 = (doublereal *)w4_;


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
    w4_dim2 = *n2 / *m2 + 8;
    w4_offset = 1 + 2 * (1 + w4_dim2);
    w4 -= w4_offset;
    w3_dim2 = *m2 + 8;
    w3_offset = 1 + 2 * (1 + w3_dim2);
    w3 -= w3_offset;
    w2_dim2 = *n2 / *m2 + 8;
    w2_offset = 1 + 2 * (1 + w2_dim2);
    w2 -= w2_offset;
    w1_dim2 = *m2 + 8;
    w1_offset = 1 + 2 * (1 + w1_dim2);
    w1 -= w1_offset;

    /* Function Body */
    pi2 = atan(1.) * 8.;
    px = -pi2 / (doublereal) (*n1 * *n2);
/* $OMP PARALLEL */
/* $OMP DO */
    i__1 = *m1;
    for (j = 1; j <= i__1; ++j) {
/* DIR$ VECTOR ALIGNED */
	i__2 = *m2;
	for (k = 1; k <= i__2; ++k) {
	    w1[(k + j * w1_dim2 << 1) + 1] = cos(px * (doublereal) (j - 1) * (
		    doublereal) (k - 1));
	    w1[(k + j * w1_dim2 << 1) + 2] = sin(px * (doublereal) (j - 1) * (
		    doublereal) (k - 1));
/* L10: */
	}
/* DIR$ VECTOR ALIGNED */
	i__2 = *n2 / *m2;
	for (is = 1; is <= i__2; ++is) {
	    w2[(is + j * w2_dim2 << 1) + 1] = cos(px * (doublereal) (j - 1) * 
		    (doublereal) (is - 1) * (doublereal) (*m2));
	    w2[(is + j * w2_dim2 << 1) + 2] = sin(px * (doublereal) (j - 1) * 
		    (doublereal) (is - 1) * (doublereal) (*m2));
/* L20: */
	}
/* L30: */
    }
/* $OMP DO */
    i__1 = *n1 / *m1;
    for (ir = 1; ir <= i__1; ++ir) {
/* DIR$ VECTOR ALIGNED */
	i__2 = *m2;
	for (k = 1; k <= i__2; ++k) {
	    w3[(k + ir * w3_dim2 << 1) + 1] = cos(px * (doublereal) (ir - 1) *
		     (doublereal) (*m1) * (doublereal) (k - 1));
	    w3[(k + ir * w3_dim2 << 1) + 2] = sin(px * (doublereal) (ir - 1) *
		     (doublereal) (*m1) * (doublereal) (k - 1));
/* L40: */
	}
/* DIR$ VECTOR ALIGNED */
	i__2 = *n2 / *m2;
	for (is = 1; is <= i__2; ++is) {
	    w4[(is + ir * w4_dim2 << 1) + 1] = cos(px * (doublereal) (ir - 1) 
		    * (doublereal) (*m1) * (doublereal) (is - 1) * (
		    doublereal) (*m2));
	    w4[(is + ir * w4_dim2 << 1) + 2] = sin(px * (doublereal) (ir - 1) 
		    * (doublereal) (*m1) * (doublereal) (is - 1) * (
		    doublereal) (*m2));
/* L50: */
	}
/* L60: */
    }
/* $OMP END PARALLEL */
    return 0;
} /* settbls_ */

