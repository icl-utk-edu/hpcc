/* kernel2.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "hpccfft.h"


/*     FFTE: A FAST FOURIER TRANSFORM PACKAGE */

/*     (C) COPYRIGHT SOFTWARE, 2000-2004, ALL RIGHTS RESERVED */
/*                BY */
/*         DAISUKE TAKAHASHI */
/*         INSTITURE OF INFORMATION SCIENCES AND ELECTRONICS */
/*         UNIVERSITY OF TSUKUBA */
/*         1-1-1 TENNODAI, TSUKUBA, IBARAKI 305-8573, JAPAN */
/*         E-MAIL: daisuke@is.tsukuba.ac.jp */


/*     RADIX-2, 3, 4, 5 AND 8 FFT KERNEL ROUTINE */

/*     FORTRAN77 SOURCE PROGRAM */

/*     WRITTEN BY DAISUKE TAKAHASHI */

/* Subroutine */ int fft2_(doublereal *a, doublereal *b, integer *m)
{
    /* System generated locals */
    integer a_dim2, a_offset, b_dim2, b_offset, i__1;

    /* Local variables */
    integer i__;
    doublereal x0, y0, x1, y1;


/* DIR$ VECTOR ALIGNED */
    /* Parameter adjustments */
    b_dim2 = *m;
    b_offset = 1 + 2 * (1 + b_dim2);
    b -= b_offset;
    a_dim2 = *m;
    a_offset = 1 + 2 * (1 + a_dim2);
    a -= a_offset;

    /* Function Body */
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x0 = a[(i__ + a_dim2 << 1) + 1];
	y0 = a[(i__ + a_dim2 << 1) + 2];
	x1 = a[(i__ + (a_dim2 << 1) << 1) + 1];
	y1 = a[(i__ + (a_dim2 << 1) << 1) + 2];
	b[(i__ + b_dim2 << 1) + 1] = x0 + x1;
	b[(i__ + b_dim2 << 1) + 2] = y0 + y1;
	b[(i__ + (b_dim2 << 1) << 1) + 1] = x0 - x1;
	b[(i__ + (b_dim2 << 1) << 1) + 2] = y0 - y1;
/* L10: */
    }
    return 0;
} /* fft2_ */

/* Subroutine */ int fft3a_(doublereal *a, doublereal *b, doublereal *w, 
	integer *l)
{
    /* Initialized data */

    static doublereal c31 = .86602540378443865;
    static doublereal c32 = .5;

    /* System generated locals */
    integer a_dim2, a_offset, i__1;

    /* Local variables */
    integer j;
    doublereal x0, y0, x1, y1, x2, y2, wi1, wi2, wr1, wr2;

    /* Parameter adjustments */
    b -= 9;
    w -= 3;
    a_dim2 = *l;
    a_offset = 1 + 2 * (1 + a_dim2);
    a -= a_offset;

    /* Function Body */

/* DIR$ VECTOR ALIGNED */
    i__1 = *l;
    for (j = 1; j <= i__1; ++j) {
	wr1 = w[(j << 1) + 1];
	wi1 = w[(j << 1) + 2];
	wr2 = wr1 * wr1 - wi1 * wi1;
	wi2 = wr1 * wi1 + wr1 * wi1;
	x0 = a[(j + (a_dim2 << 1) << 1) + 1] + a[(j + a_dim2 * 3 << 1) + 1];
	y0 = a[(j + (a_dim2 << 1) << 1) + 2] + a[(j + a_dim2 * 3 << 1) + 2];
	x1 = a[(j + a_dim2 << 1) + 1] - c32 * x0;
	y1 = a[(j + a_dim2 << 1) + 2] - c32 * y0;
	x2 = c31 * (a[(j + (a_dim2 << 1) << 1) + 2] - a[(j + a_dim2 * 3 << 1) 
		+ 2]);
	y2 = c31 * (a[(j + a_dim2 * 3 << 1) + 1] - a[(j + (a_dim2 << 1) << 1) 
		+ 1]);
	b[(j * 3 + 1 << 1) + 1] = a[(j + a_dim2 << 1) + 1] + x0;
	b[(j * 3 + 1 << 1) + 2] = a[(j + a_dim2 << 1) + 2] + y0;
	b[(j * 3 + 2 << 1) + 1] = wr1 * (x1 + x2) - wi1 * (y1 + y2);
	b[(j * 3 + 2 << 1) + 2] = wr1 * (y1 + y2) + wi1 * (x1 + x2);
	b[(j * 3 + 3 << 1) + 1] = wr2 * (x1 - x2) - wi2 * (y1 - y2);
	b[(j * 3 + 3 << 1) + 2] = wr2 * (y1 - y2) + wi2 * (x1 - x2);
/* L10: */
    }
    return 0;
} /* fft3a_ */

/* Subroutine */ int fft3b_(doublereal *a, doublereal *b, doublereal *w, 
	integer *m, integer *l)
{
    /* Initialized data */

    static doublereal c31 = .86602540378443865;
    static doublereal c32 = .5;

    /* System generated locals */
    integer a_dim2, a_dim3, a_offset, b_dim2, b_offset, i__1, i__2;

    /* Local variables */
    integer i__, j;
    doublereal x0, y0, x1, y1, x2, y2, wi1, wi2, wr1, wr2;

    /* Parameter adjustments */
    w -= 3;
    b_dim2 = *m;
    b_offset = 1 + 2 * (1 + (b_dim2 << 2));
    b -= b_offset;
    a_dim2 = *m;
    a_dim3 = *l;
    a_offset = 1 + 2 * (1 + a_dim2 * (1 + a_dim3));
    a -= a_offset;

    /* Function Body */

/* DIR$ VECTOR ALIGNED */
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x0 = a[(i__ + ((a_dim3 << 1) + 1) * a_dim2 << 1) + 1] + a[(i__ + (
		a_dim3 * 3 + 1) * a_dim2 << 1) + 1];
	y0 = a[(i__ + ((a_dim3 << 1) + 1) * a_dim2 << 1) + 2] + a[(i__ + (
		a_dim3 * 3 + 1) * a_dim2 << 1) + 2];
	x1 = a[(i__ + (a_dim3 + 1) * a_dim2 << 1) + 1] - c32 * x0;
	y1 = a[(i__ + (a_dim3 + 1) * a_dim2 << 1) + 2] - c32 * y0;
	x2 = c31 * (a[(i__ + ((a_dim3 << 1) + 1) * a_dim2 << 1) + 2] - a[(i__ 
		+ (a_dim3 * 3 + 1) * a_dim2 << 1) + 2]);
	y2 = c31 * (a[(i__ + (a_dim3 * 3 + 1) * a_dim2 << 1) + 1] - a[(i__ + (
		(a_dim3 << 1) + 1) * a_dim2 << 1) + 1]);
	b[(i__ + (b_dim2 << 2) << 1) + 1] = a[(i__ + (a_dim3 + 1) * a_dim2 << 
		1) + 1] + x0;
	b[(i__ + (b_dim2 << 2) << 1) + 2] = a[(i__ + (a_dim3 + 1) * a_dim2 << 
		1) + 2] + y0;
	b[(i__ + b_dim2 * 5 << 1) + 1] = x1 + x2;
	b[(i__ + b_dim2 * 5 << 1) + 2] = y1 + y2;
	b[(i__ + b_dim2 * 6 << 1) + 1] = x1 - x2;
	b[(i__ + b_dim2 * 6 << 1) + 2] = y1 - y2;
/* L10: */
    }
    i__1 = *l;
    for (j = 2; j <= i__1; ++j) {
	wr1 = w[(j << 1) + 1];
	wi1 = w[(j << 1) + 2];
	wr2 = wr1 * wr1 - wi1 * wi1;
	wi2 = wr1 * wi1 + wr1 * wi1;
/* DIR$ VECTOR ALIGNED */
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    x0 = a[(i__ + (j + (a_dim3 << 1)) * a_dim2 << 1) + 1] + a[(i__ + (
		    j + a_dim3 * 3) * a_dim2 << 1) + 1];
	    y0 = a[(i__ + (j + (a_dim3 << 1)) * a_dim2 << 1) + 2] + a[(i__ + (
		    j + a_dim3 * 3) * a_dim2 << 1) + 2];
	    x1 = a[(i__ + (j + a_dim3) * a_dim2 << 1) + 1] - c32 * x0;
	    y1 = a[(i__ + (j + a_dim3) * a_dim2 << 1) + 2] - c32 * y0;
	    x2 = c31 * (a[(i__ + (j + (a_dim3 << 1)) * a_dim2 << 1) + 2] - a[(
		    i__ + (j + a_dim3 * 3) * a_dim2 << 1) + 2]);
	    y2 = c31 * (a[(i__ + (j + a_dim3 * 3) * a_dim2 << 1) + 1] - a[(
		    i__ + (j + (a_dim3 << 1)) * a_dim2 << 1) + 1]);
	    b[(i__ + (j * 3 + 1) * b_dim2 << 1) + 1] = a[(i__ + (j + a_dim3) *
		     a_dim2 << 1) + 1] + x0;
	    b[(i__ + (j * 3 + 1) * b_dim2 << 1) + 2] = a[(i__ + (j + a_dim3) *
		     a_dim2 << 1) + 2] + y0;
	    b[(i__ + (j * 3 + 2) * b_dim2 << 1) + 1] = wr1 * (x1 + x2) - wi1 *
		     (y1 + y2);
	    b[(i__ + (j * 3 + 2) * b_dim2 << 1) + 2] = wr1 * (y1 + y2) + wi1 *
		     (x1 + x2);
	    b[(i__ + (j * 3 + 3) * b_dim2 << 1) + 1] = wr2 * (x1 - x2) - wi2 *
		     (y1 - y2);
	    b[(i__ + (j * 3 + 3) * b_dim2 << 1) + 2] = wr2 * (y1 - y2) + wi2 *
		     (x1 - x2);
/* L20: */
	}
/* L30: */
    }
    return 0;
} /* fft3b_ */

/* Subroutine */ int fft4a_(doublereal *a, doublereal *b, doublereal *w, 
	integer *l)
{
    /* System generated locals */
    integer a_dim2, a_offset, i__1;

    /* Local variables */
    integer j;
    doublereal x0, y0, x1, y1, x2, y2, x3, y3, wi1, wi2, wi3, wr1, wr2, wr3;


/* DIR$ VECTOR ALIGNED */
    /* Parameter adjustments */
    b -= 11;
    w -= 3;
    a_dim2 = *l;
    a_offset = 1 + 2 * (1 + a_dim2);
    a -= a_offset;

    /* Function Body */
    i__1 = *l;
    for (j = 1; j <= i__1; ++j) {
	wr1 = w[(j << 1) + 1];
	wi1 = w[(j << 1) + 2];
	wr2 = wr1 * wr1 - wi1 * wi1;
	wi2 = wr1 * wi1 + wr1 * wi1;
	wr3 = wr1 * wr2 - wi1 * wi2;
	wi3 = wr1 * wi2 + wi1 * wr2;
	x0 = a[(j + a_dim2 << 1) + 1] + a[(j + a_dim2 * 3 << 1) + 1];
	y0 = a[(j + a_dim2 << 1) + 2] + a[(j + a_dim2 * 3 << 1) + 2];
	x1 = a[(j + a_dim2 << 1) + 1] - a[(j + a_dim2 * 3 << 1) + 1];
	y1 = a[(j + a_dim2 << 1) + 2] - a[(j + a_dim2 * 3 << 1) + 2];
	x2 = a[(j + (a_dim2 << 1) << 1) + 1] + a[(j + (a_dim2 << 2) << 1) + 1]
		;
	y2 = a[(j + (a_dim2 << 1) << 1) + 2] + a[(j + (a_dim2 << 2) << 1) + 2]
		;
	x3 = a[(j + (a_dim2 << 1) << 1) + 2] - a[(j + (a_dim2 << 2) << 1) + 2]
		;
	y3 = a[(j + (a_dim2 << 2) << 1) + 1] - a[(j + (a_dim2 << 1) << 1) + 1]
		;
	b[((j << 2) + 1 << 1) + 1] = x0 + x2;
	b[((j << 2) + 1 << 1) + 2] = y0 + y2;
	b[((j << 2) + 3 << 1) + 1] = wr2 * (x0 - x2) - wi2 * (y0 - y2);
	b[((j << 2) + 3 << 1) + 2] = wr2 * (y0 - y2) + wi2 * (x0 - x2);
	b[((j << 2) + 2 << 1) + 1] = wr1 * (x1 + x3) - wi1 * (y1 + y3);
	b[((j << 2) + 2 << 1) + 2] = wr1 * (y1 + y3) + wi1 * (x1 + x3);
	b[((j << 2) + 4 << 1) + 1] = wr3 * (x1 - x3) - wi3 * (y1 - y3);
	b[((j << 2) + 4 << 1) + 2] = wr3 * (y1 - y3) + wi3 * (x1 - x3);
/* L10: */
    }
    return 0;
} /* fft4a_ */

/* Subroutine */ int fft4b_(doublereal *a, doublereal *b, doublereal *w, 
	integer *m, integer *l)
{
    /* System generated locals */
    integer a_dim2, a_dim3, a_offset, b_dim2, b_offset, i__1, i__2;

    /* Local variables */
    integer i__, j;
    doublereal x0, y0, x1, y1, x2, y2, x3, y3, wi1, wi2, wi3, wr1, wr2, wr3;


/* DIR$ VECTOR ALIGNED */
    /* Parameter adjustments */
    w -= 3;
    b_dim2 = *m;
    b_offset = 1 + 2 * (1 + b_dim2 * 5);
    b -= b_offset;
    a_dim2 = *m;
    a_dim3 = *l;
    a_offset = 1 + 2 * (1 + a_dim2 * (1 + a_dim3));
    a -= a_offset;

    /* Function Body */
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x0 = a[(i__ + (a_dim3 + 1) * a_dim2 << 1) + 1] + a[(i__ + (a_dim3 * 3 
		+ 1) * a_dim2 << 1) + 1];
	y0 = a[(i__ + (a_dim3 + 1) * a_dim2 << 1) + 2] + a[(i__ + (a_dim3 * 3 
		+ 1) * a_dim2 << 1) + 2];
	x1 = a[(i__ + (a_dim3 + 1) * a_dim2 << 1) + 1] - a[(i__ + (a_dim3 * 3 
		+ 1) * a_dim2 << 1) + 1];
	y1 = a[(i__ + (a_dim3 + 1) * a_dim2 << 1) + 2] - a[(i__ + (a_dim3 * 3 
		+ 1) * a_dim2 << 1) + 2];
	x2 = a[(i__ + ((a_dim3 << 1) + 1) * a_dim2 << 1) + 1] + a[(i__ + ((
		a_dim3 << 2) + 1) * a_dim2 << 1) + 1];
	y2 = a[(i__ + ((a_dim3 << 1) + 1) * a_dim2 << 1) + 2] + a[(i__ + ((
		a_dim3 << 2) + 1) * a_dim2 << 1) + 2];
	x3 = a[(i__ + ((a_dim3 << 1) + 1) * a_dim2 << 1) + 2] - a[(i__ + ((
		a_dim3 << 2) + 1) * a_dim2 << 1) + 2];
	y3 = a[(i__ + ((a_dim3 << 2) + 1) * a_dim2 << 1) + 1] - a[(i__ + ((
		a_dim3 << 1) + 1) * a_dim2 << 1) + 1];
	b[(i__ + b_dim2 * 5 << 1) + 1] = x0 + x2;
	b[(i__ + b_dim2 * 5 << 1) + 2] = y0 + y2;
	b[(i__ + b_dim2 * 7 << 1) + 1] = x0 - x2;
	b[(i__ + b_dim2 * 7 << 1) + 2] = y0 - y2;
	b[(i__ + b_dim2 * 6 << 1) + 1] = x1 + x3;
	b[(i__ + b_dim2 * 6 << 1) + 2] = y1 + y3;
	b[(i__ + (b_dim2 << 3) << 1) + 1] = x1 - x3;
	b[(i__ + (b_dim2 << 3) << 1) + 2] = y1 - y3;
/* L10: */
    }
    i__1 = *l;
    for (j = 2; j <= i__1; ++j) {
	wr1 = w[(j << 1) + 1];
	wi1 = w[(j << 1) + 2];
	wr2 = wr1 * wr1 - wi1 * wi1;
	wi2 = wr1 * wi1 + wr1 * wi1;
	wr3 = wr1 * wr2 - wi1 * wi2;
	wi3 = wr1 * wi2 + wi1 * wr2;
/* DIR$ VECTOR ALIGNED */
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    x0 = a[(i__ + (j + a_dim3) * a_dim2 << 1) + 1] + a[(i__ + (j + 
		    a_dim3 * 3) * a_dim2 << 1) + 1];
	    y0 = a[(i__ + (j + a_dim3) * a_dim2 << 1) + 2] + a[(i__ + (j + 
		    a_dim3 * 3) * a_dim2 << 1) + 2];
	    x1 = a[(i__ + (j + a_dim3) * a_dim2 << 1) + 1] - a[(i__ + (j + 
		    a_dim3 * 3) * a_dim2 << 1) + 1];
	    y1 = a[(i__ + (j + a_dim3) * a_dim2 << 1) + 2] - a[(i__ + (j + 
		    a_dim3 * 3) * a_dim2 << 1) + 2];
	    x2 = a[(i__ + (j + (a_dim3 << 1)) * a_dim2 << 1) + 1] + a[(i__ + (
		    j + (a_dim3 << 2)) * a_dim2 << 1) + 1];
	    y2 = a[(i__ + (j + (a_dim3 << 1)) * a_dim2 << 1) + 2] + a[(i__ + (
		    j + (a_dim3 << 2)) * a_dim2 << 1) + 2];
	    x3 = a[(i__ + (j + (a_dim3 << 1)) * a_dim2 << 1) + 2] - a[(i__ + (
		    j + (a_dim3 << 2)) * a_dim2 << 1) + 2];
	    y3 = a[(i__ + (j + (a_dim3 << 2)) * a_dim2 << 1) + 1] - a[(i__ + (
		    j + (a_dim3 << 1)) * a_dim2 << 1) + 1];
	    b[(i__ + ((j << 2) + 1) * b_dim2 << 1) + 1] = x0 + x2;
	    b[(i__ + ((j << 2) + 1) * b_dim2 << 1) + 2] = y0 + y2;
	    b[(i__ + ((j << 2) + 3) * b_dim2 << 1) + 1] = wr2 * (x0 - x2) - 
		    wi2 * (y0 - y2);
	    b[(i__ + ((j << 2) + 3) * b_dim2 << 1) + 2] = wr2 * (y0 - y2) + 
		    wi2 * (x0 - x2);
	    b[(i__ + ((j << 2) + 2) * b_dim2 << 1) + 1] = wr1 * (x1 + x3) - 
		    wi1 * (y1 + y3);
	    b[(i__ + ((j << 2) + 2) * b_dim2 << 1) + 2] = wr1 * (y1 + y3) + 
		    wi1 * (x1 + x3);
	    b[(i__ + ((j << 2) + 4) * b_dim2 << 1) + 1] = wr3 * (x1 - x3) - 
		    wi3 * (y1 - y3);
	    b[(i__ + ((j << 2) + 4) * b_dim2 << 1) + 2] = wr3 * (y1 - y3) + 
		    wi3 * (x1 - x3);
/* L20: */
	}
/* L30: */
    }
    return 0;
} /* fft4b_ */

/* Subroutine */ int fft5a_(doublereal *a, doublereal *b, doublereal *w, 
	integer *l)
{
    /* Initialized data */

    static doublereal c51 = .95105651629515357;
    static doublereal c52 = .61803398874989485;
    static doublereal c53 = .55901699437494742;
    static doublereal c54 = .25;

    /* System generated locals */
    integer a_dim2, a_offset, i__1;

    /* Local variables */
    integer j;
    doublereal x0, y0, x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, x7, y7,
	     x8, y8, x9, y9, x10, y10, wi1, wi2, wi3, wi4, wr1, wr2, wr3, wr4;

    /* Parameter adjustments */
    b -= 13;
    w -= 3;
    a_dim2 = *l;
    a_offset = 1 + 2 * (1 + a_dim2);
    a -= a_offset;

    /* Function Body */

/* DIR$ VECTOR ALIGNED */
    i__1 = *l;
    for (j = 1; j <= i__1; ++j) {
	wr1 = w[(j << 1) + 1];
	wi1 = w[(j << 1) + 2];
	wr2 = wr1 * wr1 - wi1 * wi1;
	wi2 = wr1 * wi1 + wr1 * wi1;
	wr3 = wr1 * wr2 - wi1 * wi2;
	wi3 = wr1 * wi2 + wi1 * wr2;
	wr4 = wr2 * wr2 - wi2 * wi2;
	wi4 = wr2 * wi2 + wr2 * wi2;
	x0 = a[(j + (a_dim2 << 1) << 1) + 1] + a[(j + a_dim2 * 5 << 1) + 1];
	y0 = a[(j + (a_dim2 << 1) << 1) + 2] + a[(j + a_dim2 * 5 << 1) + 2];
	x1 = a[(j + a_dim2 * 3 << 1) + 1] + a[(j + (a_dim2 << 2) << 1) + 1];
	y1 = a[(j + a_dim2 * 3 << 1) + 2] + a[(j + (a_dim2 << 2) << 1) + 2];
	x2 = c51 * (a[(j + (a_dim2 << 1) << 1) + 1] - a[(j + a_dim2 * 5 << 1) 
		+ 1]);
	y2 = c51 * (a[(j + (a_dim2 << 1) << 1) + 2] - a[(j + a_dim2 * 5 << 1) 
		+ 2]);
	x3 = c51 * (a[(j + a_dim2 * 3 << 1) + 1] - a[(j + (a_dim2 << 2) << 1) 
		+ 1]);
	y3 = c51 * (a[(j + a_dim2 * 3 << 1) + 2] - a[(j + (a_dim2 << 2) << 1) 
		+ 2]);
	x4 = x0 + x1;
	y4 = y0 + y1;
	x5 = c53 * (x0 - x1);
	y5 = c53 * (y0 - y1);
	x6 = a[(j + a_dim2 << 1) + 1] - c54 * x4;
	y6 = a[(j + a_dim2 << 1) + 2] - c54 * y4;
	x7 = x6 + x5;
	y7 = y6 + y5;
	x8 = x6 - x5;
	y8 = y6 - y5;
	x9 = y2 + c52 * y3;
	y9 = -x2 - c52 * x3;
	x10 = c52 * y2 - y3;
	y10 = x3 - c52 * x2;
	b[(j * 5 + 1 << 1) + 1] = a[(j + a_dim2 << 1) + 1] + x4;
	b[(j * 5 + 1 << 1) + 2] = a[(j + a_dim2 << 1) + 2] + y4;
	b[(j * 5 + 2 << 1) + 1] = wr1 * (x7 + x9) - wi1 * (y7 + y9);
	b[(j * 5 + 2 << 1) + 2] = wr1 * (y7 + y9) + wi1 * (x7 + x9);
	b[(j * 5 + 3 << 1) + 1] = wr2 * (x8 + x10) - wi2 * (y8 + y10);
	b[(j * 5 + 3 << 1) + 2] = wr2 * (y8 + y10) + wi2 * (x8 + x10);
	b[(j * 5 + 4 << 1) + 1] = wr3 * (x8 - x10) - wi3 * (y8 - y10);
	b[(j * 5 + 4 << 1) + 2] = wr3 * (y8 - y10) + wi3 * (x8 - x10);
	b[(j * 5 + 5 << 1) + 1] = wr4 * (x7 - x9) - wi4 * (y7 - y9);
	b[(j * 5 + 5 << 1) + 2] = wr4 * (y7 - y9) + wi4 * (x7 - x9);
/* L10: */
    }
    return 0;
} /* fft5a_ */

/* Subroutine */ int fft5b_(doublereal *a, doublereal *b, doublereal *w, 
	integer *m, integer *l)
{
    /* Initialized data */

    static doublereal c52 = .61803398874989485;
    static doublereal c53 = .55901699437494742;
    static doublereal c54 = .25;
    static doublereal c51 = .95105651629515357;

    /* System generated locals */
    integer a_dim2, a_dim3, a_offset, b_dim2, b_offset, i__1, i__2;

    /* Local variables */
    integer i__, j;
    doublereal x0, y0, x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, x7, y7,
	     x8, y8, x9, y9, x10, y10, wi1, wi2, wi3, wi4, wr1, wr2, wr3, wr4;

    /* Parameter adjustments */
    w -= 3;
    b_dim2 = *m;
    b_offset = 1 + 2 * (1 + b_dim2 * 6);
    b -= b_offset;
    a_dim2 = *m;
    a_dim3 = *l;
    a_offset = 1 + 2 * (1 + a_dim2 * (1 + a_dim3));
    a -= a_offset;

    /* Function Body */

/* DIR$ VECTOR ALIGNED */
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x0 = a[(i__ + ((a_dim3 << 1) + 1) * a_dim2 << 1) + 1] + a[(i__ + (
		a_dim3 * 5 + 1) * a_dim2 << 1) + 1];
	y0 = a[(i__ + ((a_dim3 << 1) + 1) * a_dim2 << 1) + 2] + a[(i__ + (
		a_dim3 * 5 + 1) * a_dim2 << 1) + 2];
	x1 = a[(i__ + (a_dim3 * 3 + 1) * a_dim2 << 1) + 1] + a[(i__ + ((
		a_dim3 << 2) + 1) * a_dim2 << 1) + 1];
	y1 = a[(i__ + (a_dim3 * 3 + 1) * a_dim2 << 1) + 2] + a[(i__ + ((
		a_dim3 << 2) + 1) * a_dim2 << 1) + 2];
	x2 = c51 * (a[(i__ + ((a_dim3 << 1) + 1) * a_dim2 << 1) + 1] - a[(i__ 
		+ (a_dim3 * 5 + 1) * a_dim2 << 1) + 1]);
	y2 = c51 * (a[(i__ + ((a_dim3 << 1) + 1) * a_dim2 << 1) + 2] - a[(i__ 
		+ (a_dim3 * 5 + 1) * a_dim2 << 1) + 2]);
	x3 = c51 * (a[(i__ + (a_dim3 * 3 + 1) * a_dim2 << 1) + 1] - a[(i__ + (
		(a_dim3 << 2) + 1) * a_dim2 << 1) + 1]);
	y3 = c51 * (a[(i__ + (a_dim3 * 3 + 1) * a_dim2 << 1) + 2] - a[(i__ + (
		(a_dim3 << 2) + 1) * a_dim2 << 1) + 2]);
	x4 = x0 + x1;
	y4 = y0 + y1;
	x5 = c53 * (x0 - x1);
	y5 = c53 * (y0 - y1);
	x6 = a[(i__ + (a_dim3 + 1) * a_dim2 << 1) + 1] - c54 * x4;
	y6 = a[(i__ + (a_dim3 + 1) * a_dim2 << 1) + 2] - c54 * y4;
	x7 = x6 + x5;
	y7 = y6 + y5;
	x8 = x6 - x5;
	y8 = y6 - y5;
	x9 = y2 + c52 * y3;
	y9 = -x2 - c52 * x3;
	x10 = c52 * y2 - y3;
	y10 = x3 - c52 * x2;
	b[(i__ + b_dim2 * 6 << 1) + 1] = a[(i__ + (a_dim3 + 1) * a_dim2 << 1) 
		+ 1] + x4;
	b[(i__ + b_dim2 * 6 << 1) + 2] = a[(i__ + (a_dim3 + 1) * a_dim2 << 1) 
		+ 2] + y4;
	b[(i__ + b_dim2 * 7 << 1) + 1] = x7 + x9;
	b[(i__ + b_dim2 * 7 << 1) + 2] = y7 + y9;
	b[(i__ + (b_dim2 << 3) << 1) + 1] = x8 + x10;
	b[(i__ + (b_dim2 << 3) << 1) + 2] = y8 + y10;
	b[(i__ + b_dim2 * 9 << 1) + 1] = x8 - x10;
	b[(i__ + b_dim2 * 9 << 1) + 2] = y8 - y10;
	b[(i__ + b_dim2 * 10 << 1) + 1] = x7 - x9;
	b[(i__ + b_dim2 * 10 << 1) + 2] = y7 - y9;
/* L10: */
    }
    i__1 = *l;
    for (j = 2; j <= i__1; ++j) {
	wr1 = w[(j << 1) + 1];
	wi1 = w[(j << 1) + 2];
	wr2 = wr1 * wr1 - wi1 * wi1;
	wi2 = wr1 * wi1 + wr1 * wi1;
	wr3 = wr1 * wr2 - wi1 * wi2;
	wi3 = wr1 * wi2 + wi1 * wr2;
	wr4 = wr2 * wr2 - wi2 * wi2;
	wi4 = wr2 * wi2 + wr2 * wi2;
/* DIR$ VECTOR ALIGNED */
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    x0 = a[(i__ + (j + (a_dim3 << 1)) * a_dim2 << 1) + 1] + a[(i__ + (
		    j + a_dim3 * 5) * a_dim2 << 1) + 1];
	    y0 = a[(i__ + (j + (a_dim3 << 1)) * a_dim2 << 1) + 2] + a[(i__ + (
		    j + a_dim3 * 5) * a_dim2 << 1) + 2];
	    x1 = a[(i__ + (j + a_dim3 * 3) * a_dim2 << 1) + 1] + a[(i__ + (j 
		    + (a_dim3 << 2)) * a_dim2 << 1) + 1];
	    y1 = a[(i__ + (j + a_dim3 * 3) * a_dim2 << 1) + 2] + a[(i__ + (j 
		    + (a_dim3 << 2)) * a_dim2 << 1) + 2];
	    x2 = c51 * (a[(i__ + (j + (a_dim3 << 1)) * a_dim2 << 1) + 1] - a[(
		    i__ + (j + a_dim3 * 5) * a_dim2 << 1) + 1]);
	    y2 = c51 * (a[(i__ + (j + (a_dim3 << 1)) * a_dim2 << 1) + 2] - a[(
		    i__ + (j + a_dim3 * 5) * a_dim2 << 1) + 2]);
	    x3 = c51 * (a[(i__ + (j + a_dim3 * 3) * a_dim2 << 1) + 1] - a[(
		    i__ + (j + (a_dim3 << 2)) * a_dim2 << 1) + 1]);
	    y3 = c51 * (a[(i__ + (j + a_dim3 * 3) * a_dim2 << 1) + 2] - a[(
		    i__ + (j + (a_dim3 << 2)) * a_dim2 << 1) + 2]);
	    x4 = x0 + x1;
	    y4 = y0 + y1;
	    x5 = c53 * (x0 - x1);
	    y5 = c53 * (y0 - y1);
	    x6 = a[(i__ + (j + a_dim3) * a_dim2 << 1) + 1] - c54 * x4;
	    y6 = a[(i__ + (j + a_dim3) * a_dim2 << 1) + 2] - c54 * y4;
	    x7 = x6 + x5;
	    y7 = y6 + y5;
	    x8 = x6 - x5;
	    y8 = y6 - y5;
	    x9 = y2 + c52 * y3;
	    y9 = -x2 - c52 * x3;
	    x10 = c52 * y2 - y3;
	    y10 = x3 - c52 * x2;
	    b[(i__ + (j * 5 + 1) * b_dim2 << 1) + 1] = a[(i__ + (j + a_dim3) *
		     a_dim2 << 1) + 1] + x4;
	    b[(i__ + (j * 5 + 1) * b_dim2 << 1) + 2] = a[(i__ + (j + a_dim3) *
		     a_dim2 << 1) + 2] + y4;
	    b[(i__ + (j * 5 + 2) * b_dim2 << 1) + 1] = wr1 * (x7 + x9) - wi1 *
		     (y7 + y9);
	    b[(i__ + (j * 5 + 2) * b_dim2 << 1) + 2] = wr1 * (y7 + y9) + wi1 *
		     (x7 + x9);
	    b[(i__ + (j * 5 + 3) * b_dim2 << 1) + 1] = wr2 * (x8 + x10) - wi2 
		    * (y8 + y10);
	    b[(i__ + (j * 5 + 3) * b_dim2 << 1) + 2] = wr2 * (y8 + y10) + wi2 
		    * (x8 + x10);
	    b[(i__ + (j * 5 + 4) * b_dim2 << 1) + 1] = wr3 * (x8 - x10) - wi3 
		    * (y8 - y10);
	    b[(i__ + (j * 5 + 4) * b_dim2 << 1) + 2] = wr3 * (y8 - y10) + wi3 
		    * (x8 - x10);
	    b[(i__ + (j * 5 + 5) * b_dim2 << 1) + 1] = wr4 * (x7 - x9) - wi4 *
		     (y7 - y9);
	    b[(i__ + (j * 5 + 5) * b_dim2 << 1) + 2] = wr4 * (y7 - y9) + wi4 *
		     (x7 - x9);
/* L20: */
	}
/* L30: */
    }
    return 0;
} /* fft5b_ */

/* Subroutine */ int fft8a_(doublereal *a, doublereal *b, doublereal *w, 
	integer *l)
{
    /* Initialized data */

    static doublereal c81 = .70710678118654752;

    /* System generated locals */
    integer a_dim2, a_offset, i__1;

    /* Local variables */
    integer j;
    doublereal u0, v0, u1, x0, y0, x1, y1, x2, y2, x3, y3, v1, x4, y4, x5, y5,
	     x6, y6, x7, y7, u2, v2, u3, v3, wi1, wi2, wi3, wi4, wi5, wi6, 
	    wi7, wr1, wr2, wr3, wr4, wr5, wr6, wr7;

    /* Parameter adjustments */
    b -= 19;
    w -= 3;
    a_dim2 = *l;
    a_offset = 1 + 2 * (1 + a_dim2);
    a -= a_offset;

    /* Function Body */

/* DIR$ VECTOR ALIGNED */
    i__1 = *l;
    for (j = 1; j <= i__1; ++j) {
	wr1 = w[(j << 1) + 1];
	wi1 = w[(j << 1) + 2];
	wr2 = wr1 * wr1 - wi1 * wi1;
	wi2 = wr1 * wi1 + wr1 * wi1;
	wr3 = wr1 * wr2 - wi1 * wi2;
	wi3 = wr1 * wi2 + wi1 * wr2;
	wr4 = wr2 * wr2 - wi2 * wi2;
	wi4 = wr2 * wi2 + wr2 * wi2;
	wr5 = wr2 * wr3 - wi2 * wi3;
	wi5 = wr2 * wi3 + wi2 * wr3;
	wr6 = wr3 * wr3 - wi3 * wi3;
	wi6 = wr3 * wi3 + wr3 * wi3;
	wr7 = wr3 * wr4 - wi3 * wi4;
	wi7 = wr3 * wi4 + wi3 * wr4;
	x0 = a[(j + a_dim2 << 1) + 1] + a[(j + a_dim2 * 5 << 1) + 1];
	y0 = a[(j + a_dim2 << 1) + 2] + a[(j + a_dim2 * 5 << 1) + 2];
	x1 = a[(j + a_dim2 << 1) + 1] - a[(j + a_dim2 * 5 << 1) + 1];
	y1 = a[(j + a_dim2 << 1) + 2] - a[(j + a_dim2 * 5 << 1) + 2];
	x2 = a[(j + a_dim2 * 3 << 1) + 1] + a[(j + a_dim2 * 7 << 1) + 1];
	y2 = a[(j + a_dim2 * 3 << 1) + 2] + a[(j + a_dim2 * 7 << 1) + 2];
	x3 = a[(j + a_dim2 * 3 << 1) + 2] - a[(j + a_dim2 * 7 << 1) + 2];
	y3 = a[(j + a_dim2 * 7 << 1) + 1] - a[(j + a_dim2 * 3 << 1) + 1];
	u0 = x0 + x2;
	v0 = y0 + y2;
	u1 = x0 - x2;
	v1 = y0 - y2;
	x4 = a[(j + (a_dim2 << 1) << 1) + 1] + a[(j + a_dim2 * 6 << 1) + 1];
	y4 = a[(j + (a_dim2 << 1) << 1) + 2] + a[(j + a_dim2 * 6 << 1) + 2];
	x5 = a[(j + (a_dim2 << 1) << 1) + 1] - a[(j + a_dim2 * 6 << 1) + 1];
	y5 = a[(j + (a_dim2 << 1) << 1) + 2] - a[(j + a_dim2 * 6 << 1) + 2];
	x6 = a[(j + (a_dim2 << 2) << 1) + 1] + a[(j + (a_dim2 << 3) << 1) + 1]
		;
	y6 = a[(j + (a_dim2 << 2) << 1) + 2] + a[(j + (a_dim2 << 3) << 1) + 2]
		;
	x7 = a[(j + (a_dim2 << 2) << 1) + 1] - a[(j + (a_dim2 << 3) << 1) + 1]
		;
	y7 = a[(j + (a_dim2 << 2) << 1) + 2] - a[(j + (a_dim2 << 3) << 1) + 2]
		;
	u2 = x4 + x6;
	v2 = y4 + y6;
	u3 = y4 - y6;
	v3 = x6 - x4;
	b[((j << 3) + 1 << 1) + 1] = u0 + u2;
	b[((j << 3) + 1 << 1) + 2] = v0 + v2;
	b[((j << 3) + 5 << 1) + 1] = wr4 * (u0 - u2) - wi4 * (v0 - v2);
	b[((j << 3) + 5 << 1) + 2] = wr4 * (v0 - v2) + wi4 * (u0 - u2);
	b[((j << 3) + 3 << 1) + 1] = wr2 * (u1 + u3) - wi2 * (v1 + v3);
	b[((j << 3) + 3 << 1) + 2] = wr2 * (v1 + v3) + wi2 * (u1 + u3);
	b[((j << 3) + 7 << 1) + 1] = wr6 * (u1 - u3) - wi6 * (v1 - v3);
	b[((j << 3) + 7 << 1) + 2] = wr6 * (v1 - v3) + wi6 * (u1 - u3);
	u0 = x1 + c81 * (x5 - x7);
	v0 = y1 + c81 * (y5 - y7);
	u1 = x1 - c81 * (x5 - x7);
	v1 = y1 - c81 * (y5 - y7);
	u2 = x3 + c81 * (y5 + y7);
	v2 = y3 - c81 * (x5 + x7);
	u3 = x3 - c81 * (y5 + y7);
	v3 = y3 + c81 * (x5 + x7);
	b[((j << 3) + 2 << 1) + 1] = wr1 * (u0 + u2) - wi1 * (v0 + v2);
	b[((j << 3) + 2 << 1) + 2] = wr1 * (v0 + v2) + wi1 * (u0 + u2);
	b[((j << 3) + 6 << 1) + 1] = wr5 * (u1 + u3) - wi5 * (v1 + v3);
	b[((j << 3) + 6 << 1) + 2] = wr5 * (v1 + v3) + wi5 * (u1 + u3);
	b[((j << 3) + 4 << 1) + 1] = wr3 * (u1 - u3) - wi3 * (v1 - v3);
	b[((j << 3) + 4 << 1) + 2] = wr3 * (v1 - v3) + wi3 * (u1 - u3);
	b[((j << 3) + 8 << 1) + 1] = wr7 * (u0 - u2) - wi7 * (v0 - v2);
	b[((j << 3) + 8 << 1) + 2] = wr7 * (v0 - v2) + wi7 * (u0 - u2);
/* L10: */
    }
    return 0;
} /* fft8a_ */

/* Subroutine */ int fft8b_(doublereal *a, doublereal *b, doublereal *w, 
	integer *m, integer *l)
{
    /* Initialized data */

    static doublereal c81 = .70710678118654752;

    /* System generated locals */
    integer a_dim2, a_dim3, a_offset, b_dim2, b_offset, i__1, i__2;

    /* Local variables */
    integer i__, j;
    doublereal u0, v0, u1, x0, y0, x1, y1, x2, y2, x3, y3, v1, x4, y4, x5, y5,
	     x6, y6, x7, y7, u2, v2, u3, v3, wi1, wi2, wi3, wi4, wi5, wi6, 
	    wi7, wr1, wr2, wr3, wr4, wr5, wr6, wr7;

    /* Parameter adjustments */
    w -= 3;
    b_dim2 = *m;
    b_offset = 1 + 2 * (1 + b_dim2 * 9);
    b -= b_offset;
    a_dim2 = *m;
    a_dim3 = *l;
    a_offset = 1 + 2 * (1 + a_dim2 * (1 + a_dim3));
    a -= a_offset;

    /* Function Body */

/* DIR$ VECTOR ALIGNED */
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x0 = a[(i__ + (a_dim3 + 1) * a_dim2 << 1) + 1] + a[(i__ + (a_dim3 * 5 
		+ 1) * a_dim2 << 1) + 1];
	y0 = a[(i__ + (a_dim3 + 1) * a_dim2 << 1) + 2] + a[(i__ + (a_dim3 * 5 
		+ 1) * a_dim2 << 1) + 2];
	x1 = a[(i__ + (a_dim3 + 1) * a_dim2 << 1) + 1] - a[(i__ + (a_dim3 * 5 
		+ 1) * a_dim2 << 1) + 1];
	y1 = a[(i__ + (a_dim3 + 1) * a_dim2 << 1) + 2] - a[(i__ + (a_dim3 * 5 
		+ 1) * a_dim2 << 1) + 2];
	x2 = a[(i__ + (a_dim3 * 3 + 1) * a_dim2 << 1) + 1] + a[(i__ + (a_dim3 
		* 7 + 1) * a_dim2 << 1) + 1];
	y2 = a[(i__ + (a_dim3 * 3 + 1) * a_dim2 << 1) + 2] + a[(i__ + (a_dim3 
		* 7 + 1) * a_dim2 << 1) + 2];
	x3 = a[(i__ + (a_dim3 * 3 + 1) * a_dim2 << 1) + 2] - a[(i__ + (a_dim3 
		* 7 + 1) * a_dim2 << 1) + 2];
	y3 = a[(i__ + (a_dim3 * 7 + 1) * a_dim2 << 1) + 1] - a[(i__ + (a_dim3 
		* 3 + 1) * a_dim2 << 1) + 1];
	u0 = x0 + x2;
	v0 = y0 + y2;
	u1 = x0 - x2;
	v1 = y0 - y2;
	x4 = a[(i__ + ((a_dim3 << 1) + 1) * a_dim2 << 1) + 1] + a[(i__ + (
		a_dim3 * 6 + 1) * a_dim2 << 1) + 1];
	y4 = a[(i__ + ((a_dim3 << 1) + 1) * a_dim2 << 1) + 2] + a[(i__ + (
		a_dim3 * 6 + 1) * a_dim2 << 1) + 2];
	x5 = a[(i__ + ((a_dim3 << 1) + 1) * a_dim2 << 1) + 1] - a[(i__ + (
		a_dim3 * 6 + 1) * a_dim2 << 1) + 1];
	y5 = a[(i__ + ((a_dim3 << 1) + 1) * a_dim2 << 1) + 2] - a[(i__ + (
		a_dim3 * 6 + 1) * a_dim2 << 1) + 2];
	x6 = a[(i__ + ((a_dim3 << 2) + 1) * a_dim2 << 1) + 1] + a[(i__ + ((
		a_dim3 << 3) + 1) * a_dim2 << 1) + 1];
	y6 = a[(i__ + ((a_dim3 << 2) + 1) * a_dim2 << 1) + 2] + a[(i__ + ((
		a_dim3 << 3) + 1) * a_dim2 << 1) + 2];
	x7 = a[(i__ + ((a_dim3 << 2) + 1) * a_dim2 << 1) + 1] - a[(i__ + ((
		a_dim3 << 3) + 1) * a_dim2 << 1) + 1];
	y7 = a[(i__ + ((a_dim3 << 2) + 1) * a_dim2 << 1) + 2] - a[(i__ + ((
		a_dim3 << 3) + 1) * a_dim2 << 1) + 2];
	u2 = x4 + x6;
	v2 = y4 + y6;
	u3 = y4 - y6;
	v3 = x6 - x4;
	b[(i__ + b_dim2 * 9 << 1) + 1] = u0 + u2;
	b[(i__ + b_dim2 * 9 << 1) + 2] = v0 + v2;
	b[(i__ + b_dim2 * 13 << 1) + 1] = u0 - u2;
	b[(i__ + b_dim2 * 13 << 1) + 2] = v0 - v2;
	b[(i__ + b_dim2 * 11 << 1) + 1] = u1 + u3;
	b[(i__ + b_dim2 * 11 << 1) + 2] = v1 + v3;
	b[(i__ + b_dim2 * 15 << 1) + 1] = u1 - u3;
	b[(i__ + b_dim2 * 15 << 1) + 2] = v1 - v3;
	u0 = x1 + c81 * (x5 - x7);
	v0 = y1 + c81 * (y5 - y7);
	u1 = x1 - c81 * (x5 - x7);
	v1 = y1 - c81 * (y5 - y7);
	u2 = x3 + c81 * (y5 + y7);
	v2 = y3 - c81 * (x5 + x7);
	u3 = x3 - c81 * (y5 + y7);
	v3 = y3 + c81 * (x5 + x7);
	b[(i__ + b_dim2 * 10 << 1) + 1] = u0 + u2;
	b[(i__ + b_dim2 * 10 << 1) + 2] = v0 + v2;
	b[(i__ + b_dim2 * 14 << 1) + 1] = u1 + u3;
	b[(i__ + b_dim2 * 14 << 1) + 2] = v1 + v3;
	b[(i__ + b_dim2 * 12 << 1) + 1] = u1 - u3;
	b[(i__ + b_dim2 * 12 << 1) + 2] = v1 - v3;
	b[(i__ + (b_dim2 << 4) << 1) + 1] = u0 - u2;
	b[(i__ + (b_dim2 << 4) << 1) + 2] = v0 - v2;
/* L10: */
    }
    i__1 = *l;
    for (j = 2; j <= i__1; ++j) {
	wr1 = w[(j << 1) + 1];
	wi1 = w[(j << 1) + 2];
	wr2 = wr1 * wr1 - wi1 * wi1;
	wi2 = wr1 * wi1 + wr1 * wi1;
	wr3 = wr1 * wr2 - wi1 * wi2;
	wi3 = wr1 * wi2 + wi1 * wr2;
	wr4 = wr2 * wr2 - wi2 * wi2;
	wi4 = wr2 * wi2 + wr2 * wi2;
	wr5 = wr2 * wr3 - wi2 * wi3;
	wi5 = wr2 * wi3 + wi2 * wr3;
	wr6 = wr3 * wr3 - wi3 * wi3;
	wi6 = wr3 * wi3 + wr3 * wi3;
	wr7 = wr3 * wr4 - wi3 * wi4;
	wi7 = wr3 * wi4 + wi3 * wr4;
/* DIR$ VECTOR ALIGNED */
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    x0 = a[(i__ + (j + a_dim3) * a_dim2 << 1) + 1] + a[(i__ + (j + 
		    a_dim3 * 5) * a_dim2 << 1) + 1];
	    y0 = a[(i__ + (j + a_dim3) * a_dim2 << 1) + 2] + a[(i__ + (j + 
		    a_dim3 * 5) * a_dim2 << 1) + 2];
	    x1 = a[(i__ + (j + a_dim3) * a_dim2 << 1) + 1] - a[(i__ + (j + 
		    a_dim3 * 5) * a_dim2 << 1) + 1];
	    y1 = a[(i__ + (j + a_dim3) * a_dim2 << 1) + 2] - a[(i__ + (j + 
		    a_dim3 * 5) * a_dim2 << 1) + 2];
	    x2 = a[(i__ + (j + a_dim3 * 3) * a_dim2 << 1) + 1] + a[(i__ + (j 
		    + a_dim3 * 7) * a_dim2 << 1) + 1];
	    y2 = a[(i__ + (j + a_dim3 * 3) * a_dim2 << 1) + 2] + a[(i__ + (j 
		    + a_dim3 * 7) * a_dim2 << 1) + 2];
	    x3 = a[(i__ + (j + a_dim3 * 3) * a_dim2 << 1) + 2] - a[(i__ + (j 
		    + a_dim3 * 7) * a_dim2 << 1) + 2];
	    y3 = a[(i__ + (j + a_dim3 * 7) * a_dim2 << 1) + 1] - a[(i__ + (j 
		    + a_dim3 * 3) * a_dim2 << 1) + 1];
	    u0 = x0 + x2;
	    v0 = y0 + y2;
	    u1 = x0 - x2;
	    v1 = y0 - y2;
	    x4 = a[(i__ + (j + (a_dim3 << 1)) * a_dim2 << 1) + 1] + a[(i__ + (
		    j + a_dim3 * 6) * a_dim2 << 1) + 1];
	    y4 = a[(i__ + (j + (a_dim3 << 1)) * a_dim2 << 1) + 2] + a[(i__ + (
		    j + a_dim3 * 6) * a_dim2 << 1) + 2];
	    x5 = a[(i__ + (j + (a_dim3 << 1)) * a_dim2 << 1) + 1] - a[(i__ + (
		    j + a_dim3 * 6) * a_dim2 << 1) + 1];
	    y5 = a[(i__ + (j + (a_dim3 << 1)) * a_dim2 << 1) + 2] - a[(i__ + (
		    j + a_dim3 * 6) * a_dim2 << 1) + 2];
	    x6 = a[(i__ + (j + (a_dim3 << 2)) * a_dim2 << 1) + 1] + a[(i__ + (
		    j + (a_dim3 << 3)) * a_dim2 << 1) + 1];
	    y6 = a[(i__ + (j + (a_dim3 << 2)) * a_dim2 << 1) + 2] + a[(i__ + (
		    j + (a_dim3 << 3)) * a_dim2 << 1) + 2];
	    x7 = a[(i__ + (j + (a_dim3 << 2)) * a_dim2 << 1) + 1] - a[(i__ + (
		    j + (a_dim3 << 3)) * a_dim2 << 1) + 1];
	    y7 = a[(i__ + (j + (a_dim3 << 2)) * a_dim2 << 1) + 2] - a[(i__ + (
		    j + (a_dim3 << 3)) * a_dim2 << 1) + 2];
	    u2 = x4 + x6;
	    v2 = y4 + y6;
	    u3 = y4 - y6;
	    v3 = x6 - x4;
	    b[(i__ + ((j << 3) + 1) * b_dim2 << 1) + 1] = u0 + u2;
	    b[(i__ + ((j << 3) + 1) * b_dim2 << 1) + 2] = v0 + v2;
	    b[(i__ + ((j << 3) + 5) * b_dim2 << 1) + 1] = wr4 * (u0 - u2) - 
		    wi4 * (v0 - v2);
	    b[(i__ + ((j << 3) + 5) * b_dim2 << 1) + 2] = wr4 * (v0 - v2) + 
		    wi4 * (u0 - u2);
	    b[(i__ + ((j << 3) + 3) * b_dim2 << 1) + 1] = wr2 * (u1 + u3) - 
		    wi2 * (v1 + v3);
	    b[(i__ + ((j << 3) + 3) * b_dim2 << 1) + 2] = wr2 * (v1 + v3) + 
		    wi2 * (u1 + u3);
	    b[(i__ + ((j << 3) + 7) * b_dim2 << 1) + 1] = wr6 * (u1 - u3) - 
		    wi6 * (v1 - v3);
	    b[(i__ + ((j << 3) + 7) * b_dim2 << 1) + 2] = wr6 * (v1 - v3) + 
		    wi6 * (u1 - u3);
	    u0 = x1 + c81 * (x5 - x7);
	    v0 = y1 + c81 * (y5 - y7);
	    u1 = x1 - c81 * (x5 - x7);
	    v1 = y1 - c81 * (y5 - y7);
	    u2 = x3 + c81 * (y5 + y7);
	    v2 = y3 - c81 * (x5 + x7);
	    u3 = x3 - c81 * (y5 + y7);
	    v3 = y3 + c81 * (x5 + x7);
	    b[(i__ + ((j << 3) + 2) * b_dim2 << 1) + 1] = wr1 * (u0 + u2) - 
		    wi1 * (v0 + v2);
	    b[(i__ + ((j << 3) + 2) * b_dim2 << 1) + 2] = wr1 * (v0 + v2) + 
		    wi1 * (u0 + u2);
	    b[(i__ + ((j << 3) + 6) * b_dim2 << 1) + 1] = wr5 * (u1 + u3) - 
		    wi5 * (v1 + v3);
	    b[(i__ + ((j << 3) + 6) * b_dim2 << 1) + 2] = wr5 * (v1 + v3) + 
		    wi5 * (u1 + u3);
	    b[(i__ + ((j << 3) + 4) * b_dim2 << 1) + 1] = wr3 * (u1 - u3) - 
		    wi3 * (v1 - v3);
	    b[(i__ + ((j << 3) + 4) * b_dim2 << 1) + 2] = wr3 * (v1 - v3) + 
		    wi3 * (u1 - u3);
	    b[(i__ + ((j << 3) + 8) * b_dim2 << 1) + 1] = wr7 * (u0 - u2) - 
		    wi7 * (v0 - v2);
	    b[(i__ + ((j << 3) + 8) * b_dim2 << 1) + 2] = wr7 * (v0 - v2) + 
		    wi7 * (u0 - u2);
/* L20: */
	}
/* L30: */
    }
    return 0;
} /* fft8b_ */

