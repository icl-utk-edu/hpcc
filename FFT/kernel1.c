/* kernel1.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "hpccfft.h"

/* Table of constant values */

static integer c__8 = 8;
static integer c__5 = 5;
static integer c__4 = 4;
static integer c__3 = 3;


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

/* Subroutine */ int fft235a_(doublecomplex *a, doublecomplex *b, 
	doublecomplex *w, integer *n, integer *ip)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer j, k, l, m, kp4, kp8, key;
    extern /* Subroutine */ int fft2_(doublecomplex *, doublecomplex *, 
	    integer *), fft3_(doublecomplex *, doublecomplex *, doublecomplex 
	    *, integer *, integer *), fft4_(doublecomplex *, doublecomplex *, 
	    doublecomplex *, integer *, integer *), fft5_(doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, integer *), fft8_(
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    integer *);


    /* Parameter adjustments */
    --ip;
    --w;
    --b;
    --a;

    /* Function Body */
    if (ip[1] != 1) {
	kp4 = 2 - (ip[1] + 2) % 3;
	kp8 = (ip[1] - kp4) / 3;
    } else {
	kp4 = 0;
	kp8 = 0;
    }

    key = 1;
    j = 1;
    l = *n;
    m = 1;
    i__1 = kp8;
    for (k = 1; k <= i__1; ++k) {
	l /= 8;
	if (l >= 2) {
	    if (key >= 0) {
		fft8_(&a[1], &b[1], &w[j], &m, &l);
	    } else {
		fft8_(&b[1], &a[1], &w[j], &m, &l);
	    }
	    key = -key;
	} else {
	    if (key >= 0) {
		fft8_(&a[1], &a[1], &w[j], &m, &l);
	    } else {
		fft8_(&b[1], &a[1], &w[j], &m, &l);
	    }
	}
	m <<= 3;
	j += l;
/* L10: */
    }
    i__1 = ip[3];
    for (k = 1; k <= i__1; ++k) {
	l /= 5;
	if (l >= 2) {
	    if (key >= 0) {
		fft5_(&a[1], &b[1], &w[j], &m, &l);
	    } else {
		fft5_(&b[1], &a[1], &w[j], &m, &l);
	    }
	    key = -key;
	} else {
	    if (key >= 0) {
		fft5_(&a[1], &a[1], &w[j], &m, &l);
	    } else {
		fft5_(&b[1], &a[1], &w[j], &m, &l);
	    }
	}
	m *= 5;
	j += l;
/* L20: */
    }
    i__1 = kp4;
    for (k = 1; k <= i__1; ++k) {
	l /= 4;
	if (l >= 2) {
	    if (key >= 0) {
		fft4_(&a[1], &b[1], &w[j], &m, &l);
	    } else {
		fft4_(&b[1], &a[1], &w[j], &m, &l);
	    }
	    key = -key;
	} else {
	    if (key >= 0) {
		fft4_(&a[1], &a[1], &w[j], &m, &l);
	    } else {
		fft4_(&b[1], &a[1], &w[j], &m, &l);
	    }
	}
	m <<= 2;
	j += l;
/* L30: */
    }
    i__1 = ip[2];
    for (k = 1; k <= i__1; ++k) {
	l /= 3;
	if (l >= 2) {
	    if (key >= 0) {
		fft3_(&a[1], &b[1], &w[j], &m, &l);
	    } else {
		fft3_(&b[1], &a[1], &w[j], &m, &l);
	    }
	    key = -key;
	} else {
	    if (key >= 0) {
		fft3_(&a[1], &a[1], &w[j], &m, &l);
	    } else {
		fft3_(&b[1], &a[1], &w[j], &m, &l);
	    }
	}
	m *= 3;
	j += l;
/* L40: */
    }
    if (ip[1] == 1) {
	if (key >= 0) {
	    fft2_(&a[1], &a[1], &m);
	} else {
	    fft2_(&b[1], &a[1], &m);
	}
    }
    return 0;
} /* fft235a_ */

/* Subroutine */ int fft235b_(doublecomplex *a, doublecomplex *b, 
	doublecomplex *w, integer *n, integer *ip)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer j, k, l, m, kp4, kp8, key;
    extern /* Subroutine */ int fft2_(doublecomplex *, doublecomplex *, 
	    integer *), fft3_(doublecomplex *, doublecomplex *, doublecomplex 
	    *, integer *, integer *), fft4_(doublecomplex *, doublecomplex *, 
	    doublecomplex *, integer *, integer *), fft5_(doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, integer *), fft8_(
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    integer *);


    /* Parameter adjustments */
    --ip;
    --w;
    --b;
    --a;

    /* Function Body */
    if (ip[1] != 1) {
	kp4 = 2 - (ip[1] + 2) % 3;
	kp8 = (ip[1] - kp4) / 3;
    } else {
	kp4 = 0;
	kp8 = 0;
    }

    key = 1;
    j = 1;
    l = *n;
    m = 1;
    i__1 = kp8;
    for (k = 1; k <= i__1; ++k) {
	l /= 8;
	if (l >= 2) {
	    if (key >= 0) {
		fft8_(&a[1], &b[1], &w[j], &m, &l);
	    } else {
		fft8_(&b[1], &a[1], &w[j], &m, &l);
	    }
	    key = -key;
	} else {
	    if (key >= 0) {
		fft8_(&a[1], &b[1], &w[j], &m, &l);
	    } else {
		fft8_(&b[1], &b[1], &w[j], &m, &l);
	    }
	}
	m <<= 3;
	j += l;
/* L10: */
    }
    i__1 = ip[3];
    for (k = 1; k <= i__1; ++k) {
	l /= 5;
	if (l >= 2) {
	    if (key >= 0) {
		fft5_(&a[1], &b[1], &w[j], &m, &l);
	    } else {
		fft5_(&b[1], &a[1], &w[j], &m, &l);
	    }
	    key = -key;
	} else {
	    if (key >= 0) {
		fft5_(&a[1], &b[1], &w[j], &m, &l);
	    } else {
		fft5_(&b[1], &b[1], &w[j], &m, &l);
	    }
	}
	m *= 5;
	j += l;
/* L20: */
    }
    i__1 = kp4;
    for (k = 1; k <= i__1; ++k) {
	l /= 4;
	if (l >= 2) {
	    if (key >= 0) {
		fft4_(&a[1], &b[1], &w[j], &m, &l);
	    } else {
		fft4_(&b[1], &a[1], &w[j], &m, &l);
	    }
	    key = -key;
	} else {
	    if (key >= 0) {
		fft4_(&a[1], &b[1], &w[j], &m, &l);
	    } else {
		fft4_(&b[1], &b[1], &w[j], &m, &l);
	    }
	}
	m <<= 2;
	j += l;
/* L30: */
    }
    i__1 = ip[2];
    for (k = 1; k <= i__1; ++k) {
	l /= 3;
	if (l >= 2) {
	    if (key >= 0) {
		fft3_(&a[1], &b[1], &w[j], &m, &l);
	    } else {
		fft3_(&b[1], &a[1], &w[j], &m, &l);
	    }
	    key = -key;
	} else {
	    if (key >= 0) {
		fft3_(&a[1], &b[1], &w[j], &m, &l);
	    } else {
		fft3_(&b[1], &b[1], &w[j], &m, &l);
	    }
	}
	m *= 3;
	j += l;
/* L40: */
    }
    if (ip[1] == 1) {
	if (key >= 0) {
	    fft2_(&a[1], &b[1], &m);
	} else {
	    fft2_(&b[1], &b[1], &m);
	}
    }
    return 0;
} /* fft235b_ */

/* Subroutine */ int fft3_(doublecomplex *a, doublecomplex *b, doublecomplex *
	w, integer *m, integer *l)
{
    extern /* Subroutine */ int fft3a_(doublecomplex *, doublecomplex *, 
	    doublecomplex *, integer *), fft3b_(doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, integer *);


    /* Parameter adjustments */
    --w;
    --b;
    --a;

    /* Function Body */
    if (*m == 1) {
	fft3a_(&a[1], &b[1], &w[1], l);
    } else {
	fft3b_(&a[1], &b[1], &w[1], m, l);
    }
    return 0;
} /* fft3_ */

/* Subroutine */ int fft4_(doublecomplex *a, doublecomplex *b, doublecomplex *
	w, integer *m, integer *l)
{
    extern /* Subroutine */ int fft4a_(doublecomplex *, doublecomplex *, 
	    doublecomplex *, integer *), fft4b_(doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, integer *);


    /* Parameter adjustments */
    --w;
    --b;
    --a;

    /* Function Body */
    if (*m == 1) {
	fft4a_(&a[1], &b[1], &w[1], l);
    } else {
	fft4b_(&a[1], &b[1], &w[1], m, l);
    }
    return 0;
} /* fft4_ */

/* Subroutine */ int fft5_(doublecomplex *a, doublecomplex *b, doublecomplex *
	w, integer *m, integer *l)
{
    extern /* Subroutine */ int fft5a_(doublecomplex *, doublecomplex *, 
	    doublecomplex *, integer *), fft5b_(doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, integer *);


    /* Parameter adjustments */
    --w;
    --b;
    --a;

    /* Function Body */
    if (*m == 1) {
	fft5a_(&a[1], &b[1], &w[1], l);
    } else {
	fft5b_(&a[1], &b[1], &w[1], m, l);
    }
    return 0;
} /* fft5_ */

/* Subroutine */ int fft8_(doublecomplex *a, doublecomplex *b, doublecomplex *
	w, integer *m, integer *l)
{
    extern /* Subroutine */ int fft8a_(doublecomplex *, doublecomplex *, 
	    doublecomplex *, integer *), fft8b_(doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, integer *);


    /* Parameter adjustments */
    --w;
    --b;
    --a;

    /* Function Body */
    if (*m == 1) {
	fft8a_(&a[1], &b[1], &w[1], l);
    } else {
	fft8b_(&a[1], &b[1], &w[1], m, l);
    }
    return 0;
} /* fft8_ */

/* Subroutine */ int settbl_(doublecomplex *w, integer *n)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer j, k, l, ip[3], kp4, kp8;
    extern /* Subroutine */ int factor_(integer *, integer *), settbl2_(
	    void *, integer *, integer *);


    /* Parameter adjustments */
    --w;

    /* Function Body */
    factor_(n, ip);

    if (ip[0] != 1) {
	kp4 = 2 - (ip[0] + 2) % 3;
	kp8 = (ip[0] - kp4) / 3;
    } else {
	kp4 = 0;
	kp8 = 0;
    }

    j = 1;
    l = *n;
    i__1 = kp8;
    for (k = 1; k <= i__1; ++k) {
	l /= 8;
	settbl2_(&w[j], &l, &c__8);
	j += l;
/* L10: */
    }
    i__1 = ip[2];
    for (k = 1; k <= i__1; ++k) {
	l /= 5;
	settbl2_(&w[j], &l, &c__5);
	j += l;
/* L20: */
    }
    i__1 = kp4;
    for (k = 1; k <= i__1; ++k) {
	l /= 4;
	settbl2_(&w[j], &l, &c__4);
	j += l;
/* L30: */
    }
    i__1 = ip[1];
    for (k = 1; k <= i__1; ++k) {
	l /= 3;
	settbl2_(&w[j], &l, &c__3);
	j += l;
/* L40: */
    }
    return 0;
} /* settbl_ */

/* Subroutine */ int settbl2_(void *w_, integer *l, integer *m)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__;
    doublereal px, pi2;
    doublereal *w = (doublereal *)w_;


    /* Parameter adjustments */
    w -= 3;

    /* Function Body */
    pi2 = atan(1.) * 8.;
    px = -pi2 / (doublereal) (*l * *m);
/* DIR$ VECTOR ALIGNED */
    i__1 = *l;
    for (i__ = 1; i__ <= i__1; ++i__) {
	w[(i__ << 1) + 1] = cos(px * (doublereal) (i__ - 1));
	w[(i__ << 1) + 2] = sin(px * (doublereal) (i__ - 1));
/* L10: */
    }
    return 0;
} /* settbl2_ */

/* Subroutine */ int factor_(integer *n, integer *ip)
{
    integer n2;


    /* Parameter adjustments */
    --ip;

    /* Function Body */
    ip[1] = 0;
    ip[2] = 0;
    ip[3] = 0;
    n2 = *n;
    if (*n % 2 != 0 && *n % 3 != 0 && *n % 5 != 0) {
	return 0;
    }
L10:
    if (n2 <= 1) {
	return 0;
    }
    if (n2 % 2 == 0) {
	++ip[1];
	n2 /= 2;
	goto L10;
    } else if (n2 % 3 == 0) {
	++ip[2];
	n2 /= 3;
	goto L10;
    } else if (n2 % 5 == 0) {
	++ip[3];
	n2 /= 5;
	goto L10;
    }
    return 0;
} /* factor_ */

/* Subroutine */ int factor8_(s64Int_t *n, integer *ip)
{
    s64Int_t n2;

    ip[0] = 0;
    ip[1] = 0;
    ip[2] = 0;
    n2 = *n;
    if (n2 % 2 != 0 && n2 % 3 != 0 && n2 % 5 != 0) {
	return 0;
    }
L10:
    if (n2 <= 1) {
	return 0;
    }
    if (n2 % 2 == 0) {
	++ip[0];
	n2 /= 2;
	goto L10;
    } else if (n2 % 3 == 0) {
	++ip[1];
	n2 /= 3;
	goto L10;
    } else if (n2 % 5 == 0) {
	++ip[2];
	n2 /= 5;
	goto L10;
    }
    return 0;
} /* factor8_ */
