/* ztrans.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "hpccfft.h"
#define min(a,b) ((a)<=(b)?(a):(b))


/*     FFTE: A FAST FOURIER TRANSFORM PACKAGE */

/*     (C) COPYRIGHT SOFTWARE, 2000-2004, ALL RIGHTS RESERVED */
/*                BY */
/*         DAISUKE TAKAHASHI */
/*         INSTITURE OF INFORMATION SCIENCES AND ELECTRONICS */
/*         UNIVERSITY OF TSUKUBA */
/*         1-1-1 TENNODAI, TSUKUBA, IBARAKI 305-8573, JAPAN */
/*         E-MAIL: daisuke@is.tsukuba.ac.jp */


/*     TRANSPOSE ROUTINE */

/*     FORTRAN77 SOURCE PROGRAM */

/*     CALL ZTRANS(A,B,N1,N2) */

/*     WRITTEN BY DAISUKE TAKAHASHI */

/* Subroutine */ int ztrans_(doublecomplex *a, doublecomplex *b, integer *n1, 
	integer *n2)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3, i__4, i__5, 
	    i__6;

    /* Local variables */
    integer i__, j, ii, jj;


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
    a_dim1 = *n1;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *n2;
    b_offset = 1 + b_dim1;
    b -= b_offset;

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
		    i__5 = j + i__ * b_dim1;
		    i__6 = i__ + j * a_dim1;
		    b[i__5].r = a[i__6].r, b[i__5].i = a[i__6].i;
/* L10: */
		}
/* L20: */
	    }
/* L30: */
	}
/* L40: */
    }
    return 0;
} /* ztrans_ */

