/* -*- mode: C; tab-width: 2; indent-tabs-mode: nil; -*- */

#include <hpcc.h>

/* Common Block Declarations */

struct {
    int irand[2], ias[2], ics[2];
} rancom_;

#define rancom_1 rancom_

static int
ladd_(int *j, int *k, int *i__) {

/*  -- ScaLAPACK routine (version 1.7) -- */
/*     University of Tennessee, Knoxville, Oak Ridge National Laboratory, */
/*     and University of California, Berkeley. */
/*     May 1, 1997 */

    /* Parameter adjustments */
    --i__;
    --k;
    --j;

    /* Function Body */
    i__[1] = (k[1] + j[1]) % 65536;
    i__[2] = ((k[1] + j[1]) / 65536 + k[2] + j[2]) % 32768;

    return 0;
} /* ladd_ */


static int
lmul_(int *k, int *j, int *i__) {
    int kt, lt;

/*  -- ScaLAPACK routine (version 1.7) -- */
/*     University of Tennessee, Knoxville, Oak Ridge National Laboratory, */
/*     and University of California, Berkeley. */
/*     May 1, 1997 */

    /* Parameter adjustments */
    --i__;
    --j;
    --k;

    /* Function Body */
    kt = k[1] * j[1];
    if (kt < 0) {
	kt += -2147483647;
	kt += -1;
    }
    i__[1] = kt % 65536;
    lt = k[1] * j[2] + k[2] * j[1];
    if (lt < 0) {
	lt += -2147483647;
	lt += -1;
    }
    kt = kt / 65536 + lt;
    if (kt < 0) {
	kt += -2147483647;
	kt += -1;
    }
    i__[2] = kt % 32768;

    return 0;
} /* lmul_ */


int
xjumpm_(int *jumpm, int *mult, int *iadd, 
	int *irann, int *iranm, int *iam, int *icm) {
    /* System generated locals */
    int i__1;

    int i__, j[2];

/*  -- ScaLAPACK routine (version 1.7) -- */
/*     University of Tennessee, Knoxville, Oak Ridge National Laboratory, */
/*     and University of California, Berkeley. */
/*     May 1, 1997 */

    /* Parameter adjustments */
    --icm;
    --iam;
    --iranm;
    --irann;
    --iadd;
    --mult;

    if (*jumpm > 0) {
	for (i__ = 1; i__ <= 2; ++i__) {
	    iam[i__] = mult[i__];
	    icm[i__] = iadd[i__];
/* L10: */
	}
	i__1 = *jumpm - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    lmul_(&iam[1], &mult[1], j);
	    iam[1] = j[0];
	    iam[2] = j[1];
	    lmul_(&icm[1], &mult[1], j);
	    ladd_(&iadd[1], j, &icm[1]);
/* L20: */
	}
	lmul_(&irann[1], &iam[1], j);
	ladd_(j, &icm[1], &iranm[1]);
    } else {
	iranm[1] = irann[1];
	iranm[2] = irann[2];
    }

    return 0;
} /* xjumpm_ */


int setran_(int *iran, int *ia, int *ic) {
    int i__;
/*  -- ScaLAPACK routine (version 1.7) -- */
/*     University of Tennessee, Knoxville, Oak Ridge National Laboratory, */
/*     and University of California, Berkeley. */
/*     May 1, 1997 */

    /* Parameter adjustments */
    --ic;
    --ia;
    --iran;

    for (i__ = 1; i__ <= 2; ++i__) {
	rancom_1.irand[i__ - 1] = iran[i__];
	rancom_1.ias[i__ - 1] = ia[i__];
	rancom_1.ics[i__ - 1] = ic[i__];
/* L10: */
    }

    return 0;
} /* setran_ */


int jumpit_(int *mult, int *iadd, int *irann, int *iranm) {
    int j[2];

/*  -- ScaLAPACK routine (version 1.7) -- */
/*     University of Tennessee, Knoxville, Oak Ridge National Laboratory, */
/*     and University of California, Berkeley. */
/*     May 1, 1997 */

    /* Parameter adjustments */
    --iranm;
    --irann;
    --iadd;
    --mult;

    lmul_(&irann[1], &mult[1], j);
    ladd_(j, &iadd[1], &iranm[1]);

    rancom_1.irand[0] = iranm[1];
    rancom_1.irand[1] = iranm[2];

    return 0;
} /* jumpit_ */

double
pdrand() {
    /* System generated locals */
    double ret_val;

    /* Local variables */
    static int j[2];
/*  -- ScaLAPACK routine (version 1.7) -- */
/*     University of Tennessee, Knoxville, Oak Ridge National Laboratory, */
/*     and University of California, Berkeley. */
/*     May 1, 1997 */

    ret_val = ((double) rancom_1.irand[0] + (double) rancom_1.irand[1]
	     * 65536.) / 2147483648.;

    lmul_(rancom_1.irand, rancom_1.ias, j);
    ladd_(j, rancom_1.ics, rancom_1.irand);

    return ret_val;
} /* pdrand */
