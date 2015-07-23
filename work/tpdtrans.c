
#include <stdio.h>
#include <stdlib.h>

int Cblacs_nprow, Cblacs_npcol, Cblacs_myrow, Cblacs_mycol;

void
Cblacs_gridinfo(int ctxt, int *nprow, int *npcol, int *myrow, int *mycol) {
  *nprow = Cblacs_nprow;
  *npcol = Cblacs_npcol;
  *myrow = Cblacs_myrow;
  *mycol = Cblacs_mycol;
}

void
Cblacs_dSendrecv(int ctxt, int mSrc, int nSrc, double *Asrc, int ldaSrc, int rdest, int cdest,
  int mDest, int nDest, double *Adest, int ldaDest, int rsrc, int csrc) {
  int src, dest, dataIsContiguousSrc, dataIsContiguousDest, countSrc, countDest, npcol, nprow, myrow, mycol;
  /*
  MPI_Datatype typeSrc, typeDest;
  MPI_Status stat;
  */


  if (mSrc == ldaSrc || 1 == nSrc) {
    dataIsContiguousSrc = 1;
    countSrc = mSrc * nSrc;
    /*
    typeSrc = MPI_DOUBLE;
    */
  } else {
    dataIsContiguousSrc = 0;
    countSrc = 1;
    /*
    MPI_Type_vector( nSrc, mSrc, ldaSrc, MPI_DOUBLE, &typeSrc );
    MPI_Type_commit( &typeSrc );
    */
  }
  if (mDest == ldaDest || 1 == nDest) {
    dataIsContiguousDest = 1;
    countDest = mDest * nDest;
    /*
    typeDest = MPI_DOUBLE;
    */
  } else {
    dataIsContiguousDest = 0;
    countDest = 1;
    /*
    MPI_Type_vector( nDest, mDest, ldaDest, MPI_DOUBLE, &typeDest );
    MPI_Type_commit( &typeDest );
    */
  }

  /*
  rowComm = CblacsGetRowComm( ctxt );
  MPI_Comm_size( rowComm, &npcol );
  */

  Cblacs_gridinfo( 1, &nprow, &npcol, &myrow, &mycol);

  dest = cdest + rdest * npcol;
  src  = csrc + rsrc * npcol;

  /*
  MPI_Sendrecv( Asrc, countSrc, typeSrc, dest, 0, Adest, countDest, typeDest, src, 0, comm, &stat );
   */
  printf( "%s(%d): %d %d %d %d %d %d\n", __FILE__, __LINE__, mSrc, nSrc, dataIsContiguousSrc, dataIsContiguousDest, countSrc, countDest );

  /*
  if (! dataIsContiguousSrc) MPI_Type_free( &typeSrc );
  if (! dataIsContiguousDest) MPI_Type_free( &typeDest );
  */
}

void
HPL_dscal(int n, double alpha, double *x,  int incx){
  if (n < 1 || alpha != alpha || ! x || incx == 0)
    abort();
}

void
dtr2mx_(double *A, int *lda, double *beta, double *B, int *ldb, int *m, int *n, int *mb, int *nb, int *r, int *c) {
  if (!A || ! lda || ! beta || ! B || ! ldb || ! m || ! n || ! mb || ! nb || ! r || ! c)
    abort();
}

void
dtr2bf_(double *A, int *lda, double *B, int *ldb, int *m, int *n, int *mb, int *nb, int *r, int *c) {
  if (!A || ! lda || ! B || ! ldb || ! m || ! n || ! mb || ! nb || ! r || ! c)
    abort();
}

void
dmv2mx_(double *A, int *lda, double *beta, double *B, int *ldb, int *m, int *n, int *mb, int *nb, int *r, int *c) {
  if (!A || ! lda || ! beta || ! B || ! ldb || ! m || ! n || ! mb || ! nb || ! r || ! c)
    abort();
}

struct {
    int iaz, jaz, itz, jtz;
} commtrb_;

#define commtrb_1 commtrb_

static void
pxerbla(int *ctxt, const char *rname, int *info) {
  fprintf( stderr, "%d %s %d\n", *ctxt, rname, *info );
}

int ilcm_(int *m, int *n)
{
  /* System generated locals */
  int ret_val;

  /* Local variables */
  int ia, iq, ir;

  if (*m >= *n) {
    ia = *m;
    ret_val = *n;
  } else {
    ia = *n;
    ret_val = *m;
  }

  for (;;) {
    iq = ia / ret_val;
    ir = ia - iq * ret_val;
    if (ir == 0) {
      ret_val = *m * *n / ret_val;
      return ret_val;
    }
    ia = ret_val;
    ret_val = ir;
  }
} /* ilcm_ */

int
iceil_(int *inum, int *idenom) {
    /* System generated locals */
    int ret_val;

    ret_val = (*inum + *idenom - 1) / *idenom;

    return ret_val;
} /* iceil_ */

static int
numroc_(int *n, int *nb, int *iproc, int *isrcproc, int *nprocs) {
    /* System generated locals */
    int ret_val;

    /* Local variables */
    int extrablks, mydist, nblocks;


/*     Figure PROC's distance from source process */

    mydist = (*nprocs + *iproc - *isrcproc) % *nprocs;

/*     Figure the total number of whole NB blocks N is split up into */

    nblocks = *n / *nb;

/*     Figure the minimum number of rows/cols a process can have */

    ret_val = nblocks / *nprocs * *nb;

/*     See if there are any extra blocks */

    extrablks = nblocks % *nprocs;

/*     If I have an extra block */

    if (mydist < extrablks) {
        ret_val += *nb;

/*         If I have last block, it may be a partial block */

    } else if (mydist == extrablks) {
        ret_val += *n % *nb;
    }

    return ret_val;
} /* numroc_ */
        

static int c__0 = 0;

int
pdtrans(char *trans, int *m, int *n, int * mb, int *nb, double *a, int *lda, double *beta,
	double *c__, int *ldc, int *imrow, int *imcol, double *work, int *iwork) {
  int ctxt = 1;
    /* System generated locals */
    long a_dim1, a_offset, c_dim1, c_offset;
    int i__1, i__2, i__3, i__4;

    /* Local variables */
    int j1, k1, k2, ml, nl, mp, mq, np, nq, mb0, mb1, mb2, nb0,
	    nb1, nb2, kia, kja, kic, kjc, lbm, lbn, lcm, ldt, lbm0, lbm1,
	     lbm2, lbn0, lbn1, lbn2, igcd;
    long ipt;
    int mcol, info, lcmp, lcmq, item, ncol, kmod1, kmod2;
    double tbeta;
    int kpcol, mpcol, npcol, mrcol, mycol, kprow, mprow, nprow, mrrow, myrow;

/*     Get grid parameters */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;
    --iwork;

    /* Function Body */
    Cblacs_gridinfo( 1, &nprow, &npcol, &myrow, &mycol);

/*     Test for the input parameters. */

    info = 0;
    if (*trans != 'T' && *trans != 'C') {
	info = 1;
    } else if (*m < 0) {
	info = 2;
    } else if (*n < 0) {
	info = 3;
    } else if (*mb < 1) {
	info = 4;
    } else if (*nb < 1) {
	info = 5;
    } else if (*lda < 1) {
	info = 7;
    } else if (*ldc < 1) {
	info = 10;
    } else if (*imrow < 0 || *imrow >= nprow) {
	info = 11;
    } else if (*imcol < 0 || *imcol >= npcol) {
	info = 12;
    }

L10:
    if (info != 0) {
	pxerbla( &ctxt, "PDTRANS", &info );
	return 0;
    }

/*     Initialize parameters */

    mprow = nprow + myrow;
    mpcol = npcol + mycol;
    mrrow = (mprow - *imrow) % nprow;
    mrcol = (mpcol - *imcol) % npcol;

    lcm = ilcm_(&nprow, &npcol);
    lcmp = lcm / nprow;
    lcmq = lcm / npcol;
    igcd = nprow / lcmq;

    mp = numroc_(m, mb, &mrrow, &c__0, &nprow);
    mq = numroc_(m, mb, &mrcol, &c__0, &npcol);
    np = numroc_(n, nb, &mrrow, &c__0, &nprow);
    nq = numroc_(n, nb, &mrcol, &c__0, &npcol);

    i__1 = iceil_(m, mb);
    lbm = iceil_(&i__1, &lcm);
    i__1 = iceil_(n, nb);
    lbn = iceil_(&i__1, &lcm);

/*     Test for the input parameters again with local parameters */

    if (*lda < mp) {
	info = 7;
    } else if (*ldc < np) {
	info = 10;
    }
    if (info != 0) {
	goto L10;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0) {
	return 0;
    }

/*     At first, scale C with beta if beta != 0.0 & beta != 1.0 */

    tbeta = *beta;
    if (*beta != 0. && *beta != 1.) {
	i__1 = mq;
	for (j1 = 1; j1 <= i__1; ++j1) {
	    HPL_dscal( np, *beta, &c__[j1 * c_dim1 + 1], 1 );
/* L20: */
	}
	tbeta = 1.;
    }

    commtrb_1.iaz = lcmp * *mb;
    commtrb_1.jaz = lcmq * *nb;
    commtrb_1.itz = lcmp * *nb;
    commtrb_1.jtz = lcmq * *mb;

    ml = lbm * *mb;
    nl = lbn * *nb;
    ipt = (long)ml * (long)nl + 1;
    ldt = nl;
    kprow = mrrow + nprow;
    kpcol = mrcol + npcol;

/*     Initialize Parameters -- Compute the positions of subblocks */

    i__1 = npcol - 1;
    for (k1 = 0; k1 <= i__1; ++k1) {
	ncol = (kpcol - k1) % npcol;
	i__2 = lcmq - 1;
	for (j1 = 0; j1 <= i__2; ++j1) {
	    item = npcol * j1 + ncol;
	    if (item % nprow == mrrow) {
		iwork[ncol * 3 + 1] = item / nprow;
	    }
/* L30: */
	}
    }

    i__2 = lcmq - 1;
    for (j1 = 0; j1 <= i__2; ++j1) {
	item = (npcol * j1 + mrcol) % nprow;
	iwork[item * 3 + 2] = j1;
	iwork[item * 3 + 3] = j1;
	i__1 = igcd - 1;
	for (k1 = 1; k1 <= i__1; ++k1) {
	    iwork[(item + nprow - k1) % nprow * 3 + 2] = j1;
	    iwork[(item + k1) % nprow * 3 + 3] = j1;
/* L40: */
	}
    }

/*     Set parameters for efficient copying */

    lbm0 = lbm;
    lbm1 = lbm;
    lbm2 = lbm;
    lbn0 = lbn;
    lbn1 = lbn;
    lbn2 = lbn;
    mb0 = *mb;
    mb1 = *mb;
    mb2 = *mb;
    nb0 = *nb;
    nb1 = *nb;
    nb2 = *nb;

    if (nprow == npcol) {
	lbm0 = 1;
	lbn0 = 1;
	mb0 = mp;
	nb0 = nq;
    }
    if (nprow == lcm) {
	lbm1 = 1;
	lbn2 = 1;
	mb1 = mp;
	nb2 = np;
    }
    if (npcol == lcm) {
	lbn1 = 1;
	lbm2 = 1;
	nb1 = nq;
	mb2 = mq;
    }

/*     For each K2 loop (rowwise), Copy A' to WORK & Send it to KTPROC */
/*                                 then, Receive WORK and Copy WORK to C */

    kmod1 = (nprow + mrcol - mrrow) % igcd;
    kmod2 = (igcd - kmod1) % igcd;

    i__1 = lcmp - 1;
    for (k2 = 0; k2 <= i__1; ++k2) {

/*        Copy A' to WORK in the appropriate order & Send it */

	k1 = k2 * igcd + kmod1;
	mcol = (kpcol - k1) % npcol;
	kia = iwork[mcol * 3 + 1] * *mb;
	mcol = (mcol + *imcol) % npcol;
	ncol = (mrcol + k2 * igcd + kmod2) % npcol;
	kic = iwork[ncol * 3 + 1] * *nb;
	ncol = (ncol + *imcol) % npcol;

	i__2 = lcmq - 1;
	for (j1 = 0; j1 <= i__2; ++j1) {
	    kja = iwork[(mrrow + igcd * j1) % nprow * 3 + 2] * *nb;

	    if (myrow == (myrow + igcd * j1 + kmod1) % nprow && mycol == mcol)
		     {
		kjc = iwork[(kprow - igcd * j1) % nprow * 3 + 3] * *mb;
		i__3 = mp - kia;
		i__4 = nq - kja;
		dtr2mx_(&a[kia + 1 + (kja + 1) * a_dim1], lda, &tbeta, &c__[
			kic + 1 + (kjc + 1) * c_dim1], ldc, &lbm0, &lbn0, &
			mb0, &nb0, &i__3, &i__4);

	    } else {
		i__3 = mp - kia;
		i__4 = nq - kja;
		dtr2bf_(&a[kia + 1 + (kja + 1) * a_dim1], lda, &work[1], &ldt,
			 &lbm1, &lbn1, &mb1, &nb1, &i__3, &i__4);

		if (nprow == npcol && *beta == 0. && *ldc == ldt) {
		    i__3 = (myrow + igcd * j1 + kmod1) % nprow;
		    i__4 = (mprow - igcd * j1 - kmod2) % nprow;
		    kjc = iwork[(kprow - igcd * j1) % nprow * 3 + 3] * *mb;
#if 0
		    Cdgesd2d(context_1.ictxt,nl,ml,&work[1],nl,i__3,mcol);
		    Cdgerv2d(context_1.ictxt,nl,ml,&c__[(kjc + 1) * c_dim1 + 1],*ldc,i__4,ncol);
#else
		    Cblacs_dSendrecv( ctxt,
                          nl, ml, &work[1], nl, i__3, mcol,
                          nl, ml, &c__[(kjc + 1) * c_dim1 + 1], *ldc, i__4, ncol );
#endif

		} else {
		    i__3 = (myrow + igcd * j1 + kmod1) % nprow;
		    i__4 = (mprow - igcd * j1 - kmod2) % nprow;
#if 0
		    Cdgesd2d(context_1.ictxt,nl,ml,&work[1],nl,i__3,mcol);
		    Cdgerv2d(context_1.ictxt,nl,ml,&work[ipt],nl, i__4,ncol);
#else
        Cblacs_dSendrecv( ctxt,
                          nl, ml, &work[1],   nl, i__3, mcol,
                          nl, ml, &work[ipt], nl, i__4, ncol );
#endif

		    kjc = iwork[(kprow - igcd * j1) % nprow * 3 + 3] * *mb;
		    i__3 = np - kic;
		    i__4 = mq - kjc;
		    dmv2mx_(&work[ipt], &ldt, &tbeta, &c__[kic + 1 + (kjc + 1)
			     * c_dim1], ldc, &lbn2, &lbm2, &nb2, &mb2, &i__3,
			    &i__4);
		}
	    }
	}
    }

    return 0;
} /* pdtrans_ */

int
main(int argc, char *argv[]) {
  int n = 55000, nb=80, row=0, col=1, *iwork, nprow, npcol;
  double beta=1.0, A[1], C[1], work[1];

  if (argc <= 1 || sscanf( argv[1], "%d", &nprow ) != 1 || nprow < 1) nprow = 2;
  if (argc <= 2 || sscanf( argv[2], "%d", &npcol ) != 1 || npcol < 1) npcol = 2;
  if (argc <= 3 || sscanf( argv[3], "%d", &n ) != 1  || n < 1) n = 55000;
  if (argc <= 4 || sscanf( argv[4], "%d", &nb ) != 1 || nb < 1) nb = 80;

  Cblacs_nprow = nprow;
  Cblacs_npcol = npcol;
  Cblacs_myrow = 0;
  Cblacs_mycol = 0;

  printf( "%dx%d %d %d\n", nprow, npcol, n, nb );

  iwork = (int *)malloc( (sizeof *iwork) * 3 * nprow * npcol );

  pdtrans( "T", &n, &n, &nb, &nb, A, &n, &beta, C, &n, &row, &col, work, iwork );

  free( iwork );

  return 0;
}
