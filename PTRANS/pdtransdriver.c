/* -*- mode: C; tab-width: 2; indent-tabs-mode: nil; -*- */
/*
pdtransdriver.c
*/

#include <hpcc.h>

/* Common Block Declarations */

struct {
    int ictxt;
} context_;

#define context_1 context_

/* Table of constant values */

static int c__1 = 1;
static int c__0 = 0;

int
PTRANS(HPCC_Params *params) {
  /* calculation of passed/failed/skipped tests assumes that MPI rank 0 is 0x0 in CBLACS */
  int ktests = 0;
  int kpass = 0;
  int kfail = 0;
  int kskip = 0;

  int i__, j, m, n;
  int mb, nb, ii, mg, ng, mp, mq, np, nq;
  int mp0, mq0, np0, nq0, lda, ldc, iam, lcm;
  double eps, *mem;
  int *imem;
  long ipa, ipc, ipw, ipiw, isw;
  int nmat, *mval, ierr[1], *nval;
    int nbmat, *mbval, imcol, *nbval;
    double ctime[2], resid, resid0;
    int npcol, *npval, mycol, *nqval;
    double wtime[2];
    int imrow, nprow, myrow, iaseed;
    char *passed;
    int ngrids;
    double thresh;
    int nprocs;


/*  -- PUMMA Package routine (version 2.1) -- */
/*     Jaeyoung Choi, Oak Ridge National Laboratory. */
/*     Jack Dongarra, Univ. of Tennessee, Oak Ridge National Laboratory. */
/*     David Walker,  Oak Ridge National Laboratory. */
/*     March 26, 1995. */

/*  Purpose: Driver routine for testing the full matrix transpose. */

/*    The user should modify TOTMEM to indicate the maximum amount of memory */
/*  in bytes his system has.  Remember to leave room in memory for operating */
/*  system, the buffer of the communication package, etc ... */

/*    The constants INTGSZ and DBLESZ indicate the length in bytes on the */
/*  given platform for an integer and a double precision real. For example, */
/*  on a system with 8 MB of memory, the parameters we use are */
/*  TOTMEM=6200000 (leaving 1.8 MB for OS, code, communication buffer, etc). */
/*  However, for PVM, we usually set TOTMEM = 2000000. */
/*  The length of a double precision real is 8, and an integer takes up 4 bytes. */

/*    Some playing around to discover what the maximum value you can set */
/*  TOTMEM to may be required. All arrays used by the factorization and */
/*  check are allocated out of the array called MEM. The integer IPA, */
/*  for example, indicates the element of MEM that the answer vector(s) */
/*  A begin(s) on. */

  FILE *outFile;
  double curGBs, cpuGBs, *GBs;
  int AllocSuccessful;
  int icseed = 200;
  double d_One = 1.0;
  long dMemSize, li;

  GBs = &params->PTRANSrdata.GBs;
  *GBs = curGBs = 0.0;

  Cblacs_pinfo(&iam, &nprocs);

  if (0 == iam) {
    outFile = fopen( params->outFname, "a" );
    if (! outFile) outFile = stderr;
  }

  nmat = params->PTRANSns;
  mval = params->PTRANSnval;
  nval = params->PTRANSnval;

  nbmat = params->PTRANSnbs;
  mbval = params->PTRANSnbval;
  nbval = params->PTRANSnbval;

  ngrids = params->PTRANSnpqs;
  npval = params->PTRANSpval;
  nqval = params->PTRANSqval;

  thresh = params->test.thrsh;
  eps = params->test.epsil;

  iaseed = 100;
  imrow = imcol = 0;

  /* calculate and allocate memory */
  AllocSuccessful = 0;
  MaxMem( nprocs, imrow, imcol, nmat, mval, nval, nbmat, mbval, nbval, ngrids, npval, nqval,
          &dMemSize );
  mem = NULL; imem = NULL;
  if (dMemSize > 0) {
    mem = XMALLOC( double, dMemSize );
    imem = XMALLOC( int, (3 * nprocs) );
    if (mem && imem) AllocSuccessful = 1;
  }

  MPI_Allreduce( &AllocSuccessful, ierr, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD );
  if (ierr[0] < 1) {
    if (imem) free(imem);
    if (mem) free(mem);
    if (0 == iam) fprintf( outFile, "Failed to allocate %ld doubles\n", dMemSize );
    goto mem_failure;
  }

  /* initialize working arrays; it is necessary because on some systems it will contain NaNs
   * (Not a Number) and NaTs (Not a Thing) and this makes pdmatgen() work incorrectly
   * (0.0 * NaN may cause exception) */
  for (li = 0; li < dMemSize; li++) mem[li] = 0.0;
  for (j = 0; j < 3 * nprocs; j++) imem[j] = 0;

  /* Print headings */
  if (iam == 0) {
    /* matrix sizes */
    fprintf( outFile, "M:" );
    for (j = 0; j < nmat; j++) fprintf( outFile, " %d", mval[j] ); fprintf( outFile, "\n" );
    fprintf( outFile, "N:" );
    for (j = 0; j < nmat; j++) fprintf( outFile, " %d", nval[j] ); fprintf( outFile, "\n" );
    /* block sizes */
    fprintf( outFile, "MB:" );
    for (j = 0; j < nbmat; j++) fprintf( outFile, " %d", mbval[j] ); fprintf( outFile, "\n" );
    fprintf( outFile, "NB:" );
    for (j = 0; j < nbmat; j++) fprintf( outFile, " %d", nbval[j] ); fprintf( outFile, "\n" );
    /* process grids */
    fprintf( outFile, "P:" );
    for (j = 0; j < ngrids; j++) fprintf( outFile, " %d", npval[j] ); fprintf( outFile, "\n" );
    fprintf( outFile, "Q:" );
    for (j = 0; j < ngrids; j++) fprintf( outFile, " %d", nqval[j] ); fprintf( outFile, "\n" );

    fprintf( outFile,
             "TIME   M     N    MB  NB  P   Q     TIME   CHECK   GB/s   RESID\n"
             "---- ----- ----- --- --- --- --- -------- ------ -------- -----\n" );
    fflush( outFile );
  }

/*     Loop over different process grids */

  for (j = 0; j < ngrids; ++j) {
    nprow = npval[j];
    npcol = nqval[j];

/*        Make sure grid information is correct */

    ierr[0] = 0;
    if (nprow < 1) {
      if (iam == 0) {
        fprintf( outFile, "ILLEGAL %s: %s = %d; It should be at least 1", "GRID", "nprow", nprow );
      }
     ierr[0] = 1;
    } else if (npcol < 1) {
      if (iam == 0) {
        fprintf( outFile, "ILLEGAL %s: %s = %d; It should be at least 1", "GRID", "npcol", npcol );
      }
      ierr[0] = 1;
    } else if (nprow * npcol > nprocs) {
      if (iam == 0) {
        fprintf( outFile, "ILLEGAL GRID: nprow*npcol = %d. It can be at most %d.",
                 nprow * npcol, nprocs );
      }
      ierr[0] = 1;
    }

    if (ierr[0] > 0) {
      if (iam == 0) {
        fprintf( outFile, "Bad %s parameters: going on to next test case.", "grid" );
      }
      ++kskip;
      goto L30;
    }

/*        Define process grid */

    Cblacs_get(-1, 0, &context_1.ictxt);
    Cblacs_gridinit(&context_1.ictxt, "Row-major", nprow, npcol);
    Cblacs_gridinfo(context_1.ictxt, &nprow, &npcol, &myrow, &mycol);

/*        Go to bottom of process grid loop if this case doesn't use my process */

    if (myrow >= nprow || mycol >= npcol) {
      goto L30;
    }

    for (i__ = 0; i__ < nmat; ++i__) {
      m = mval[i__];
      n = nval[i__];

/*           Make sure matrix information is correct */

      ierr[0] = 0;
      if (m < 1) {
        if (iam == 0) {
          fprintf( outFile, "ILLEGAL %s: %s = %d; It should be at least 1", "MATRIX", "M", m );
        }
        ierr[0] = 1;
      } else if (n < 1) {
        if (iam == 0) {
          fprintf( outFile, "ILLEGAL %s: %s = %d; It should be at least 1", "MATRIX", "N", n );
        }
        ierr[0] = 1;
      }

/*           Make sure no one had error */

      Cigsum2d(context_1.ictxt,"a","h",1,1,ierr,1,-1,0);

      if (ierr[0] > 0) {
        if (iam == 0) {
          fprintf( outFile, "Bad %s parameters: going on to next test case.", "matrix" );
        }
        ++kskip;
        goto L20;
      }

/*           Loop over different block sizes */

      for (ii = 1; ii <= nbmat; ++ii) {

        mb = mbval[ii - 1];
        nb = nbval[ii - 1];

/*              Make sure blocking sizes are legal */

        ierr[0] = 0;
        if (mb < 1) {
          ierr[0] = 1;
          if (iam == 0) {
            fprintf( outFile, "ILLEGAL %s: %s = %d; It should be at least 1", "MB", "MB", mb );
          }
        } else if (nb < 1) {
          ierr[0] = 1;
          if (iam == 0) {
            fprintf( outFile, "ILLEGAL %s: %s = %d; It should be at least 1", "NB", "NB", nb );
          }
        }

/*              Make sure no one had error */

        Cigsum2d(context_1.ictxt,"a","h",1,1,ierr, 1,-1,0);

        if (ierr[0] > 0) {
          if (iam == 0) {
            fprintf( outFile, "Bad %s parameters: going on to next test case.", "NB" );
          }
          ++kskip;
          goto L10;
        }

        mp = numroc_(&m, &mb, &myrow, &imrow, &nprow);
        mq = numroc_(&m, &mb, &mycol, &imcol, &npcol);
        np = numroc_(&n, &nb, &myrow, &imrow, &nprow);
        nq = numroc_(&n, &nb, &mycol, &imcol, &npcol);

        mg = iceil_(&m, &mb);
        ng = iceil_(&n, &nb);

        mp0 = iceil_(&mg, &nprow) * mb;
        mq0 = iceil_(&mg, &npcol) * mb;
        np0 = iceil_(&ng, &nprow) * nb;
        nq0 = iceil_(&ng, &npcol) * nb;

        lcm = ilcm_(&nprow, &npcol);
        ipc = 1;
        ipa = ipc + (long)np0 * (long)mq0;
        ipiw = (long)mp0 * (long)nq0 + ipa;
        ipw = ipiw;
        isw = ipw + (long)(iceil_(&mg, &lcm) << 1) * (long)mb * (long)iceil_(&ng, &lcm) * (long)nb;

/*              Make sure have enough memory to handle problem */

        if (isw > dMemSize) {
          if (iam == 0) {
            fprintf( outFile, "Unable to perform %s: need memory of at least %ld doubles\n",
                     "PTRANS", isw );
          }
          ierr[0] = 1;
        }

/*              Make sure no one had error */

        Cigsum2d(context_1.ictxt,"a","h",1,1,ierr, 1,-1,0);

        if (ierr[0] > 0) {
          if (iam == 0) {
            fprintf( outFile, "Bad %s parameters: going on to next test case.", "MEMORY" );
          }
          ++kskip;
          goto L10;
        }

/*              Generate matrix A */

        lda = MAX(1,mp);
        /* A = rand(m, n, iaseed) */
        pdmatgen(&context_1.ictxt, "N", "N", &m, &n, &mb, &nb, &mem[ipa - 1], &lda, &imrow, &imcol,
                 &iaseed, &c__0, &mp, &c__0, &nq, &myrow, &mycol, &nprow, &npcol, 0.0);
        /* C = rand(n, m, icseed) */
        pdmatgen(&context_1.ictxt, "T", "N", &n, &m, &nb, &mb, &mem[ipc - 1], &lda, &imrow, &imcol,
                 &icseed, &c__0, &np, &c__0, &mq, &myrow, &mycol, &nprow, &npcol, 0.0);

        slboot_();
        Cblacs_barrier(context_1.ictxt, "All");
        sltimer_(&c__1);

/*              Perform the matrix transpose */

        ldc = MAX(1,np);
        /* C := A' + d_One * C */
        pdtrans( "T", &m, &n, &mb, &nb, &mem[ipa - 1], &lda, &d_One, &mem[ipc - 1], &ldc, &imrow,
                 &imcol, &mem[ipw - 1], imem );

        sltimer_(&c__1);

        if (thresh > 0.0) {

/*                  Regenerate matrix A in transpose form (A') */

          lda = MAX(1,np);
          /* A = rand(n, m, icseed) */
          pdmatgen(&context_1.ictxt, "T", "N", &n, &m, &nb, &mb, &
                   mem[ipa - 1], &lda, &imrow, &imcol, &icseed, &
                   c__0, &np, &c__0, &mq, &myrow, &mycol, &nprow, &npcol, 0.0);
          /* A += rand(m, n, iaseed) */
          pdmatgen(&context_1.ictxt, "T", "N", &m, &n, &mb, &nb, &mem[ipa - 1], &lda, &imrow,
                   &imcol, &iaseed, &c__0, &mp, &c__0, &nq, &myrow, &mycol, &nprow, &npcol, 1.0);

/*                  Compare A' to C */

          pdmatcmp(&context_1.ictxt, &np, &mq, &mem[ipa - 1], &lda, &mem[ipc - 1], &ldc, &resid);
          resid0 = resid;

          resid /= eps * MAX( m, n );
          if (resid <= thresh && resid - resid == 0.0) { /* if `resid' is small and is not NaN */
            ++kpass;
            passed = "PASSED";
          } else {
            ++kfail;
            passed = "FAILED";
          }
        } else {

/*                  Don't perform the checking, only the timing operation */

          ++kpass;
          resid -= resid;
          passed = "BYPASS";
        }

/*              Gather maximum of all CPU and WALL clock timings */

        slcombine_(&context_1.ictxt, "All", ">", "W", &c__1, &c__1, wtime);
        slcombine_(&context_1.ictxt, "All", ">", "C", &c__1, &c__1, ctime);

/*              Print results */

        if (iam == 0) {

/*                  Print WALL time if machine supports it */

          if (wtime[0] > 0.0) {
            curGBs = 1e-9 / wtime[0] * m * n * sizeof(double);
            if (curGBs > *GBs) {
              *GBs = curGBs;
              params->PTRANSrdata.time = wtime[0];
              params->PTRANSrdata.residual = resid0;
              params->PTRANSrdata.n = n;
              params->PTRANSrdata.nb = nb;
              params->PTRANSrdata.nprow = nprow;
              params->PTRANSrdata.npcol = npcol;
            }
            fprintf( outFile, "WALL %5d %5d %3d %3d %3d %3d %8.2f %s %8.3f %5.2f\n",
                     m, n, mb, nb, nprow, npcol, wtime[0], passed, curGBs, resid );
          }

/*                  Print CPU time if machine supports it */

          if (ctime[0] > 0.0) {
            cpuGBs = 1e-9 / ctime[0] * m * n * sizeof(double);
            fprintf( outFile, "CPU  %5d %5d %3d %3d %3d %3d %8.2f %s %8.3f %5.2f\n",
                     m, n, mb, nb, nprow, npcol, ctime[0], passed, cpuGBs, resid );
          }
        }
      L10:
        ;
      }
  L20:
      ;
  }

  Cblacs_gridexit(context_1.ictxt);

  L30:
  ;
  }

  if (imem) free( imem );
  if (mem) free( mem );

  mem_failure:

/*     Print out ending messages and close output file */

  if (iam == 0) {
    ktests = kpass + kfail + kskip;

    fprintf( outFile, "\nFinished %4d tests, with the following results:\n", ktests );

    if (thresh > 0.0) {
      fprintf( outFile, "%5d tests completed and passed residual checks.\n", kpass );
      fprintf( outFile, "%5d tests completed and failed residual checks.\n", kfail );
    } else {
      fprintf( outFile, "%5d tests completed without checking.\n", kpass );
    }
    fprintf( outFile, "%5d tests skipped because of illegal input values.\n", kskip );


    fprintf( outFile, "\nEND OF TESTS.\n" );

    if (outFile != stdout && outFile != stderr) fclose( outFile );
  }

  Cblacs_exit(1);

  /* if at least one test failed or was skipped then it's a total failure */
  MPI_Reduce( &kfail, &ktests, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD );
  if (ktests) params->Failure = 1;
  MPI_Reduce( &kskip, &ktests, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD );
  if (ktests) params->Failure = 1;

  return 0;
}
