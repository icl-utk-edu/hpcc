/* -*- mode: C; tab-width: 2; indent-tabs-mode: nil; -*- */

/* 
 * This code has been contributed by the DARPA HPCS program.  Contact 
 * David Koester (MITRE) or Bob Lucas (USC/ISI) if you have questions.
 *
 * Information on the rules to run RandomAccess or modify it to optimize
 * performance can be found at http://icl.cs.utk.edu/hpcc/
 *
 * The benchmark is successfully completed if at least 99% of the table
 * is correctly built; it *is* OK, for example, to have a few errors
 * enter into the table due to rare race conditions between processors
 * updating the same memory location.
 *
 * It is expected that in a multi-processor version, one will start the
 * random number generator at a location along its cycle proportional t
 * the local table size.  Each processor then steps the generator through
 * its section of the cycle.  (The starts routine provided can be used to 
 * start the random number generator at any specified point.)
 *
 * The benchmark objective is to measure interprocessor and local memory 
 * bandwidth.  It is ok to bucket sort locally and use a message passing 
 * protocol between processors if this is faster than a global shared 
 * memory approach.  More sophisticated optimizations which attempt to 
 * reorder the loop so that all updates are local are not considered "in 
 * the spirit" of the benchmark.
 *
 * This code has been designed to run on any number of processors.  However, 
 * due to the reduced complexity of operations when there are a power of two
 * number of processors - the processor number can be determined from the
 * proper bits in the random stream versus requiring an integer division - 
 * the two cases are examined separately.
 *
 */

#include <assert.h>
#include <hpcc.h>

#include "RandomAccess.h"

/* Allocate main table (in global memory) */
/*unsigned long *Table;*/
u64Int *Table;

#define TRUE 1
#define FALSE 0
#define DONE 0

#define FINISHED_TAG 1
#define UPDATE_TAG   2
#define USE_NONBLOCKING_SEND 1
              
void
Sum64(void *invec, void *inoutvec, int *len, MPI_Datatype *datatype) {
  int i, n = *len; s64Int *invec64 = invec, *inoutvec64 = inoutvec;
  for (i = n; i; i--, invec64++, inoutvec64++) *inoutvec64 += *invec64;
}

void
AnyNodesMPIRandomAccessUpdate(u64Int logTableSize,
                              u64Int TableSize,
                              u64Int LocalTableSize,
                              u64Int MinLocalTableSize,
                              u64Int GlobalStartMyProc,
                              u64Int Top,
                              int logNumProcs,
                              int NumProcs,
                              int Remainder,
                              int MyProc,
                              MPI_Datatype INT64_DT)
{
  s64Int i, k, m;
  int proc_count;

  s64Int SendCnt;
  u64Int Ran, RanTemp;
  s64Int lcRecvDone =  FALSE;
  s64Int WhichPe;
  u64Int GlobalOffset, LocalOffset;
  int NumberReceiving = NumProcs - 1;
  u64Int inmsg;
  MPI_Request inreq, outreq = MPI_REQUEST_NULL;
  MPI_Status status;
  int have_done;

  /* Initialize main table */
  for (i=0; i<LocalTableSize; i++)
    Table[i] = i + GlobalStartMyProc;

  /* Perform updates to main table.  The scalar equivalent is:
   *
   *     u64Int Ran;
   *     Ran = 1;
   *     for (i=0; i<NUPDATE; i++) {
   *       Ran = (Ran << 1) ^ (((s64Int) Ran < 0) ? POLY : 0);
   *       Table[Ran & (TABSIZE-1)] ^= Ran;
   *     }
   */
   SendCnt = (4 * LocalTableSize);
   Ran = starts (4 * GlobalStartMyProc);

   MPI_Irecv(&inmsg, 1, INT64_DT, MPI_ANY_SOURCE, MPI_ANY_TAG,
             MPI_COMM_WORLD, &inreq);

   for (i = 0 ; i < SendCnt ;) {
     if (outreq == MPI_REQUEST_NULL) {
       Ran = (Ran << 1) ^ ((s64Int) Ran < ZERO64B ? POLY : ZERO64B);
       GlobalOffset = Ran & (TableSize-1);
       if ( GlobalOffset < Top) 
         WhichPe = ( GlobalOffset / (MinLocalTableSize + 1) );
       else
         WhichPe = ( (GlobalOffset - Remainder) / MinLocalTableSize );

       if (WhichPe == MyProc) {
         LocalOffset = (Ran & (TableSize - 1)) - GlobalStartMyProc;
         Table[LocalOffset] ^= Ran;
         i++;
         continue;
       } else {
#if USE_NONBLOCKING_SEND
         MPI_Isend(&Ran, 1, INT64_DT, (int) WhichPe, UPDATE_TAG, 
                   MPI_COMM_WORLD, &outreq);
#else
         MPI_Send(&Ran, 1, INT64_DT, (int) WhichPe, UPDATE_TAG, 
                  MPI_COMM_WORLD);
         i++;
#endif
       }
     }
#if USE_NONBLOCKING_SEND
     /* no way we can get here unless outreq is an active request */
     MPI_Test(&outreq, &have_done, MPI_STATUS_IGNORE);
     if (have_done) {
       outreq = MPI_REQUEST_NULL;
       i++;
     }
#endif

     /* receive messages */
     do {
       MPI_Test(&inreq, &have_done, &status);
       if (have_done) {
         if (status.MPI_TAG == UPDATE_TAG) {
           LocalOffset = (inmsg & (TableSize - 1)) - GlobalStartMyProc;
           Table[LocalOffset] ^= inmsg;
         } else if (status.MPI_TAG == FINISHED_TAG) {
           /* we got a done message.  Thanks for playing... */
           NumberReceiving--;
         } else {
           abort();
         }
         MPI_Irecv(&inmsg, 1, INT64_DT, MPI_ANY_SOURCE, MPI_ANY_TAG,
                   MPI_COMM_WORLD, &inreq);
       }
     } while (have_done && NumberReceiving > 0);
   }

   /* send our done messages */
   for (proc_count = 0 ; proc_count < NumProcs ; ++proc_count) {
     if (proc_count == MyProc) continue;
     /* send garbage - who cares, no one will look at it */
     MPI_Send(&Ran, 1, INT64_DT, proc_count, FINISHED_TAG, 
              MPI_COMM_WORLD);
   }

   /* Finish everyone else up... */
   while (NumberReceiving > 0) {
     MPI_Wait(&inreq, &status);
     if (status.MPI_TAG == UPDATE_TAG) {
       LocalOffset = (inmsg & (TableSize - 1)) - GlobalStartMyProc;
       Table[LocalOffset] ^= inmsg;
     } else if (status.MPI_TAG == FINISHED_TAG) {
       /* we got a done message.  Thanks for playing... */
       NumberReceiving--;
     } else {
       abort();
     }
     MPI_Irecv(&inmsg, 1, INT64_DT, MPI_ANY_SOURCE, MPI_ANY_TAG,
               MPI_COMM_WORLD, &inreq);
   }

#if 0
   {
     char foo[200];
     snprintf(foo, 200, "mpi_table.%d", MyProc);
     FILE *fp;
     fp = fopen(foo, "w");
     for (proc_count = 0 ; proc_count < LocalTableSize ; ++proc_count) {
       fprintf(fp, "%lld\n", Table[proc_count]);
     }
     fclose(fp);
   }
#endif

   /* Be nice and clean up after ourselves */
   MPI_Cancel(&inreq);

/* end multiprocessor code */
}


void
Power2NodesMPIRandomAccessUpdate(u64Int logTableSize,
                                 u64Int TableSize,
                                 u64Int LocalTableSize,
                                 u64Int MinLocalTableSize,
                                 u64Int GlobalStartMyProc,
                                 u64Int Top,
                                 int logNumProcs,
                                 int NumProcs,
                                 int Remainder,
                                 int MyProc,
                                 MPI_Datatype INT64_DT)
{
  s64Int i, k, m;
  int proc_count;

  s64Int SendCnt;
  u64Int Ran, RanTemp;
  s64Int lcRecvDone =  FALSE;
  s64Int WhichPe;
  u64Int GlobalOffset, LocalOffset;
  int NumberReceiving = NumProcs - 1;
  u64Int inmsg;
  MPI_Request inreq, outreq = MPI_REQUEST_NULL;
  MPI_Status status;
  int have_done;

  /* Initialize main table */
  for (i=0; i<LocalTableSize; i++)
    Table[i] = i + GlobalStartMyProc;

  /* Perform updates to main table.  The scalar equivalent is:
   *
   *     u64Int Ran;
   *     Ran = 1;
   *     for (i=0; i<NUPDATE; i++) {
   *       Ran = (Ran << 1) ^ (((s64Int) Ran < 0) ? POLY : 0);
   *       Table[Ran & (TABSIZE-1)] ^= Ran;
   *     }
   */
   SendCnt = (4 * LocalTableSize);
   Ran = starts (4 * GlobalStartMyProc);

   MPI_Irecv(&inmsg, 1, INT64_DT, MPI_ANY_SOURCE, MPI_ANY_TAG,
             MPI_COMM_WORLD, &inreq);

   for (i = 0 ; i < SendCnt ;) {
     if (outreq == MPI_REQUEST_NULL) {
       Ran = (Ran << 1) ^ ((s64Int) Ran < ZERO64B ? POLY : ZERO64B);
       WhichPe = (Ran >> (logTableSize - logNumProcs)) & (NumProcs - 1);

       if (WhichPe == MyProc) {
         LocalOffset = (Ran & (TableSize - 1)) - GlobalStartMyProc;
         Table[LocalOffset] ^= Ran;
         i++;
         continue;
       } else {
#if USE_NONBLOCKING_SEND
         MPI_Isend(&Ran, 1, INT64_DT, (int) WhichPe, UPDATE_TAG, 
                   MPI_COMM_WORLD, &outreq);
#else
         MPI_Send(&Ran, 1, INT64_DT, (int) WhichPe, UPDATE_TAG, 
                  MPI_COMM_WORLD);
         i++;
#endif
       }
     }
#if USE_NONBLOCKING_SEND
     /* no way we can get here unless outreq is an active request */
     MPI_Test(&outreq, &have_done, MPI_STATUS_IGNORE);
     if (have_done) {
       outreq = MPI_REQUEST_NULL;
       i++;
     }
#endif

     /* receive messages */
     do {
       MPI_Test(&inreq, &have_done, &status);
       if (have_done) {
         if (status.MPI_TAG == UPDATE_TAG) {
           Table[(inmsg & ((1<<(logTableSize - logNumProcs))-1))] ^= inmsg;
         } else if (status.MPI_TAG == FINISHED_TAG) {
           /* we got a done message.  Thanks for playing... */
           NumberReceiving--;
         } else {
           abort();
         }
         MPI_Irecv(&inmsg, 1, INT64_DT, MPI_ANY_SOURCE, MPI_ANY_TAG,
                   MPI_COMM_WORLD, &inreq);
       }
     } while (have_done && NumberReceiving > 0);
   }

   /* send our done messages */
   for (proc_count = 0 ; proc_count < NumProcs ; ++proc_count) {
     if (proc_count == MyProc) continue;
     /* send garbage - who cares, no one will look at it */
     MPI_Send(&Ran, 1, INT64_DT, proc_count, FINISHED_TAG, 
              MPI_COMM_WORLD);
   }

   /* Finish everyone else up... */
   while (NumberReceiving > 0) {
     MPI_Wait(&inreq, &status);
     if (status.MPI_TAG == UPDATE_TAG) {
       Table[(inmsg & ((1<<(logTableSize - logNumProcs))-1))] ^= inmsg;
     } else if (status.MPI_TAG == FINISHED_TAG) {
       /* we got a done message.  Thanks for playing... */
       NumberReceiving--;
     } else {
       abort();
     }
     MPI_Irecv(&inmsg, 1, INT64_DT, MPI_ANY_SOURCE, MPI_ANY_TAG,
               MPI_COMM_WORLD, &inreq);
   }

#if 0
   {
     char foo[200];
     snprintf(foo, 200, "mpi_table.%d", MyProc);
     FILE *fp;
     fp = fopen(foo, "w");
     for (proc_count = 0 ; proc_count < LocalTableSize ; ++proc_count) {
       fprintf(fp, "%lld\n", Table[proc_count]);
     }
     fclose(fp);
   }
#endif

   /* Be nice and clean up after ourselves */
   MPI_Cancel(&inreq);

/* end multiprocessor code */
}


int
MPIRandomAccess(HPCC_Params *params) {
  s64Int i;
  s64Int SendCnt,GlbSendCnt;

  int NumProcs, logNumProcs, MyProc;
  s64Int WhichPe;
  u64Int GlobalOffset, GlobalStartMyProc, LocalOffset;

  u64Int logTableSize, TableSize;  
  u64Int RanTemp;               /* Current random number */ 
  double CPUTime;               /* CPU  time to update table */
  double RealTime;              /* Real time to update table */
  double TotalMem;
  int sAbort, rAbort;
  int PowerofTwo;

  u64Int Nupdate;           /* Number of updates to table (suggested: 4x number of table entries) */
  int Remainder;            /* Number of processors with (LocalTableSize + 1) entries */
  u64Int Top;               /* Number of table entries in top of Table */
  u64Int LocalTableSize;    /* Local table width */
  u64Int MinLocalTableSize; /* Integer ratio TableSie/NumProcs */

  FILE *outFile;
  MPI_Op sum64;
  double *GUPs;

  MPI_Datatype INT64_DT;

#ifdef LONG_IS_64BITS
  INT64_DT = MPI_LONG;
#else
  INT64_DT = MPI_LONG_LONG_INT;
#endif

  GUPs = &params->MPIGUPs;

  MPI_Comm_size( MPI_COMM_WORLD, &NumProcs );
  MPI_Comm_rank( MPI_COMM_WORLD, &MyProc );

  if (0 == MyProc) {
    outFile = fopen( params->outFname, "a" );
    if (! outFile) outFile = stderr;
  }

  sAbort = 0; if (sizeof(u64Int) < 8) sAbort = 1;

  MPI_Allreduce( &sAbort, &rAbort, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
  if (rAbort > 0) {
    if (MyProc == 0) fprintf( outFile, "No 64-bit integer type available.\n" );
    goto comp_end;
  }

  TotalMem = params->HPLMaxProcMem; /* max single node memory */
  /*  TotalMem *= NumProcs; */             /* max memory in NumProcs nodes */
  TotalMem /= sizeof(u64Int);

  /* calculate TableSize --- the size of update array (must be a power of 2) */
  for (TotalMem *= 0.5, logTableSize = 0, TableSize = 1;
       TotalMem >= 1.0;
       TotalMem *= 0.5, logTableSize++, TableSize <<= 1)
    ; /* EMPTY */


  /* determine whether the number of processors is a power of 2 */

  for (i = 1, logNumProcs = 0; ; logNumProcs++, i <<= 1) {
    if (i == NumProcs) { 
      PowerofTwo = TRUE;
      Remainder = 0;
      Top = 0;
      MinLocalTableSize = (TableSize / NumProcs);
      LocalTableSize = MinLocalTableSize;
      GlobalStartMyProc = (MinLocalTableSize * MyProc);
      break;

    /* number of processes is not a power 2 (too many shifts may introduce negative values or 0) */
    } 
    else if (i > NumProcs || i <= 0) {
      PowerofTwo = FALSE;
/* Minimum local table size --- some processors have an additional entry */
      MinLocalTableSize = (TableSize / NumProcs);
/* Number of processors with (LocalTableSize + 1) entries */
      Remainder = TableSize  - (MinLocalTableSize * NumProcs);
/* Number of table entries in top of Table */
      Top = (MinLocalTableSize + 1) * Remainder;
/* Local table size */
      if (MyProc < Remainder) {
          LocalTableSize = (MinLocalTableSize + 1);
          GlobalStartMyProc = ( (MinLocalTableSize + 1) * MyProc);
        }
        else {
          LocalTableSize = MinLocalTableSize;
          GlobalStartMyProc = ( (MinLocalTableSize * MyProc) + Remainder );
        }
      break;

    } /* end else if */
  } /* end for i */


  Table = XMALLOC( u64Int, LocalTableSize);
  sAbort = 0; if (! Table) sAbort = 1;

  MPI_Allreduce( &sAbort, &rAbort, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
  if (rAbort > 0) {
    if (MyProc == 0) fprintf(outFile, "Failed to allocate memory for the main table.\n");
    goto failed_table;
  }

  /* Number of global updates to table - 4x number of table entries */
  Nupdate = 4 * TableSize;

  if (MyProc == 0) {
    fprintf( outFile, "Running on %d processors%s\n", NumProcs, PowerofTwo ? " (PowerofTwo)" : "");
    fprintf( outFile, "Total Main table size = 2^" FSTR64 " = " FSTR64 " words\n",
             logTableSize, TableSize );
    if (PowerofTwo)
        fprintf( outFile, "PE Main table size = 2^" FSTR64 " = " FSTR64 " words/PE\n",
                 (logTableSize - logNumProcs), TableSize/NumProcs );
      else
        fprintf( outFile, "PE Main table size = (2^" FSTR64 ")/%d  = " FSTR64 " words/PE MAX\n",
                 logTableSize, NumProcs, LocalTableSize);
   
    fprintf( outFile, "Number of updates = " FSTR64 "\n", Nupdate);
  }

  MPI_Barrier( MPI_COMM_WORLD );

  /* Begin timing here */
  CPUTime = -CPUSEC();
  RealTime = -RTSEC();

  if (PowerofTwo) {
    Power2NodesMPIRandomAccessUpdate(logTableSize, TableSize, LocalTableSize, 
                                     MinLocalTableSize, GlobalStartMyProc, Top, 
                                     logNumProcs, NumProcs, Remainder, 
                                     MyProc, INT64_DT);
  } else {
    AnyNodesMPIRandomAccessUpdate(logTableSize, TableSize, LocalTableSize, 
                                  MinLocalTableSize, GlobalStartMyProc, Top, 
                                  logNumProcs, NumProcs, Remainder, 
                                  MyProc, INT64_DT);
  }

  /* End timed section */
  CPUTime += CPUSEC();
  RealTime += RTSEC();

  /* Print timing results */
  if (MyProc == 0){
    *GUPs = 1e-9*Nupdate / RealTime;
    fprintf( outFile, "CPU time used = %.6f seconds\n", CPUTime );
    fprintf( outFile, "Real time used = %.6f seconds\n", RealTime );
    fprintf( outFile, "%.9f Billion(10^9) Updates    per second [GUP/s]\n", 
             *GUPs );
    fprintf( outFile, "%.9f Billion(10^9) Updates/PE per second [GUP/s]\n",
             *GUPs / NumProcs );
    *GUPs /= NumProcs;
  }
  /* distribute result to all nodes */
  MPI_Bcast( GUPs, 1, MPI_INT, 0, MPI_COMM_WORLD ); 

  /* Verification of results (in serial or "safe" mode; optional) */
  CPUTime  = -CPUSEC();
  RealTime = -RTSEC();
  RanTemp = 0x1;
  for (i=0; i<Nupdate; i++) {
    RanTemp = (RanTemp << 1) ^ (((s64Int) RanTemp < ZERO64B) ? POLY : ZERO64B);
    if (PowerofTwo) {
      WhichPe = (RanTemp >> (logTableSize - logNumProcs)) & (NumProcs - 1);
    }
    else {
      GlobalOffset = RanTemp & (TableSize-1);
      if ( GlobalOffset < Top)
        WhichPe = ( GlobalOffset / (MinLocalTableSize + 1) );
      else
        WhichPe = ( (GlobalOffset - Remainder) / MinLocalTableSize );
    }

    if(WhichPe == MyProc) {
      if (PowerofTwo) {
        Table[(RanTemp & ((1<<(logTableSize - logNumProcs))-1))]
           ^= RanTemp;
        }
        else {
          GlobalOffset = RanTemp & (TableSize - 1);
          LocalOffset = (GlobalOffset - GlobalStartMyProc);
          Table[LocalOffset] ^= RanTemp;
        }
      }

  }

  SendCnt = 0;
  for (i=0; i<TableSize/NumProcs; i++)
    if (Table[i] != i + GlobalStartMyProc)
      SendCnt++;

#ifdef LONG_IS_64BITS
  MPI_Reduce( &SendCnt, &GlbSendCnt, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD );
#else
  /* MPI 1.1 standard (obsolete at this point) doesn't define MPI_SUM
    to work on `long long':
    http://www.mpi-forum.org/docs/mpi-11-html/node78.html and
    therefore LAM 6.5.6 chooses not to implement it (even though there
    is code for it in LAM and for other reductions work OK,
    e.g. MPI_MAX). MPICH 1.2.5 doesn't complain about MPI_SUM but it
    doesn't have MPI_UNSIGNED_LONG_LONG (but has MPI_LONG_LONG_INT):
    http://www.mpi-forum.org/docs/mpi-20-html/node84.htm So I need to
    create a trivial summation operation. */
  MPI_Op_create( Sum64, 1, &sum64 );
  MPI_Reduce( &SendCnt, &GlbSendCnt, 1, INT64_DT, sum64, 0, MPI_COMM_WORLD );
  MPI_Op_free( &sum64 );
#endif

  CPUTime += CPUSEC();
  RealTime += RTSEC();
  if(MyProc == 0){
    fprintf( outFile, "Check  CPU time used = %.6f seconds\n", CPUTime);
    fprintf( outFile, "Check Real time used = %.6f seconds\n", RealTime);
    fprintf( outFile, "Found " FSTR64 " errors in " FSTR64 " locations (%s).\n",
	    GlbSendCnt, TableSize, (GlbSendCnt <= 0.01*TableSize) ?
	    "passed" : "failed");
    if (GlbSendCnt > 0.01*TableSize) params->Failure = 1;
  }

  /* Deallocate memory (in reverse order of allocation which should
     help fragmentation) */
  free( Table );

  failed_table:
  comp_end:

  if (0 == MyProc) if (outFile != stderr) fclose( outFile );

  MPI_Barrier( MPI_COMM_WORLD );

  return 0;
}
