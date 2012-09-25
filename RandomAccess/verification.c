#include <hpcc.h>
#include "RandomAccess.h"

#include "verification.h"

/* Verification phase: local buckets to sort into */
#define BUCKET_SIZE 1024
#define SLOT_CNT 1
#define FIRST_SLOT 2

void
HPCC_Power2NodesMPIRandomAccessCheck(u64Int logTableSize,
                                u64Int TableSize,
                                s64Int LocalTableSize,
                                u64Int GlobalStartMyProc,
                                int logNumProcs,
                                int NumProcs,
                                int MyProc,
				u64Int ProcNumUpdates,
                                MPI_Datatype UINT64_DT,
                                s64Int *NumErrors) {
  u64Int Ran, RanTmp;
  s64Int NextSlot, WhichPe, PeBucketBase, SendCnt, errors, *PeCheckDone, il;
  int i, j, n;
  int LocalAllDone =  HPCC_FALSE;
  int sAbort, rAbort;

  u64Int *LocalBuckets, *GlobalBuckets; /* buckets used in verification phase */


  LocalBuckets = XMALLOC( u64Int, (NumProcs*(BUCKET_SIZE+FIRST_SLOT)));
  sAbort = 0; if (! LocalBuckets) sAbort = 1;
  MPI_Allreduce( &sAbort, &rAbort, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
  if (rAbort > 0) {
    if (MyProc == 0) fprintf(stderr, "Failed to allocate memory for local buckets.\n");
    goto failed_localbuckets;
  }
  GlobalBuckets = XMALLOC( u64Int, (NumProcs*(BUCKET_SIZE+FIRST_SLOT)));
  sAbort = 0; if (! GlobalBuckets) sAbort = 1;
  MPI_Allreduce( &sAbort, &rAbort, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
  if (rAbort > 0) {
    if (MyProc == 0) fprintf(stderr, "Failed to allocate memory for global buckets.\n");
    goto failed_globalbuckets;
  }


  SendCnt = ProcNumUpdates; /*  SendCnt = 4 * LocalTableSize; */
  Ran = HPCC_starts (4 * GlobalStartMyProc);

  PeCheckDone = XMALLOC ( s64Int, NumProcs);

  for (i=0; i<NumProcs; i++)
    PeCheckDone[i] = HPCC_FALSE;

  while(LocalAllDone == HPCC_FALSE){
    if (SendCnt > 0) {
      /* Initalize local buckets */
      for (i=0; i<NumProcs; i++){
        PeBucketBase = i * (BUCKET_SIZE+FIRST_SLOT);
        LocalBuckets[PeBucketBase+SLOT_CNT] = FIRST_SLOT;
        LocalBuckets[PeBucketBase+HPCC_DONE] = HPCC_FALSE;
      }

      /* Fill local buckets until one is full or out of data */
      NextSlot = FIRST_SLOT;
      while(NextSlot != (BUCKET_SIZE+FIRST_SLOT) && SendCnt>0 ) {
        Ran = (Ran << 1) ^ ((s64Int) Ran < ZERO64B ? POLY : ZERO64B);
        WhichPe = (Ran >> (logTableSize - logNumProcs)) & (NumProcs - 1);
        PeBucketBase = WhichPe * (BUCKET_SIZE+FIRST_SLOT);
        NextSlot = LocalBuckets[PeBucketBase+SLOT_CNT];
        LocalBuckets[PeBucketBase+NextSlot] = Ran;
        LocalBuckets[PeBucketBase+SLOT_CNT] = ++NextSlot;
        SendCnt--;
      }

      if (SendCnt == 0)
        for (i=0; i<NumProcs; i++)
          LocalBuckets[i*(BUCKET_SIZE+FIRST_SLOT)+HPCC_DONE] = HPCC_TRUE;

    } /* End of sending loop */

    MPI_Barrier(MPI_COMM_WORLD);

    LocalAllDone = HPCC_TRUE;

    /* Now move all the buckets to the appropriate pe */
    MPI_Alltoall(LocalBuckets, (BUCKET_SIZE+FIRST_SLOT), UINT64_DT,
                 GlobalBuckets, (BUCKET_SIZE+FIRST_SLOT), UINT64_DT,
                 MPI_COMM_WORLD);

    for (i = 0; i < NumProcs; i ++) {
      if(PeCheckDone[i] == HPCC_FALSE) {
        PeBucketBase = i * (BUCKET_SIZE+FIRST_SLOT);
        PeCheckDone[i] = GlobalBuckets[PeBucketBase+HPCC_DONE];
        n = (int)(GlobalBuckets[PeBucketBase+SLOT_CNT]);
        for (j = FIRST_SLOT; j < n; ++j) {
          RanTmp = GlobalBuckets[PeBucketBase+j];
          HPCC_Table[RanTmp & (LocalTableSize-1)] ^= RanTmp;
        }
        LocalAllDone &= PeCheckDone[i];
      }
    }
  }

  errors = 0;
  for (il=0; il < LocalTableSize; il++)
    if (HPCC_Table[il] != il + GlobalStartMyProc)
      errors++;

  *NumErrors = errors;

  free( PeCheckDone );

  free( GlobalBuckets );

  failed_globalbuckets:

  free( LocalBuckets );

  failed_localbuckets:
  return;
}

void
HPCC_AnyNodesMPIRandomAccessCheck(u64Int logTableSize,
                             u64Int TableSize,
                             s64Int LocalTableSize,
                             u64Int MinLocalTableSize,
                             u64Int GlobalStartMyProc,
                             u64Int Top,
                             int logNumProcs,
                             int NumProcs,
                             int Remainder,
                             int MyProc,
			     u64Int ProcNumUpdates,
                             MPI_Datatype UINT64_DT,
                             s64Int *NumErrors) {
  u64Int Ran, RanTmp;
  s64Int WhichPe, LocalOffset, NextSlot, PeBucketBase, SendCnt, errors, *PeCheckDone, il;
  u64Int GlobalOffset;
  int i, j, n;
  int LocalAllDone =  HPCC_FALSE;
  int sAbort, rAbort;

  u64Int *LocalBuckets, *GlobalBuckets; /* buckets used in verification phase */

  LocalBuckets = XMALLOC( u64Int, (NumProcs*(BUCKET_SIZE+FIRST_SLOT)));
  sAbort = 0; if (! LocalBuckets) sAbort = 1;
  MPI_Allreduce( &sAbort, &rAbort, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
  if (rAbort > 0) {
    if (MyProc == 0) fprintf(stderr, "Failed to allocate memory for local buckets.\n");
    goto failed_localbuckets;
  }
  GlobalBuckets = XMALLOC( u64Int, (NumProcs*(BUCKET_SIZE+FIRST_SLOT)));
  sAbort = 0; if (! GlobalBuckets) sAbort = 1;
  MPI_Allreduce( &sAbort, &rAbort, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
  if (rAbort > 0) {
    if (MyProc == 0) fprintf(stderr, "Failed to allocate memory for global buckets.\n");
    goto failed_globalbuckets;
  }


  SendCnt = ProcNumUpdates; /* SendCnt = 4 * LocalTableSize; */
  Ran = HPCC_starts (4 * GlobalStartMyProc);

  PeCheckDone = XMALLOC (s64Int, NumProcs);
  for (i=0; i<NumProcs; i++)
    PeCheckDone[i] = HPCC_FALSE;

  while(LocalAllDone == HPCC_FALSE){
    if (SendCnt > 0) {
      /* Initalize local buckets */
      for (i=0; i<NumProcs; i++){
        PeBucketBase = i * (BUCKET_SIZE+FIRST_SLOT);
        LocalBuckets[PeBucketBase+SLOT_CNT] = FIRST_SLOT;
        LocalBuckets[PeBucketBase+HPCC_DONE] = HPCC_FALSE;
      }

      /* Fill local buckets until one is full or out of data */
      NextSlot = FIRST_SLOT;
      while(NextSlot != (BUCKET_SIZE+FIRST_SLOT) && SendCnt>0 ) {
        Ran = (Ran << 1) ^ ((s64Int) Ran < ZERO64B ? POLY : ZERO64B);
        GlobalOffset = Ran & (TableSize-1);
        if ( GlobalOffset < Top)
          WhichPe = ( GlobalOffset / (MinLocalTableSize + 1) );
        else
          WhichPe = ( (GlobalOffset - Remainder) / MinLocalTableSize );
        PeBucketBase = WhichPe * (BUCKET_SIZE+FIRST_SLOT);
        NextSlot = LocalBuckets[PeBucketBase+SLOT_CNT];
        LocalBuckets[PeBucketBase+NextSlot] = Ran;
        LocalBuckets[PeBucketBase+SLOT_CNT] = ++NextSlot;
        SendCnt--;
      }

      if (SendCnt == 0)
        for (i=0; i<NumProcs; i++)
          LocalBuckets[i*(BUCKET_SIZE+FIRST_SLOT)+HPCC_DONE] = HPCC_TRUE;

    } /* End of sending loop */

    MPI_Barrier(MPI_COMM_WORLD);

    LocalAllDone = HPCC_TRUE;

    /* Now move all the buckets to the appropriate pe*/
    MPI_Alltoall(LocalBuckets, (BUCKET_SIZE+FIRST_SLOT), UINT64_DT,
                 GlobalBuckets, (BUCKET_SIZE+FIRST_SLOT), UINT64_DT,
                 MPI_COMM_WORLD);

    for (i = 0; i < NumProcs; i ++) {
      if(PeCheckDone[i] == HPCC_FALSE) {
        PeBucketBase = i * (BUCKET_SIZE+FIRST_SLOT);
        PeCheckDone[i] = GlobalBuckets[PeBucketBase+HPCC_DONE];
        n = (int)(GlobalBuckets[PeBucketBase+SLOT_CNT]);
        for (j = FIRST_SLOT; j < n; ++j) {
          RanTmp = GlobalBuckets[PeBucketBase+j];
          GlobalOffset = RanTmp & (TableSize - 1);
          LocalOffset = GlobalOffset - GlobalStartMyProc;
          HPCC_Table[LocalOffset] ^= RanTmp;
        }
        LocalAllDone &= PeCheckDone[i];
      }
    }

  } /* no more local data */

  errors  = 0;
  for (il=0; il < LocalTableSize; il++)
    if (HPCC_Table[il] != il + GlobalStartMyProc)
      errors++;

  *NumErrors = errors;

  free( PeCheckDone );

  free( GlobalBuckets );

  failed_globalbuckets:

  free( LocalBuckets );

  failed_localbuckets:

  return;
}
