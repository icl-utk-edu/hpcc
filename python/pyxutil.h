/* -*- C -*- */

#define pyxmpi_alloc(t,n) ((MPI##t *)malloc(sizeof(MPI##t)*(n)))

typedef unsigned long long u64Int;
typedef long long s64Int;
#define ZERO64B 0LL
#define POLY 7ULL
#define PERIOD 1317624576693539401LL

#define HPCC_FALSE 0
#define HPCC_TRUE 1
#define BUCKET_SIZE 1024
#define HPCC_DONE 0
#define SLOT_CNT 1
#define FIRST_SLOT 2

#define XMALLOC(t,s) malloc(sizeof(t)*(s))

u64Int
HPCC_starts(s64Int n)
{
  int i, j;
  u64Int m2[64];
  u64Int temp, ran;

  while (n < 0) n += PERIOD;
  while (n > PERIOD) n -= PERIOD;
  if (n == 0) return 0x1;

  temp = 0x1;
  for (i=0; i<64; i++) {
    m2[i] = temp;
    temp = (temp << 1) ^ ((s64Int) temp < 0 ? POLY : 0);
    temp = (temp << 1) ^ ((s64Int) temp < 0 ? POLY : 0);
  }

  for (i=62; i>=0; i--)
    if ((n >> i) & 1)
      break;

  ran = 0x2;
  while (i > 0) {
    temp = 0;
    for (j=0; j<64; j++)
      if ((ran >> j) & 1)
        temp ^= m2[j];
    ran = temp;
    i -= 1;
    if ((n >> i) & 1)
      ran = (ran << 1) ^ ((s64Int) ran < 0 ? POLY : 0);
  }

  return ran;
}
void
HPCC_AnyNodesMPIRandomAccessCheck(u64Int logTableSize,
                             u64Int TableSize,
                             u64Int LocalTableSize,
                             u64Int MinLocalTableSize,
                             u64Int GlobalStartMyProc,
                             u64Int Top,
                             int logNumProcs,
                             int NumProcs,
                             int Remainder,
                             int MyProc,
			     u64Int ProcNumUpdates, 
                             MPI_Datatype UINT64_DT,
                             s64Int *NumErrors,
                             u64Int *HPCC_Table)
{
  u64Int Ran;
  u64Int RanTmp;
  s64Int WhichPe;
  s64Int LocalOffset;
  s64Int GlobalOffset;
  s64Int NextSlot;
  s64Int PeBucketBase;
  s64Int SendCnt;
  s64Int errors;
  int i;
  int j;
  s64Int *PeCheckDone;
  int LocalAllDone =  HPCC_FALSE;
  int sAbort, rAbort;

  u64Int *LocalBuckets;     /* buckets used in verification phase */
  u64Int *GlobalBuckets;    /* buckets used in verification phase */

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
        for (j = FIRST_SLOT; j < GlobalBuckets[PeBucketBase+SLOT_CNT]; j ++) {
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
  for (i=0; i<LocalTableSize; i++)
    if (HPCC_Table[i] != i + GlobalStartMyProc)
      errors++;
  
  *NumErrors = errors;

  printf( "%d %d %d\n", (int)errors, (int)(4 * GlobalStartMyProc), (int)HPCC_starts(4 * GlobalStartMyProc) );

  failed_localbuckets:
  failed_globalbuckets:

  
  free (PeCheckDone);  
  free(GlobalBuckets);
  free(LocalBuckets);

}
#define UINT64_DT MPI_LONG_LONG
int Cverify(int logTableSize, int TableSize, int LocalTableSize, int MinLocalTableSize, int GlobalStartMyProc, int Top, int Remainder, int ProcNumUpdates, void *Table) {
  u64Int NumErrors;
  int MyProc, NumProcs;

  MPI_Comm_rank( MPI_COMM_WORLD, &MyProc );
  MPI_Comm_size( MPI_COMM_WORLD, &NumProcs );

  HPCC_AnyNodesMPIRandomAccessCheck( logTableSize, TableSize, LocalTableSize, MinLocalTableSize, GlobalStartMyProc, Top, 0, NumProcs, Remainder, MyProc, ProcNumUpdates, UINT64_DT, &NumErrors, Table);
  return NumErrors;
}
