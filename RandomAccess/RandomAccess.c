/* -*- mode: C; tab-width: 2; indent-tabs-mode: nil; -*- 
 *
 * See RandomAccess.h for a comprehensive description of this test and
 * its goals.
 *
 * This file contains the computational core of the single cpu version
 * of GUPS.  The inner loop should easily be vectorized by compilers
 * with such support.
 *
 * This core is used by both the single_cpu and star_single_cpu tests.
 */

#include <hpcc.h>
#include "RandomAccess.h"

/* Number of updates to table (suggested: 4x number of table entries) */
#define NUPDATE (4 * TableSize)

/* Allocate main table (in global memory) */
u64Int *Table;

void
RandomAccessUpdate(u64Int TableSize) {
  s64Int i;
  u64Int ran[128];              /* Current random numbers */
  int j;

  /* Initialize main table */
  for (i=0; i<TableSize; i++) Table[i] = i;

  /* Perform updates to main table.  The scalar equivalent is:
   *
   *     u64Int ran;
   *     ran = 1;
   *     for (i=0; i<NUPDATE; i++) {
   *       ran = (ran << 1) ^ (((s64Int) ran < 0) ? POLY : 0);
   *       table[ran & (TableSize-1)] ^= ran;
   *     }
   */
  for (j=0; j<128; j++)
    ran[j] = starts ((NUPDATE/128) * j);

  for (i=0; i<NUPDATE/128; i++) {
/* #pragma ivdep */
    for (j=0; j<128; j++) {
      ran[j] = (ran[j] << 1) ^ ((s64Int) ran[j] < 0 ? POLY : 0);
      Table[ran[j] & (TableSize-1)] ^= ran[j];
    }
  }

#if 0
  {
    char foo[200];
    snprintf(foo, 200, "table");
    FILE *fp;
    fp = fopen(foo, "w");
    for (j = 0 ; j < TableSize ; ++j) {
      fprintf(fp, "%lld\n", Table[j]);
    }
    fclose(fp);
  }
#endif
}

int
RandomAccess(HPCC_Params *params, int doIO, double *GUPs, int *failure) {
  s64Int i;
  u64Int temp;
  double cputime;               /* CPU time to update table */
  double realtime;              /* Real time to update table */
  double totalMem;
  u64Int logTableSize, TableSize;
  FILE *outFile;

  if (doIO) {
    outFile = fopen( params->outFname, "a" );
    if (! outFile) {
      outFile = stderr;
      fprintf( outFile, "Cannot open output file.\n" );
      return 1;
    }
  }

  /* calculate local memory per node for the update table */
  totalMem = params->HPLMaxProcMem;
  totalMem /= sizeof(u64Int);

  /* calculate the size of update array (must be a power of 2) */
  for (totalMem *= 0.5, logTableSize = 0, TableSize = 1;
       totalMem >= 1.0;
       totalMem *= 0.5, logTableSize++, TableSize <<= 1)
    ; /* EMPTY */

  Table = XMALLOC( u64Int, TableSize );
  if (! Table) {
    if (doIO) {
      fprintf( outFile, "Failed to allocate memory for the update table (" FSTR64 ").\n", TableSize);
      fclose( outFile );
    }
    return 1;
  }

  /* Print parameters for run */
  if (doIO) {
  fprintf( outFile, "Main table size   = 2^" FSTR64 " = " FSTR64 " words\n", logTableSize,TableSize);
  fprintf( outFile, "Number of updates = " FSTR64 "\n", NUPDATE);
  }

  /* Begin timing here */
  cputime = -CPUSEC();
  realtime = -RTSEC();

  RandomAccessUpdate( TableSize );

  /* End timed section */
  cputime += CPUSEC();
  realtime += RTSEC();

  /* make sure no division by zero */
  *GUPs = (realtime > 0.0 ? 1.0 / realtime : -1.0);
  *GUPs *= 1e-9*NUPDATE;
  /* Print timing results */
  if (doIO) {
  fprintf( outFile, "CPU time used  = %.6f seconds\n", cputime);
  fprintf( outFile, "Real time used = %.6f seconds\n", realtime);
  fprintf( outFile, "%.9f Billion(10^9) Updates    per second [GUP/s]\n", *GUPs );
  }

  /* Verification of results (in serial or "safe" mode; optional) */
#if 1
  temp = 0x1;
  for (i=0; i<NUPDATE; i++) {
    temp = (temp << 1) ^ (((s64Int) temp < 0) ? POLY : 0);
    Table[temp & (TableSize-1)] ^= temp;
  }

  temp = 0;
  for (i=0; i<TableSize; i++)
    if (Table[i] != i)
      temp++;

  if (doIO) {
  fprintf( outFile, "Found " FSTR64 " errors in " FSTR64 " locations (%s).\n",
	  temp, TableSize, (temp <= 0.01*TableSize) ? "passed" : "failed");
  }
  if (temp <= 0.01*TableSize) *failure = 0;
  else *failure = 1;
#endif

  free( Table );

  if (doIO) {
    fflush( outFile );
    fclose( outFile );
  }

  return 0;
}

/* Utility routine to start random number generator at Nth step */
u64Int
starts(s64Int n)
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
