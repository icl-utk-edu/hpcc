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
*/

#include <hpcc.h>

#include "RandomAccess.h"

/* Log size of main table (suggested: half of global memory) */
#define LTABSIZE 25            
#define TABSIZE (1L << LTABSIZE)

/* Log size of substitution table (suggested: half of primary cache) */
#define LSTSIZE 9
#define STSIZE (1 << LSTSIZE)

/* Number of updates to table (suggested: 4x number of table entries) */
#define NUPDATE (4 * TableSize)

/* Allocate main table (in global memory) */
u64Int *Table;

void
RandomAccessUpdate(u64Int TableSize, u64Int *stable) {
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
   *       table[ran & (TableSize-1)] ^= stable[ran >> (64-LSTSIZE)];
   *     }
   */
  for (j=0; j<128; j++)
    ran[j] = starts ((NUPDATE/128) * j);
  for (i=0; i<NUPDATE/128; i++) {
/* #pragma ivdep */
    for (j=0; j<128; j++) {
      ran[j] = (ran[j] << 1) ^ ((s64Int) ran[j] < 0 ? POLY : 0);
      Table[ran[j] & (TableSize-1)] ^= stable[ran[j] >> (64-LSTSIZE)];
    }
  }
}

int
RandomAccess(HPCC_Params *params, int doIO, double *GUPs, int *failure) {
  s64Int i;
  u64Int stable[STSIZE];        /* Substitution table */
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
  fprintf( outFile, "Subst table size  = 2^%d = %d words\n", LSTSIZE, STSIZE);
  fprintf( outFile, "Number of updates = " FSTR64 "\n", NUPDATE);
  }

  /* Initialize substitution table (not time critical) */
  stable[0] = 0;
  for (i=1; i<STSIZE; i++)
    stable[i] = stable[i-1] +
#ifdef LONG_IS_64BITS
      0x0123456789abcdefL;
#else
      0x0123456789abcdefLL;
#endif

  /* Begin timing here */
  cputime = -CPUSEC();
  realtime = -RTSEC();

  RandomAccessUpdate( TableSize, stable );

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
    Table[temp & (TableSize-1)] ^= stable[temp >> (64-LSTSIZE)];
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
starts(s64Int n) {
/* s64Int i, j; */
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
