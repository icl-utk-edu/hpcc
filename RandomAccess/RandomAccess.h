/* -*- mode: C; tab-width: 2; indent-tabs-mode: nil; -*- */

#if 0
static FILE *LogFile;
static void LogBegin(int p) {char fname[100]; sprintf(fname, "%d.log", p); LogFile = fopen(fname, "a");}
static void LogEnd() {if(LogFile) fclose(LogFile);}
#define DLOG(i,v) do{if(LogFile){fprintf(LogFile,__FILE__ "(%d)@%d:" #v "=%g\n",__LINE__,i,(double)(v));fflush(LogFile);}}while(0)
#endif

/* Types used by program (should be 64 bits) */
#ifdef LONG_IS_64BITS
typedef unsigned long u64Int;
typedef long s64Int;
#define FSTR64 "%ld"
#define ZERO64B 0L
#else
typedef unsigned long long u64Int;
typedef long long s64Int;
#define FSTR64 "%lld"
#define ZERO64B 0LL
#endif

/* Random number generator */
#ifdef LONG_IS_64BITS
#define POLY 0x0000000000000007UL
#define PERIOD 1317624576693539401L
#else
#define POLY 0x0000000000000007ULL
#define PERIOD 1317624576693539401LL
#endif

/* Macros for timing */
#define CPUSEC() ((double)clock()/CLOCKS_PER_SEC)
#define RTSEC() (MPI_Wtime())

extern u64Int starts (s64Int);
