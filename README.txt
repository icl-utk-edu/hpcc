DARPA/DOE HPC Challange Benchmark

Piotr Luszczek

   UTK, ICL
   University of Tennessee Knoxville, Innovative Computing Laboratory
     _________________________________________________________________

   Table of Contents

   [1]Introduction
   [2]Compiling
   [3]Configuration
   [4]Running

Introduction

   This is a suite of benchmarks that measure performance of CPU, memory
   subsytem and the interconnect. For details refer to the [5]HPC
   Challenge web site.

   In essence, HPC Challange consists of a number of subbenchmarks each
   of which tests different aspect of the system.

   If you familiar with the HPL benchmark code (see the [6]HPL web site)
   then you can reuse the configuration file (input for make(1)) and the
   input file that you already have for HPL. The HPC Challange benchmark
   includes HPL and uses its configuration and input files with only
   slight modifications.

Compiling

   The first step is to create a configuration file that reflects
   characteristics of your machine. The configuration file should be
   created in the hpl directory. This directory contains instructions
   (the files README and INSTALL) on how to create the configuration
   file. The directory hpl/setup contains many examples of configuration
   files. A good approach is to copy one of them to the hpl directory and
   if it doesn't work then change it. This file is reused by all the
   components of the HPC Challange suite.

   When configuration is done, a file should exist in the hpl directory
   whose name starts with Make. and ends with the name for the system
   used for tests. For example, if the name of the system is Unix, the
   file should be named Make.Unix.

   To build the benchmark executable (for the system named Unix) type:
   make arch=Unix. This command should be run in the top directory (not
   in the hpl directory). It will look in the hpl directory for the
   configuration file and use it to build the benchmark executable.

Configuration

   The HPC Challange is driven by a short input file named hpccinf.txt
   that is almost the same as the input file for HPL (customarily called
   HPL.dat). Refer to the file hpl/www/tuning.html for details about the
   input file for HPL. A sample input file is included with the HPC
   Challange distribution.

   The differences between HPL input file and HPC Challange input file
   can be summarized as follows:
     * Lines 3 and 4 are ignored. The output always goes to the file
       named hpccoutf.txt.
     * There are additional lines (starting with line 33) that may (but
       do not have to) be used to customize the HPC Challenge benchmark.
       They are described below.

   The additional lines in the HPC Challenge input file (compared to the
   HPL input file) are:
     * Lines 33 and 34 describe additional matrix sizes to be used for
       running the PTRANS benchmark (one of the components of the HPC
       Challange benchmark).
     * Lines 35 and 36 describe additional blocking factors to be used
       for running PTRANS benchmark.

   Just for completeness, here is the list of lines of the HPC
   Challange's input file with brief descriptions of their meaning:
     * Line 1: ignored
     * Line 2: ignored
     * Line 3: ignored
     * Line 4: ignored
     * Line 5: number of matrix sizes for HPL (and PTRANS)
     * Line 6: matrix sizes for HPL (and PTRANS)
     * Line 7: number of blocking factors for HPL (and PTRANS)
     * Line 8: blocking factors for HPL (and PTRANS)
     * Line 9: type of process ordering for HPL
     * Line 10: number of process grids for HPL (and PTRANS)
     * Line 11: numbers of process rows of each process grid for HPL (and
       PTRANS)
     * Line 12: numbers of process columns of each process grid for HPL
       (and PTRANS)
     * Line 13: threshold value not to be exceeded by scaled residual for
       HPL (and PTRANS)
     * Line 14: number of panel factorization methods for HPL
     * Line 15: panel factorization methods for HPL
     * Line 16: number of recursive stopping criteria for HPL
     * Line 17: recursive stopping criteria for HPL
     * Line 18: number of recursion panel counts for HPL
     * Line 19: recursion panel counts for HPL
     * Line 20: number of recursive panel factorization methods for HPL
     * Line 21: recursive panel factorization methods for HPL
     * Line 22: number of broadcast methods for HPL
     * Line 23: broadcast methods for HPL
     * Line 24: number of look-ahead depths for HPL
     * Line 25: look-ahead depths for HPL
     * Line 26: swap methods for HPL
     * Line 27: swapping threshold for HPL
     * Line 28: form of L1 for HPL
     * Line 29: form of U for HPL
     * Line 30: value that specifies whether equilibration should be used
       by HPL
     * Line 31: memory alignment for HPL
     * Line 32: ignored
     * Line 33: number of additional problem sizes for PTRANS
     * Line 34: additional problem sizes for PTRANS
     * Line 35: number of additional blocking factors for PTRANS
     * Line 36: additional blocking factors for PTRANS

Running

   It is hard to describe all the possible ways in which the HPC
   Challange benchmark could be run on various systems. An example
   command to run the benchmark could like like this: mpirun -np 4 hpcc.
   The meaning of the command's components is as follows:
     * mpirun is the command that starts execution of an MPI code.
       Depending on the system, it might also be aprun, mpiexec, mprun,
       poe, or something appropriate for your computer.
     * -np 4 is the argument that specifies that 4 MPI processes should
       be started. The number of MPI processes should be large enough to
       accomodate all the process grids specified in the hpccinf.txt
       file.
     * hpcc is the name of the HPC Challange executable to run.

   After the run, a file called hpccoutf.txt is created which contains
   results of the benchmark. This file should be uploaded through the web
   form at the HPC Challenge website.

References

   1. README.html#id2442246
   2. README.html#id2442312
   3. README.html#id2479724
   4. README.html#id2433842
   5. http://icl.cs.utk.edu/hpcc/
   6. http://www.netlib.org/benchmark/hpl/
