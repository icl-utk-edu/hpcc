/*
C
C     FFTE: A FAST FOURIER TRANSFORM PACKAGE
C
C     (C) COPYRIGHT SOFTWARE, 2000-2004, ALL RIGHTS RESERVED
C                BY
C         DAISUKE TAKAHASHI
C         INSTITURE OF INFORMATION SCIENCES AND ELECTRONICS
C         UNIVERSITY OF TSUKUBA
C         1-1-1 TENNODAI, TSUKUBA, IBARAKI 305-8573, JAPAN
C         E-MAIL: daisuke@is.tsukuba.ac.jp
C
C
C     GLOBAL TRANSPOSE ROUTINE (MPI VERSION)
C
C     FORTRAN77 SOURCE PROGRAM
C
C     WRITTEN BY DAISUKE TAKAHASHI
C
*/

#include "hpccfft.h"
#include "wrapmpifftw.h"

int
pztrans_(fftw_complex *a, fftw_complex *b, s64Int_t *nn, int *npu, hpcc_fftw_mpi_plan p) {
  int i, nn2;

  nn2 = *nn / *npu;

  if (1 == *npu)
    for (i = 0; i < nn2; i++) b[i] = a[i];
  else
    MPI_Alltoall( a, nn2, p->cmplx, b, nn2, p->cmplx, p->comm );

  return 0;
}
