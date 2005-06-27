#! /usr/bin/env python
# -*- mode: Python; tab-width: 4; indent-tabs-mode: nil; fill-column: 79;  coding: iso-latin-1-unix -*-
#
# HPC Challenge Benchmark
# Piotr Luszczek
# University of Tennessee Knoxville
# Innovative Computing Laboratories                                 
# (C) Copyright 2005 All Rights Reserved                       
#

import sys

import numarray as N
import numarray.ufunc as U

import mpi

SIZEOF_INT = 4
SIZEOF_DOUBLE = 8
SIZEOF_UINT64 = 8

LONG_IS_64BITS = False

class HPCC:
    def __init__(self, in_f, out_f):
        self.HPL_max_proc_mem = 1024**2
        self.comm = mpi.COMM_WORLD
        self.rank = self.comm.rank()
        self.size = self.comm.size()
        if 0 == self.rank:
            self.log_file = open(out_f, "a")

        self.results = {}

        if LONG_IS_64BITS:
            self.mpi_int64dt = mpi.LONG
        else:
            self.mpi_int64dt = mpi.LONG_LONG_INT

    def __del__(self):
        if 0 == self.rank:
            self.log_file.close()

    def log(self, s):
        if 0 == self.rank:
            self.log_file.write(s)
            self.log_file.write("\n")
            self.log_file.flush()

    def local_vector_size(self, vec_cnt, size, pow2):
        max_int_bits2 = SIZEOF_INT * 8 - 2

        flg2 = 1
        while self.HPL_max_proc_mem / size / vec_cnt >> flg2:
            flg2 += 1

        if flg2 <= max_int_bits2:
            if pow2:
                return 1 << flg2
            return self.HPL_max_proc_mem / size / vec_cnt

        return 1 << max_int_bits2

    def stream(self):
        NTIMES = 10

        times = N.array(type=N.Float64, shape=(4, NTIMES))

        VectorSize = self.local_vector_size( 3, SIZEOF_DOUBLE, False )

        scalar = 3.0

        a = N.array(type=N.Float64, shape=VectorSize)
        b = N.array(type=N.Float64, shape=VectorSize)
        c = N.array(type=N.Float64, shape=VectorSize)

        a[:] = 1.0
        b[:] = 2.0
        c[:] = 0.0

        for i in range(NTIMES):
            times[0, i] = -mpi.Wtime()
            # COPY
            a[:] = c
            times[0, i] += mpi.Wtime()

            times[1, i] = -mpi.Wtime()
            # SCALE: b = scalar*c
            U.multiply(scalar, c, b)
            times[1, i] += mpi.Wtime()

            times[2, i] = -mpi.Wtime()
            # ADD: c = a+b
            U.add(a, b, c)
            times[2, i] += mpi.Wtime()

            times[3, i] = -mpi.Wtime()
            # TRIAD: a = b+scalar*c
            U.add(b, scalar*c, a)
            times[3, i] += mpi.Wtime()

        return 1e-9 * N.array([2, 2, 3, 3]) * SIZEOF_DOUBLE / N.array(map(min, times[:, 1:])) * VectorSize

    def ep_stream(self):
        lr = self.stream()
        gr = N.array(shape=len(lr), type=N.Float64)
        self.comm.Allreduce(lr, gr, mpi.SUM)

        i = 0
        for s in ("ep_stream_copy", "ep_stream_scale", "ep_stream_add", "ep_stream_triad"):
            self.results[s] = gr[i] / self.size
            i += 1

    def mpira_any(self, logTableSize, TableSize, LocalTableSize, MinLocalTableSize, GlobalStartMyProc, Top,
                  Remainder, ProcNumUpdates):
        pass

    def mpira(self):
        total_mem = self.HPL_max_proc_mem * self.size / SIZEOF_UINT64

        TableSize = 1
        logTableSize = 0
        total_mem /= 2
        while total_mem >= 1:
            total_mem /= 2
            logTableSize += 1
            TableSize <<= 1

        MinLocalTableSize = TableSize / NumProcs # Minimum local table size - some processors have an additional entry
        Remainder = TableSize  - MinLocalTableSize * NumProcs # Number of processors with (LocalTableSize + 1) entries
        Top = (MinLocalTableSize + 1) * Remainder # Number of table entries in top of Table

        # Local table size
        if MyProc < Remainder:
            LocalTableSize = MinLocalTableSize + 1
            GlobalStartMyProc = (MinLocalTableSize + 1) * MyProc
        else:
            LocalTableSize = MinLocalTableSize
            GlobalStartMyProc = (MinLocalTableSize * MyProc) + Remainder

        self.Table = N.array(type=N.UInt64, shape=LocalTableSize)

        NumUpdates_Default = 4 * TableSize # Default number of global updates to table: 4x number of table entries

        ProcNumUpdates = 4*LocalTableSize
        NumUpdates = NumUpdates_Default

        mpi.Barrier(self.comm)

        t = -mpi.Wtime()

        self.mpira_any(logTableSize, TableSize, LocalTableSize, MinLocalTableSize, GlobalStartMyProc, Top,
                       Remainder, ProcNumUpdates)

        mpi.Barrier(self.comm)
        t += mpi.Wtime()

        del self.Table

def main(argv):
    mpi.Init(argv)

    hpcc = HPCC("hpccinf.txt", "hpccoutf.txt")

    hpcc.ep_stream()

    keys = hpcc.results.keys()
    keys.sort()
    for k in keys:
        hpcc.log( "%s=%f" % (k, hpcc.results[k]) )

    mpi.Finalize()
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
