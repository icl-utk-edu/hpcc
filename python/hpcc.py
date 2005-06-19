#! /usr/bin/env python
# -*- mode: Python; tab-width: 4; indent-tabs-mode: nil; fill-column: 79;  coding: iso-latin-1-unix -*-
#

import sys
from time import time

import numarray as N
import numarray.ufunc as U

SIZEOF_INT = 4
SIZEOF_DOUBLE = 8

class HPCC:
    def __init__(self, in_f, out_f):
        pass

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

        times = N.array(type=Float64, shape=(4, NTIMES))

        VectorSize = self.local_vector_size( 3, SIZEOF_DOUBLE, False )

        scalar = 3.0

        a = N.array(type=Float64, shape=VectorSize)
        b = N.array(type=Float64, shape=VectorSize)
        c = N.array(type=Float64, shape=VectorSize)

        a[:] = 1.0
        b[:] = 2.0
        c[:] = 0.0

        for i in range(NTIMES):
            times[0, i] = -time()
            # COPY
            a[:] = c
            times[0, i] += time()

            times[1, i] = -time()
            # SCALE: b = scalar*c
            U.multiply(scalar, c, b)
            times[1, i] += time()

            times[2, i] = -time()
            # ADD: c = a+b
            U.add(a, b, c)
            times[2, i] += time()

            times[3, i] = -time()
            # TRIAD: a = b+scalar*c
            U.add(b, scalar*c, a)
            times[3, i] += time()

def main(argv):
    mpi.Init(argv)

    hpcc = HPCC("hpccinf.txt", "hpccoutf.txt")

    mpi.Finalize()
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
