#! /usr/bin/env python
# -*- mode: Python; tab-width: 4; indent-tabs-mode: nil; fill-column: 79;  coding: iso-latin-1-unix -*-
#

import sys

from time import time

import numarray as N
import numarray.ufunc as U

POLY, PERIOD = 7, 1317624576693539401L

def starts(n):
    n = N.array([n], N.Int64)
    m2 = N.zeros([64], N.UInt64)

    while n[0] < 0:      n += PERIOD
    while n[0] > PERIOD: n -= PERIOD
    if n[0] == 0: return 1

    temp = N.array([1], N.UInt64)
    ival = N.array([0], N.UInt64)
    for i in range(64):
        m2[i] = temp[0]
        temp = (temp << 1) ^ ((temp.astype(N.Int64)<0)[0] and POLY or 0)
        temp = (temp << 1) ^ ((temp.astype(N.Int64)<0)[0] and POLY or 0)
    for i in range(62, -1, -1):
        if ((n>>i) & 1)[0]: break

    ran = N.array([2], N.UInt64)
    while (i > 0):
        temp[0] = 0
        for j in range(64):
            if ((ran>>j) & 1)[0]:
                temp ^= m2[j:j+1]
        ran[0] = temp[0]
        i -= 1
        if ((n>>i) & 1)[0]:
            ran = (ran << 1) ^ ((ran.astype(N.Int64)<0)[0] and POLY or 0)
    return ran[0]

def RandomAccess(nupdate, n, table):
    table[:] = N.arange(n, type=N.UInt64)
    temp = N.array([1], N.UInt64) << 63

    ran = N.array(type=N.UInt64, shape=128)

    for j in range(128):
        ran[j] = starts(nupdate / 128 * j)

    for i in range(nupdate / 128):
        ran = (ran << 1) ^ (((ran & temp) >> 63) * 7).astype(N.UInt64)
        #
        # This loop cannot be "vectorized" with numarray because ^= doesn't combine results.
        # As a simple example consider:
        # a=N.zeros(1);a[N.zeros(2)]+=1
        # a[1] is 1 because a[index_array] makes a copy then performs an operation and then copies back (or something to that effect).

        for j in range(128):
            table[ran[j] & (n-1)] ^= ran[j]

        #N.put(table, ran & (n-1), N.take(table, ran & (n-1)) ^ ran)
        #N.put(table, ran & (n-1), table[ran & (n-1)] ^ ran)
        #table[:] = N.take(table, ran & (n-1)) ^ ran

    return

    for i in range(nupdate / 128):
        ran = (ran << 1) ^ ((ran.astype(N.Int64)<0)*7).astype(N.UInt64)
        print ran & (n-1), ran.type()
        table[ran & (n-1)] ^= ran

    return

    for i in range(nupdate / 128):
        for j in range(128):
            temp[:] = ran[j:j+1]
            temp = (temp << 1) ^ ((temp.astype(N.Int64)<0)[0] and POLY or 0)
            ran[j] = temp[0]
            table[ran[j] & (n-1)] ^= ran[j]


def RandomAccessCheck(nupdate, n, table):
    temp = N.array([1], N.UInt64)
    for i in range(nupdate):
        temp = (temp << 1) ^ ((temp.astype(N.Int64)<0)[0] and POLY or 0)
        table[temp & (n-1)] ^= temp

    errors = n - (table-N.arange(n, type=N.UInt64) == 0).sum()
    ferr = float(errors)/float(n) * 100.0
    return errors, ferr, ferr <= 0.01

def main(argv):
    n = 128 / 4
    if len(argv) > 1:
        try:
            n = int(argv[1])
        except: pass

    nupdate = 4 * n
    table = N.array(type=N.UInt64, shape=n)
    RandomAccess(nupdate, n, table)
    print table
    errors, ferr, success = RandomAccessCheck(nupdate, n, table)
    print errors, ferr, success

if __name__ == "__main__":
    sys.exit(main(sys.argv))
