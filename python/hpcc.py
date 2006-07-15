#! /usr/bin/env python
# -*- mode: Python; tab-width: 4; indent-tabs-mode: nil; fill-column: 79;  coding: iso-latin-1-unix -*-
#
# HPC Challenge Benchmark
# Piotr Luszczek
# University of Tennessee Knoxville
# Innovative Computing Laboratory
# (C) Copyright 2005 All Rights Reserved
#

import sys

import numarray as N
import numarray.ufunc as U

import mpi

G_SIZEOF_INT = 4
G_SIZEOF_DOUBLE = 8
G_SIZEOF_UINT64 = 8

G_LONG_IS_64BITS = False

G_MAX_TOTAL_PENDING_UPDATES = 1024
G_LOCAL_BUFFER_SIZE = G_MAX_TOTAL_PENDING_UPDATES
G_USE_MULTIPLE_RECV = True or  False
G_MAX_RECV = G_USE_MULTIPLE_RECV and 16 or 1

G_INT64_DT = mpi.data_type(N.array(0, type=N.UInt64))

G_FINISHED_TAG = 1
G_UPDATE_TAG = 2

G_POLY, G_PERIOD = 7, 1317624576693539401L

def starts(n):

    n = N.array([n], N.Int64)
    m2 = N.zeros([64], N.UInt64)

    while n[0] < 0:        n += G_PERIOD
    while n[0] > G_PERIOD: n -= G_PERIOD
    if n[0] == 0: return 1

    temp = N.array([1], N.UInt64)
    ival = N.array([0], N.UInt64)
    for i in range(64):
        m2[i] = temp[0]
        temp = (temp << 1) ^ ((temp.astype(N.Int64)<0)[0] and G_POLY or 0)
        temp = (temp << 1) ^ ((temp.astype(N.Int64)<0)[0] and G_POLY or 0)
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
            ran = (ran << 1) ^ ((ran.astype(N.Int64)<0)[0] and G_POLY or 0)
    return ran[0]

class Buckets:
    def __init__(self, n_keys, n_vals):
        self.d = {}

    def InsertUpdate(self, v, key):
        d = self.d
        if not d.has_key(key):
            d[key] = [v]
        else:
            d[key].append(v)

    def GetUpdates(self, buf): # get an entry with
        d = self.d

        lmax = 0
        for k in d:
            l = len(d[k])
            if l > lmax:
                lmax = l
                key = k

        #print lmax, d[key], len(d[key]), buf.type()
        #buf[:lmax] = d[key]
        for i in range(lmax):
            buf[i:i+1] = d[key][i]
        del d[key]

        return key, lmax

class HPCC:
    def __init__(self, in_f, out_f):
        self.HPL_max_proc_mem = 1024**2
        self.comm = mpi.COMM_WORLD
        self.rank = self.comm.rank()
        self.size = self.comm.size()
        if 0 == self.rank:
            self.log_file = open(out_f, "a")

        self.results = {"mpi_comm_size": self.size}

        if G_LONG_IS_64BITS:
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
        max_int_bits2 = G_SIZEOF_INT * 8 - 2

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

        VectorSize = self.local_vector_size( 3, G_SIZEOF_DOUBLE, False )

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

        return 1e-9 * N.array([2, 2, 3, 3]) * G_SIZEOF_DOUBLE / N.array(map(min, times[:, 1:])) * VectorSize

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

        NumProcs, MyProc = self.size, self.rank
        NumberReceiving = NumProcs - 1

        Table = N.arange(LocalTableSize, type=N.UInt64)
        Table += GlobalStartMyProc

        if G_USE_MULTIPLE_RECV:
            inreq = mpi.Request(G_MAX_RECV)
        else:
            inreq = mpi.Request()
        outreq = mpi.Request(); outreq[:] = mpi.REQUEST_NULL

        finish_statuses = mpi.Status(NumProcs)
        finish_req = mpi.Request(NumProcs)

        status = mpi.Status()

        recv_size = G_LOCAL_BUFFER_SIZE
        LocalSendBuffer = N.array(type=N.UInt64, shape=recv_size)
        LocalRecvBuffer = N.array(type=N.UInt64, shape=G_MAX_RECV*recv_size)

        buckets = Buckets(NumProcs, recv_size)

        SendCnt = ProcNumUpdates
        Ran = N.array([starts(4 * GlobalStartMyProc)], type=N.UInt64)

        pendingUpdates = 0
        i = 0

        if G_USE_MULTIPLE_RECV:
            NumRecvs = NumProcs > 4 and min(4,G_MAX_RECV) or 1
            for j in range(NumRecvs):
                self.comm.Irecv(LocalRecvBuffer[j*recv_size:(j+1)*recv_size], mpi.ANY_SOURCE, mpi.ANY_TAG, inreq[j])
        else:
            self.comm.Irecv(LocalRecvBuffer, mpi.ANY_SOURCE, mpi.ANY_TAG, inreq)

        while i < SendCnt:
            while True: # receive messages
                if G_USE_MULTIPLE_RECV:
                    index, have_done = inreq[:NumRecvs].Testany(status)
                else:
                    have_done = inreq.Test(status)

                if have_done:

                    if status.MPI_TAG == G_UPDATE_TAG:
                        recvUpdates = status.Get_count(G_INT64_DT)

                        if G_USE_MULTIPLE_RECV:
                            bufferBase = index*recv_size
                        else:
                            bufferBase = 0

                        for j in range(recvUpdates):
                            inmsg = N.array(LocalRecvBuffer[bufferBase+j], type=N.UInt64)
                            LocalOffset = (inmsg & (TableSize - 1)) - GlobalStartMyProc
                            Table[LocalOffset] ^= inmsg

                    elif status.MPI_TAG == G_FINISHED_TAG: # we got a done message.  Thanks for playing...
                        NumberReceiving -= 1
                    else:
                        self.comm.Abort()

                    if G_USE_MULTIPLE_RECV:
                        self.comm.Irecv(LocalRecvBuffer[index*recv_size:(index+1)*recv_size], mpi.ANY_SOURCE, mpi.ANY_TAG, inreq[index])
                    else:
                        self.comm.Irecv(LocalRecvBuffer, mpi.ANY_SOURCE, mpi.ANY_TAG, inreq)

                if have_done and NumberReceiving > 0:
                    continue
                break

            if pendingUpdates < recv_size:
                Ran = (Ran << 1) ^ ((Ran.astype(N.Int64) < 0)[0] and G_POLY or 0)
                GlobalOffset = (Ran & (TableSize-1))[0]

                if GlobalOffset < Top:
                    WhichPe = GlobalOffset / (MinLocalTableSize + 1)
                else:
                    WhichPe = (GlobalOffset - Remainder) / MinLocalTableSize

                if WhichPe == MyProc:
                    LocalOffset = GlobalOffset - GlobalStartMyProc
                    Table[LocalOffset] ^= Ran[0]
                else:
                    buckets.InsertUpdate(Ran.copy(), WhichPe)
                    pendingUpdates += 1

                i += 1

            else:
                have_done = outreq.Test(mpi.STATUS_IGNORE)
                if have_done:
                    outreq[:] = mpi.REQUEST_NULL
                    pe, peUpdates = buckets.GetUpdates(LocalSendBuffer)
                    self.comm.Isend(LocalSendBuffer[:peUpdates], pe, G_UPDATE_TAG, outreq)
                    pendingUpdates -= peUpdates

        while pendingUpdates > 0: # send remaining updates in buckets

            while True: # receive messages
                if G_USE_MULTIPLE_RECV:
                    index, have_done = inreq[:NumRecvs].Testany(status)
                else:
                    have_done = inreq.Test(status)

                if have_done:
                    if status.MPI_TAG == G_UPDATE_TAG:
                        recvUpdates = status.Get_count(G_INT64_DT)
                        if G_USE_MULTIPLE_RECV:
                            bufferBase = index*recv_size
                        else:
                            bufferBase = 0

                        for j in range(recvUpdates):
                            inmsg = N.array(LocalRecvBuffer[bufferBase+j], type=N.UInt64)
                            LocalOffset = (inmsg & (TableSize - 1)) - GlobalStartMyProc
                            Table[LocalOffset] ^= inmsg

                    elif status.MPI_TAG == G_FINISHED_TAG: # we got a done message.  Thanks for playing...
                        NumberReceiving -= 1
                    else:
                        self.comm.Abort()

                    if G_USE_MULTIPLE_RECV:
                        self.comm.Irecv(LocalRecvBuffer[index*recv_size:(index+1)*recv_size], mpi.ANY_SOURCE, mpi.ANY_TAG, inreq[index])
                    else:
                        self.comm.Irecv(LocalRecvBuffer, mpi.ANY_SOURCE, mpi.ANY_TAG, inreq)

                if have_done and NumberReceiving > 0:
                    continue
                break

            have_done = outreq.Test(mpi.STATUS_IGNORE)
            if have_done:
                outreq[:] = mpi.REQUEST_NULL
                pe, peUpdates = buckets.GetUpdates(LocalSendBuffer)
                self.comm.Isend(LocalSendBuffer[:peUpdates], pe, G_UPDATE_TAG, outreq)
                pendingUpdates -= peUpdates

        for proc_count in range(NumProcs): # send our done messages
            if proc_count == MyProc:
                finish_req[MyProc] = mpi.REQUEST_NULL
                continue

            # send garbage - who cares, no one will look at it
            self.comm.Isend(Ran, proc_count, G_FINISHED_TAG, finish_req[proc_count])

        while NumberReceiving > 0: # Finish everyone else up...
            if G_USE_MULTIPLE_RECV:
                index = inreq[:NumRecvs].Waitany(status)
            else:
                inreq.Wait(status)

            if status.MPI_TAG == G_UPDATE_TAG:
                recvUpdates = status.Get_count(G_INT64_DT)

                if G_USE_MULTIPLE_RECV:
                    bufferBase = index * recv_size
                else:
                    bufferBase = 0

                for j in range(recvUpdates):
                    inmsg = N.array(LocalRecvBuffer[bufferBase+j], type=N.UInt64)
                    LocalOffset = (inmsg & (TableSize - 1)) - GlobalStartMyProc
                    Table[LocalOffset] ^= inmsg

            elif status.MPI_TAG == G_FINISHED_TAG: # we got a done message.  Thanks for playing...
                NumberReceiving -= 1
            else:
                self.comm.Abort()

            if G_USE_MULTIPLE_RECV:
                self.comm.Irecv(LocalRecvBuffer[index*recv_size:(index+1)*recv_size], mpi.ANY_SOURCE, mpi.ANY_TAG, inreq[index])
            else:
                self.comm.Irecv(LocalRecvBuffer, mpi.ANY_SOURCE, mpi.ANY_TAG, inreq)

        finish_req.Waitall(finish_statuses)

        if G_USE_MULTIPLE_RECV:
            for j in range(NumRecvs):
                inreq[j].Cancel()
                inreq[j].Wait(mpi.STATUS_IGNORE)
        else:
            inreq.Cancel()
            inreq.Wait(mpi.STATUS_IGNORE)

        return Table

    def mpira_verify(self, logTableSize, TableSize, LocalTableSize, MinLocalTableSize, GlobalStartMyProc, Top,
                  Remainder, ProcNumUpdates, Table):
        FALSE = 0
        TRUE = 1
        BUCKET_SIZE = 1024
        DONE = 0

        SLOT_CNT = 1
        FIRST_SLOT = 2

        NumProcs, MyProc = self.size, self.rank

        LocalAllDone = FALSE

        LocalBuckets =  N.array(shape=NumProcs*(BUCKET_SIZE+FIRST_SLOT), type=N.UInt64)
        GlobalBuckets = N.array(shape=NumProcs*(BUCKET_SIZE+FIRST_SLOT), type=N.UInt64)

        SendCnt = ProcNumUpdates
        Ran = N.array([starts(4 * GlobalStartMyProc)], type=N.UInt64)

        PeCheckDone = N.array(shape=NumProcs, type=N.Int64)
        PeCheckDone[:] = FALSE

        while LocalAllDone == FALSE:
            if SendCnt > 0:
                for i in range(NumProcs): # Initalize local buckets
                    PeBucketBase = i * (BUCKET_SIZE+FIRST_SLOT)
                    LocalBuckets[PeBucketBase+SLOT_CNT] = FIRST_SLOT
                    LocalBuckets[PeBucketBase+DONE] = FALSE

                # Fill local buckets until one is full or out of data
                NextSlot = FIRST_SLOT
                while NextSlot != (BUCKET_SIZE+FIRST_SLOT) and SendCnt>0:
                    Ran = (Ran << 1) ^ ((Ran.astype(N.Int64) < 0)[0] and G_POLY or 0)
                    GlobalOffset = (Ran & (TableSize-1))[0]
                    if GlobalOffset < Top:
                        WhichPe = GlobalOffset / (MinLocalTableSize + 1)
                    else:
                        WhichPe = (GlobalOffset - Remainder) / MinLocalTableSize
                    PeBucketBase = WhichPe * (BUCKET_SIZE+FIRST_SLOT)
                    NextSlot = LocalBuckets[PeBucketBase+SLOT_CNT]
                    LocalBuckets[PeBucketBase+NextSlot] = Ran[0]
                    NextSlot += 1
                    LocalBuckets[PeBucketBase+SLOT_CNT] = NextSlot
                    SendCnt -= 1

                if SendCnt == 0:
                    LocalBuckets[DONE:DONE+NumProcs*(BUCKET_SIZE+FIRST_SLOT):BUCKET_SIZE+FIRST_SLOT] = TRUE

            # End of sending loop

            self.comm.Barrier()

            LocalAllDone = TRUE

            # Now move all the buckets to the appropriate pe
            self.comm.Alltoall(LocalBuckets[:BUCKET_SIZE+FIRST_SLOT], GlobalBuckets[:BUCKET_SIZE+FIRST_SLOT])

            for i in range(NumProcs):
                if PeCheckDone[i] == FALSE:
                    PeBucketBase = i * (BUCKET_SIZE+FIRST_SLOT)
                    PeCheckDone[i] = GlobalBuckets[PeBucketBase+DONE]
                    for j in range (FIRST_SLOT, GlobalBuckets[PeBucketBase+SLOT_CNT]):
                        RanTmp = N.array([GlobalBuckets[PeBucketBase+j]], type=N.UInt64)
                        GlobalOffset = RanTmp & (TableSize - 1)
                        LocalOffset = GlobalOffset - GlobalStartMyProc
                        Table[LocalOffset] ^= RanTmp
                    LocalAllDone &= PeCheckDone[i]
        # no more local data

        # FIXME: make it faster
        errors = 0
        for i in range(LocalTableSize):
            if Table[i] != i + GlobalStartMyProc:
                errors += 1

        return errors

    def mpira(self):
        total_mem = self.HPL_max_proc_mem * self.size / G_SIZEOF_UINT64

        TableSize = 1
        logTableSize = 0
        total_mem /= 2
        while total_mem >= 1:
            total_mem /= 2
            logTableSize += 1
            TableSize <<= 1

        NumProcs, MyProc = self.size, self.rank

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

        NumUpdates_Default = 4 * TableSize # Default number of global updates to table: 4x number of table entries

        ProcNumUpdates = 4*LocalTableSize
        NumUpdates = NumUpdates_Default

        self.comm.Barrier()

        t = -mpi.Wtime()

        Table = self.mpira_any(logTableSize, TableSize, LocalTableSize, MinLocalTableSize, GlobalStartMyProc, Top,
                       Remainder, ProcNumUpdates)

        self.comm.Barrier()

        t += mpi.Wtime()

        self.results["mpira_gups"] = 1e-9 * ProcNumUpdates / t
        self.results["mpira_updates"] = ProcNumUpdates
        self.results["mpira_time"] = t

        t = -mpi.Wtime()
        errors = self.mpira_verify(logTableSize, TableSize, LocalTableSize, MinLocalTableSize, GlobalStartMyProc,
        #errors = mpi.verify(logTableSize, TableSize, LocalTableSize, MinLocalTableSize, GlobalStartMyProc,
                          Top, Remainder, ProcNumUpdates, Table)
        t += mpi.Wtime()

        self.results["mpira_vtime"] = t
        self.results["mpira_errors"] = errors

def main(argv):
    mpi.Init(argv)

    hpcc = HPCC("hpccinf.txt", "hpccoutf.txt")

    hpcc.ep_stream()
    hpcc.mpira()

    keys = hpcc.results.keys()
    keys.sort()
    for k in keys:
        hpcc.log( "%s=%f" % (k, hpcc.results[k]) )

    mpi.Finalize()
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
