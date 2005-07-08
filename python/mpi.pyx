# -*- Pyrex -*-
#
# Pyrex functions/methods receive C-typed argments (other than primitive C types), only generic Python types or Pyrex C types (except of C-only functions defined with `cdef')
#

cdef extern from "numarray/libnumarray.h":
    ctypedef int maybelong

    cdef struct PyArray_Descr:
        int type_num # PyArray_TYPES
        int elsize   # bytes for 1 element
        char type    # One of "cb1silfdFD "  Object array not supported
        # function pointers omitted

    ctypedef class numarray._numarray._numarray [object PyArrayObject]:
        cdef char *data
        cdef int nd
        cdef maybelong *dimensions
        cdef maybelong *strides
        cdef object base
        cdef PyArray_Descr *descr
        cdef int flags

    ctypedef double Float64

    cdef struct s_Complex64:
        Float64 r
        Float64 i

    ctypedef s_Complex64 Complex64

cdef extern from "Python.h":
    char *PyString_AsString(object string)

cdef extern from "mpi.h":
    ctypedef struct MPI_Comm:
        pass
    ctypedef struct MPI_Datatype:
        pass
    ctypedef struct MPI_Op:
        pass
    ctypedef struct MPI_Request:
        pass
    ctypedef struct MPI_Status:
        int MPI_TAG

    MPI_Comm MPI_COMM_WORLD
    MPI_Datatype MPI_DOUBLE
    MPI_Datatype MPI_INT
    MPI_Datatype MPI_LONG
    MPI_Datatype MPI_LONG_LONG_INT
    MPI_Op MPI_SUM
    MPI_Request MPI_REQUEST_NULL
    MPI_Status *MPI_STATUS_IGNORE
    MPI_Status *MPI_STATUSES_IGNORE

    int MPI_SUCCESS
    int MPI_ANY_SOURCE
    int MPI_ANY_TAG
    int MPI_Abort(MPI_Comm comm, int errcode)
    int MPI_Allreduce(void *sbuf, void *rbuf, int count, MPI_Datatype dtype, MPI_Op op, MPI_Comm comm)
    int MPI_Alltoall(void *sbuf, int scount, MPI_Datatype sdtype, void* rbuf, int rcount, MPI_Datatype rdtype, MPI_Comm comm)
    int MPI_Barrier(MPI_Comm comm)
    int MPI_Cancel(MPI_Request *preq)
    int MPI_Comm_rank(MPI_Comm comm, int *rank)
    int MPI_Comm_size(MPI_Comm comm, int *rank)
    int MPI_Finalize()
    int MPI_Get_count(MPI_Status *stat, MPI_Datatype dtype, int *count)
    int MPI_Init(int *argc, char ***argv)
    int MPI_Irecv(void *buf, int count, MPI_Datatype dtype, int src, int tag, MPI_Comm comm, MPI_Request *req)
    int MPI_Isend(void *buf, int count, MPI_Datatype dtype, int dest, int tag, MPI_Comm comm, MPI_Request *req)
    int MPI_Send(void *buf, int count, MPI_Datatype dtype, int dest, int tag, MPI_Comm comm)
    int MPI_Test(MPI_Request *req, int *flag, MPI_Status *stat)
    int MPI_Testany(int count, MPI_Request *reqs, int *index, int *flag, MPI_Status *stat)
    int MPI_Wait(MPI_Request *preq, MPI_Status *stat)
    int MPI_Waitall(int count, MPI_Request *reqs, MPI_Status *stats)
    int MPI_Waitany(int count, MPI_Request *reqs, int *index, MPI_Status *stat)
    double MPI_Wtime()
    double MPI_Wtick()

cdef extern from "pyxutil.h":
    cdef enum _Type:
        _Request
        _Status
    void *pyxmpi_alloc(_Type t, int n) # allocator for MPI data types
    int Cverify(int logTableSize, int TableSize, int LocalTableSize, int MinLocalTableSize, int GlobalStartMyProc, int Top, int Remainder, int ProcNumUpdates, void *Table) 

def verify(int logTableSize, int TableSize, int LocalTableSize, int MinLocalTableSize, int GlobalStartMyProc, int Top, int Remainder, int ProcNumUpdates, _numarray Table):
    return Cverify(logTableSize, TableSize, LocalTableSize, MinLocalTableSize, GlobalStartMyProc, Top, Remainder, ProcNumUpdates, Table.data)

cdef MPI_Datatype GetTC(_numarray a, int *c):
    c[0] = a.dimensions[0]
    typ = a.descr.type
    #print c[0], typ, a.descr.elsize
    if 105 == typ: # i (Int)
        return MPI_INT
    elif 100 == typ: # d (Float64)
        return MPI_DOUBLE
    elif typ in (78, 85): # U (UInt64), N (Int64)
        return MPI_LONG_LONG_INT

cdef class Op:
    cdef MPI_Op op
    def __init__(self, op):
        if "sum" == op:
            self.op = MPI_SUM

cdef class Datatype:
    cdef MPI_Datatype dtype
    def __init__(self, t):
        if "LONG" == t:
            self.dtype = MPI_LONG
        else:
            self.dtype = MPI_LONG_LONG_INT

cdef class Status:
    cdef MPI_Status *status
    cdef int n

    def __init__(self, int n=1):
        assert n > 0
        self.status = <MPI_Status *>pyxmpi_alloc(_Status, n)
        self.n = n

    property MPI_TAG:
        def __get__(self):
            assert self.n == 1
            return self.status[0].MPI_TAG

    def Get_count(self, Datatype dtype):
        cdef int count

        assert self.n == 1

        #FIXME: return value of MPI_Get_count
        MPI_Get_count(self.status, dtype.dtype, &count)
        return count

cdef class Request:
    cdef MPI_Request *req
    cdef int n

    def __init__(self, int n=1, int view=0, object base=None, int idx=0):
        #print n, view, idx
        cdef Request base_
        if view:
            base_ = base
            self.req = base_.req + idx
            self.n = n
            return

        assert n > 0
        self.req = <MPI_Request *>pyxmpi_alloc(_Request, n)
        self.n = n

    def __getitem__(self, object x):
        cdef int start, stop #, step
        if type(x) is type(slice(None)):
            start = 0
            if x.start is not None:
                start = x.start
                assert 0 <= start < n
            stop = self.n
            if x.stop is not None:
                stop = x.stop
                assert start < stop <= self.n
            #print x, start, stop
            return Request(stop-start, 1, self, start)
        elif type(x) is type(0):
            return Request(1, 1, self, x)
        return None
        assert 0 <= i < self.n
        cdef Request r

        r = Request(i, 1, self)
        #r.req[0] = self.req[i]
        return r

    def __setitem__(self, object x, object y):
        cdef int i
        assert y is None
        #print x, y, type(x)
        if type(x) is type(slice(None)):
            for i from 0 <= i < self.n:
                self.req[i] = MPI_REQUEST_NULL
        elif type(x) is type(0):
            i = x
            self.req[i] = MPI_REQUEST_NULL

    def Cancel(self):
        assert self.n == 1
        MPI_Cancel(self.req)

    def Testany(self, Status status):
        cdef int index, flag
        cdef MPI_Status *lstat

        lstat = MPI_STATUS_IGNORE
        if status is not None:
            lstat = status.status

        #FIXME: return value of MPI_Testany() ?
        MPI_Testany(self.n, self.req, &index, &flag, lstat)
        return index, flag

    def Test(self, Status status):
        cdef int flag
        cdef MPI_Status *lstat

        lstat = MPI_STATUS_IGNORE
        if status is not None:
            lstat = status.status

        #FIXME: return value of MPI_Test() ?
        MPI_Test(self.req, &flag, lstat)
        return flag

    def Wait(self, Status status):
        cdef MPI_Status *lstat

        assert self.n == 1

        lstat = MPI_STATUSES_IGNORE
        if status is not None:
            lstat = status.status

        return MPI_Wait(self.req, lstat)

    def Waitall(self, Status status):
        cdef MPI_Status *lstat

        lstat = MPI_STATUSES_IGNORE
        if status is not None:
            lstat = status.status

        return MPI_Waitall(self.n, self.req, lstat)

    def Waitany(self, Status status):
        cdef int index
        cdef MPI_Status *lstat

        lstat = MPI_STATUS_IGNORE
        if status is not None:
            lstat = status.status

        #FIXME: return value of MPI_Waitany() ?
        MPI_Waitany(self.n, self.req, &index, lstat)
        return index

cdef class Comm:
    cdef MPI_Comm comm
    def _set_world(self):
        self.comm = MPI_COMM_WORLD

    def rank(self):
        cdef int r
        MPI_Comm_rank(self.comm, &r)
        return r

    def size(self):
        cdef int r
        MPI_Comm_size(self.comm, &r)
        return r

    def Abort(self, errcode):
        cdef int rv
        rv = MPI_Abort(self.comm, errcode)
        return rv

    def Alltoall(self, _numarray sbuf, _numarray rbuf):
        cdef int scount, rcount
        cdef MPI_Datatype sdtype, rdtype

        sdtype = GetTC(sbuf, &scount)
        rdtype = GetTC(rbuf, &rcount)

        return MPI_Alltoall(sbuf.data, scount, sdtype, rbuf.data, rcount, rdtype, self.comm)

    def Allreduce(self, _numarray src, _numarray dst, Op op):
        cdef int count
        cdef MPI_Datatype dtype

        dtype = GetTC(src, &count)
        return MPI_Allreduce(src.data, dst.data, count, dtype, op.op, self.comm)

    def Barrier(self): MPI_Barrier(self.comm)

    def Irecv(self, _numarray buf, int src, int tag, Request req):
        cdef int count
        cdef MPI_Datatype dtype

        assert req.n == 1

        dtype = GetTC(buf, &count)
        return MPI_Irecv(buf.data, count, dtype, src, tag, self.comm, req.req)

    def Isend(self, _numarray buf, int dest, int tag, Request req):
        cdef int count
        cdef MPI_Datatype dtype

        assert req.n == 1

        dtype = GetTC(buf, &count)
        return MPI_Isend(buf.data, count, dtype, dest, tag, self.comm, req.req)

    def Send(self, _numarray a, int dest, int tag):
        cdef int count
        cdef MPI_Datatype dtype

        dtype = GetTC(a, &count)
        return MPI_Send(a.data, count, dtype, dest, tag, self.comm)

def data_type(_numarray a):
    cdef Datatype dtype
    cdef int count

    dtype = Datatype("")
    dtype.dtype = GetTC(a, &count)

    return dtype

def Init(argv):
    cdef char *cargv[101]
    cdef char **wargv[1]
    cdef int i, l

    l = len(argv)
    if l > 100:
        l = 100

    for i from 0 <= i < l:
        cargv[i] = PyString_AsString(argv[i])

    cargv[l] = NULL
    wargv[0] = cargv

    i = MPI_Init( &i, wargv )

    return i

def Finalize():
    cdef int i
    i = MPI_Finalize()
    return i

def Comm_rank(Comm comm):
    return comm.rank()

def Comm_size(Comm comm):
    return comm.size()

def Wtime():
    return MPI_Wtime()

cdef Wtick():
    return MPI_Wtick()

global SUM
SUM = Op("sum")

global LONG
LONG = Datatype("LONG")

global LONG_LONG_INT
LONG_LONG_INT = Datatype("LONG_LONG_INT")

global COMM_WORLD
COMM_WORLD = Comm()
COMM_WORLD._set_world()

global ANY_SOURCE
ANY_SOURCE = MPI_ANY_SOURCE

global ANY_TAG
ANY_TAG = MPI_ANY_TAG

global STATUS_IGNORE
STATUS_IGNORE = None

global REQUEST_NULL
REQUEST_NULL = None
