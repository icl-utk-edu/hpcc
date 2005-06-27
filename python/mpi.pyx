# -*- Pyrex -*-

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

    MPI_Comm MPI_COMM_WORLD
    MPI_Datatype MPI_DOUBLE
    MPI_Datatype MPI_INT
    MPI_Datatype MPI_LONG
    MPI_Datatype MPI_LONG_LONG_INT
    MPI_Op MPI_SUM

    int MPI_SUCCESS
    int MPI_Abort(MPI_Comm comm, int errcode)
    int MPI_Allreduce(void *sbuf, void *rbuf, int count, MPI_Datatype dtype, MPI_Op op, MPI_Comm comm)
    int MPI_Init(int *argc, char ***argv)
    int MPI_Finalize()
    int MPI_Comm_rank(MPI_Comm comm, int *rank)
    int MPI_Comm_size(MPI_Comm comm, int *rank)
    int MPI_Send(void *buf, int count, MPI_Datatype dtype, int dest, int tag, MPI_Comm comm)
    double MPI_Wtime()
    double MPI_Wtick()

cdef MPI_Datatype GetTC(_numarray a, int *c):
    c[0] = a.dimensions[0]
    return MPI_DOUBLE

cdef class Op:
    cdef MPI_Op op
    def __init__(self, op):
        if "sum" == op:
            self.op = MPI_SUM

cdef class Datatype:
    cdef MPI_Datatype tpe
    def __init__(self, t):
        if "LONG" == t:
            self.tpe = MPI_LONG
        else:
            self.tpe = MPI_LONG_LONG_INT

cdef class Comm:
    cdef MPI_Comm comm
    def __new__(self):
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

    def Allreduce(self, _numarray src, _numarray dst, Op op):
        cdef int count
        cdef MPI_Datatype typ

        typ = GetTC(src, &count)
        return MPI_Allreduce(src.data, dst.data, count, typ, op.op, self.comm)

    def Send(self, _numarray a, int dest, int tag):
        cdef int count
        cdef MPI_Datatype typ

        typ = GetTC(a, &count)
        return MPI_Send(a.data, count, typ, dest, tag, self.comm)

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

def Wtick():
    return MPI_Wtick()

global SUM
SUM = Op("sum")

global LONG
LONG = Datatype("LONG")

global LONG_LONG_INT
LONG_LONG_INT = Datatype("LONG_LONG_INT")

global COMM_WORLD
COMM_WORLD = Comm()
