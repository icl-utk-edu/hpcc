#! /usr/bin/env python
# -*- mode: Python; tab-width: 4; indent-tabs-mode: nil; fill-column: 79;  coding: iso-latin-1-unix -*-
#
# Use this file to generate ../hpl/Makefile.hpcc
#

import os, string, sys

hauxil = "../../../include/hpl_auxil.h"
hblas = "../../../include/hpl_blas.h"
hcomm = "../../../include/hpl_comm.h"
hgesv = "../../../include/hpl_gesv.h"
hgrid = "../../../include/hpl_grid.h"
hmatgen = "../../../include/hpl_matgen.h"
hmisc = "../../../include/hpl_misc.h"
hpanel = "../../../include/hpl_panel.h"
hpauxil = "../../../include/hpl_pauxil.h"
hpfact = "../../../include/hpl_pfact.h"
hpgesv = "../../../include/hpl_pgesv.h"
hpmatgen = "../../../include/hpl_pmatgen.h"
hpmisc = "../../../include/hpl_pmisc.h"
hptest = "../../../include/hpl_ptest.h"
hptimer = "../../../include/hpl_ptimer.h"
htest = "../../../include/hpl_test.h" # this file seems to contain unsued functions
htimer = "../../../include/hpl_timer.h"
hhpl = "../../../include/hpl.h"
hhpcc = "../../../../include/hpcc.h"
hhpccv = "../../../../include/hpccver.h"
hhpccm = "../../../include/hpccmema.h"

allDeps = (
    ("src/auxil/HPL_", (hmisc, hblas, hauxil, htest),
     ("dlacpy", "dlatcpy", "fprintf", "warn", "abort", "dlaprnt", "dlange"), ""),

    ("src/auxil/HPL_", (hmisc, hblas, hauxil),
     ("dlamch", ), "$(CCNOOPT)"),

    ("src/blas/HPL_", (hmisc, hblas),
     ("dcopy", "daxpy", "dscal", "idamax", "dgemv", "dtrsv", "dger", "dgemm", "dtrsm"), ""),

    ("src/comm/HPL_", (hmisc, hpmisc, hgrid, hpanel, hpgesv),
     ("1ring", "1rinM", "2ring", "2rinM", "blong", "blonM", "packL", "copyL", "binit", "bcast", "bwait", "send", "recv", "sdrv"), ""),

    ("src/grid/HPL_", (hmisc, hpmisc, hgrid),
     ("grid_init", "pnum", "grid_info", "grid_exit", "broadcast", "reduce", "all_reduce", "barrier", "min", "max", "sum"), ""),

    ("src/panel/HPL_", (hmisc, hblas, hauxil, hpmisc, hgrid, hcomm, hpauxil, hpanel, hpfact, hpgesv),
     ("pdpanel_new", "pdpanel_init", "pdpanel_disp", "pdpanel_free"), ""),

    ("src/pauxil/HPL_", (hmisc, hblas, hauxil, hpmisc, hgrid, hpauxil),
     ("indxg2l", "indxg2lp", "indxg2p", "indxl2g", "infog2l", "numroc", "numrocI", "dlaswp00N", "dlaswp10N", "dlaswp01N", "dlaswp01T", "dlaswp02N", "dlaswp03N", "dlaswp03T", "dlaswp04N", "dlaswp04T", "dlaswp05N", "dlaswp05T", "dlaswp06N", "dlaswp06T", "pwarn", "pabort", "pdlaprnt", "pdlamch", "pdlange"), ""),

    ("src/pfact/HPL_", (hmisc, hblas, hauxil, hpmisc, hpauxil, hpfact),
     ("dlocmax", "dlocswpN", "dlocswpT", "pdmxswp", "pdpancrN", "pdpancrT", "pdpanllN", "pdpanllT", "pdpanrlN", "pdpanrlT", "pdrpanllN", "pdrpanllT", "pdrpancrN", "pdrpancrT", "pdrpanrlN", "pdrpanrlT", "pdfact"), ""),

    ("src/pgesv/HPL_", (hmisc, hblas, hauxil, hpmisc, hgrid, hcomm, hpauxil, hpanel, hpfact, hpgesv),
     ("pipid", "plindx0", "pdlaswp00N", "pdlaswp00T", "perm", "logsort", "plindx10", "plindx1", "spreadN", "spreadT", "rollN", "rollT", "equil", "pdlaswp01N", "pdlaswp01T", "pdupdateNN", "pdupdateNT", "pdupdateTN", "pdupdateTT", "pdtrsv", "pdgesv0", "pdgesvK1", "pdgesvK2", "pdgesv"), ""),

    ("testing/matgen/HPL_", (hmisc, hblas, hauxil, hmatgen),
     ("dmatgen", "ladd", "lmul", "xjumpm", "jumpit", "rand", "setran"), ""),

    ("testing/timer/HPL_", (hpmisc, htimer),
     ("timer", "timer_cputime", "timer_walltime"), ""),

    ("testing/pmatgen/HPL_", (hpmisc, hmatgen, hpmisc, hpauxil, hpmatgen),
     ("pdmatgen",), ""),

    ("testing/ptimer/HPL_", (hpmisc, hptimer),
     ("ptimer", "ptimer_cputime", "ptimer_walltime"), ""),

    ("testing/ptest/HPL_", (hmisc, hblas, hauxil, hgesv, hpmisc, hpauxil, hpanel, hpmatgen, hpgesv, hptimer, hptest),
     ("pddriver", "pdinfo", "pdtest"), ""),

    ("../RandomAccess/", (hhpcc, hhpl, "../../../../RandomAccess/RandomAccess.h",
                          "../../../../RandomAccess/buckets.h", "../../../../RandomAccess/heap.h",
                          "../../../../RandomAccess/pool.h", "../../../../RandomAccess/time_bound.h"),
     ("MPIRandomAccess", "buckets", "core_single_cpu_lcg", "core_single_cpu", "heap", "pool", "single_cpu_lcg",
      "single_cpu", "star_single_cpu_lcg", "star_single_cpu", "time_bound", "time_bound_lcg", "utility", "verification_lcg", "verification",
      "MPIRandomAccess_vanilla", "MPIRandomAccess_opt", "MPIRandomAccessLCG", "MPIRandomAccessLCG_vanilla",
      "MPIRandomAccessLCG_opt"), "-I../../../../include $(CCFLAGS)"),

    ("../STREAM/", (hhpcc, hhpl),
     ("onecpu", "stream"), "-I../../../../include $(CCFLAGS)"),

    ("../PTRANS/", (hhpcc, hhpl, "../../../../PTRANS/cblacslt.h"),
     ("pmatgeninc", "pdmatgen", "pdtransdriver", "pdmatcmp", "pdtrans", "sclapack", "cblacslt", "mem"), "-I../../../../include $(CCFLAGS)"),

    ("../src/", (hhpcc, hhpccv, hhpccm, hhpl),
     ("bench_lat_bw_1.5.2",  "hpcc",  "io", "extinit", "extfinalize"), "-I../../../../include $(CCFLAGS)"),

    ("../src/", (hhpcc, hhpl),
     ("HPL_slamch", "noopt"), "-I../../../../include $(CCNOOPT)"),

    ("../DGEMM/", (hhpcc, hhpl),
     ("tstdgemm",  "onecpu"), "-I../../../../include $(CCFLAGS)"),

    ("../FFT/", (hhpcc, hhpl, "../../../../FFT/hpccfft.h", "../../../../FFT/wrapfftw.h", "../../../../FFT/wrapmpifftw.h"),
     ("bcnrand", "fft235", "zfft1d", "pzfft1d", "onecpu", "tstfft", "wrapfftw", "wrapmpifftw", "mpifft"), "-I../../../../include $(CCFLAGS)")

    )

def Gen(deps):
    ilines = ("# -*- Makefile -*-", "", "arch = UNKNOWN", "include ../../../Make.$(arch)", "")
    for l in ilines:
        print l

    allobj = []

    i = 0
    for d in deps:
        i = i + 1

        prfx, hfiles, files, flags = d
        hvar = "HDEP" + str(i)
        ovar = "OBJS" + str(i)

        prfx = "../../../" + prfx

        allobj.append("$(%s)" % ovar)

        print hvar, "=", reduce(lambda x, y: x + " " + y, hfiles)

        print ovar, "=",
        for ff in files:
            f = prfx + ff
            fo = f + ".o"
            print fo,
        print

    exe = "../../../../hpcc"
    sallobj = reduce(lambda x, y: x + " " + y, allobj)
    print
    print exe, ": $(HPLlib)"
    print "\t$(LINKER) $(LINKFLAGS) -o", exe, "$(HPL_LIBS)"
    print
    print "$(HPLlib):", sallobj
    print "\t$(ARCHIVER) $(ARFLAGS) $(HPLlib)", sallobj
    print "\t$(RANLIB) $(HPLlib)"
    print
    print "clean:"
    print "\t$(RM)", exe, "$(HPLlib)"
    print "\t$(RM)", sallobj
    print

    i = 0
    for d in deps:
        i = i + 1

        prfx, hfiles, files, flags = d
        hvar = "HDEP" + str(i)

        prfx = "../../../" + prfx

        if not flags:
            flags = "$(CCFLAGS)"

        for ff in files:
            f = prfx + ff
            fo = f + ".o"
            fc = f + ".c"
            print fo, ":", fc, "$(%s)" % hvar
            print "\t$(CC) -o", fo, "-c", fc, flags
            print

def DirVisitor(newItems, dirname, items):
    if dirname[-3:] == "CVS": # skip CVS directories
        return

    for name in items:
        path = os.path.join(dirname, name)
        if name[-1] == "~": continue
        elif name in ("semantic.cache", "CVS"): continue
        elif os.path.isfile(path):
            newItems.append(path)

def TraverseDirs(prfx, items):
    newItems = []
    for i in items:
        name = os.path.join(prfx, i)

        if os.path.isdir(name):
            os.path.walk(name, DirVisitor, newItems)
        else:
            newItems.append(name)

    l = len(prfx) + 1
    for i in range(len(newItems)):
        newItems[i] = newItems[i][l:]

    return newItems

def Dist(deps, prfx="hpcc"):
    allPrfx = prfx + "/"

    addItems = ["Makefile", "README.html", "README.txt", "_hpccinf.txt",
                "hpl/Make.UNKNOWN", "hpl/BUGS", "hpl/COPYRIGHT", "hpl/HISTORY",
                "hpl/INSTALL", "hpl/Make.top", "hpl/Makefile",
                "hpl/README", "hpl/TODO", "hpl/TUNING", "hpl/lib/arch/build/Makefile.hpcc",
                "hpl/makes", "hpl/man", "hpl/setup", "hpl/www"]

    addItems = TraverseDirs(prfx, addItems)
    #print string.join(addItems, "\n")

    mf = os.path.join(prfx, "README.tex")
    mt = os.path.getmtime(mf)
    for f in ("README.html", "README.txt"):
        t = os.path.getmtime(os.path.join(prfx, f))
        if mt >= t:
            raise RuntimeError, "File " + f + " is older than " + mf


    allFiles = []
    hDict = {}
    for d in deps:
        prfx, hfiles, files, flags = d

        for f in files:
            if prfx[0] != ".": # if it's an HPL file
                ff = allPrfx + "hpl/" + prfx + f + ".c"
            else: # if it's not an HPL file
                ff = allPrfx + prfx[3:] + f + ".c"
            allFiles.append(ff)

        for h in hfiles:
            if h[9:11] == "..": # if it's not an HPL header
                hprfx = allPrfx
                idx = 12
            else:
                hprfx = allPrfx + "hpl/"
                idx = 9
            hh = hprfx + h[idx:]
            hDict[hh] = 1

    allItems = map(lambda x, p=allPrfx: p + x, addItems) + allFiles + hDict.keys()
    allItems.sort()
    # check existence and type of all files
    for fname in allItems:
        if not os.path.isfile(fname):
            raise RuntimeError, "File " + fname + " doesn't exist."

    if sys.platform in ("linux", "linux2"):
        cmnd = "tar"
        group = "root"
    elif sys.platform == "darwin":
        cmnd = "gnutar"
        group = "wheel"

    print "%s --group=%s --owner=root -cvohf " % (cmnd, group) + allPrfx[:-1] + ".tar", string.join( allItems, " " )


def main(argv):
    global allDeps

    if len(argv) > 1 and argv[1] == "dist":
        if len(argv) > 2: # use custom directory prefix
          prefix = argv[2]
        else:
          prefix = os.path.split(os.path.dirname(argv[0]))[-2]
        Dist(allDeps, prefix)
    else:
        Gen(allDeps)

    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
