#! /usr/bin/env python
# -*- mode: Python; tab-width: 4; indent-tabs-mode: nil; fill-column: 79;  coding: iso-latin-1-unix -*-
#
# Use this file to generate ../hpl/Makefile.hpcc
#

import sys

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
htimer = "../../../include/hpl_timer.h"
hhpl = "../../../include/hpl.h"
hhpcc = "../../../../include/hpcc.h"
hhpccv = "../../../../include/hpccver.h"

deps = (
    ("src/auxil/HPL_", (hmisc, hblas, hauxil),
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

    ("../RandomAccess/", (hhpcc, hhpl, "../../../../RandomAccess/RandomAccess.h"),
     ("MPIRandomAccess", "RandomAccess", "onecpu"), "-I../../../../include $(CCFLAGS)"),

    ("../STREAM/", (hhpcc, hhpl),
     ("onecpu", "stream"), "-I../../../../include $(CCFLAGS)"),

    ("../PTRANS/", (hhpcc, hhpl),
     ("pmatgeninc", "pdmatgen", "pdtransdriver", "pdmatcmp", "pdtrans", "sclapack", "cblacslt", "mem"), "-I../../../../include $(CCFLAGS)"),

    ("../src/", (hhpcc, hhpccv, hhpl),
     ("bench_lat_bw_1.5.1",  "hpcc",  "io"), "-I../../../../include $(CCFLAGS)"),

    ("../DGEMM/", (hhpcc, hhpl),
     ("tstdgemm",  "onecpu"), "-I../../../../include $(CCFLAGS)")

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
    print "\t$(LINKER) $(LINKFLAGS) -o", exe, sallobj, "$(HPL_LIBS)"
    print
    print "$(HPLlib):", sallobj
    print "\t$(ARCHIVER) $(ARFLAGS) $(HPLlib)", sallobj
    print "\t$(RANLIB) $(HPLlib)"
    print
    print "clean:"
    print "\t$(RM)", exe
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

def main(argv):
    global deps
    Gen(deps)
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
