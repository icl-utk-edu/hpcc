#! /usr/bin/env python
# -*- mode: Python; tab-width: 4; indent-tabs-mode: nil; fill-column: 79;  coding: iso-latin-1-unix -*-
#

import sys

hauxil = "hpl/include/hpl_auxil.h"
hblas = "hpl/include/hpl_blas.h"
hcomm = "hpl/include/hpl_comm.h"
hgesv = "hpl/include/hpl_gesv.h"
hgrid = "hpl/include/hpl_grid.h"
hmatgen = "hpl/include/hpl_matgen.h"
hmisc = "hpl/include/hpl_misc.h"
hpanel = "hpl/include/hpl_panel.h"
hpauxil = "hpl/include/hpl_pauxil.h"
hpfact = "hpl/include/hpl_pfact.h"
hpgesv = "hpl/include/hpl_pgesv.h"
hpmatgen = "hpl/include/hpl_pmatgen.h"
hpmisc = "hpl/include/hpl_pmisc.h"
hptest = "hpl/include/hpl_ptest.h"
hptimer = "hpl/include/hpl_ptimer.h"
htimer = "hpl/include/hpl_timer.h"

deps = (
    ("hpl/src/auxil/HPL_", (hmisc, hblas, hauxil),
     ("dlacpy", "dlatcpy", "fprintf", "warn", "abort", "dlaprnt", "dlange"), ""),

    ("hpl/src/auxil/HPL_", (hmisc, hblas, hauxil),
     ("dlamch", ), "$(CCNOOPT)"),

    ("hpl/src/blas/HPL_", (hmisc, hblas),
     ("dcopy", "daxpy", "dscal", "idamax", "dgemv", "dtrsv", "dger", "dgemm", "dtrsm"), ""),

    ("hpl/src/comm/HPL_", (hmisc, hpmisc, hgrid, hpanel, hpgesv),
     ("1ring", "1rinM", "2ring", "2rinM", "blong", "blonM", "packL", "copyL", "binit", "bcast", "bwait", "send", "recv", "sdrv"), ""),

    ("hpl/src/grid/HPL_", (hmisc, hpmisc, hgrid),
     ("grid_init", "pnum", "grid_info", "grid_exit", "broadcast", "reduce", "all_reduce", "barrier", "min", "max", "sum"), ""),

    ("hpl/src/panel/HPL_", (hmisc, hblas, hauxil, hpmisc, hgrid, hcomm, hpauxil, hpanel, hpfact, hpgesv),
     ("pdpanel_new", "pdpanel_init", "pdpanel_disp", "pdpanel_free"), ""),

    ("hpl/src/pauxil/HPL_", (hmisc, hblas, hauxil, hpmisc, hgrid, hpauxil),
     ("indxg2l", "indxg2lp", "indxg2p", "indxl2g", "infog2l", "numroc", "numrocI", "dlaswp00N", "dlaswp10N", "dlaswp01N", "dlaswp01T", "dlaswp02N", "dlaswp03N", "dlaswp03T", "dlaswp04N", "dlaswp04T", "dlaswp05N", "dlaswp05T", "dlaswp06N", "dlaswp06T", "pwarn", "pabort", "pdlaprnt", "pdlamch", "pdlange"), ""),

    ("hpl/src/pfact/HPL_", (hmisc, hblas, hauxil, hpmisc, hpauxil, hpfact),
     ("dlocmax", "dlocswpN", "dlocswpT", "pdmxswp", "pdpancrN", "pdpancrT", "pdpanllN", "pdpanllT", "pdpanrlN", "pdpanrlT", "pdrpanllN", "pdrpanllT", "pdrpancrN", "pdrpancrT", "pdrpanrlN", "pdrpanrlT", "pdfact"), ""),

    ("hpl/src/pgesv/HPL_", (hmisc, hblas, hauxil, hpmisc, hgrid, hcomm, hpauxil, hpanel, hpfact, hpgesv),
     ("pipid", "plindx0", "pdlaswp00N", "pdlaswp00T", "perm", "logsort", "plindx10", "plindx1", "spreadN", "spreadT", "rollN", "rollT", "equil", "pdlaswp01N", "pdlaswp01T", "pdupdateNN", "pdupdateNT", "pdupdateTN", "pdupdateTT", "pdtrsv", "pdgesv0", "pdgesvK1", "pdgesvK2", "pdgesv"), ""),

    ("hpl/testing/matgen/HPL_", (hmisc, hblas, hauxil, hmatgen),
     ("dmatgen", "ladd", "lmul", "xjumpm", "jumpit", "rand", "setran"), ""),

    ("hpl/testing/timer/HPL_", (hpmisc, htimer),
     ("timer", "timer_cputime", "timer_walltime"), ""),

    ("hpl/testing/pmatgen/HPL_", (hpmisc, hmatgen, hpmisc, hpauxil, hpmatgen),
     ("pdmatgen",), ""),

    ("hpl/testing/ptimer/HPL_", (hpmisc, hptimer),
     ("ptimer", "ptimer_cputime", "ptimer_walltime"), ""),

    ("hpl/src/pgesv/HPL_", (hmisc, hblas, hauxil, hgesv, hpmisc, hpauxil, hpanel, hpmatgen, hpgesv, hptimer, hptest),
     ("pddriver", "pdinfo", "pdtest"), "")
    )

def Gen(deps):
    print "# -*- Makefile -*-"
    print

    i = 0
    for d in deps:
        i = i + 1

        prfx, hfiles, files, flags = d
        hvar = "HDEP" + str(i)

        print hvar, "=", reduce(lambda x, y: x + " " + y, hfiles)

    print

    i = 0
    for d in deps:
        i = i + 1

        prfx, hfiles, files, flags = d
        hvar = "HDEP" + str(i)

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
