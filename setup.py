
import os
import sys

def err_message(msg):
  sys.stderr.write(msg)
  sys.stderr.write("\n")
  sys.stderr.flush()

def out_message(msg):
  sys.stdout.write(msg)
  sys.stdout.write("\n")
  sys.stdout.flush()

def check_pwd(show=1):
  for fname in ("hpl/setup/Make.UNKNOWN.in", "_hpccinf.txt"):
    if not os.path.exists(fname):
      if show:
        err_message("Not in HPCC top-level directory!")

        return -1

  return 0

def create_generic_makefile():
  generic_makefile = os.path.join("hpl", "setup", "Make.UNKNOWN.in")

  hpcc_arch=sys.platform
  makefile = os.path.join("hpl", "Make." + hpcc_arch)

  #if os.path.exists(generic_makefile):
    #shutil.copy(generic_makefile, makefile)

  return hpcc_arch, makefile, open(generic_makefile).read()

def get_makefile_defaults():
  def dict_from_args(**args):
    return args
  return dict_from_args(
    SHELL="sh",
    CD="cd",
    CP="cp",
    LN_S="ln -s",
    MKDIR="mkdir",
    RM="rm",
    TOUCH="touch",
    MPDIR="/mnt/scratch/sw/openmpi-1.3.1-gcc",
    MPINC="",
    MPLIB="-Wl,-rpath,$(MPdir)/lib",
    LADIR="/mnt/scratch/sw/atlas-gcc",
    LAINC="",
    #LALIB="-L$(LAdir)/lib -lf77blas -latlas -lm",
    LALIB="-L$(LAdir)/lib -lcblas -latlas -L/mnt/scratch/sw/intel/11.1.069/lib/intel64 -Wl,-rpath,/mnt/scratch/sw/intel/11.1.069/lib/intel64 -lirc -lm",
    #F2CDEFS="-DAdd_",
    F2CDEFS="-DHPL_CALL_CBLAS -DHPL_DETAILED_TIMING",
    CCNOOPT="",
    CCFLAGS="-O2 -DRA_SANDIA_NOPT=1",
    LINKER="mpicc",
    LINKFLAGS="",
    ARCHIVER="ar",
    ARFLAGS="rc",
    RANLIB="ranlib",
    CC="mpicc"
  )

def main(argv):

  if check_pwd(1):
    return 127

  hpcc_arch, makefile, makefile_text = create_generic_makefile()

  makefile_defaults = get_makefile_defaults()

  for key in makefile_defaults.keys():
    makefile_text = makefile_text.replace("@" + key + "@", makefile_defaults[key])

  open(makefile, "w").write(makefile_text)

  out_message("Created: " + makefile)
  out_message("Type: make arch=" + hpcc_arch)

  return 0

if "__main__" == __name__:
  sys.exit(main(sys.argv))
