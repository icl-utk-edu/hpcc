# -*- Makefile -*-

arch = UNKNOWN
include hpl/Make.$(arch)

all:
	- $(MKDIR) hpl/lib/$(arch)
	( $(CD) hpl/lib/arch/build ; $(MAKE) -f ../../../Makefile.hpcc )

clean:
	- $(MKDIR) hpl/lib/$(arch)
	( $(CD) hpl/lib/arch/build ; $(MAKE) -f ../../../Makefile.hpcc clean )
