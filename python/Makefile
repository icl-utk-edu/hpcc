# -*- Makefile -*-
#

PYTHON = python
RM = rm -rf
PYREX = $(HOME)/build/Pyrex-0.9.3/pyrexc.py 

all: mpi.so

mpi.c: mpi.pyx
	$(PYTHON) $(PYREX) mpi.pyx

mpi.so: mpi.c
	$(PYTHON) setup.py build_ext --inplace

.PHONY: clean dist dist-clean

clean:
	$(RM) mpi.c mpi.so *.o *.a *.pxi

dist:
	cd ../.. ; tar cvf hpcc.tar hpcc/python/mpi.pyx
	gzip --best ../hpcc.tar

dist-clean: clean
	$(RM) build *.pyc
