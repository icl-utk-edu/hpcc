#! /usr/bin/env python
# -*- mode: Python; tab-width: 4; indent-tabs-mode: nil; fill-column: 79;  coding: iso-latin-1-unix -*-
#

import os, sys

from distutils.core import setup, Extension

module = Extension("mpi",
                   libraries = ["mpi", "lam"],
                   sources = ["mpi.c"])
setup (name = "mpi",
       version = "0.1",
       description = "MPI binding",
       author = "Piotr Luszczek",
       author_email = "luszczek __at__ cs __dot__ utk __dot__ edu",
       url = "http://icl.cs.utk.edu/hpcc/",
       long_description = """MPI Python binding using numarray.""",
       ext_modules = [module])
