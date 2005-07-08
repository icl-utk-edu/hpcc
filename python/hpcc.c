/* -*- mode: C; tab-width: 2; indent-tabs-mode: nil; -*- */

#include <Python.h>

DL_EXPORT(void) init_netlib(void); /*proto*/

int
main(int argc, char *argv[]) {
  PyObject *pName, *pModule, *pDict, *pFunc;

  Py_Initialize();

  if (argc < 2) {
    fprintf( stderr, "%s file.py\n", argv[0] );
    return 0;
  }

  PyRun_SimpleString( "import sys; sys.argv = ['hpcc.py']" );
  initmpi();

  /*
  pName = PyString_FromString("mpi");
  pModule = PyImport_Import(pName);
  if (!pModule) {
    PyErr_Print();
    return 0;
  }
  */

  PyRun_SimpleString( "import sys; sys.path.append('.')" );
  PyRun_SimpleString( "import hpcc; hpcc.main('hpcc.py')" );
  /*
  PyRun_SimpleString( "import sys; print sys.path" );
  PyRun_SimpleString( "execfile(\"server.py\")" );
  PyRun_SimpleFile( stderr, argv[1] );
  */

  Py_Finalize();
  return 0;
}
