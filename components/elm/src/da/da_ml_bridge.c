#include <stdio.h>
#include <Python.h>
#include "da.h"
#include <sys/time.h>

int
ml_init_c()
{
    int i;
    PyImport_AppendInittab("da", PyInit_da);
    Py_Initialize();
    PyImport_ImportModule("da");

    python_init();

    return 0;
}
