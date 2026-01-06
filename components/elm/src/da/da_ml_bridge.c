#include <stdio.h>
#include <Python.h>
#include "da.h"
#include <sys/time.h>

int
ml_init_c()
{
    int i;
    //PyImport_AppendInittab("ml_reservoir_release", PyInit_ml_reservoir_release);
    Py_Initialize();
    //PyImport_ImportModule("ml_reservoir_release");

    python_init();

    return 0;
}
