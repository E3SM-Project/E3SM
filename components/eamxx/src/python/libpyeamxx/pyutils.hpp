#ifndef PYUTILS_HPP
#define PYUTILS_HPP

#include <nanobind/nanobind.h>
#include <mpi4py/mpi4py.h>
#include <mpi.h>

#include <typeinfo>

MPI_Comm get_c_comm (nb::object py_comm) {
  if (import_mpi4py() < 0) {
    throw nb::python_error();
  }
  auto py_src = py_comm.ptr();
  if (not PyObject_TypeCheck(py_src, &PyMPIComm_Type)) {
    throw std::bad_cast();
  }

  auto comm_ptr = PyMPIComm_Get(py_src);
  return *comm_ptr;
}

#endif // PYUTILS_HPP
