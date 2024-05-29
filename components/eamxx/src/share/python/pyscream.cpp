#include "scream_session.hpp"
#include "pyfield.hpp"
#include "pygrid.hpp"
#include "pyatmproc.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace scream {

void initialize () {
  initialize_scream_session(true);
}
void finalize () {
  finalize_scream_session();
}

PYBIND11_MODULE (pyscream,m) {

  m.doc() = "Python interfaces to certain EAMxx infrastructure code";

  // Scream Session
  m.def("init",&initialize);
  m.def("finalize",&finalize);

  // Field class
  py::class_<PyField>(m,"Field")
    .def("get",&PyField::get)
    .def("print",&PyField::print);

  // Grid
  py::class_<PyGrid>(m,"Grid")
    .def(py::init<const std::string&,int,int>());

  // Atm process
  py::class_<PyAtmProc>(m,"AtmProc")
    .def(py::init<const PyGrid&>())
    .def("get_arr",&PyAtmProc::get_arr)
    .def("initialize",py::overload_cast<>(&PyAtmProc::initialize))
    .def("initialize",py::overload_cast<const std::string&>(&PyAtmProc::initialize));
}

} // namespace scream
