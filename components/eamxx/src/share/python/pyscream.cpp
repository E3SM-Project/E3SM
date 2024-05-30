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
  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);
}
void finalize () {
  scorpio::finalize_subsystem();
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
    .def("sync_to_host",&PyField::sync_to_host)
    .def("print",&PyField::print);

  // Grid
  py::class_<PyGrid>(m,"Grid")
    .def(py::init<const std::string&,int,int>());

  // Atm process
  py::class_<PyAtmProc>(m,"AtmProc")
    .def(py::init<const PyGrid&>())
    .def("get_arr",&PyAtmProc::get_arr)
    .def("initialize",&PyAtmProc::initialize)
    .def("read_ic",&PyAtmProc::read_ic,py::arg("ic_filename") = "",py::arg("default_init") = true);
}

} // namespace scream
