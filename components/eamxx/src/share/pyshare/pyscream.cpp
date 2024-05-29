#include "scream_session.hpp"
#include "pyfield.hpp"
#include "pygrid.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace scream {

void initialize () {
  initialize_scream_session(true);
}
void finalize () {
  finalize_scream_session();
}

namespace py = pybind11;

PYBIND11_MODULE (pyscream,m) {

  m.doc() = "Python interfaces to certain EAMxx infrastructure code";

  // Scream Session
  m.def("init",&initialize);
  m.def("finalize",&finalize);

  // Field class
  py::class_<PyField>(m,"Field")
    .def(py::init<const std::string&,
                  const std::vector<int>&>())
    .def(py::init<const std::string&,
                  py::array_t<double>
                 >())
    .def("print",&PyField::print)
    .def("cleanup",&PyField::cleanup);

  // PointGrid
  py::class_<PyPointGrid>(m,"PointGrid")
    .def(py::init<const std::string&,int,int>())
    .def("cleanup",&PyPointGrid::cleanup);
}

} // namespace scream
