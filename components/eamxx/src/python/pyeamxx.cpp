#include <share/scream_session.hpp>
#include "pyfield.hpp"
#include "pygrid.hpp"
#include "pyatmproc.hpp"
#include "pyparamlist.hpp"
#include "pyutils.hpp"

#include <ekat/mpi/ekat_comm.hpp>

#include <pybind11/pybind11.h>

#include <mpi.h>

namespace py = pybind11;

namespace scream {

void initialize (MPI_Comm mpi_comm) {
  ekat::Comm comm(mpi_comm);
  initialize_scream_session(comm.am_i_root());
  scorpio::init_subsystem(comm);
}

void initialize () {
  initialize(MPI_COMM_WORLD);
}
void initialize (pybind11::object py_comm) {
  initialize(get_c_comm(py_comm));
}
void finalize () {
  scorpio::finalize_subsystem();
  finalize_scream_session();
}

PYBIND11_MODULE (pyeamxx,m) {

  m.doc() = "Python interfaces to certain EAMxx infrastructure code";

  // Scream Session
  m.def("init",py::overload_cast<>(&initialize));
  m.def("init",py::overload_cast<pybind11::object>(&initialize));
  m.def("finalize",&finalize);

  // Call all other headers' registration routines
  pybind_pyparamlist(m);
  pybind_pyfield(m);
  pybind_pygrid(m);
  pybind_pyatmproc(m);
}

} // namespace scream
