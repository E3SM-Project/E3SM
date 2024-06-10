#include "scream_session.hpp"
#include "pyfield.hpp"
#include "pygrid.hpp"
#include "pyatmproc.hpp"
#include "pyparamlist.hpp"

#include <pybind11/pybind11.h>

namespace py = pybind11;
namespace scream {

void initialize () {
  ekat::Comm comm(MPI_COMM_WORLD);
  initialize_scream_session(comm.am_i_root());
  scorpio::init_subsystem(comm);
}
void finalize () {
  scorpio::finalize_subsystem();
  finalize_scream_session();
}

PYBIND11_MODULE (pyeamxx,m) {

  m.doc() = "Python interfaces to certain EAMxx infrastructure code";

  // Scream Session
  m.def("init",&initialize);
  m.def("finalize",&finalize);

  // Call all other headers' registration routines
  pybind_pyparamlist(m);
  pybind_pyfield(m);
  pybind_pygrid(m);
  pybind_pyatmproc(m);
}

} // namespace scream
