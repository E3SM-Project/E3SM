#include <share/eamxx_session.hpp>
#include "pyfield.hpp"
#include "pygrid.hpp"
#include "pyatmproc.hpp"
#include "pyparamlist.hpp"
#include "pyutils.hpp"

#include <ekat/mpi/ekat_comm.hpp>

#include <nanobind/nanobind.h>

#include <mpi.h>

namespace nb = nanobind;

namespace scream {

void initialize (MPI_Comm mpi_comm) {
  ekat::Comm comm(mpi_comm);
  initialize_eamxx_session(comm.am_i_root());
  scorpio::init_subsystem(comm);

  // Register everything in the eamxx factories
  register_physics();
  register_dynamics();
  register_diagnostics();

  auto& s = PySession::get();
  s.comm = comm;
  s.inited = true;
}

void initialize () {
  initialize(MPI_COMM_WORLD);
}
void initialize (nb::object py_comm) {
  initialize(get_c_comm(py_comm));
}
void finalize () {
  auto& s = PySession::get();
  s.gm = nullptr;
  s.inited = false;

  scorpio::finalize_subsystem();
  finalize_eamxx_session();
}

NB_MODULE (pyeamxx_ext,m) {

  m.doc() = "Python interfaces to certain EAMxx infrastructure code";

  // Scream Session
  m.def("init",nb::overload_cast<>(&initialize));
  m.def("init",nb::overload_cast<nb::object>(&initialize));
  m.def("finalize",&finalize);

  // Call all other headers' registration routines
  nb_pyparamlist(m);
  nb_pyfield(m);
  nb_pygrid(m);
  nb_pyatmproc(m);
}

} // namespace scream
