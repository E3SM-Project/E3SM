#ifndef PYSCREAM_HPP
#define PYSCREAM_HPP

#include "physics/register_physics.hpp"
#include "diagnostics/register_diagnostics.hpp"
#include "dynamics/register_dynamics.hpp"

#include "share/grid/grids_manager.hpp"

#include <ekat/mpi/ekat_comm.hpp>

namespace scream {

struct PySession {
  static PySession& get () {
    static PySession s;
    return s;
  }

  ekat::Comm comm;
  std::shared_ptr<GridsManager> gm;
  bool inited = false;
private:

  PySession () = default;
};

} // namespace scream

#endif // PYSCREAM_HPP
