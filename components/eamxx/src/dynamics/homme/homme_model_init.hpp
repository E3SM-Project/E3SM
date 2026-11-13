#ifndef EAMXX_HOMME_MODEL_INIT_HPP
#define EAMXX_HOMME_MODEL_INIT_HPP

#include "share/data_managers/model_init.hpp"

namespace scream
{

// A ModelInit for runs where Homme is the dynamics package.
//
// On top of the generic init of the atm input fields (done by ModelInit), this class takes
// care of initializing Homme itself, so that, by the time run() returns, both the atm fields
// and Homme's internal state are fully valid:
//  - complete Homme's init1 phase, and set the views of the dyn-grid fields into Homme's C++
//    data structures (see HommeInstance);
//  - for an initial run, init Homme's dyn state from the GLL initial conditions; for a
//    restart run, finalize the restart of the dyn state (read from the restart file by the
//    base class) and restore Homme's time step counter;
//  - complete Homme's model init;
//  - if the physics grid is a FV one (pg2), remap Homme's state to it, to init the physics
//    state (for a restart, the pg2 state is in the restart file, so this is only used to
//    compute the T and uv at the end of the last Homme step).
//
// The initial conditions files are always on the GLL grid. Therefore, for a pg2 run, the
// state variables that Homme owns on the physics grid (horiz_winds, T_mid, ps, phis, and the
// tracers) are not inited on the pg2 grid. Instead, we create temporary copies of them on the
// physics_gll grid (in a separate field manager), init them there as if the physics grid were
// GLL (constants, IC file, topography, IOP, perturbation), and use them to init Homme's dyn
// state. The pg2 fields are then obtained from the dyn state, via Homme's own remap.
// This way, we never go GLL->pg2->GLL, and the pg2 fields are valid once run() is over.
//
// This class requires the HommeInstance to have been setup by the HommeDynamics atm process
// (during its create_requests and set_computed_group_impl calls), and that the atm processes
// buffers have been allocated.
class HommeModelInit : public ModelInit {
public:
  using ModelInit::ModelInit;

  void run (const std::shared_ptr<FieldManager>& fm,
            const util::TimeStamp& t0,
            const RunType run_type) override;

protected:

  // For a pg2 run: nothing on physics_gll is part of a restart file; on the pg2 grid, the state
  // is inited from the staging field manager. Anything else on pg2 is inited only if it has
  // a constant, or it is a topography field (the IC file has no data on the pg2 grid).
  std::vector<Field>
  get_fields (const std::shared_ptr<FieldManager>& fm,
              const std::string& group_name,
              const std::string& grid_name) override;

  // For a pg2 run, phis is computed by remapping the dyn state, not read
  strmap_t<std::string>
  get_topography_file_names (const std::string& grid_name) const override;

private:

  // Create a field manager with temporary copies, on the physics_gll grid, of the state fields
  // that Homme owns on the pg2 grid.
  std::shared_ptr<FieldManager>
  create_gll_ic_field_manager (const std::shared_ptr<FieldManager>& fm) const;

  // N in pgN, and the name of the pgN grid. N<0 if the physics grid is not a FV one.
  int         m_pgN = -1;
  std::string m_pg_grid;
};

} // namespace scream

#endif // EAMXX_HOMME_MODEL_INIT_HPP
