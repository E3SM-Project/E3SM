#ifndef SCREAM_DUMMY_ATM_SETUP_HPP
#define SCREAM_DUMMY_ATM_SETUP_HPP

#include "share/atm_process/atmosphere_process.hpp"
#include "share/grid/user_provided_grids_manager.hpp"
#include "share/grid/point_grid.hpp"

#include "dummy_atm_proc.hpp"

namespace scream {

inline void dummy_atm_init (const int num_cols, const int nvl, const ekat::Comm& comm) {
  using namespace scream;

  // Need to register products in the factory *before* we create any AtmosphereProcessGroup,
  // which rely on factory for process creation. The initialize method of the AD does that.
  // While we're at it, check that the case insensitive key of the factory works.
  auto& proc_factory = AtmosphereProcessFactory::instance();
  proc_factory.register_product("Dummy",&create_atmosphere_process<DummyProcess>);

  // Need to register grids managers before we create the driver
  auto& gm_factory = GridsManagerFactory::instance();
  gm_factory.register_product("User Provided",create_user_provided_grids_manager);

  // Set the dummy grid in the UserProvidedGridManager
  // Recall that this class stores *static* members, so whatever
  // we set here, will be reflected in the GM built by the factory.
  UserProvidedGridsManager upgm;
  auto dummy_grid_a = create_point_grid("Point Grid",num_cols,nvl,comm);

  upgm.set_grid(dummy_grid_a);
  upgm.set_reference_grid("Point Grid");
}

inline void dummy_atm_cleanup () {
  UserProvidedGridsManager upgm;
  upgm.clean_up();

  AtmosphereProcessFactory::instance().clean_up();
  GridsManagerFactory::instance().clean_up();
}

} // namespace scream

#endif // SCREAM_DUMMY_ATM_SETUP_HPP
