#ifndef SCREAM_DUMMY_ATM_SETUP_HPP
#define SCREAM_DUMMY_ATM_SETUP_HPP

#include "share/atm_process/atmosphere_process.hpp"
#include "share/grid/user_provided_grids_manager.hpp"

#include "dummy_atm_proc.hpp"
#include "dummy_grid.hpp"

namespace scream {

void dummy_atm_init (const int num_cols, const int nvl, const ekat::Comm& comm) {
  using namespace scream;
  using device_type = DefaultDevice;

  // Need to register products in the factory *before* we create any AtmosphereProcessGroup,
  // which rely on factory for process creation. The initialize method of the AD does that.
  // While we're at it, check that the case insensitive key of the factory works.
  auto& proc_factory = AtmosphereProcessFactory::instance();
  proc_factory.register_product("physics_fwd",&create_atmosphere_process<DummyProcess<device_type,2,true>>);
  proc_factory.register_product("physics_bwd",&create_atmosphere_process<DummyProcess<device_type,4,false>>);

  // Need to register grids managers before we create the driver
  auto& gm_factory = GridsManagerFactory::instance();
  gm_factory.register_product("User Provided",create_user_provided_grids_manager);

  // Set the dummy grid in the UserProvidedGridManager
  // Recall that this class stores *static* members, so whatever
  // we set here, will be reflected in the GM built by the factory.
  UserProvidedGridsManager upgm;
  auto dummy_grid_fwd = std::make_shared<SimpleGrid>("Physics_fwd",num_cols,nvl,comm);
  auto dummy_grid_bwd = std::make_shared<SimpleGrid>("Physics_bwd",num_cols,nvl,comm);

  upgm.set_grid(dummy_grid_fwd);
  upgm.set_grid(dummy_grid_bwd);
  upgm.set_reference_grid("Physics_fwd");
  using remapper_type = DummyPhysicsGridRemapper<Real,device_type>;
  upgm.set_remapper(std::make_shared<remapper_type>(dummy_grid_fwd,dummy_grid_bwd));
  upgm.set_remapper(std::make_shared<remapper_type>(dummy_grid_bwd,dummy_grid_fwd));
}

void dummy_atm_cleanup () {
  UserProvidedGridsManager upgm;
  upgm.clean_up();

  AtmosphereProcessFactory::instance().clean_up();
  GridsManagerFactory::instance().clean_up();
}

} // namespace scream

#endif // SCREAM_DUMMY_ATM_SETUP_HPP
