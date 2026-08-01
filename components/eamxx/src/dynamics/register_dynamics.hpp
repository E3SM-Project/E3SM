#ifndef SCREAM_REGISTER_DYNAMICS_PROCESS_HPP
#define SCREAM_REGISTER_DYNAMICS_PROCESS_HPP

#include "share/atm_process/atmosphere_process.hpp"
#include "share/data_managers/model_init.hpp"

#ifdef EAMXX_HAS_HOMME
#include "homme/eamxx_homme_process_interface.hpp"
#include "homme/homme_grids_manager.hpp"
#endif

namespace scream {

inline void register_dynamics () {
  auto& proc_factory = AtmosphereProcessFactory::instance();
  auto& gm_factory   = GridsManagerFactory::instance();

#ifdef EAMXX_HAS_HOMME
  proc_factory.register_product("homme",&create_atmosphere_process<HommeDynamics>);

  gm_factory.register_product("homme",&create_homme_grids_manager);

  auto& mi_factory = ModelInitFactory::instance();
  mi_factory.register_product("homme",&create_model_init);
#endif
  (void) proc_factory;
  (void) gm_factory;
}

} // namespace scream

#endif // SCREAM_REGISTER_DYNAMICS_PROCESS_HPP
