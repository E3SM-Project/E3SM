#ifndef SCREAM_REGISTER_SURFACE_COUPLING_PROCESSES_HPP
#define SCREAM_REGISTER_SURFACE_COUPLING_PROCESSES_HPP

#include "share/atm_process/atmosphere_process.hpp"

// TODO: in the future, you may want to guard these headers,
//       in case we add more options for each parametrization,
//       and we want to only register the ones built,
//       without hardcoding all of them.

#include "control/atmosphere_surface_coupling_importer.hpp"
#include "control/atmosphere_surface_coupling_exporter.hpp"

namespace scream {

inline void register_surface_coupling () {
  auto& proc_factory = AtmosphereProcessFactory::instance();

  proc_factory.register_product("sc_import",&create_atmosphere_process<SurfaceCouplingImporter>);
  proc_factory.register_product("sc_export",&create_atmosphere_process<SurfaceCouplingExporter>);
}

} // namespace scream

#endif // SCREAM_REGISTER_SURFACE_COUPLING_PROCESSES_HPP
