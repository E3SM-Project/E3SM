#ifndef SCREAM_ATMOSPHERE_PROCESS_UTILS_HPP
#define SCREAM_ATMOSPHERE_PROCESS_UTILS_HPP

#include "ekat/ekat_assert.hpp"

#include <string>

namespace scream {

enum class CheckFailHandling {
  Fatal,
  Warning
};

enum class AtmosphereProcessType {
  Dynamics,                // Process responsible of handling the dynamics
  Physics,                 // Process handling a physics parametrization
  SurfaceCouplingImporter, // Process handling the transfers from surface models to atm
  SurfaceCouplingExporter, // Process handling the transfers from atm to surface models
  Group,                   // Process that groups a bunch of processes (so they look as a single process)
  Diagnostic               // Process that handles a diagnostic output
};

inline std::string e2str (const AtmosphereProcessType ap_type) {
  switch (ap_type) {
    case AtmosphereProcessType::Dynamics:                return "Atmosphere Dynamics";
    case AtmosphereProcessType::Physics:                 return "Atmosphere Physics Parametrization";
    case AtmosphereProcessType::SurfaceCouplingImporter: return "Surface Coupling Importer";
    case AtmosphereProcessType::SurfaceCouplingExporter: return "Surface Coupling Exporter";
    case AtmosphereProcessType::Group:                   return "Atmosphere Process Group";
    case AtmosphereProcessType::Diagnostic:              return "Atmosphere Diagnostic";
    default:
      ekat::error::runtime_abort("Error! Unrecognized atmosphere process type.\n");
  }
  return "INVALID";
}

// This enum is mostly used by AtmosphereProcessGroup to establish whether
// its atm procs are to be run concurrently or sequentially.
// We put the enum here so other files can easily access it.
enum class ScheduleType {
  Sequential,
  Parallel
};

// Enum used for disinguishing between pre/postcondition
// property checks for output.
enum PropertyCheckCategory {
  Precondition,
  Postcondition
};

} // namespace scream

#endif // SCREAM_ATMOSPHERE_PROCESS_UTILS_HPP
