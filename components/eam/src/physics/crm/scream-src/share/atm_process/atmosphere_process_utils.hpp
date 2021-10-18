#ifndef SCREAM_ATMOSPHERE_PROCESS_UTILS_HPP
#define SCREAM_ATMOSPHERE_PROCESS_UTILS_HPP

#include "ekat/ekat_assert.hpp"

#include <string>

namespace scream {

enum class AtmosphereProcessType {
  Coupling,   // Process responsible of interfacing with the component coupler
  Dynamics,   // Process responsible of handling the dynamics
  Physics,    // Process handling a physics parametrization
  Group       // Process that groups a bunch of processes (so they look as a single process)
};

inline std::string e2str (const AtmosphereProcessType ap_type) {
  switch (ap_type) {
    case AtmosphereProcessType::Coupling:  return "Surface Coupling";
    case AtmosphereProcessType::Dynamics:  return "Atmosphere Dynamics";
    case AtmosphereProcessType::Physics:   return "Atmosphere Physics Parametrization";
    case AtmosphereProcessType::Group:     return "Atmosphere Process Group";
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

} // namespace scream

#endif // SCREAM_ATMOSPHERE_PROCESS_UTILS_HPP
