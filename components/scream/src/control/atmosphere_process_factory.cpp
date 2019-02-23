#include "control/atmosphere_process_factory.hpp"

#include "control/surface_coupling.hpp"
#include "control/atmosphere_process_group.hpp"
#include "dynamics/atmosphere_dynamics.hpp"

namespace scream {

AtmosphereProcess*
AtmosphereProcessFactory::
create (const std::string& name, const ParameterList& params) {
  AtmosphereProcess* ptr = nullptr;
  if (name=="Surface Coupling") {
    ptr = new SurfaceCoupling(params);
  } else if (name=="Process Group") {
    ptr = new AtmosphereProcessGroup(params);
  } else if (name=="Dynamics") {
    ptr = new AtmosphereDynamics(params);
  } else {
    error::runtime_abort("Error! Unknown atmosphere process name '" + name + "'.\n");
  }

  return ptr;
}

} // namespace scream
