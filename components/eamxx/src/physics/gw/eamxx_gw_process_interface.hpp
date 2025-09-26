#ifndef SCREAM_GW_MICROPHYSICS_HPP
#define SCREAM_GW_MICROPHYSICS_HPP

#include "share/atm_process/atmosphere_process.hpp"
#include "physics/gw/gw_functions.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"

#include <ekat_parameter_list.hpp>

#include <string>

namespace scream
{
/*
 * The class responsible to handle the gravity wave physics
 *
 * The AD should store exactly ONE instance of this class stored
 * in its list of subcomponents (the AD should make sure of this).
 *
 * This is currently just a placeholder
*/

class GWMicrophysics : public AtmosphereProcess
{
 public:
  // Constructors
  GWMicrophysics (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The type of subcomponent
  AtmosphereProcessType type () const { return AtmosphereProcessType::Physics; }

  // The name of the subcomponent
  std::string name () const { return "gw"; }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

 protected:
  std::shared_ptr<const AbstractGrid> m_grid;
}; // class GWMicrophysics

} // namespace scream

#endif // SCREAM_GW_MICROPHYSICS_HPP
