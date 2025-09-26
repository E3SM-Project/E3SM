#include "gw_functions.hpp"
#include "eamxx_gw_process_interface.hpp"

#include <ekat_assert.hpp>
#include <ekat_units.hpp>

#include <array>

namespace scream
{

// =========================================================================================
GWMicrophysics::GWMicrophysics(const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  // Nothing to do here
}

// =========================================================================================
void GWMicrophysics::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  m_grid = grids_manager->get_grid("physics");
}

// =========================================================================================
} // namespace scream
