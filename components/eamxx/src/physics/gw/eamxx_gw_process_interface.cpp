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
void GWMicrophysics::create_requests()
{
  m_grid = m_grids_manager->get_grid("physics");
}

// =========================================================================================
} // namespace scream
