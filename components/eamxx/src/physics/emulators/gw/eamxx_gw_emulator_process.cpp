#include "eamxx_gw_emulator_process.hpp"

#include <ekat_assert.hpp>
#include <ekat_units.hpp>

namespace scream
{

GWEmulator::
GWEmulator(const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  // Nothing to do here
}

void GWEmulator::
set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  m_grid = grids_manager->get_grid("physics");

  auto 
}

void GWEmulator::initialize_impl (const RunType /* run_type */)
{

}

void GWEmulator::run_impl (const double dt)
{

}

void GWEmulator::finalize_impl ()
{

}

} // namespace scream
