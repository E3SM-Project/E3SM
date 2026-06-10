#include "eamxx_water_tracers_process_interface.hpp"

#include <ekat_assert.hpp>

namespace scream
{

// =========================================================================================
WaterTracers::WaterTracers(const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
  , m_tracer_count(0)
{
  // Read tracer count from parameter list (default 0)
  m_tracer_count = m_params.get<int>("tracer_count", 0);
}

// =========================================================================================
void WaterTracers::create_requests()
{
  // Get the grid for this process
  m_grid = m_grids_manager->get_grid("physics");
  m_num_cols = m_grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = m_grid->get_num_vertical_levels();  // Number of levels per column

  // TODO (spec 002): Define field requests for water tracer arrays
  // For now, this process requires no fields and computes no fields
  // Field definitions will be added in spec 002
}

// =========================================================================================
void WaterTracers::initialize_impl(const RunType /* run_type */)
{
  // TODO (spec 002+): Initialize tracer field arrays and any precomputed data
  // For now, this is a no-op since no fields are defined yet
}

// =========================================================================================
void WaterTracers::run_impl(const double /* dt */)
{
  // TODO (spec 002+): Implement tracer transport and physics
  // This stub implementation performs no operations
  // Tracer physics will be added in subsequent specs of the water isotope campaign
}

// =========================================================================================
void WaterTracers::finalize_impl()
{
  // Nothing to finalize at this stage
}

} // namespace scream
