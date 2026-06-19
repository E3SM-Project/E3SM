#include "eamxx_water_isotopes_process_interface.hpp"

namespace scream
{

// =========================================================================================
WaterIsotopes::WaterIsotopes(const ekat::Comm& comm, const ekat::ParameterList& params)
  : WaterTracers(comm, params)
{
  // Water isotopes will inherit all tracer handling from WaterTracers
  // No additional initialization needed at this stage
}

// =========================================================================================
void WaterIsotopes::run_impl(const double dt)
{
  // Call base class tracer physics (currently a no-op)
  WaterTracers::run_impl(dt);

  // TODO: Add fractionation physics here
}

} // namespace scream
