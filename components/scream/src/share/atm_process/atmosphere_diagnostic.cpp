#include "share/atm_process/atmosphere_diagnostic.hpp"

namespace scream
{

AtmosphereDiagnostic::
AtmosphereDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm,params)
{
  // Nothing to do here
}

} // namespace scream

