#include "share/atm_process/atmosphere_diagnostic.hpp"

namespace scream
{

AtmosphereDiagnostic::
AtmosphereDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm,params)
{
  // Nothing to do here
}

// Function to retrieve the diagnostic output which is stored in m_diagnostic_output
Field AtmosphereDiagnostic::get_diagnostic (const Real dt) {
  run(dt);
  return m_diagnostic_output.get_const();
}

void AtmosphereDiagnostic::set_computed_field (const Field& /* f */) {
  EKAT_ERROR_MSG("Error! Diagnostics are not allowed to compute fields. See " + name() + ".\n");
}

void AtmosphereDiagnostic::set_computed_group (const FieldGroup& /* group */) {
  EKAT_ERROR_MSG("Error! Diagnostics are not allowed to compute field groups. See " + name() + ".\n");
}

} // namespace scream

