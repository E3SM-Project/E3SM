#include "share/atm_process/atmosphere_diagnostic.hpp"

namespace scream
{

AtmosphereDiagnostic::
AtmosphereDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm,params)
{
  // Nothing to do here
}


void AtmosphereDiagnostic::set_computed_field_impl (const Field& f) {
  EKAT_ERROR_MSG("Error! Diagnostics are not allowed to compute fields. See " + name() + ".\n");
}

void AtmosphereDiagnostic::set_computed_group_impl (const FieldGroup& group) {
  EKAT_ERROR_MSG("Error! Diagnostics are not allowed to compute field groups. See " + name() + ".\n");
}

} // namespace scream

