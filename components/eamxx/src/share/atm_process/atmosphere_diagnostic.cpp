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
Field AtmosphereDiagnostic::get_diagnostic () const {
  EKAT_REQUIRE_MSG (m_diagnostic_output.is_allocated(),
      "Error! Getting a diagnostic field before it is allocated is suspicious at best.\n"
      "       We chose to throw an error, but if this is a legit use, please, contact developers.\n");
  return m_diagnostic_output.get_const();
}

void AtmosphereDiagnostic::compute_diagnostic (const double dt) {
  // Some diagnostics need the timestep, store in case.
  m_dt = dt;

  compute_diagnostic_impl ();

  // Set the timestamp of the diagnostic to the most
  // recent timestamp among the inputs
  const auto& inputs = get_fields_in();
  util::TimeStamp ts;
  for (const auto& f : inputs) {
    const auto& fts = f.get_header().get_tracking().get_time_stamp();
    if (not ts.is_valid() || ts<fts) {
      ts = fts;
    }
  }

  // If all inputs have invalid timestamps, we have a problem.
  EKAT_REQUIRE_MSG (ts.is_valid(),
      "Error! All inputs to diagnostic have invalid timestamp.\n"
      "  - Diag name: " + name() + "\n");

  m_diagnostic_output.get_header().get_tracking().update_time_stamp(ts);
}


void AtmosphereDiagnostic::run_impl (const double dt) {
  compute_diagnostic(dt);
}
void AtmosphereDiagnostic::set_computed_field (const Field& /* f */) {
  EKAT_ERROR_MSG("Error! Diagnostics are not allowed to compute fields. See " + name() + ".\n");
}

void AtmosphereDiagnostic::set_computed_group (const FieldGroup& /* group */) {
  EKAT_ERROR_MSG("Error! Diagnostics are not allowed to compute field groups. See " + name() + ".\n");
}

} // namespace scream

