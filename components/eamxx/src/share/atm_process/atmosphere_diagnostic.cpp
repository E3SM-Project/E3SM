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
  return m_diagnostic_output;
}

void AtmosphereDiagnostic::compute_diagnostic (const double dt) {
  // Some diagnostics need the timestep, store in case.
  m_dt = dt;

  // Set the timestamp of the diagnostic to the most recent timestamp among the inputs
  // NOTE: it's a corner case, but it can happen that the diag has NO input fields.
  //       This can happen for BinaryOpsDiag, in the case where both args are physics
  //       constants. In that case, we assume the diagnostic can be pre-computed
  //       during initialiation, so that every call to compute_diagnostic sees
  //       ts==m_last_eval_ts, and does nothing.
  const auto& inputs = get_fields_in();
  util::TimeStamp ts = m_last_eval_ts;
  if (inputs.size()>0) {
    for (const auto& f : inputs) {
      const auto& fts = f.get_header().get_tracking().get_time_stamp();
      if (not ts.is_valid() || ts<fts) {
        ts = fts;
      }
    }

    // If all inputs have invalid timestamps, we have a problem.
    auto fname = [](const Field& f) { return f.name(); };
    EKAT_REQUIRE_MSG (ts.is_valid(),
        "Error! All inputs to diagnostic have invalid timestamp.\n"
        "  - Diag name: " + name() + "\n"
        "  - Diag field name: " + m_diagnostic_output.name() + "\n"
        "  - Diag inputs names: " + ekat::join(inputs,fname,",") + "\n");

  }

  if (ts==m_last_eval_ts) {
    // No need to compute the diag again
    return;
  }

  m_diagnostic_output.get_header().get_tracking().update_time_stamp(ts);

  // Note: call the impl method *after* setting the diag time stamp.
  // Some derived classes may "refuse" to compute the diag, due to some
  // inconsistency of data. In that case, they can reset the diag time stamp
  // to something invalid, which can be used by downstream classes to determine
  // if the diag has been successfully computed or not.
  compute_diagnostic_impl ();

  m_last_eval_ts = ts;
}

void AtmosphereDiagnostic::run_impl (const double dt) {
  compute_diagnostic(dt);
}

void AtmosphereDiagnostic::set_field_impl (const Field& f)
{
  // Check that the field is not writable, as we diags do not "compute" fields
  // in the FieldRequest sense..
  // Check that the field has the pack size that was requested
  // TODO: I don't think diagnostics should "request" a pack size.
  //       Diags should work with whatever the AD is storing.
  //       That's b/c the field is already allocated by the time
  //       we create any diagnostic.
  //       While we fix all diags, this method will at least
  //       throw an error if the pack size that the diag "requested"
  //       is not compatible with the field alloc props.
  for (const auto& r : get_field_requests()) {
    if (r.fid.name()==f.name()) {
      EKAT_REQUIRE_MSG(not (r.usage & Computed),
          "Error! Diagnostics are not allowed to compute fields.\n"
          " - diag name: " + name() + ".\n");
      const auto& fap = f.get_header().get_alloc_properties();
      EKAT_REQUIRE_MSG (fap.get_largest_pack_size()>=r.pack_size,
          "Error! Diagnostic input field cannot accommodate the needed pack size.\n"
          "  - diag name: " + name() + "\n"
          "  - input field: " + f.name() + "\n"
          "  - requested pack size: " + std::to_string(r.pack_size) + "\n"
          "  - field max pack size: " + std::to_string(fap.get_largest_pack_size()) + "\n");
      break;
    }
  }
}

void AtmosphereDiagnostic::set_group_impl (const FieldGroup& group)
{
  // Check this is NOT a computed group
  for (const auto& r : get_group_requests()) {
    if (r.name==group.m_info->m_group_name and r.grid==group.grid_name()) {
      EKAT_REQUIRE_MSG(not (r.usage & Computed),
          "Error! Diagnostics are not allowed to compute fields.\n"
          " - diag name: " + name() + ".\n");
      break;
    }
  }
}

} // namespace scream

