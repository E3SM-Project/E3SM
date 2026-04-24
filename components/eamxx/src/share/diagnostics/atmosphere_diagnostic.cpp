#include "share/diagnostics/atmosphere_diagnostic.hpp"

#include <ekat_string_utils.hpp>

namespace scream
{

AtmosphereDiagnostic::
AtmosphereDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params)
  : m_comm(comm)
  , m_params(params)
  , m_dt(0)
{
  // Nothing to do here
}

void AtmosphereDiagnostic::set_grid (const std::shared_ptr<const AbstractGrid>& grid)
{
  EKAT_REQUIRE_MSG (grid!=nullptr,
      "Error! Input grid pointer is null.\n");
  m_grid = grid;
  create_requests();
}

void AtmosphereDiagnostic::initialize ()
{
  initialize_impl();
  m_is_initialized = true;
}

void AtmosphereDiagnostic::set_required_field (const Field& f)
{
  // Validate pack size if this field was requested with a specific pack size
  for (const auto& r : m_field_requests) {
    if (r.fid.name()==f.name() || (r.incomplete && r.fid.name()==f.name())) {
      const auto& fap = f.get_header().get_alloc_properties();
      EKAT_REQUIRE_MSG (fap.get_largest_pack_size()>=r.pack_size,
          "Error! Diagnostic input field cannot accommodate the needed pack size.\n"
          "  - diag name: " + name() + "\n"
          "  - diag field: " + m_diagnostic_output.name() + "\n"
          "  - input field: " + f.name() + "\n"
          "  - requested pack size: " + std::to_string(r.pack_size) + "\n"
          "  - field max pack size: " + std::to_string(fap.get_largest_pack_size()) + "\n");
      break;
    }
  }

  m_fields_in.push_back(f);
  m_fields_in_pointers[f.name()] = &m_fields_in.back();
}

const Field& AtmosphereDiagnostic::
get_field_in (const std::string& field_name) const
{
  auto it = m_fields_in_pointers.find(field_name);
  EKAT_REQUIRE_MSG (it!=m_fields_in_pointers.end(),
      "Error! Field '" + field_name + "' not found in diagnostic inputs.\n"
      "  - Diag name: " + name() + "\n");
  return *it->second;
}

Field& AtmosphereDiagnostic::
get_field_in (const std::string& field_name)
{
  auto it = m_fields_in_pointers.find(field_name);
  EKAT_REQUIRE_MSG (it!=m_fields_in_pointers.end(),
      "Error! Field '" + field_name + "' not found in diagnostic inputs.\n"
      "  - Diag name: " + name() + "\n");
  return *it->second;
}

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
  //       during initialization, so that every call to compute_diagnostic sees
  //       ts==m_last_eval_ts, and does nothing.
  const auto& inputs = m_fields_in;
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

} // namespace scream
