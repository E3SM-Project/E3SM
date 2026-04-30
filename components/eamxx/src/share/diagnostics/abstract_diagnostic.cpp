#include "share/diagnostics/abstract_diagnostic.hpp"

#include <ekat_std_utils.hpp>

namespace scream
{

AbstractDiagnostic::
AbstractDiagnostic (const ekat::Comm& comm,
                    const ekat::ParameterList& params,
                    const std::shared_ptr<const AbstractGrid>& grid)
 : m_comm(comm)
 , m_params(params)
 , m_grid(grid)
{
  EKAT_REQUIRE_MSG (grid, "[AbstractDiagnostic] Error! Invalid grid pointer.\n");
}

void AbstractDiagnostic::initialize ()
{
  initialize_impl();
  m_is_initialized = true;
}

void AbstractDiagnostic::set_input_field (const Field& f)
{
  // Safety check
  EKAT_REQUIRE_MSG(ekat::contains(m_field_in_names,f.name()),
      "Error! Setting a field in the diagnostic that was not requested.\n"
      " - diag name: " + name() + "\n"
      " - field name: " + f.name() + "\n");
  m_fields_in[f.name()] = f;
}

Field AbstractDiagnostic::get () const
{
  EKAT_REQUIRE_MSG (m_diagnostic_output.is_allocated(),
      "Error! Getting a diagnostic field before it is allocated is suspicious at best.\n"
      "       We chose to throw an error, but if this is a legit use, please, contact developers.\n"
      " - diag field name: " + m_diagnostic_output.name() + "\n");
  return m_diagnostic_output;
}

void AbstractDiagnostic::compute (const util::TimeStamp& ts)
{
  // Compute a hash of ts with all the timestamps of the input fields
  bfbhash::HashType tsh = 0;
  for (auto it : m_fields_in) {
    const auto& fts = it.second.get_header().get_tracking().get_time_stamp();
    util::hash(fts,tsh);
  }
  util::hash(ts,tsh);

  // If the hash matches the last evaluation hash, then nothing has really
  // changed, so the stored diagnostic field does not have to be recomputed
  if (tsh==m_last_eval_ts_hash) {
    return;
  }

  compute_impl ();

  // Update timestamp info
  m_diagnostic_output.get_header().get_tracking().update_time_stamp(ts);
  m_last_eval_ts = ts;
  m_last_eval_ts_hash = tsh;
}

} // namespace scream
