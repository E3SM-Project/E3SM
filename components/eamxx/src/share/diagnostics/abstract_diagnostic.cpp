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

Field AbstractDiagnostic::get_diagnostic () const
{
  EKAT_REQUIRE_MSG (m_diagnostic_output.is_allocated(),
      "Error! Getting a diagnostic field before it is allocated is suspicious at best.\n"
      "       We chose to throw an error, but if this is a legit use, please, contact developers.\n"
      " - diag field name: " + m_diagnostic_output.name() + "\n");
  return m_diagnostic_output;
}

void AbstractDiagnostic::compute_diagnostic (const util::TimeStamp& ts)
{
  if (m_last_eval_ts.is_valid() and ts==m_last_eval_ts) {
    // No need to compute the diag again
    return;
  }

  m_last_eval_ts = ts;
  m_diagnostic_output.get_header().get_tracking().update_time_stamp(ts);

  // Note: call the impl method *after* setting the diag time stamp.
  // Some derived classes may "refuse" to compute the diag, due to some
  // inconsistency of data. In that case, they can reset the diag time stamp
  // to something invalid, which can be used by downstream classes to determine
  // if the diag has been successfully computed or not.
  compute_diagnostic_impl ();
}

} // namespace scream
