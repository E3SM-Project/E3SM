#include "share/diagnostics/abstract_diagnostic.hpp"

#include "share/util/eamxx_timing.hpp"

#include <ekat_std_utils.hpp>

#include <algorithm>

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
  EKAT_REQUIRE_MSG (not m_is_initialized,
      "Error! Diagnostic was already initialized.\n"
      " - diag class: " + name() + "\n");

  // A diag with several outputs declares them all in its constructor (see add_output)
  const auto num_ctor_outputs = m_outputs.size();

  initialize_impl();

  if (m_outputs.empty()) {
    EKAT_REQUIRE_MSG (m_diagnostic_output.is_allocated(),
        "Error! A diagnostic did not allocate its output field during initialization.\n"
        " - diag class: " + name() + "\n"
        " - diag output: " + m_diagnostic_output.name() + "\n");
  } else {
    EKAT_REQUIRE_MSG (not m_diagnostic_output.is_allocated(),
        "Error! A diagnostic must either set m_diagnostic_output or call add_output, not both.\n"
        " - diag class: " + name() + "\n"
        " - diag outputs: " + outputs_names() + "\n");
    EKAT_REQUIRE_MSG (m_outputs.size()==1 or m_outputs.size()==num_ctor_outputs,
        "Error! A diagnostic with several outputs must declare them all in its constructor,\n"
        "       so that customers can share it among requests of any of its outputs.\n"
        " - diag class: " + name() + "\n"
        " - diag outputs: " + outputs_names() + "\n");
  }

  // A diag names its output field from its params, e.g. BinaryOp builds
  // "<arg1>_<op>_<arg2>". That is not always what the customer asked for, and
  // customers look the field up by name: an expression diag is looked up by its
  // expression. So honor an explicit name when given one.
  // NOTE: an alias shares tracking, alloc props and extra data, so masks,
  //       timestamps and avg-cnt bookkeeping carry over.
  if (m_params.isParameter("output_field_name")) {
    const auto& name = m_params.get<std::string>("output_field_name");
    EKAT_REQUIRE_MSG (m_outputs.size()<=1,
        "Error! Cannot rename the outputs of a diagnostic with several outputs.\n"
        " - diag class: " + this->name() + "\n"
        " - diag outputs: " + outputs_names() + "\n"
        " - requested name: " + name + "\n"
        "Request an output by its own name. To write it under another name, use\n"
        "  <name> := <output>\n");
    auto& f = m_outputs.empty() ? m_diagnostic_output : m_outputs.front();
    if (name!=f.name()) {
      f = f.alias(name);
    }
  }

  if (m_outputs.empty()) {
    m_outputs.push_back(m_diagnostic_output);
  }

  // Outputs whose masks have the same name share an avg count in output streams,
  // and the horizontal remapper remaps such masks once: they must be the same mask.
  for (size_t i=0; i<m_outputs.size(); ++i) {
    for (size_t j=i+1; j<m_outputs.size(); ++j) {
      const auto& fi = m_outputs[i];
      const auto& fj = m_outputs[j];
      if (not fi.has_valid_mask() or not fj.has_valid_mask() or
          fi.get_valid_mask().name()!=fj.get_valid_mask().name()) {
        continue;
      }
      EKAT_REQUIRE_MSG (fi.get_valid_mask().is_aliasing(fj.get_valid_mask()) and
                        fi.get_header().get_identifier().get_layout()==
                        fj.get_header().get_identifier().get_layout(),
          "Error! Outputs with masks of the same name must share the mask, and have the same layout.\n"
          " - diag class: " + name() + "\n"
          " - outputs: " + fi.name() + ", " + fj.name() + "\n"
          " - mask name: " + fi.get_valid_mask().name() + "\n");
    }
  }

  m_is_initialized = true;
}

void AbstractDiagnostic::add_output (const Field& f)
{
  EKAT_REQUIRE_MSG (not m_is_initialized,
      "Error! Cannot add an output to a diagnostic after it is initialized.\n"
      " - diag class: " + name() + "\n"
      " - field name: " + f.name() + "\n");
  EKAT_REQUIRE_MSG (f.is_allocated(),
      "Error! The output of a diagnostic must be allocated before it is added.\n"
      " - diag class: " + name() + "\n"
      " - field name: " + f.name() + "\n");
  EKAT_REQUIRE_MSG (f.get_header().get_identifier().get_grid_name()==m_grid->name(),
      "Error! The output of a diagnostic must be on the diagnostic grid.\n"
      " - diag class: " + name() + "\n"
      " - field name: " + f.name() + "\n"
      " - field grid: " + f.get_header().get_identifier().get_grid_name() + "\n"
      " - diag grid: " + m_grid->name() + "\n");
  for (const auto& o : m_outputs) {
    EKAT_REQUIRE_MSG (not ekat::caseInsensitiveEqualString(o.name(),f.name()),
        "Error! A diagnostic cannot have two outputs with the same name.\n"
        " - diag class: " + name() + "\n"
        " - field name: " + f.name() + "\n"
        " - diag outputs: " + outputs_names() + "\n");
  }
  m_outputs.push_back(f);
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
  EKAT_REQUIRE_MSG (m_outputs.size()<=1,
      "Error! Cannot get the output of a diagnostic with several outputs.\n"
      " - diag class: " + name() + "\n"
      " - diag outputs: " + outputs_names() + "\n"
      "Use get(name) to get one of them.\n");
  // Before initialize(), a diag with one output may not have added it to m_outputs yet
  const auto& f = m_outputs.empty() ? m_diagnostic_output : m_outputs.front();
  EKAT_REQUIRE_MSG (f.is_allocated(),
      "Error! Getting a diagnostic field before it is allocated is suspicious at best.\n"
      "       We chose to throw an error, but if this is a legit use, please, contact developers.\n"
      " - diag class: " + name() + "\n"
      " - diag field name: " + f.name() + "\n");
  return f;
}

Field AbstractDiagnostic::get (const std::string& name) const
{
  // Before initialize(), a diag with one output may not have added it to m_outputs yet
  if (m_outputs.empty() and ekat::caseInsensitiveEqualString(m_diagnostic_output.name(),name)) {
    return get();
  }
  auto it = std::find_if(m_outputs.begin(),m_outputs.end(),
                         [&](const Field& f) { return ekat::caseInsensitiveEqualString(f.name(),name); });
  EKAT_REQUIRE_MSG (it!=m_outputs.end(),
      "Error! The diagnostic has no output with the requested name.\n"
      " - diag class: " + this->name() + "\n"
      " - requested output: " + name + "\n"
      " - available outputs: " + outputs_names() + "\n");
  return *it;
}

std::string AbstractDiagnostic::outputs_names () const
{
  if (m_outputs.empty()) {
    return m_diagnostic_output.is_allocated() ? m_diagnostic_output.name() : "(none yet)";
  }
  return ekat::join(m_outputs,[](const Field& f) { return f.name(); },", ");
}

void AbstractDiagnostic::compute (const util::TimeStamp& ts)
{
  EKAT_REQUIRE_MSG (m_is_initialized,
      "Error! Cannot compute a diagnostic before it is initialized.\n"
      " - diag class: " + name() + "\n");

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

  // One timer per diagnostic class (not per instance)
  const std::string timer_name = "EAMxx::diag::" + name();
  start_timer(timer_name);
  compute_impl ();
  stop_timer(timer_name);

  // Update timestamp info
  for (auto& f : m_outputs) {
    f.get_header().get_tracking().update_time_stamp(ts);
  }
  m_last_eval_ts = ts;
  m_last_eval_ts_hash = tsh;
}

} // namespace scream
