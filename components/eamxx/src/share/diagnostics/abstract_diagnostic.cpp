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

void AbstractDiagnostic::add_output (const Field& f)
{
  EKAT_REQUIRE_MSG (f.is_allocated(),
      "Error! Declaring a diagnostic output that is not allocated.\n"
      " - diag name : " + name() + "\n"
      " - field name: " + f.name() + "\n");
  EKAT_REQUIRE_MSG (not has_output(f.name()),
      "Error! Declaring two diagnostic outputs with the same name.\n"
      " - diag name : " + name() + "\n"
      " - field name: " + f.name() + "\n");
  m_extra_outputs.push_back(f);
}

const Field* AbstractDiagnostic::find_output (const std::string& n) const
{
  // An unset m_diagnostic_output is a default-constructed Field, whose name is
  // a placeholder rather than something a caller could be asking for.
  if (m_diagnostic_output.is_allocated() and m_diagnostic_output.name()==n) {
    return &m_diagnostic_output;
  }
  for (const auto& f : m_extra_outputs) {
    if (f.name()==n) {
      return &f;
    }
  }
  return nullptr;
}

bool AbstractDiagnostic::has_output (const std::string& n) const
{
  return find_output(n)!=nullptr;
}

std::vector<Field> AbstractDiagnostic::get_outputs () const
{
  std::vector<Field> outputs;
  if (m_diagnostic_output.is_allocated()) {
    outputs.push_back(m_diagnostic_output);
  }
  outputs.insert(outputs.end(),m_extra_outputs.begin(),m_extra_outputs.end());
  return outputs;
}

std::vector<std::string> AbstractDiagnostic::get_output_names () const
{
  std::vector<std::string> names;
  for (const auto& f : get_outputs()) {
    names.push_back(f.name());
  }
  return names;
}

Field AbstractDiagnostic::get () const
{
  // A diag computing several fields need not single one out as the primary, in
  // which case the first one it declared stands in for it.
  if (not m_diagnostic_output.is_allocated() and not m_extra_outputs.empty()) {
    return m_extra_outputs.front();
  }

  EKAT_REQUIRE_MSG (m_diagnostic_output.is_allocated(),
      "Error! Getting a diagnostic field before it is allocated is suspicious at best.\n"
      "       We chose to throw an error, but if this is a legit use, please, contact developers.\n"
      " - diag field name: " + m_diagnostic_output.name() + "\n");
  return m_diagnostic_output;
}

Field AbstractDiagnostic::get (const std::string& n) const
{
  const auto* f = find_output(n);
  EKAT_REQUIRE_MSG (f!=nullptr,
      "Error! This diagnostic does not compute a field by this name.\n"
      " - diag name    : " + name() + "\n"
      " - field name   : " + n + "\n"
      " - diag computes: " + ekat::join(get_output_names(),", ") + "\n");
  return *f;
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

  // Update timestamp info. One compute_impl() call produces all the outputs,
  // so they are all as fresh as the timestamp we were given.
  for (auto f : get_outputs()) {
    f.get_header().get_tracking().update_time_stamp(ts);
  }
  m_last_eval_ts = ts;
  m_last_eval_ts_hash = tsh;
}

DiagOutputsRegistry& DiagOutputsRegistry::instance ()
{
  static DiagOutputsRegistry r;
  return r;
}

void DiagOutputsRegistry::
register_outputs (const std::string& factory_key,
                  const std::vector<std::string>& outputs)
{
  for (const auto& o : outputs) {
    auto it = m_provider.find(o);
    EKAT_REQUIRE_MSG (it==m_provider.end() or it->second==factory_key,
        "Error! Two diagnostics claim to compute the same output field.\n"
        " - field name: " + o + "\n"
        " - diag 1    : " + it->second + "\n"
        " - diag 2    : " + factory_key + "\n");
    m_provider[o] = factory_key;
  }
}

std::string DiagOutputsRegistry::
provider_of (const std::string& output_name) const
{
  auto it = m_provider.find(output_name);
  return it==m_provider.end() ? "" : it->second;
}

} // namespace scream
