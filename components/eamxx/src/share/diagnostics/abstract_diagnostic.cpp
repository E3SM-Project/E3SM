#include "share/diagnostics/abstract_diagnostic.hpp"

#include "share/util/eamxx_timing.hpp"

#include <ekat_std_utils.hpp>

#include <algorithm>
#include <limits>
#include <sstream>

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

  // A diag with several outputs declares them all in its constructor, so that
  // its outputs are known right after construction (see the class documentation)
  const auto num_ctor_outputs = m_outputs.size();

  initialize_impl();

  EKAT_REQUIRE_MSG (not m_outputs.empty(),
      "Error! A diagnostic did not declare its output field (see add_output).\n"
      " - diag class: " + name() + "\n");
  for (const auto& [n,f] : m_outputs) {
    EKAT_REQUIRE_MSG (f.is_allocated(),
        "Error! A diagnostic did not allocate its output field during initialization.\n"
        " - diag class: " + name() + "\n"
        " - diag output: " + f.name() + "\n");
  }
  EKAT_REQUIRE_MSG (m_outputs.size()==1 or m_outputs.size()==num_ctor_outputs,
      "Error! A diagnostic with several outputs must declare them all in its constructor,\n"
      "       so that customers can share it among requests of any of its outputs.\n"
      " - diag class: " + name() + "\n"
      " - diag outputs: " + outputs_names() + "\n");

  // A diag names its output field from its params, e.g. BinaryOp builds
  // "<arg1>_<op>_<arg2>". That is not always what the customer asked for, and
  // customers look the field up by name: an expression diag is looked up by its
  // expression. So honor an explicit name when given one.
  // NOTE: an alias shares tracking, alloc props and extra data, so masks,
  //       timestamps and valid-count bookkeeping carry over.
  if (m_params.isParameter("output_field_name")) {
    const auto& name = m_params.get<std::string>("output_field_name");
    EKAT_REQUIRE_MSG (m_outputs.size()<=1,
        "Error! Cannot rename the outputs of a diagnostic with several outputs.\n"
        " - diag class: " + this->name() + "\n"
        " - diag outputs: " + outputs_names() + "\n"
        " - requested name: " + name + "\n"
        "Request an output by its own name. To write it under another name, use\n"
        "  <name> := <output>\n");
    auto f = m_outputs.begin()->second;
    if (name!=f.name()) {
      m_outputs.clear();
      m_outputs.emplace(name,f.alias(name));
    }
  }

  m_is_initialized = true;
}

Field& AbstractDiagnostic::add_output (const Field& f)
{
  EKAT_REQUIRE_MSG (not m_is_initialized,
      "Error! Cannot add an output to a diagnostic after it is initialized.\n"
      " - diag class: " + name() + "\n"
      " - field name: " + f.name() + "\n");
  EKAT_REQUIRE_MSG (f.get_header().get_identifier().get_grid_name()==m_grid->name(),
      "Error! The output of a diagnostic must be on the diagnostic grid.\n"
      " - diag class: " + name() + "\n"
      " - field name: " + f.name() + "\n"
      " - field grid: " + f.get_header().get_identifier().get_grid_name() + "\n"
      " - diag grid: " + m_grid->name() + "\n");
  auto [it,inserted] = m_outputs.emplace(f.name(),f);
  EKAT_REQUIRE_MSG (inserted,
      "Error! A diagnostic cannot have two outputs with the same name.\n"
      " - diag class: " + name() + "\n"
      " - field name: " + f.name() + "\n"
      " - diag outputs: " + outputs_names() + "\n");
  return it->second;
}

Field& AbstractDiagnostic::add_output (const FieldIdentifier& fid)
{
  return add_output(Field(fid,true));
}

Field& AbstractDiagnostic::add_output (const std::string& name, const FieldLayout& layout,
                                       const ekat::units::Units& units, const DataType dt)
{
  return add_output(FieldIdentifier(name,layout,units,m_grid->name(),dt));
}

std::string AbstractDiagnostic::hash () const
{
  std::stringstream ss;
  ss << name() << "\n";

  // The first line printed is the name of the list
  std::stringstream params;
  // The default precision cannot tell apart two floating point params that are close
  params.precision(std::numeric_limits<double>::max_digits10);
  m_params.print(params);
  const auto str = params.str();
  // NOTE: print shows only some types. Equal prints must not hide different params,
  //       so refuse a diag with a param that print cannot show
  EKAT_REQUIRE_MSG (str.find("Cannot print type")==std::string::npos,
      "Error! Cannot hash a diagnostic with a parameter that cannot be printed.\n"
      " - diag: " + name() + "\n" + str);
  ss << str.substr(str.find('\n')+1);

  for (const auto& n : m_field_in_names) {
    const auto it = m_fields_in.find(n);
    EKAT_REQUIRE_MSG (it!=m_fields_in.end(),
        "Error! Cannot hash a diagnostic before all its input fields are set.\n"
        " - diag: " + name() + "\n"
        " - missing input: " + n + "\n");
    // Not the field identifier (name, layout, ...), which two different fields may share (e.g., an
    // alias with the same name in two streams, bound to different fields). The address of the
    // alloc props is the same for aliases, whatever their names, and different for relatives:
    // a subfield has its own (e.g., horiz_winds, and its components U and V).
    ss << n << "=" << static_cast<const void*>(&it->second.get_header().get_alloc_properties()) << "\n";
  }
  return ss.str();
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
  EKAT_REQUIRE_MSG (m_outputs.size()==1,
      "Error! The diagnostic has not declared its output yet (see add_output).\n"
      " - diag class: " + name() + "\n");
  return m_outputs.begin()->second;
}

Field AbstractDiagnostic::get (const std::string& name) const
{
  auto it = m_outputs.find(name);
  EKAT_REQUIRE_MSG (it!=m_outputs.end(),
      "Error! The diagnostic has no output with the requested name.\n"
      " - diag class: " + this->name() + "\n"
      " - requested output: " + name + "\n"
      " - available outputs: " + outputs_names() + "\n");
  return it->second;
}

std::string AbstractDiagnostic::outputs_names () const
{
  if (m_outputs.empty()) {
    return "(none yet)";
  }
  return ekat::join(m_outputs,[](const auto& it) { return it.second.name(); },", ");
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
  for (auto& [n,f] : m_outputs) {
    f.get_header().get_tracking().update_time_stamp(ts);
  }
  m_last_eval_ts = ts;
  m_last_eval_ts_hash = tsh;
}

} // namespace scream
