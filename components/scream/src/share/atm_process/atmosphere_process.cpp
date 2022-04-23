#include "share/atm_process/atmosphere_process.hpp"

#include "ekat/ekat_assert.hpp"

#include <set>
#include <stdexcept>
#include <string>

namespace scream
{

AtmosphereProcess::
AtmosphereProcess (const ekat::Comm& comm, const ekat::ParameterList& params)
  : m_comm       (comm)
  , m_params     (params)
{
  if (m_params.isParameter("Logger")) {
    m_atm_logger = m_params.get<std::shared_ptr<logger_t>>("Logger");
  } else {
    // Create a console-only logger, that logs all ranks
    using namespace ekat::logger;
    using logger_impl_t = Logger<LogNoFile,LogAllRanks>;
    m_atm_logger = std::make_shared<logger_impl_t>("",LogLevel::trace,m_comm);
  }

  if (m_params.isParameter("Number of Subcycles")) {
    m_num_subcycles = m_params.get<int>("Number of Subcycles");
  }
  EKAT_REQUIRE_MSG (m_num_subcycles>0,
      "Error! Invalid number of subcycles in param list " + m_params.name() + ".\n"
      "  - Num subcycles: " + std::to_string(m_num_subcycles) + "\n");
}

void AtmosphereProcess::initialize (const TimeStamp& t0, const RunType run_type) {
  set_fields_and_groups_pointers();
  m_time_stamp = t0;
  initialize_impl(run_type);
}

void AtmosphereProcess::run (const int dt) {
  if (m_params.get("Enable Precondition Checks", true)) {
    // Run 'pre-condition' property checks stored in this AP
    run_precondition_checks();
  }

  EKAT_REQUIRE_MSG ( (dt % m_num_subcycles)==0,
      "Error! The number of subcycle iterations does not exactly divide the time step.\n"
      "  - Atm proc name: " + this->name() + "\n"
      "  - Num subcycles: " + std::to_string(m_num_subcycles) + "\n"
      "  - Time step    : " + std::to_string(dt) + "\n");

  // Let the derived class do the actual run
  auto dt_sub = dt / m_num_subcycles;
  for (m_subcycle_iter=0; m_subcycle_iter<m_num_subcycles; ++m_subcycle_iter) {
    run_impl(dt_sub);
  }

  if (m_params.get("Enable Postcondition Checks", true)) {
    // Run 'post-condition' property checks stored in this AP
    run_postcondition_checks();
  }

  m_time_stamp += dt;
  if (m_update_time_stamps) {
    // Update all output fields time stamps
    update_time_stamps ();
  }
}

void AtmosphereProcess::finalize (/* what inputs? */) {
  finalize_impl(/* what inputs? */);
}

void AtmosphereProcess::set_required_field (const Field& f) {
  // Sanity check
  EKAT_REQUIRE_MSG (has_required_field(f.get_header().get_identifier()),
    "Error! Input field is not required by this atm process.\n"
    "    field id: " + f.get_header().get_identifier().get_id_string() + "\n"
    "    atm process: " + this->name() + "\n"
    "Something is wrong up the call stack. Please, contact developers.\n");

  if (not ekat::contains(m_fields_in,f)) {
    m_fields_in.emplace_back(f);
  }

  // AtmosphereProcessGroup is just a "container" of *real* atm processes,
  // so don't add me as customer if I'm an atm proc group.
  if (this->type()!=AtmosphereProcessType::Group) {
    // Add myself as customer to the field
    add_me_as_customer(f);
  }

  set_required_field_impl (f);
}

void AtmosphereProcess::set_computed_field (const Field& f) {
  // Sanity check
  EKAT_REQUIRE_MSG (has_computed_field(f.get_header().get_identifier()),
    "Error! Input field is not computed by this atm process.\n"
    "   field id: " + f.get_header().get_identifier().get_id_string() + "\n"
    "   atm process: " + this->name() + "\n"
    "Something is wrong up the call stack. Please, contact developers.\n");

  if (not ekat::contains(m_fields_out,f)) {
    m_fields_out.emplace_back(f);
  }

  // AtmosphereProcessGroup is just a "container" of *real* atm processes,
  // so don't add me as provider if I'm an atm proc group.
  if (this->type()!=AtmosphereProcessType::Group) {
    // Add myself as provider for the field
    add_me_as_provider(f);
  }

  set_computed_field_impl (f);
}

void AtmosphereProcess::set_required_group (const FieldGroup& group) {
  // Sanity check
  EKAT_REQUIRE_MSG (has_required_group(group.m_info->m_group_name,group.grid_name()),
    "Error! This atmosphere process does not require the input group.\n"
    "   group name: " + group.m_info->m_group_name + "\n"
    "   grid name : " + group.grid_name() + "\n"
    "   atm process: " + this->name() + "\n"
    "Something is wrong up the call stack. Please, contact developers.\n");

  if (not ekat::contains(m_groups_in,group)) {
    m_groups_in.emplace_back(group);
    // AtmosphereProcessGroup is just a "container" of *real* atm processes,
    // so don't add me as customer if I'm an atm proc group.
    if (this->type()!=AtmosphereProcessType::Group) {
      if (group.m_bundle) {
        add_me_as_customer(*group.m_bundle);
      } else {
        for (auto& it : group.m_fields) {
          add_me_as_customer(*it.second);
        }
      }
    }
  }

  set_required_group_impl(group);
}

void AtmosphereProcess::set_computed_group (const FieldGroup& group) {
  // Sanity check
  EKAT_REQUIRE_MSG (has_computed_group(group.m_info->m_group_name,group.grid_name()),
    "Error! This atmosphere process does not compute the input group.\n"
    "   group name: " + group.m_info->m_group_name + "\n"
    "   grid name : " + group.grid_name() + "\n"
    "   atm process: " + this->name() + "\n"
    "Something is wrong up the call stack. Please, contact developers.\n");

  if (not ekat::contains(m_groups_out,group)) {
    m_groups_out.emplace_back(group);
    // AtmosphereProcessGroup is just a "container" of *real* atm processes,
    // so don't add me as provider if I'm an atm proc group.
    if (this->type()!=AtmosphereProcessType::Group) {
      if (group.m_bundle) {
        add_me_as_provider(*group.m_bundle);
      } else {
        for (auto& it : group.m_fields) {
          add_me_as_provider(*it.second);
        }
      }
    }
  }

  set_computed_group_impl(group);
}

void AtmosphereProcess::run_precondition_checks () const {
  // Run all pre-condition property checks
  for (const auto& it : m_precondition_checks) {
    const auto& pc = it.second;

    auto check_result = pc->check();
    if (check_result.pass) {
      continue;
    }

    // Failed check.
    if (pc->can_repair()) {
      // Ok, just fix it
      pc->repair();
    } else if (it.first==CheckFailHandling::Warning) {
      // Still ok, but warn the user
      log (LogLevel::warn,
        "WARNING: Pre-condition property check failed.\n"
        "  - Property check name: " + pc->name() + "\n"
        "  - Atmosphere process name: " + name() + "\n"
        "  - Atmosphere process MPI Rank: " + std::to_string(m_comm.rank()) + "\n"
        "  - Message: " + check_result.msg);
    } else {
      // No hope. Crash.
      EKAT_ERROR_MSG(
          "Error! Failed pre-condition check (cannot be repaired).\n"
          "  - Atm process name: " + name() + "\n"
          "  - Property check name: " + pc->name() + "\n"
          "  - Atmosphere process MPI Rank: " + std::to_string(m_comm.rank()) + "\n"
          "  - Message: " + check_result.msg);
    }
  }
}

void AtmosphereProcess::run_postcondition_checks () const {
  // Run all post-condition property checks
  for (const auto& it : m_postcondition_checks) {
    const auto& pc = it.second;

    auto check_result = pc->check();
    if (check_result.pass) {
      continue;
    }

    // Failed check.
    if (pc->can_repair()) {
      // Ok, just fix it
      pc->repair();
      std::cout << "WARNING: Post-condition property check failed and repaired.\n"
        "  - Property check name: " + pc->name() + "\n"
        "  - Atmosphere process name: " + name() + "\n"
        "  - Atmosphere process MPI Rank: " + std::to_string(m_comm.rank()) + "\n";
    } else if (it.first==CheckFailHandling::Warning) {
      // Still ok, but warn the user
      log (LogLevel::warn,
        "WARNING: Post-condition property check failed.\n"
        "  - Property check name: " + pc->name() + "\n"
        "  - Atmosphere process name: " + name() + "\n"
        "  - Atmosphere process MPI Rank: " + std::to_string(m_comm.rank()) + "\n"
        "  - Error message: " + check_result.msg);
    } else {
      // No hope. Crash.
      EKAT_ERROR_MSG(
          "Error! Failed post-condition check (cannot be repaired).\n"
          "  - Atm process name: " + name() + "\n"
          "  - Property check name: " + pc->name() + "\n"
          "  - Atmosphere process MPI Rank: " + std::to_string(m_comm.rank()) + "\n"
          "  - Error message: " + check_result.msg);
    }
  }
}

bool AtmosphereProcess::has_required_field (const FieldIdentifier& id) const {
  for (const auto& it : m_required_field_requests) {
    if (it.fid==id) {
      return true;
    }
  }
  return false;
}

bool AtmosphereProcess::has_computed_field (const FieldIdentifier& id) const {
  for (const auto& it : m_computed_field_requests) {
    if (it.fid==id) {
      return true;
    }
  }
  return false;
}

bool AtmosphereProcess::has_required_group (const std::string& name, const std::string& grid) const {
  for (const auto& it : m_required_group_requests) {
    if (it.name==name && it.grid==grid) {
      return true;
    }
  }
  return false;
}

bool AtmosphereProcess::has_computed_group (const std::string& name, const std::string& grid) const {
  for (const auto& it : m_computed_group_requests) {
    if (it.name==name && it.grid==grid) {
      return true;
    }
  }
  return false;
}

void AtmosphereProcess::log (const LogLevel lev, const std::string& msg) const {
  m_atm_logger->log(lev,msg);
}

void AtmosphereProcess::set_update_time_stamps (const bool do_update) {
  m_update_time_stamps = do_update;
}

void AtmosphereProcess::update_time_stamps () {
  const auto& t = timestamp();

  // Update *all* output fields/groups, regardless of whether
  // they were touched at all during this time step.
  // TODO: this might have to be changed
  for (auto& f : m_fields_out) {
    f.get_header().get_tracking().update_time_stamp(t);
  }
  for (auto& g : m_groups_out) {
    if (g.m_bundle) {
      g.m_bundle->get_header().get_tracking().update_time_stamp(t);
    } else {
      for (auto& f : g.m_fields) {
        f.second->get_header().get_tracking().update_time_stamp(t);
      }
    }
  }
}

void AtmosphereProcess::add_me_as_provider (const Field& f) {
  f.get_header_ptr()->get_tracking().add_provider(weak_from_this());
}

void AtmosphereProcess::add_me_as_customer (const Field& f) {
  f.get_header_ptr()->get_tracking().add_customer(weak_from_this());
}

void AtmosphereProcess::add_internal_field (const Field& f) {
  m_internal_fields.push_back(f);
}

const Field& AtmosphereProcess::
get_field_in(const std::string& field_name, const std::string& grid_name) const {
  return get_field_in_impl(field_name,grid_name);
}

Field& AtmosphereProcess::
get_field_in(const std::string& field_name, const std::string& grid_name) {
  return get_field_in_impl(field_name,grid_name);
}

const Field& AtmosphereProcess::
get_field_in(const std::string& field_name) const {
  return get_field_in_impl(field_name);
}

Field& AtmosphereProcess::
get_field_in(const std::string& field_name) {
  return get_field_in_impl(field_name);
}

const Field& AtmosphereProcess::
get_field_out(const std::string& field_name, const std::string& grid_name) const {
  return get_field_out_impl(field_name,grid_name);
}

Field& AtmosphereProcess::
get_field_out(const std::string& field_name, const std::string& grid_name) {
  return get_field_out_impl(field_name,grid_name);
}

const Field& AtmosphereProcess::
get_field_out(const std::string& field_name) const {
  return get_field_out_impl (field_name);
}

Field& AtmosphereProcess::
get_field_out(const std::string& field_name) {
  return get_field_out_impl (field_name);
}

const FieldGroup& AtmosphereProcess::
get_group_in(const std::string& group_name, const std::string& grid_name) const {
  return get_group_in_impl (group_name,grid_name);
}

FieldGroup& AtmosphereProcess::
get_group_in(const std::string& group_name, const std::string& grid_name) {
  return get_group_in_impl (group_name,grid_name);
}

const FieldGroup& AtmosphereProcess::
get_group_in(const std::string& group_name) const {
  return get_group_in_impl(group_name);
}

FieldGroup& AtmosphereProcess::
get_group_in(const std::string& group_name) {
  return get_group_in_impl(group_name);
}

const FieldGroup& AtmosphereProcess::
get_group_out(const std::string& group_name, const std::string& grid_name) const {
  return get_group_out_impl(group_name,grid_name);
}

FieldGroup& AtmosphereProcess::
get_group_out(const std::string& group_name, const std::string& grid_name) {
  return get_group_out_impl(group_name,grid_name);
}

const FieldGroup& AtmosphereProcess::
get_group_out(const std::string& group_name) const {
  return get_group_out_impl(group_name);
}

FieldGroup& AtmosphereProcess::
get_group_out(const std::string& group_name) {
  return get_group_out_impl(group_name);
}

const Field& AtmosphereProcess::
get_internal_field(const std::string& field_name, const std::string& grid_name) const {
  return get_internal_field_impl(field_name,grid_name);
}

Field& AtmosphereProcess::
get_internal_field(const std::string& field_name, const std::string& grid_name) {
  return get_internal_field_impl(field_name,grid_name);
}

Field& AtmosphereProcess::
get_internal_field(const std::string& field_name) {
  return get_internal_field_impl(field_name);
}

const Field& AtmosphereProcess::
get_internal_field(const std::string& field_name) const {
  return get_internal_field_impl(field_name);
}

void AtmosphereProcess::set_fields_and_groups_pointers () {
  for (auto& f : m_fields_in) {
    const auto& fid = f.get_header().get_identifier();
    m_fields_in_pointers[fid.name()][fid.get_grid_name()] = &f;
  }
  for (auto& f : m_fields_out) {
    const auto& fid = f.get_header().get_identifier();
    m_fields_out_pointers[fid.name()][fid.get_grid_name()] = &f;
  }
  for (auto& g : m_groups_in) {
    const auto& group_name = g.m_info->m_group_name;
    m_groups_in_pointers[group_name][g.grid_name()] = &g;
  }
  for (auto& g : m_groups_out) {
    const auto& group_name = g.m_info->m_group_name;
    m_groups_out_pointers[group_name][g.grid_name()] = &g;
  }
  for (auto& f : m_internal_fields) {
    const auto& fid = f.get_header().get_identifier();
    m_internal_fields_pointers[fid.name()][fid.get_grid_name()] = &f;
  }
}

void AtmosphereProcess::
alias_field_in (const std::string& field_name,
                const std::string& grid_name,
                const std::string& alias_name)
{
  try {
    m_fields_in_pointers[alias_name][grid_name] = m_fields_in_pointers.at(field_name).at(grid_name);
  } catch (const std::out_of_range&) {
    // std::out_of_range message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate input field for aliasing request.\n"
        "    - field name: " + field_name + "\n"
        "    - grid name:  " + grid_name + "\n");
  }
}

void AtmosphereProcess::
alias_field_out (const std::string& field_name,
                 const std::string& grid_name,
                 const std::string& alias_name)
{
  try {
    m_fields_out_pointers[alias_name][grid_name] = m_fields_out_pointers.at(field_name).at(grid_name);
  } catch (const std::out_of_range&) {
    // std::out_of_range message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate output field for aliasing request.\n"
        "    - field name: " + field_name + "\n"
        "    - grid name:  " + grid_name + "\n");
  }
}

void AtmosphereProcess::
alias_group_in (const std::string& group_name,
                const std::string& grid_name,
                const std::string& alias_name)
{
  try {
    m_groups_in_pointers[alias_name][grid_name] = m_groups_in_pointers.at(group_name).at(grid_name);
  } catch (const std::out_of_range&) {
    // std::out_of_range message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate input group for aliasing request.\n"
        "    - group name: " + group_name + "\n"
        "    - grid name:  " + grid_name + "\n");
  }
}

void AtmosphereProcess::
alias_group_out (const std::string& group_name,
                 const std::string& grid_name,
                 const std::string& alias_name)
{
  try {
    m_groups_out_pointers[alias_name][grid_name] = m_groups_out_pointers.at(group_name).at(grid_name);
  } catch (const std::out_of_range&) {
    // std::out_of_range message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate output group for aliasing request.\n"
        "    - group name: " + group_name + "\n"
        "    - grid name:  " + grid_name + "\n");
  }
}

Field& AtmosphereProcess::
get_field_in_impl(const std::string& field_name, const std::string& grid_name) const {
  try {
    return *m_fields_in_pointers.at(field_name).at(grid_name);
  } catch (const std::out_of_range&) {
    // std::out_of_range message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate input field in this atm proces.\n"
        "   atm proc name: " + this->name() + "\n"
        "   field name: " + field_name + "\n"
        "   grid name: " + grid_name + "\n");
  }
}

Field& AtmosphereProcess::
get_field_in_impl(const std::string& field_name) const {
  try {
    auto& copies = m_fields_in_pointers.at(field_name);
    EKAT_REQUIRE_MSG (copies.size()==1,
        "Error! Attempt to find input field providing only the name,\n"
        "       but multiple copies (on different grids) are present.\n"
        "  field name: " + field_name + "\n"
        "  atm process: " + this->name() + "\n"
        "  number of copies: " + std::to_string(copies.size()) + "\n");
    return *copies.begin()->second;
  } catch (const std::out_of_range&) {
    // std::out_of_range message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate input field in this atm proces.\n"
        "   atm proc name: " + this->name() + "\n"
        "   field name: " + field_name + "\n");
  }
}

Field& AtmosphereProcess::
get_field_out_impl(const std::string& field_name, const std::string& grid_name) const {
  try {
    return *m_fields_out_pointers.at(field_name).at(grid_name);
  } catch (const std::out_of_range&) {
    // std::out_of_range message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate output field in this atm proces.\n"
        "   atm proc name: " + this->name() + "\n"
        "   field name: " + field_name + "\n"
        "   grid name: " + grid_name + "\n");
  }
}

Field& AtmosphereProcess::
get_field_out_impl(const std::string& field_name) const {
  try {
    auto& copies = m_fields_out_pointers.at(field_name);
    EKAT_REQUIRE_MSG (copies.size()==1,
        "Error! Attempt to find output field providing only the name,\n"
        "       but multiple copies (on different grids) are present.\n"
        "  field name: " + field_name + "\n"
        "  atm process: " + this->name() + "\n"
        "  number of copies: " + std::to_string(copies.size()) + "\n");
    return *copies.begin()->second;
  } catch (const std::out_of_range&) {
    // std::out_of_range message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate output field in this atm proces.\n"
        "   atm proc name: " + this->name() + "\n"
        "   field name: " + field_name + "\n");
  }
}

FieldGroup& AtmosphereProcess::
get_group_in_impl(const std::string& group_name, const std::string& grid_name) const {
  try {
    return *m_groups_in_pointers.at(group_name).at(grid_name);
  } catch (const std::out_of_range&) {
    // std::out_of_range message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate input group in this atm proces.\n"
        "   atm proc name: " + this->name() + "\n"
        "   group name: " + group_name + "\n"
        "   grid name: " + grid_name + "\n");
  }
}

FieldGroup& AtmosphereProcess::
get_group_in_impl(const std::string& group_name) const {
  try {
    auto& copies = m_groups_in_pointers.at(group_name);
    EKAT_REQUIRE_MSG (copies.size()==1,
        "Error! Attempt to find input group providing only the name,\n"
        "       but multiple copies (on different grids) are present.\n"
        "  group name: " + group_name + "\n"
        "  atm process: " + this->name() + "\n"
        "  number of copies: " + std::to_string(copies.size()) + "\n");
    return *copies.begin()->second;
  } catch (const std::out_of_range&) {
    // std::out_of_range message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate input group in this atm proces.\n"
        "   atm proc name: " + this->name() + "\n"
        "   group name: " + group_name + "\n");
  }
}

FieldGroup& AtmosphereProcess::
get_group_out_impl(const std::string& group_name, const std::string& grid_name) const {
  try {
    return *m_groups_out_pointers.at(group_name).at(grid_name);
  } catch (const std::out_of_range&) {
    // std::out_of_range message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate output group in this atm proces.\n"
        "   atm proc name: " + this->name() + "\n"
        "   group name: " + group_name + "\n"
        "   grid name: " + grid_name + "\n");
  }
}

FieldGroup& AtmosphereProcess::
get_group_out_impl(const std::string& group_name) const {
  try {
    auto& copies = m_groups_out_pointers.at(group_name);
    EKAT_REQUIRE_MSG (copies.size()==1,
        "Error! Attempt to find output group providing only the name,\n"
        "       but multiple copies (on different grids) are present.\n"
        "  group name: " + group_name + "\n"
        "  atm process: " + this->name() + "\n"
        "  number of copies: " + std::to_string(copies.size()) + "\n");
    return *copies.begin()->second;
  } catch (const std::out_of_range&) {
    // std::out_of_range would message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate output group in this atm proces.\n"
        "   atm proc name: " + this->name() + "\n"
        "   group name: " + group_name + "\n");
  }
}

Field& AtmosphereProcess::
get_internal_field_impl(const std::string& field_name, const std::string& grid_name) const {
  try {
    return *m_internal_fields_pointers.at(field_name).at(grid_name);
  } catch (const std::out_of_range&) {
    // std::out_of_range message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate internal field in this atm proces.\n"
        "   atm proc name: " + this->name() + "\n"
        "   field name: " + field_name + "\n"
        "   grid name: " + grid_name + "\n");
  }
}

Field& AtmosphereProcess::
get_internal_field_impl(const std::string& field_name) const {
  try {
    auto& copies = m_internal_fields_pointers.at(field_name);
    EKAT_REQUIRE_MSG (copies.size()==1,
        "Error! Attempt to find internal field providing only the name,\n"
        "       but multiple copies (on different grids) are present.\n"
        "  field name: " + field_name + "\n"
        "  atm process: " + this->name() + "\n"
        "  number of copies: " + std::to_string(copies.size()) + "\n");
    return *copies.begin()->second;
  } catch (const std::out_of_range&) {
    // std::out_of_range message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate internal field in this atm proces.\n"
        "   atm proc name: " + this->name() + "\n"
        "   field name: " + field_name + "\n");
  }
}

} // namespace scream
