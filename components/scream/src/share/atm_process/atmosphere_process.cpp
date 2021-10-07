#include "share/atm_process/atmosphere_process.hpp"

#include "ekat/ekat_assert.hpp"

#include <set>
#include <stdexcept>

namespace scream
{

AtmosphereProcess::AtmosphereProcess (const ekat::Comm& comm, const ekat::ParameterList& params)
  : m_comm  (comm)
  , m_params(params)
{}

void AtmosphereProcess::initialize (const TimeStamp& t0) {
  set_fields_and_groups_pointers();
  m_time_stamp = t0;
  initialize_impl();
}

void AtmosphereProcess::run (const int dt) {
  // Make sure required fields are valid
  check_required_fields();

  // Let the derived class do the actual run
  run_impl(dt);

  // Make sure computed fields are valid
  check_computed_fields();

  // Update all output fields time stamps
  m_time_stamp += dt;
  update_time_stamps ();
}

void AtmosphereProcess::finalize (/* what inputs? */) {
  finalize_impl(/* what inputs? */);
}

void AtmosphereProcess::set_required_field (const Field<const Real>& f) {
  // Sanity check
  EKAT_REQUIRE_MSG (requires_field(f.get_header().get_identifier()),
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

void AtmosphereProcess::set_computed_field (const Field<Real>& f) {
  // Sanity check
  EKAT_REQUIRE_MSG (computes_field(f.get_header().get_identifier()),
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

void AtmosphereProcess::set_required_group (const FieldGroup<const Real>& group) {
  // Sanity check
  EKAT_REQUIRE_MSG (requires_group(group.m_info->m_group_name,group.grid_name()),
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

void AtmosphereProcess::set_computed_group (const FieldGroup<Real>& group) {
  // Sanity check
  EKAT_REQUIRE_MSG (computes_group(group.m_info->m_group_name,group.grid_name()),
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

void AtmosphereProcess::check_required_fields () const { 
  // AtmosphereProcessGroup is just a "container" of *real* atm processes,
  // so don't run checks here, and let the *real* atm process do the checks
  if (this->type()==AtmosphereProcessType::Group) {
    return;
  }

  // First run any process specific checks.
  // check_required_fields_impl();
  // Now run all field property checks on all fields
  for (const auto& field : m_fields_in) {
    EKAT_REQUIRE_MSG (field.get_header().get_tracking().get_time_stamp().is_valid(),
        "Error! Found an input field that is still not initialized.\n"
        "    field: " + field.get_header().get_identifier().name() + "\n"
        "    grid name: " + field.get_header().get_identifier().get_grid_name() + "\n"
        "    atm process: " + this->name() + "\n");
    for (const auto& pc : field.get_property_checks()) {
      EKAT_REQUIRE_MSG(pc.check(field),
         "Error: Input field field property check failed.\n"
         "   field: " + field.get_header().get_identifier().name() + "\n"
         "   grid name: " + field.get_header().get_identifier().get_grid_name() + "\n"
         "   property check: " + pc.name() + "\n"
         "   atm process: " + this->name() + "\n");
  }}
  for (const auto& group : m_groups_in) {
    if (group.m_bundle) {
      auto& field = *group.m_bundle;
      EKAT_REQUIRE_MSG (field.get_header().get_tracking().get_time_stamp().is_valid(),
          "Error! Found an input group bundled field that is still not initialized.\n"
          "    group name: " + group.m_info->m_group_name + "\n"
          "    grid name: " + group.grid_name() + "\n"
          "    atm process: " + this->name() + "\n");
      for (const auto& pc : field.get_property_checks()) {
        EKAT_REQUIRE_MSG(pc.check(field),
           "Error: Input group bundled field field property check failed.\n"
           "   group name: " + group.m_info->m_group_name + "\n"
           "   grid name: " + group.grid_name() + "\n"
           "   property check: " + pc.name() + "\n"
           "   atm process: " + this->name() + "\n");
      }
    }
    for (auto it : group.m_fields) {
      auto& field = *it.second;
      EKAT_REQUIRE_MSG (field.get_header().get_tracking().get_time_stamp().is_valid(),
          "Error! Found an input group containing a field that is still not initialized.\n"
          "    atm process: " + this->name() + "\n"
          "    group name: " + group.m_info->m_group_name + "\n"
          "    grid name: " + group.grid_name() + "\n"
          "    field name: " + field.get_header().get_identifier().name() + "\n");
      for (const auto& pc : field.get_property_checks()) {
        EKAT_REQUIRE_MSG(pc.check(field),
           "Error: Input group field field property check failed.\n"
           "   group name: " + group.m_info->m_group_name + "\n"
           "   grid name: " + group.grid_name() + "\n"
           "   field: " + field.get_header().get_identifier().name() + "\n"
           "   property check: " + pc.name() + "\n"
           "   atm process: " + this->name() + "\n");
      }
    }
  }
}

void AtmosphereProcess::check_computed_fields () {
  // AtmosphereProcessGroup is just a "container" of *real* atm processes,
  // so don't run checks here, and let the *real* atm process do the checks
  if (this->type()==AtmosphereProcessType::Group) {
    return;
  }
  // First run any process specific checks, so that derived class have a chance
  // to repair computed fields if desired/doable/appropriate.
  check_computed_fields_impl();
  // Now run all field property checks on all fields
  for (const auto& field : m_fields_out) {
    for (const auto& pc : field.get_property_checks()) {
      EKAT_REQUIRE_MSG(pc.check(field),
         "Error: Output field field property check failed.\n"
         "   field: " + field.get_header().get_identifier().name() + "\n"
         "   grid name: " + field.get_header().get_identifier().get_grid_name() + "\n"
         "   property check: " + pc.name() + "\n"
         "   atm process: " + this->name() + "\n");
  }}
  for (const auto& group : m_groups_out) {
    if (group.m_bundle) {
      auto& field = *group.m_bundle;
      for (const auto& pc : field.get_property_checks()) {
        EKAT_REQUIRE_MSG(pc.check(field),
           "Error: Output group bundled field field property check failed.\n"
           "   group name: " + group.m_info->m_group_name + "\n"
           "   grid name: " + group.grid_name() + "\n"
           "   property check: " + pc.name() + "\n"
           "   atm process: " + this->name() + "\n");
      }
    }
    for (auto it : group.m_fields) {
      auto& field = *it.second;
      for (const auto& pc : field.get_property_checks()) {
        EKAT_REQUIRE_MSG(pc.check(field),
           "Error: Output group field field property check failed.\n"
           "   group name: " + group.m_info->m_group_name + "\n"
           "   grid name: " + group.grid_name() + "\n"
           "   field: " + field.get_header().get_identifier().name() + "\n"
           "   property check: " + pc.name() + "\n"
           "   atm process: " + this->name() + "\n");
      }
    }
  }
}


bool AtmosphereProcess::requires_field (const FieldIdentifier& id) const {
  for (const auto& it : m_required_field_requests) {
    if (it.fid==id) {
      return true;
    }
  }
  return false;
}
bool AtmosphereProcess::computes_field (const FieldIdentifier& id) const {
  for (const auto& it : m_computed_field_requests) {
    if (it.fid==id) {
      return true;
    }
  }
  return false;
}

bool AtmosphereProcess::requires_group (const std::string& name, const std::string& grid) const {
  for (const auto& it : m_required_group_requests) {
    if (it.name==name && it.grid==grid) {
      return true;
    }
  }
  return false;
}
bool AtmosphereProcess::computes_group (const std::string& name, const std::string& grid) const {
  for (const auto& it : m_computed_group_requests) {
    if (it.name==name && it.grid==grid) {
      return true;
    }
  }
  return false;
}

void AtmosphereProcess::update_time_stamps () {
  const auto& t = timestamp();

  // Update *all* output fields, regardless of whether they were touched
  // at all during this time step.
  for (auto& f : m_fields_out) {
    f.get_header().get_tracking().update_time_stamp(t);
  }
}

void AtmosphereProcess::add_me_as_provider (const Field<Real>& f) {
  f.get_header_ptr()->get_tracking().add_provider(weak_from_this());
}

void AtmosphereProcess::add_me_as_customer (const Field<const Real>& f) {
  f.get_header_ptr()->get_tracking().add_customer(weak_from_this());
}

// Convenience function to retrieve input/output fields
Field<const Real>& AtmosphereProcess::
get_field_in(const std::string& field_name, const std::string& grid_name) {
  try {
    return *m_fields_in_pointers.at(field_name).at(grid_name);
  } catch (const std::out_of_range&) {
    // std::out_of_range would message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate input field in this atm proces.\n"
        "   atm proc name: " + this->name() + "\n"
        "   field name: " + field_name + "\n"
        "   grid name: " + grid_name + "\n");
  }
}

Field<const Real>& AtmosphereProcess::
get_field_in(const std::string& field_name) {
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
    // std::out_of_range would message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate input field in this atm proces.\n"
        "   atm proc name: " + this->name() + "\n"
        "   field name: " + field_name + "\n");
  }
}

Field<Real>& AtmosphereProcess::
get_field_out(const std::string& field_name, const std::string& grid_name) {
  try {
    return *m_fields_out_pointers.at(field_name).at(grid_name);
  } catch (const std::out_of_range&) {
    // std::out_of_range would message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate output field in this atm proces.\n"
        "   atm proc name: " + this->name() + "\n"
        "   field name: " + field_name + "\n"
        "   grid name: " + grid_name + "\n");
  }
}

Field<Real>& AtmosphereProcess::
get_field_out(const std::string& field_name) {
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
    // std::out_of_range would message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate output field in this atm proces.\n"
        "   atm proc name: " + this->name() + "\n"
        "   field name: " + field_name + "\n");
  }
}

FieldGroup<const Real>& AtmosphereProcess::
get_group_in(const std::string& group_name, const std::string& grid_name) {
  try {
    return *m_groups_in_pointers.at(group_name).at(grid_name);
  } catch (const std::out_of_range&) {
    // std::out_of_range would message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate input group in this atm proces.\n"
        "   atm proc name: " + this->name() + "\n"
        "   group name: " + group_name + "\n"
        "   grid name: " + grid_name + "\n");
  }
}

FieldGroup<const Real>& AtmosphereProcess::
get_group_in(const std::string& group_name) {
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
    // std::out_of_range would message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate input group in this atm proces.\n"
        "   atm proc name: " + this->name() + "\n"
        "   group name: " + group_name + "\n");
  }
}

FieldGroup<Real>& AtmosphereProcess::
get_group_out(const std::string& group_name, const std::string& grid_name) {
  try {
    return *m_groups_out_pointers.at(group_name).at(grid_name);
  } catch (const std::out_of_range&) {
    // std::out_of_range would message would not help detecting where
    // the exception originated, so print a more meaningful message.
    EKAT_ERROR_MSG (
        "Error! Could not locate output group in this atm proces.\n"
        "   atm proc name: " + this->name() + "\n"
        "   group name: " + group_name + "\n"
        "   grid name: " + grid_name + "\n");
  }
}

FieldGroup<Real>& AtmosphereProcess::
get_group_out(const std::string& group_name) {
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
}

void AtmosphereProcess::
alias_field_in (const std::string& field_name,
                const std::string& grid_name,
                const std::string& alias_name)
{
  try {
    m_fields_in_pointers[alias_name][grid_name] = m_fields_in_pointers.at(field_name).at(grid_name);
  } catch (const std::out_of_range&) {
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
    EKAT_ERROR_MSG (
        "Error! Could not locate output group for aliasing request.\n"
        "    - group name: " + group_name + "\n"
        "    - grid name:  " + grid_name + "\n");
  }
}

} // namespace scream
