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

void AtmosphereProcess::initialize (const TimeStamp& t0, const RunType run_type) {
  set_fields_and_groups_pointers();
  m_time_stamp = t0;
  initialize_impl(run_type);
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

void AtmosphereProcess::set_computed_field (const Field<Real>& f) {
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

void AtmosphereProcess::set_required_group (const FieldGroup<const Real>& group) {
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

void AtmosphereProcess::set_computed_group (const FieldGroup<Real>& group) {
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

void AtmosphereProcess::check_required_fields () const { 
  // Loop over all the field checks on input fields, and execute them
  for (const auto& fp_it : m_property_checks_in) {
    const auto& fid   = fp_it.first;
    const auto& fname = fid.name();
    const auto& gname = fid.get_grid_name();
    const auto& prop_check = *fp_it.second;
    if (has_required_field(fid)) {
      auto& f = get_field_in(fname,gname);
      EKAT_REQUIRE_MSG(prop_check.check(f),
          "Error! Input field property check failed.\n"
          "       Atm proc: " + this->name() + "\n"
          "       Field:    " + fname + "\n"
          "       Check:    " + prop_check.name() + "\n"
          " NOTE: we don't allow repairing input fields.\n");
    } else {
      auto& group = get_group_in(fname,gname);
      if (group.m_bundle) {
        EKAT_REQUIRE_MSG(prop_check.check(*group.m_bundle),
            "Error! Input field property check failed.\n"
            "       Atm proc: " + this->name() + "\n"
            "       Field:    " + group.m_info->m_group_name + "\n"
            "       Check:    " + prop_check.name() + "\n"
            " NOTE: we don't allow repairing input fields.\n");
      } else {
        for (const auto& it : group.m_fields) {
          const auto& f = *it.second;
          EKAT_REQUIRE_MSG(prop_check.check(f),
              "Error! Input field property check failed.\n"
              "       Atm proc: " + this->name() + "\n"
              "       Field:    " + fname + "\n"
              "       Check:    " + prop_check.name() + "\n"
              " NOTE: we don't allow repairing input fields.\n");
        }
      }
    }
  }
}

void AtmosphereProcess::check_computed_fields () {
  // Loop over all the field checks on input fields, and execute them
  for (const auto& fp_it : m_property_checks_out) {
    const auto& fid   = fp_it.first;
    const auto& fname = fid.name();
    const auto& gname = fid.get_grid_name();
    const auto& prop_check = *fp_it.second;
    if (has_computed_field(fid)) {
      auto& f = get_field_out(fname,gname);
      const auto check = prop_check.check(f);
      EKAT_REQUIRE_MSG(check || prop_check.can_repair(),
          "Error! Output field property check failed (and cannot be repaired).\n"
          "       Atm proc: " + this->name() + "\n"
          "       Field:    " + fname + "\n"
          "       Check:    " + prop_check.name() + "\n");
      if (not check) {
        prop_check.repair(f);
      }
    } else {
      auto& group = get_group_out(fname,gname);
      if (group.m_bundle) {
        auto& f = *group.m_bundle;
        const auto check = prop_check.check(f);
        EKAT_REQUIRE_MSG(check || prop_check.can_repair(),
            "Error! Output field property check failed (and cannot be repaired).\n"
            "       Atm proc: " + this->name() + "\n"
            "       Field:    " + group.m_info->m_group_name + "\n"
            "       Check:    " + prop_check.name() + "\n");
        if (not check) {
          prop_check.repair(f);
        }
      } else {
        for (const auto& it : group.m_fields) {
          auto& f = *it.second;
          const auto check = prop_check.check(f);
          EKAT_REQUIRE_MSG(check || prop_check.can_repair(),
              "Error! Output field property check failed (and cannot be repaired).\n"
              "       Atm proc: " + this->name() + "\n"
              "       Field:    " + it.first + "\n"
              "       Check:    " + prop_check.name() + "\n");
          if (not check) {
            prop_check.repair(f);
          }
        }
      }
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

void AtmosphereProcess::add_me_as_provider (const Field<Real>& f) {
  f.get_header_ptr()->get_tracking().add_provider(weak_from_this());
}

void AtmosphereProcess::add_me_as_customer (const Field<const Real>& f) {
  f.get_header_ptr()->get_tracking().add_customer(weak_from_this());
}

void AtmosphereProcess::add_internal_field (const Field<Real>& f) {
  m_internal_fields.push_back(f);
}

const Field<const Real>& AtmosphereProcess::
get_field_in(const std::string& field_name, const std::string& grid_name) const {
  return get_field_in_impl(field_name,grid_name);
}

Field<const Real>& AtmosphereProcess::
get_field_in(const std::string& field_name, const std::string& grid_name) {
  return get_field_in_impl(field_name,grid_name);
}

const Field<const Real>& AtmosphereProcess::
get_field_in(const std::string& field_name) const {
  return get_field_in_impl(field_name);
}

Field<const Real>& AtmosphereProcess::
get_field_in(const std::string& field_name) {
  return get_field_in_impl(field_name);
}

const Field<Real>& AtmosphereProcess::
get_field_out(const std::string& field_name, const std::string& grid_name) const {
  return get_field_out_impl(field_name,grid_name);
}

Field<Real>& AtmosphereProcess::
get_field_out(const std::string& field_name, const std::string& grid_name) {
  return get_field_out_impl(field_name,grid_name);
}

const Field<Real>& AtmosphereProcess::
get_field_out(const std::string& field_name) const {
  return get_field_out_impl (field_name);
}

Field<Real>& AtmosphereProcess::
get_field_out(const std::string& field_name) {
  return get_field_out_impl (field_name);
}

const FieldGroup<const Real>& AtmosphereProcess::
get_group_in(const std::string& group_name, const std::string& grid_name) const {
  return get_group_in_impl (group_name,grid_name);
}

FieldGroup<const Real>& AtmosphereProcess::
get_group_in(const std::string& group_name, const std::string& grid_name) {
  return get_group_in_impl (group_name,grid_name);
}

const FieldGroup<const Real>& AtmosphereProcess::
get_group_in(const std::string& group_name) const {
  return get_group_in_impl(group_name);
}

FieldGroup<const Real>& AtmosphereProcess::
get_group_in(const std::string& group_name) {
  return get_group_in_impl(group_name);
}

const FieldGroup<Real>& AtmosphereProcess::
get_group_out(const std::string& group_name, const std::string& grid_name) const {
  return get_group_out_impl(group_name,grid_name);
}

FieldGroup<Real>& AtmosphereProcess::
get_group_out(const std::string& group_name, const std::string& grid_name) {
  return get_group_out_impl(group_name,grid_name);
}

const FieldGroup<Real>& AtmosphereProcess::
get_group_out(const std::string& group_name) const {
  return get_group_out_impl(group_name);
}

FieldGroup<Real>& AtmosphereProcess::
get_group_out(const std::string& group_name) {
  return get_group_out_impl(group_name);
}

const Field<Real>& AtmosphereProcess::
get_internal_field(const std::string& field_name, const std::string& grid_name) const {
  return get_internal_field_impl(field_name,grid_name);
}

Field<Real>& AtmosphereProcess::
get_internal_field(const std::string& field_name, const std::string& grid_name) {
  return get_internal_field_impl(field_name,grid_name);
}

Field<Real>& AtmosphereProcess::
get_internal_field(const std::string& field_name) {
  return get_internal_field_impl(field_name);
}

const Field<Real>& AtmosphereProcess::
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

Field<const Real>& AtmosphereProcess::
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

Field<const Real>& AtmosphereProcess::
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

Field<Real>& AtmosphereProcess::
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

Field<Real>& AtmosphereProcess::
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

FieldGroup<const Real>& AtmosphereProcess::
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

FieldGroup<const Real>& AtmosphereProcess::
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

FieldGroup<Real>& AtmosphereProcess::
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

FieldGroup<Real>& AtmosphereProcess::
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

Field<Real>& AtmosphereProcess::
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

Field<Real>& AtmosphereProcess::
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
