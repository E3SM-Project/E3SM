#include "share/field/field_group.hpp"

namespace scream {

FieldGroup::FieldGroup (const std::string& name, const std::string& grid_name)
 : FieldGroup(std::make_shared<FieldGroupInfo>(name),grid_name)
{
  // Nothing to do here
}

FieldGroup::
FieldGroup (const std::shared_ptr<FieldGroupInfo>& info,
            const std::string& grid_name)
 : m_info (info)
 , m_grid_name(grid_name)
{
  m_individual_fields = std::make_shared<std::map<ci_string,Field>>();
}

FieldGroup
FieldGroup::get_const () const
{
  FieldGroup cgroup(m_info,m_grid_name);
  if (m_monolithic_field!=nullptr)
    cgroup.m_monolithic_field = std::make_shared<Field>(m_monolithic_field->get_const());

  for (const auto& it : *m_individual_fields)
    cgroup.set_field(it.second.get_const());

  return cgroup;
}

void FieldGroup::
set_field (const Field& f)
{
  EKAT_REQUIRE_MSG (m_individual_fields->count(f.name())==0,
      "[FieldGroup] Error! Cannot reset an individual field once set.\n"
      " - group name: " + m_info->m_group_name + "\n"
      " - grid name : " + m_grid_name + "\n"
      " - field name: " + f.name() + "\n");
  EKAT_REQUIRE_MSG (ekat::contains(m_info->m_fields_names,f.name()),
      "[FieldGroup] Error! Field was not listed as part of the group.\n"
      " - group name: " + m_info->m_group_name + "\n"
      " - grid name : " + m_grid_name + "\n"
      " - field name: " + f.name() + "\n");


  (*m_individual_fields)[f.name()] = f;
}

void FieldGroup::set_monolithic_field (const Field& f)
{
  EKAT_REQUIRE_MSG (m_monolithic_field==nullptr,
      "[FieldGroup] Error! Cannot reset the monolithic field once set.\n"
      " - group name: " + m_info->m_group_name + "\n"
      " - grid name : " + m_grid_name + "\n"
      " - field name: " + f.name() + "\n");

  EKAT_REQUIRE_MSG (m_info->m_monolithic_allocation,
      "[FieldGroup] Error! Cannot set the monolithic field if monolithic allocation was not requested.\n"
      " - group name: " + m_info->m_group_name + "\n"
      " - grid name : " + m_grid_name + "\n"
      " - field name: " + f.name() + "\n");

  m_monolithic_field = std::make_shared<Field>(f);
}

} // namespace scream
