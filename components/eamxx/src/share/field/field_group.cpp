#include "share/field/field_group.hpp"

namespace scream {

FieldGroup::FieldGroup (const std::string& name, const std::string& grid_name)
 : m_name (name)
 , m_grid_name(grid_name)
{
  m_individual_fields = std::make_shared<std::map<std::string,Field>>();
}

FieldGroup
FieldGroup::get_const () const
{
  FieldGroup cgroup(m_name,m_grid_name);
  if (m_monolithic_field!=nullptr)
    cgroup.set_monolithic_field(m_monolithic_field->get_const());
  else
    for (const auto& it : *m_individual_fields)
      cgroup.set_field(it.second.get_const());

  return cgroup;
}

void FieldGroup::
set_field (const Field& f)
{
  EKAT_REQUIRE_MSG (m_monolithic_field==nullptr,
      "[FieldGroup] Error! Cannot set an individual field once monolithic field has been set.\n"
      " - group name: " + m_name + "\n"
      " - grid name : " + m_grid_name + "\n"
      " - field name: " + f.name() + "\n");

  EKAT_REQUIRE_MSG (m_individual_fields->count(f.name())==0,
      "[FieldGroup] Error! Cannot reset an individual field once set.\n"
      " - group name: " + m_name + "\n"
      " - grid name : " + m_grid_name + "\n"
      " - field name: " + f.name() + "\n");

  (*m_individual_fields)[f.name()] = f;
}

void FieldGroup::
set_monolithic_field (const Field& f, const std::vector<std::string>& names,
                      const int subview_dim, const int subview_beg)
{
  EKAT_REQUIRE_MSG (m_monolithic_field==nullptr,
      "[FieldGroup] Error! Cannot reset the monolithic field once set.\n"
      " - group name: " + m_name + "\n"
      " - grid name : " + m_grid_name + "\n"
      " - field name: " + f.name() + "\n");
  EKAT_REQUIRE_MSG (m_individual_fields->size()==0,
      "[FieldGroup] Error! Cannot set the monolithic field once one of the fields have been set.\n"
      " - group name: " + m_name + "\n"
      " - grid name : " + m_grid_name + "\n"
      " - field name: " + f.name() + "\n");

  m_monolithic_field = std::make_shared<Field>(f);

  int offset = subview_beg;
  for (const auto& n : names) {
    (*m_individual_fields)[n] = m_monolithic_field->subfield(n,subview_dim,offset);
    ++offset;
  }
}

} // namespace scream
