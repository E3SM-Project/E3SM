#include "share/field/field_group.hpp"

namespace scream {

FieldGroup::FieldGroup (const std::string& name)
 : m_info (new FieldGroupInfo(name))
{
  // Nothing to do here
}

FieldGroup::
FieldGroup (const FieldGroupInfo& info)
 : m_info (new FieldGroupInfo(info))
{
  // Nothing to do here
}

FieldGroup
FieldGroup::get_const () const {
  FieldGroup gr(*m_info);
  if (m_info->m_monolithic_allocation) {
    gr.m_monolithic_field = std::make_shared<Field>(m_monolithic_field->get_const());
  }
  for (const auto& it : m_individual_fields) {
    gr.m_individual_fields[it.first] = std::make_shared<Field>(it.second->get_const());
  }
  return gr;
}

const std::string& FieldGroup::grid_name () const {
  EKAT_REQUIRE_MSG(m_individual_fields.size()>0 || m_monolithic_field,
      "Error! Cannot establish the group grid name until fields have been added.\n");

  if (m_monolithic_field) {
    const auto& id = m_monolithic_field->get_header().get_identifier();
    return id.get_grid_name();
  } else {
    const auto& id = m_individual_fields.begin()->second->get_header().get_identifier();
    return id.get_grid_name();
  }
}

void FieldGroup::copy_fields (const FieldGroup& src) {
  m_monolithic_field = src.m_monolithic_field;
  for (auto it : src.m_individual_fields) {
    m_individual_fields[it.first] = it.second;
  }
}

} // namespace scream
