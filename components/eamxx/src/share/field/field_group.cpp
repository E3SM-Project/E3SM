#include "share/field/field_group.hpp"

namespace scream {

FieldGroup::FieldGroup (const std::string& name, const std::string& grid_name)
 : m_name (name)
 , m_grid_name(grid_name)
{
  m_individual_fields = std::make_shared<std::map<std::string,Field>>();
  m_monolithic_field  = std::make_shared<Field>();
}

FieldGroup
FieldGroup::get_const () const
{
  FieldGroup cgroup(m_name,m_grid_name);
  cgroup.m_monolithic_field = std::make_shared<Field>(m_monolithic_field->get_const());
  for (auto [fn,f] : *m_individual_fields)
    (*cgroup.m_individual_fields)[fn] = f.get_const();

  return cgroup;
}

void FieldGroup::
set_field (const Field& f)
{
  EKAT_REQUIRE_MSG (m_individual_fields->count(f.name())==0,
      "[FieldGroup] Error! Cannot reset an individual field once set.\n"
      " - group name: " + m_name + "\n"
      " - grid name : " + m_grid_name + "\n"
      " - field name: " + f.name() + "\n");

  if (has_monolithic_field()) {
    // Run additional checks, to ensure the field is TRULY a subfield of the monolithic one
    auto fp = f.get_header().get_parent();
    EKAT_REQUIRE_MSG (fp!=nullptr,
        "[FieldGroup::set_field] Error! Input field has no parent, but the grup stores a monolithic_field.\n"
        " - group name: " + m_name + "\n"
        " - grid name : " + m_grid_name + "\n"
        " - field name: " + f.name() + "\n");

    auto mp = m_monolithic_field->get_header().get_parent();
    if (mp==nullptr)
      EKAT_REQUIRE_MSG (fp->get_identifier().name()==m_monolithic_field->name(),
          "[FieldGroup::set_field] Error! Input field does not seem to be a subfield of the monolithic field.\n"
          " - group name: " + m_name + "\n"
          " - grid name : " + m_grid_name + "\n"
          " - field name: " + f.name() + "\n");
    else {
      // Same parent field
      EKAT_REQUIRE_MSG (fp==mp,
          "[FieldGroup::set_field] Error! Input field does not seem to be a subfield of the monolithic field.\n"
          " - group name: " + m_name + "\n"
          " - grid name : " + m_grid_name + "\n"
          " - field name: " + f.name() + "\n");

      // The slice idx of f is in [m_beg,m_end] slice idx range
      auto m_beg = m_monolithic_field->get_header().get_alloc_properties().get_subview_info().slice_idx;
      auto m_end = m_monolithic_field->get_header().get_alloc_properties().get_subview_info().slice_idx_end;
      auto f_idx = f.get_header().get_alloc_properties().get_subview_info().slice_idx;

      EKAT_REQUIRE_MSG (m_beg<=f_idx and f_idx<m_end,
          "[FieldGroup::set_field] Error! Input field does not seem to be a subfield of the monolithic field.\n"
          " - group name: " + m_name + "\n"
          " - grid name : " + m_grid_name + "\n"
          " - field name: " + f.name() + "\n");
    }
  }

  EKAT_REQUIRE_MSG (f.is_allocated(),
      "[FieldGroup::set_field] Error! Input field has not yet been allocated.\n"
      " - group name: " + m_name + "\n"
      " - grid name : " + m_grid_name + "\n"
      " - field name: " + f.name() + "\n");

  (*m_individual_fields)[f.name()] = f;
}

void FieldGroup::
set_monolithic_field (const Field& f)
{
  EKAT_REQUIRE_MSG (not m_monolithic_field->is_allocated(),
      "[FieldGroup] Error! Cannot reset the monolithic field once set.\n"
      " - group name: " + m_name + "\n"
      " - grid name : " + m_grid_name + "\n"
      " - field name: " + f.name() + "\n");

  *m_monolithic_field = f;
}

} // namespace scream
