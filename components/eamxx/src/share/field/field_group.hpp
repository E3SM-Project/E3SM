#ifndef SCREAM_FIELD_GROUP_HPP
#define SCREAM_FIELD_GROUP_HPP

#include "share/field/field.hpp"

namespace scream {

/*
 * A FieldGroup is a small structure storing some info on a group of fields
 * as well as pointers to the fields.
 *
 * A group is basically a "label" attached to fields, to allow users to
 * query a FieldManager for all fields that have such label attached. A field can
 * belong to any number of groups, or no group at all.
 *
 * A FieldGroup contains:
 *
 *   - a name and a grid name (where the group lives)
 *   - a map name->field_pointer
 *   - possibly, a monolithic version of the field
 *
 * Notice that, if the group has a monolithic allocation, the monolithic field is allocated
 * with layout given by grid->get_Xd_vector_layout(), where grid is the
 * grid object where the fields are defined on, and X=2 or 3.
 * Each field is then subviewed at entry k (different for each field)
 * along dimension I (same for all field) of the monolithic field.
 *
 * E.g., say we have 3d scalar fields F1,F2,F3,F4 belonging to group MyGroup,
 * which is then allocated as a monolithic field F. F will have layout
 * given by grid->get_3d_vector_layout(). Say this layout is (COL,CMP,LEV).
 * Each field is subviewed along m_subview_dim=1, at entry 0,1,2,3 respectively.
 * Note: as of 02/2021 m_subview_dim is *always* 1, but we store this bit
 *       of info nevertheless, in case things change later on.
 */

// In order to allow downstream code to still use FieldGroup<T> during the refactor,
// the non-templated class need to use a different name. When refactor is complete,
// and all downstream code uses FieldGroup instead of FieldGroup<T>, we can remove this
// macro, and sed s/FieldGroup/FieldGroup/g all over the repo.

class FieldGroup {
public:

  FieldGroup (const std::string& name, const std::string& grid_name);
  FieldGroup (const FieldGroup&) = default;

  FieldGroup& operator= (const FieldGroup& src) = default;

  FieldGroup get_const () const;

  const std::string& name () const { return m_name; }
  const std::string& grid_name () const { return m_grid_name; }

  const std::map<std::string,Field>& individual_fields () const { return *m_individual_fields; }
        std::map<std::string,Field>& individual_fields ()       { return *m_individual_fields; }

  bool has_monolithic_field () const { return m_monolithic_field->is_allocated(); }
  const Field& monolithic_field () const { return *m_monolithic_field; }
        Field& monolithic_field ()       { return *m_monolithic_field; }

  void set_field (const Field& f);
  void set_monolithic_field (const Field& f);

  std::size_t size () const { return m_individual_fields->size(); }

private:

  // The name of this group
  std::string m_name;

  // The name of the grid where the group lives
  std::string m_grid_name;

  // Only used inside this class;
  FieldGroup () = default;

  std::shared_ptr<std::map<std::string,Field>> m_individual_fields;

  // If m_info->m_monolithic_alloc is true, this is the field
  // that all fields in m_individual_fields are a subview of.
  std::shared_ptr<Field> m_monolithic_field;
};

// We use this to find a FieldGroup in a std container.
// We do NOT allow two entries with same group name and grid name in such containers.
inline bool operator== (const FieldGroup& lhs, const FieldGroup& rhs) {
  return lhs.name() == rhs.name() &&
         lhs.grid_name() == rhs.grid_name();
}

} // namespace scream

#endif // SCREAM_FIELD_GROUP_HPP
