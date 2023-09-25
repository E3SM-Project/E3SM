#ifndef SCREAM_FIELD_GROUP_HPP
#define SCREAM_FIELD_GROUP_HPP

#include "share/field/field_group_info.hpp"
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
 *   - a FieldGroupInfo struct
 *   - a list of fields pointers
 *   - a grid name (the grid where the fields are)
 *
 * The same FieldGroupInfo can be recycled for several FieldGroup's, each living
 * on a different grid.
 *
 * Notice that, if the allocation was bundled, the big bundle is allocated
 * with layout given by grid->get_Xd_vector_layout(), where grid is the
 * grid object where the fields are defined on, and X=2 or 3.
 * Each field is then subviewed at entry k (different for each field)
 * along dimension I (same for all field) of the big bundle field.
 *
 * E.g., say we have 3d scalar fields F1,F2,F3,F4 belonging to group MyGroup,
 * which is then allocated as a bundled field F. F will have layout
 * given by grid->get_3d_vector_layout(). Say this layout is (COL,CMP,LEV).
 * Each field is subviewed along m_subview_dim=1, at entry 0,1,2,3 respectively.
 * Note: as of 02/2021 m_subview_dim is *always* 1, but we store this bit
 *       of info nevertheless, in case things change later on.
 */

// In order to allow downstream code to still use FieldGroup<T> during the refactor,
// the non-templated class need to use a different name. When refactor is complete,
// and all downstream code uses FieldGroup instead of FieldGroup<T>, we can remove this
// macro, and sed s/FieldGroup/FieldGroup/g all over the repo.

struct FieldGroup {
  using ci_string = FieldGroupInfo::ci_string;

  FieldGroup (const std::string& name);
  FieldGroup (const FieldGroupInfo& info);

  FieldGroup (const FieldGroup&) = default;

  FieldGroup get_const () const;

  FieldGroup& operator= (const FieldGroup& src) = default;

  const std::string& grid_name () const;

  // The fields in this group
  std::map<ci_string,std::shared_ptr<Field>> m_fields;

  // If m_info->m_bundled is true, this is the field that all fields
  // in m_fields are a subview of.
  std::shared_ptr<Field> m_bundle;

  // The info of this group.
  std::shared_ptr<FieldGroupInfo>  m_info;

private:

  // Only used inside this class;
  FieldGroup () = default;

  void copy_fields (const FieldGroup& src);
};

// We use this to find a FieldGroup in a std container.
// We do NOT allow two entries with same group name and grid name in such containers.
inline bool operator== (const FieldGroup& lhs, const FieldGroup& rhs) {
  return lhs.m_info->m_group_name == rhs.m_info->m_group_name &&
         lhs.grid_name() == rhs.grid_name();
}

} // namespace scream

#endif // SCREAM_FIELD_GROUP_HPP
