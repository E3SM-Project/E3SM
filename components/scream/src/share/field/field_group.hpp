#ifndef SCREAM_FIELD_GROUP_HPP
#define SCREAM_FIELD_GROUP_HPP

#include <ekat/util/ekat_string_utils.hpp>
#include <ekat/std_meta/ekat_std_utils.hpp>

#include <set>
#include <map>

namespace scream {

/*
 * A FieldGroup is a small structure storing some info on a group of fields
 * (wrapped in a FieldGroupInfo sturct), as well as pointers to the fields.
 *
 * A group is basically a "label" attached to fields, to allow users to
 * query the repo for all fields that have such label attached. A field can
 * belong to any number of groups, or no group at all.
 *
 * A FieldGroupInfo stores:
 *
 *   - a list of the field names associated to this group;
 *   - whether the field were allocated as a single "bundled" field,
 *     with each field extracted as a "subview" of the bundled one;
 *   - if the allocation was "bundled", also store for each field
 *     the index that was used to extract the corresponding subview.
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

struct FieldGroupInfo
{
  using ci_string  = ekat::CaseInsensitiveString;

  // Default initialize everything
  FieldGroupInfo (const ci_string& group_name)
    : m_group_name (group_name)
    , m_fields_names{}
    , m_bundled (false)
    , m_subview_dim(-1)
    , m_subview_idx{}
  {
    // Nothing to do here
  }

  int size() const { return m_fields_names.size(); }

  // The name of the group
  ci_string m_group_name;

  // The names of the fields in this group
  std::set<ci_string>   m_fields_names;

  // Whether the group was allocated as a bundle
  bool m_bundled;

  // If bundled, each field is subviewed along a different entry
  // along the same dimension.
  int m_subview_dim;

  // If bundled, for each field name, store the idx used
  // to subview each field from the bundle.
  std::map<ci_string,int>  m_subview_idx;
};

template<typename RealType>
class Field;

template<typename RealType>
struct FieldGroup {
  using field_type = Field<RealType>;
  using ci_string = FieldGroupInfo::ci_string;

  FieldGroup (const std::string& name, const std::string& grid)
   : m_info (new FieldGroupInfo(name))
   , m_grid_name(grid)
  {
    // Nothing to do here
  }

  FieldGroup (const FieldGroup&) = default;

  // The fields in this group
  std::map<ci_string,std::shared_ptr<field_type>> m_fields;

  // If m_info->m_bundled is true, this is the field that all fields
  // in m_fields are a subview of.
  std::shared_ptr<field_type> m_bundle;

  // The info of this group.
  std::shared_ptr<FieldGroupInfo>  m_info;

  // The grid where these field live on
  ci_string  m_grid_name;
};

} // namespace scream

#endif // SCREAM_FIELD_GROUP_HPP
