#ifndef SCREAM_FIELD_GROUP_INFO_HPP
#define SCREAM_FIELD_GROUP_INFO_HPP

#include <ekat/util/ekat_string_utils.hpp>

#include <list>
#include <map>

namespace scream {

/*
 * A FieldGroupInfo is a small structure storing some info on a group of fields.
 *
 * A group is basically a "label" attached to fields, to allow users to
 * query a FieldManager for all fields that have such label attached. A field can
 * belong to any number of groups, or no group at all.
 *
 * A FieldGroupInfo stores:
 *
 *   - a list of the field names associated to this group;
 *   - whether the field were allocated as a single "bundled" field,
 *     with each field extracted as a "subview" of the bundled one;
 *   - if the allocation was "bundled", also store for each field
 *     the index that was used to extract the corresponding subview.
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

  FieldGroupInfo (const FieldGroupInfo& src) = default;

  int size() const { return m_fields_names.size(); }

  bool empty() const { return size()==0; }

  // The name of the group
  ci_string m_group_name;

  // The names of the fields in this group
  std::list<ci_string>   m_fields_names;

  // Whether the group was allocated as a bundle
  bool m_bundled;

  // If bundled, each field is subviewed along a different entry
  // along the same dimension.
  int m_subview_dim;

  // If bundled, for each field name, store the idx used
  // to subview each field from the bundle.
  std::map<ci_string,int>  m_subview_idx;
};

inline bool operator== (const FieldGroupInfo& lhs,
                        const FieldGroupInfo& rhs)
{
  return lhs.m_group_name==rhs.m_group_name &&
         lhs.m_fields_names==rhs.m_fields_names &&
         lhs.m_bundled==rhs.m_bundled &&
         lhs.m_subview_dim==rhs.m_subview_dim &&
         lhs.m_subview_idx==rhs.m_subview_idx;
}

} // namespace scream

#endif // SCREAM_FIELD_GROUP_INFO_HPP
