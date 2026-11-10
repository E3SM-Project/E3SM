#ifndef SCREAM_FIELD_GROUP_INFO_HPP
#define SCREAM_FIELD_GROUP_INFO_HPP

#include <ekat_string_utils.hpp>

#include <list>
#include <set>
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
 *   - whether the field were allocated as a single monolithic field,
 *     with each field extracted as a "subview" of the monolithic one;
 */

struct FieldGroupInfo
{
  using ci_string  = ekat::CaseInsensitiveString;

  // Default initialize everything
  FieldGroupInfo (const ci_string& group_name)
    : m_group_name (group_name)
    , m_fields_names{}
    , m_monolithic_allocation (false)
  {
    // Nothing to do here
  }

  FieldGroupInfo (const FieldGroupInfo& src) = default;

  int size() const { return m_fields_names.size(); }

  bool empty() const { return size()==0; }

  // The name of the group
  ci_string m_group_name;

  // The names of the fields in this group
  std::set<ci_string>   m_fields_names;

  // Store the grid which registered each field
  std::map<ci_string, std::set<ci_string>> m_grid_registered;

  // Store any grid that is requested for a group.
  // This is useful in the case where we allocate
  // a monolithic field, we can add a grid that may
  // not have any registered fields, but that we want
  // the group to exist.
  std::set<ci_string> m_requested_grids;

  // Whether the group allocated a monolithic field.
  // NOTE: if monolithic, it will be monolithic across grids,
  //       and contain the same fields across all grids
  bool m_monolithic_allocation = false;
};

inline bool operator== (const FieldGroupInfo& lhs,
                        const FieldGroupInfo& rhs)
{
  return lhs.m_group_name==rhs.m_group_name &&
         lhs.m_fields_names==rhs.m_fields_names &&
         lhs.m_monolithic_allocation==rhs.m_monolithic_allocation;
}

} // namespace scream

#endif // SCREAM_FIELD_GROUP_INFO_HPP
