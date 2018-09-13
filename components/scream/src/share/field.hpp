#ifndef SCREAM_FIELD_HPP
#define SCREAM_FIELD_HPP

#include "scream_types.hpp"

#include <vector>
#include <memory>

// I have to decide where to store a field's providers and customers, that is,
// a list of parametrizations that compute or need a given field. Should it
// be done by the FieldRepository? Or should a Field store its own? I lean toward
// the latter.
// NOTE: if two parametrizations both compute/update a field, we need to ensure 
//       that the two parametrizations are run sequentially. We need to devise
//       a mechanism to ensure that, or at the very least, a test to check it.
//       Notice that some fields are needed/computed/updated only within Scream, so
//       we can store a pointer to the Parametrization, but other are computed/needed
//       in other components, so for those we can at best store a string that
//       identifies the component (at least informally).

namespace scream
{

// A small structure to hold info about a field
struct FieldHeader {
  // These could actually be retrieved from the Kokkos View, but it probably makes sens
  // to be able to retrieve info from a header, without having to query the view
  // Besides, we may want to be able to query name/rank/dims BEFORE the view is actually instantiated.
  // Finally, if you look at Field, you may notice that the way we store the view
  // in the manager may lose info about the rank/dims.
  std::string       m_name;
  int               m_rank;
  std::vector<int>  m_dims;
  // Something about output/restart?
  // Perhaps something about the timestamp of the field (when it was last updated)?
};

// A field should be composed of metadata info (the header) and a pointer to the view
template<typename ViewType>
struct Field {
  using header_type = FieldHeader;
  using view_type   = ViewType;

  // Metadata (name, rank, dims,...)
  header_type                 m_header;
  // Actual data.
  // Note: when copying a Field, alwyays perform a shallow copy, so that if someone
  //       reassigns the view, all customers will see it.
  // I'm not 100% sure we need a shared_ptr. If we use Field's always by
  // reference, we should be fine with just a ViewType. However, if parametrizations
  // want to store copies of a Field, we need the shared_ptr
  std::shared_ptr<view_type>  m_field;
};

} // namespace scream

#endif // SCREAM_FIELD_HPP
