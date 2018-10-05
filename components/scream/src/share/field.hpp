#ifndef SCREAM_FIELD_HPP
#define SCREAM_FIELD_HPP

#include "field_header.hpp"
#include "scream_types.hpp"

#include <memory>   // For std::shared_ptr

// I have to decide where to store a field's providers and customers, that is,
// a list of parametrizations that compute or need a given field. Should it
// be done by the FieldRepository? Or should a Field store its own? I lean toward
// the latter.
// NOTE: if two parametrizations both compute/update a field, we need to ensure 
//       that the two parametrizations are run sequentially. The FieldTracking
//       structure can be used (among other things) to enforce this requirement,
//       by checking the time stamps of the fields.

namespace scream
{

// ======================== FIELD ======================== //

// A field should be composed of metadata info (the header) and a pointer to the view
template<typename DataType, typename MemSpace>
class Field {
public:
  using header_type   = FieldHeader;
  using view_type     = ViewManaged<DataType,MemSpace>;
  using view_ptr_type = std::shared_ptr<view_type>;

  // Constructor(s)
  Field () = default;
  Field (const Field&) = default;
  Field (const header_type& header);

  // Assignment (defaulted)
  Field& operator= (const Field&) = default;

  // Getters
  const header_type& header () const { return m_header; }
  const view_type&   view   () const { error::runtime_check(m_view, "Error! View not set yet.", -1); return m_view; }

protected:
  // Metadata (name, rank, dims,...)
  header_type                 m_header;
  // Actual data.
  // Note: when copying a Field, alwyays perform a shallow copy, so that if someone
  //       reassigns the view, all customers will see it.
  // I'm not 100% sure we need a shared_ptr. If we use Field's always by
  // reference, we should be fine with just a ViewType. However, if parametrizations
  // want to store copies of a Field, we need the shared_ptr
  std::shared_ptr<view_type>  m_view;
};

template<typename DataType, typename MemSpace>
Field<DataType,MemSpace>::Field (const header_type& header) {
  m_header = header;
  Kokkos::LayoutRight layout;
  for (int idim=0; idim<m_header.rank(); ++idim) {
    layout.dimension[idim] = m_header.dim(idim);
  }

  m_view = std::make_shared<view_type>(m_header.name(),layout);
}

} // namespace scream

#endif // SCREAM_FIELD_HPP
