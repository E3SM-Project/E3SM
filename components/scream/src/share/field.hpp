#ifndef SCREAM_FIELD_HPP
#define SCREAM_FIELD_HPP

#include "field_header.hpp"
#include "scream_types.hpp"
#include "share/util/scream_kokkos_utils.hpp"

#include <memory>   // For std::shared_ptr

namespace scream
{

// ======================== FIELD ======================== //

// A field should be composed of metadata info (the header) and a pointer to the view
// Note: if StoreAs1D is true, the underlying data type will be a 1d array, regardless
//       of whether DataType is 1d or n-dimensional. If it is false, then the underlying
//       data type is exactly DataType. This is needed because the FieldRepository
//       stores everything as 1d arrays (otherwise would need one repository per data type).
//       On the other hand, the field's header *always* stores the *real* field dimensions.
//       This can create an incompatibility between data type in the Kokkos view, and dimensions
//       in the field header, which can cause issues at compile time. To overcome this issue,
//       we keep track of whether or not the view is 'flattened' with an extra template parameter.
template<typename DataType, typename MemSpace, bool StoreAs1D = false>
class Field {
public:

  using scalar_type      = typename util::ScalarType<DataType>::type;
  using data_type_impl   = typename std::conditional<StoreAs1D,scalar_type*,DataType>::type;
  using view_type        = ViewManaged<data_type_impl,MemSpace>;
  using traits           = typename view_type::traits;
  using const_field_type = Field<typename traits::const_data_type,MemSpace,StoreAs1D>;
  using header_type      = FieldHeader;

  // Constructor(s)
  Field () = delete;
  Field (const Field&) = default;
  Field (const std::shared_ptr<header_type> header);

  // This is to allow the following *almost copy* constructor to work without the need of extra getters
  template<typename,typename,bool> friend class Field;

  // This constructor serves two purposes:
  //   1) allows creation of a field to const by simply invoking copy constructor syntax.
  //   2) allows to create a reshaped view (of data type DataType), from/to a 1d view view
  //      (provided DT is a data type with a cv not stronger than DataType)
  template<typename SrcDataType,bool SrcStoreAs1D>
  Field (const Field<SrcDataType,MemSpace,SrcStoreAs1D>& src);

  // Assignment (defaulted)
  Field& operator= (const Field&) = default;

  // If the header was created without specifying dimensions
  // (e.g., they were not known at construction time),
  // call this method to set them.
  // NOTE: this method cannot be called when dimensions are set.
  //       Check the field header to make sure they are not.
  void set_dimensions (const std::vector<int>& dims);

  // Getters
  const header_type& header () const { return *m_header; }
        view_type    view   () const;

  const_field_type get_const () const;
  

protected:

  // Allocate the view
  void allocate ();

  // Metadata (name, rank, dims,...)
  std::shared_ptr<header_type>    m_header;
  // Actual data.
  // Note: when copying a Field, alwyays perform a shallow copy, so that if someone
  //       reassigns the view, all customers will see it.
  // I'm not 100% sure we need a shared_ptr. If we use Field's always by
  // reference, we should be fine with just a ViewType. However, if parametrizations
  // want to store copies of a Field, we need the shared_ptr
  std::shared_ptr<view_type>      m_view;
};

template<typename DataType, typename MemSpace, bool StoreAs1D>
Field<DataType,MemSpace,StoreAs1D>::
Field (std::shared_ptr<header_type> header)
 : m_header (header)
{
  if (m_header->dimensions_set()) {
    allocate ();
  }
}
template<typename T>
struct MyDebug {};

template<typename DataType, typename MemSpace, bool StoreAs1D>
template<typename SrcDataType, bool SrcStoreAs1D>
Field<DataType,MemSpace,StoreAs1D>::
Field (const Field<SrcDataType,MemSpace,SrcStoreAs1D>& src)
 : m_header(src.m_header)
{
//typename MyDebug<DataType>::toh v;
  static_assert(!std::is_const<SrcDataType>::value || std::is_const<DataType>::value, "Error! Cannot copy from const to non-const.\n");
  error::runtime_check(src.m_header->dimensions_set(), "Error! Cannot reshape a view that has not been allocated yet.\n", -1);
  if (StoreAs1D) {
    // Regardless of input view layout, we switch to 1d here
    m_view = std::make_shared<view_type>(src.m_view->data(),m_header->size());
  } else {
    // Regardless of input view layout, we switch to n-dimensional here
    Kokkos::LayoutRight layout;
    for (int idim=0; idim<m_header->rank(); ++idim) {
      layout.dimension[idim] = m_header->dim(idim);
    }
    m_view = std::make_shared<view_type>(src.m_view->data(),layout);
  }
}

template<typename DataType, typename MemSpace, bool StoreAs1D>
void Field<DataType,MemSpace,StoreAs1D>::
set_dimensions (const std::vector<int>& dims) {
  m_header->set_dimensions(dims);
  allocate ();
}

template<typename DataType, typename MemSpace, bool StoreAs1D>
typename Field<DataType,MemSpace,StoreAs1D>::view_type
Field<DataType,MemSpace,StoreAs1D>::view () const {
  error::runtime_check(static_cast<bool>(m_view), "Error! View not set yet.", -1);
  return *m_view;
}

template<typename DataType, typename MemSpace, bool StoreAs1D>
inline typename Field<DataType,MemSpace,StoreAs1D>::const_field_type 
Field<DataType,MemSpace,StoreAs1D>::get_const () const {
  return const_field_type(*this);
}

template<typename DataType, typename MemSpace, bool StoreAs1D>
void Field<DataType,MemSpace,StoreAs1D>::allocate () {
  if (StoreAs1D) {
    // We store a 1d view, regardless of the layout in the header
    m_view = std::make_shared<view_type>(m_header->name(),m_header->size());
  } else {
    // We store a view whose dat type reflects the layout in the header
    Kokkos::LayoutRight layout;
    for (int idim=0; idim<m_header->rank(); ++idim) {
      layout.dimension[idim] = m_header->dim(idim);
    }
    m_view = std::make_shared<view_type>(m_header->name(),layout);
  }
}

} // namespace scream

#endif // SCREAM_FIELD_HPP
