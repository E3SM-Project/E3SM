#ifndef SCREAM_FIELD_HPP
#define SCREAM_FIELD_HPP

#include "field_header.hpp"
#include <share/scream_types.hpp>
#include <share/util/scream_std_meta.hpp>
#include <share/util/scream_kokkos_utils.hpp>

#include <memory>   // For std::shared_ptr

namespace scream
{

// ======================== FIELD ======================== //

// A field should be composed of metadata info (the header) and a pointer to the view

template<typename DataType, typename MemSpace, typename MemManagement>
class Field {
public:
  using data_type           = DataType;
  using memory_space        = MemSpace;
  using memory_management   = MemManagement;
  using value_type          = typename util::ScalarType<data_type>::type;
  using view_type           = ViewType<data_type,MemSpace,MemManagement>;
  using const_data_type     = typename view_type::traits::const_data_type;
  using non_const_data_type = typename view_type::traits::non_const_data_type;
  using const_field_type    = Field<const_data_type,MemSpace,MemManagement>;
  using header_type         = FieldHeader;
  using identifier_type     = header_type::identifier_type;

  // Constructor(s)
  Field () = delete;
  Field (const identifier_type& id);
  Field (const std::shared_ptr<FieldHeader>& fh, const view_type& view);

  // This constructor allows const->const, nonconst->nonconst, and nonconst->const copies
  template<typename SrcDT>
  Field (const Field<SrcDT,MemSpace,MemManagement>& src);

  // Assignment: allows const->const, nonconst->nonconst, and nonconst->const copies
  template<typename SrcDT>
  Field& operator= (const Field<SrcDT,MemSpace,MemManagement>& src);

  // Getters
  const header_type& get_header () const { return *m_header; }
        header_type& get_header ()       { return *m_header; }
  const view_type&   get_view   () const { return  m_view;   }

  const std::shared_ptr<header_type>& get_header_ptr () const { return m_header; }

protected:

  // Metadata (name, rank, dims, customere/providers, time stamp, ...)
  std::shared_ptr<header_type>    m_header;

  // Actual data.
  view_type                       m_view;
};

// ================================= IMPLEMENTATION ================================== //

template<typename DataType, typename MemSpace, typename MemManagement>
Field<DataType,MemSpace,MemManagement>::
Field (const identifier_type& id)
 : m_header (new header_type(id))
{
  // Check the identifier has all the dimensions set
  error::runtime_check(m_header->get_identifier().dimensions_set(), "Error! Cannot create a field until all the field's dimensions are set.\n");

  typename view_type::traits::array_layout layout;
  for (int idim=0; idim<m_header->get_identifier().rank(); ++idim) {
    layout.dimension[idim] = m_header->get_identifier().dim(idim);
  }
  m_view = view_type(m_header->get_identifier().name(),layout);
}

template<typename DataType, typename MemSpace, typename MemManagement>
template<typename SrcDT>
Field<DataType,MemSpace,MemManagement>::
Field (const Field<SrcDT,MemSpace,MemManagement>& src)
 : m_header (src.get_header_ptr())
 , m_view   (src.get_view())
{
  // Check that underlying array layout is the same
  static_assert(std::is_same<non_const_data_type,
                             typename util::remove_all_consts<SrcDT>::type
                            >::value,
                "Error! Cannot use copy constructor if the underlying array layout is different.\n");
  // Check that destination is const or source is nonconst
  static_assert(std::is_same<data_type,const_data_type>::value ||
                std::is_same<SrcDT,non_const_data_type>::value,
                "Error! Cannot create a nonconst field from a const field.\n");
}

template<typename DataType, typename MemSpace, typename MemManagement>
Field<DataType,MemSpace,MemManagement>::
Field (const std::shared_ptr<FieldHeader>& header, const view_type& view)
 : m_header (header)
 , m_view   (view)
{
  // Perform some checks
  error::runtime_check(view.size()==0 || static_cast<bool>(m_header),
                       "Error! Input header is null, but view has non-trivial size.\n");

  error::runtime_check(static_cast<int>(view.size())==m_header->get_identifier().size(),
                       "Error! Input header and view allocation sizes do not match.\n");

  for (int idim=0; idim<view_type::Rank; ++idim) {
    error::runtime_check(view.extent_int(idim) == m_header->get_identifier().dim(idim),
                         "Error! Input header and view dimensions do not match.\n");
  }
}

template<typename DataType, typename MemSpace, typename MemManagement>
template<typename SrcDT>
Field<DataType,MemSpace,MemManagement>&
Field<DataType,MemSpace,MemManagement>::
operator= (const Field<SrcDT,MemSpace,MemManagement>& src) {
  // Check that underlying array layout is the same
  static_assert(std::is_same<non_const_data_type,
                             typename util::remove_all_consts<SrcDT>::type
                            >::value,
                "Error! Cannot use copy constructor if the underlying array layout is different.\n");
  // Check that destination is const or source is nonconst
  static_assert(std::is_same<data_type,const_data_type>::value ||
                std::is_same<SrcDT,non_const_data_type>::value,
                "Error! Cannot create a nonconst field from a const field.\n");
  if (&src!=*this) {
    m_header = src.get_header_ptr();
    m_view   = src.get_view();
  }
}

// ----------------- Free functions for field reshaping ----------------- //

// These free functions allow to create an unmanaged view of a Field
// changing the layout of the field. This is handy because the FieldRepository
// stores *only* 1d views (otherwise it would need multiple containers, one per
// data type). If a customer of the field wants to store a field with a particular
// layout (e.g., using compile time dimensions for performance), they need to 'view'
// the field with a different data type. That, by the way, is only possible
// if the resulting view is unmanaged.

template<typename DstDT, typename SrcDT, typename MemSpace, typename MemManagement>
Field<DstDT,MemSpace,MemoryUnmanaged>
reinterpret_field (const Field<SrcDT,MemSpace,MemManagement>& f,
                   std::shared_ptr<FieldHeader> header) {
  using DstField = Field<DstDT,MemSpace,MemoryUnmanaged>;
  using DstView  = typename DstField::view_type;
  error::runtime_check(header->get_identifier().dimensions_set(), "Error! Target header dimensions have not been set yet.\n");
  error::runtime_check(header->get_identifier().size()==f.get_header().get_identifier().size(), "Error! You can only reinterpret a field to a data type with the same flattened 1d length.\n");

  Kokkos::LayoutRight layout;
  for (int idim=0; idim<header->get_identifier().rank(); ++idim) {
    layout.dimension[idim] = header->get_identifier().dim(idim);
  }

  DstView dst_view (f.get_view().data(),layout);

  return DstField(header,dst_view);
}

template<typename DstDT, typename SrcDT, typename MemSpace, typename MemManagement>
Field<DstDT,MemSpace,MemoryUnmanaged>
reinterpret_field (const Field<SrcDT,MemSpace,MemManagement>& f) {
  return reinterpret_field<DstDT>(f,f.get_header_ptr());
}

} // namespace scream

#endif // SCREAM_FIELD_HPP
