#ifndef SCREAM_FIELD_HPP
#define SCREAM_FIELD_HPP

#include "field_header.hpp"
#include <share/scream_types.hpp>
#include <share/util/scream_kokkos_utils.hpp>

#include <memory>   // For std::shared_ptr

namespace scream
{

// ======================== FIELD ======================== //

// A field should be composed of metadata info (the header) and a pointer to the view

template<typename DataType, typename MemSpace, typename MemManagment>
class Field {
public:
  using data_type           = DataType;
  using memory_space        = MemSpace;
  using memory_management   = MemManagment;
  using value_type          = typename util::ScalarType<data_type>::type;
  using view_type           = ViewType<data_type,MemSpace,MemManagment>;
  using const_data_type     = typename view_type::traits::const_data_type;
  using non_const_data_type = typename view_type::traits::non_const_data_type;
  using const_field_type    = Field<const_data_type,MemSpace,MemManagment>;
  using header_type         = FieldHeader;
  using identifier_type     = header_type::identifier_type;

  // Constructor(s)
  Field () = delete;
  Field (const Field&) = default;
  Field (const identifier_type& id);

  // Assignment (defaulted)
  Field& operator= (const Field&) = default;

  // Getters
        header_type& get_header ()       { return *m_header; }
  const header_type& get_header () const { return *m_header; }
        view_type    get_view   () const { return *m_view;   }

  std::shared_ptr<const header_type> get_header_ptr () const { return m_header; }
  std::shared_ptr<const view_type>   get_view_ptr   () const { return m_view;   }

  // Get a const version of this field
  const_field_type get_const () const;

  // Check whether the kokkos view was already allocated
  bool is_allocated () const { return m_view->size()>0; }

  // Allocate the view.
  // WARNING: you can only call this method once! Use is_allocated to check if you
  //          are allowed to call it.
  void allocate ();

protected:

  // Make the nonconst data_type version of this class a friend, so they can call the protected constructor
  friend class Field<non_const_data_type,memory_space,memory_management>;

  // This protected constructor is used internally
  Field (std::shared_ptr<FieldHeader> fh, std::shared_ptr<view_type> view);

  // Metadata (name, rank, dims, customere/providers, time stamp, ...)
  std::shared_ptr<header_type>    m_header;

  // Actual data.
  std::shared_ptr<view_type>      m_view;
};

// ================================= IMPLEMENTATION ================================== //

template<typename DataType, typename MemSpace, typename MemManagment>
Field<DataType,MemSpace,MemManagment>::
Field (const identifier_type& id)
 : m_header (new header_type(id))
 , m_view   (new view_type())
{
  if (m_header->get_identifier().dimensions_set()) {
    allocate ();
  }
}

template<typename DataType, typename MemSpace, typename MemManagment>
Field<DataType,MemSpace,MemManagment>::
Field (std::shared_ptr<FieldHeader> header, std::shared_ptr<view_type> view)
 : m_header (header)
 , m_view   (view)
{
  // Nothing to be done here
}

template<typename DataType, typename MemSpace, typename MemManagment>
void Field<DataType,MemSpace,MemManagment>::allocate () {
  // Check we can indeed proceed with allocation
  error::runtime_check(!is_allocated(), "Error! Cannot call 'allocate' twice. Bad things could happen...\n");
  error::runtime_check(m_header->get_identifier().dimensions_set(), "Error! Cannot call 'allocate' before all the field's dimensions are set.\n");

  Kokkos::LayoutRight layout;
  for (int idim=0; idim<m_header->get_identifier().rank(); ++idim) {
    layout.dimension[idim] = m_header->get_identifier().dim(idim);
  }
  *m_view = view_type(m_header->get_identifier().name(),layout);
}

template<typename DataType, typename MemSpace, typename MemManagment>
typename Field<DataType,MemSpace,MemManagment>::const_field_type
Field<DataType,MemSpace,MemManagment>::get_const () const {
  auto const_view = std::make_shared<typename const_field_type::view_type>(*m_view);
  return const_field_type(m_header,const_view);
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
reinterpret_field (const Field<SrcDT,MemSpace,MemManagement>& f) {
  return reinterpret_field(f,f.get_header_ptr());
}

template<typename DstDT, typename SrcDT, typename MemSpace, typename MemManagement>
Field<DstDT,MemSpace,MemoryUnmanaged>
reinterpret_field (const Field<SrcDT,MemSpace,MemManagement>& f,
                   std::shared_ptr<FieldHeader> header) {
  using DstField = Field<DstDT,MemSpace,MemoryUnmanaged>;
  using DstView  = typename DstField::view_type;
  error::runtime_check(f.is_allocated(), "Error! Cannot reshape a field that has not been allocated yet.\n");
  error::runtime_check(header->get_identifier().dimensions_set(), "Error! Target header dimensions have not been set yet.\n");
  error::runtime_check(header->get_identifier().size()==f.get_header()->get_identifier().size(), "Error! You can only reinterpret a field to a data type with the same flattened 1d length.\n");

  Kokkos::LayoutRight layout;
  for (int idim=0; idim<header->get_identifier().rank(); ++idim) {
    layout.dimension[idim] = header->get_identifier().dim(idim);
  }

  std::shared_ptr<DstView> dst_view = std::make_shared<DstView>(header->get_identifier().name(),f.get_view().data(),layout);

  return DstField(header,dst_view);
}


} // namespace scream

#endif // SCREAM_FIELD_HPP
