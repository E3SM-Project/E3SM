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

template<typename ScalarType, typename MemSpace, typename MemManagement>
class Field {
private:
public:

  using view_type            = ViewType<ScalarType*,MemSpace,MemManagement>;
  using memory_space         = MemSpace;
  using memory_management    = MemManagement;
  using value_type           = typename view_type::traits::value_type;
  using const_value_type     = typename view_type::traits::const_value_type;
  using non_const_value_type = typename view_type::traits::non_const_value_type;
  using const_field_type     = Field<const_value_type,MemSpace,MemManagement>;
  using header_type          = FieldHeader;
  using identifier_type      = header_type::identifier_type;

  // Statically check that ScalarType is not an array.
  static_assert(view_type::Rank==1, "Error! ScalarType should not be a (multidimensional) array.\n");

  // Constructor(s)
  Field () = delete;
  Field (const identifier_type& id);
  Field (const header_type& header);

  // WARNING: this is an expert user constructor. We do perform some checks, to make
  //          sure things are consistent, but it is not to be used frequently.
  Field (const std::shared_ptr<header_type>& fh, const view_type& view);

  // This constructor allows const->const, nonconst->nonconst, and nonconst->const copies
  template<typename SrcDT>
  Field (const Field<SrcDT,MemSpace,MemManagement>& src);

  // Assignment: allows const->const, nonconst->nonconst, and nonconst->const copies
  template<typename SrcDT>
  Field& operator= (const Field<SrcDT,MemSpace,MemManagement>& src);

  // ---- Getters ---- //
  const header_type& get_header () const { return *m_header; }
        header_type& get_header ()       { return *m_header; }
  const view_type&   get_view   () const { return  m_view;   }

  const std::shared_ptr<header_type>& get_header_ptr () const { return m_header; }

  bool is_allocated () const { return m_allocated; }

  // ---- Setters ---- //

  //Allocate the actual view
  void allocate_view ();

protected:

  // Metadata (name, rank, dims, customere/providers, time stamp, ...)
  std::shared_ptr<header_type>    m_header;

  // Actual data.
  view_type                       m_view;

  // Keep track of whether the field has been allocated
  bool                            m_allocated;
};

// ================================= IMPLEMENTATION ================================== //

template<typename ScalarType, typename MemSpace, typename MemManagement>
Field<ScalarType,MemSpace,MemManagement>::
Field (const identifier_type& id)
 : m_header    (new header_type(id))
 , m_allocated (false)
{
  // At the very least, the allocation properties need to accommodate this field's value_type.
  m_header->get_alloc_properties().request_value_type_allocation<value_type>();
}

template<typename ScalarType, typename MemSpace, typename MemManagement>
template<typename SrcScalarType>
Field<ScalarType,MemSpace,MemManagement>::
Field (const Field<SrcScalarType,MemSpace,MemManagement>& src)
 : m_header    (src.get_header_ptr())
 , m_view      (src.get_view())
 , m_allocated (src.is_allocated())
{
  using src_field_type = decltype(src);

  // Check that underlying value type
  static_assert(std::is_same<non_const_value_type,
                             typename src_field_type::non_const_value_type
                            >::value,
                "Error! Cannot use copy constructor if the underlying value_type is different.\n");
  // Check that destination is const or source is nonconst
  static_assert(std::is_same<value_type,const_value_type>::value ||
                std::is_same<typename src_field_type::value_type,non_const_value_type>::value,
                "Error! Cannot create a nonconst field from a const field.\n");
}

template<typename ScalarType, typename MemSpace, typename MemManagement>
Field<ScalarType,MemSpace,MemManagement>::
Field (const std::shared_ptr<header_type>& header, const view_type& view)
 : m_header    (header)
 , m_view      (view)
 , m_allocated (true)
{
  // Perform some checks
  error::runtime_check(static_cast<bool>(header),"Error! Input header is null.\n");
  error::runtime_check(view.size()>0, "Error! Input view has 0 size. This constructor should only be used with an allocated view\n");
}

template<typename ScalarType, typename MemSpace, typename MemManagement>
template<typename SrcScalarType>
Field<ScalarType,MemSpace,MemManagement>&
Field<ScalarType,MemSpace,MemManagement>::
operator= (const Field<SrcScalarType,MemSpace,MemManagement>& src) {

  using src_field_type = decltype(src);

  // Check that underlying value type
  static_assert(std::is_same<non_const_value_type,
                             typename src_field_type::non_const_value_type
                            >::value,
                "Error! Cannot use copy constructor if the underlying value_type is different.\n");
  // Check that destination is const or source is nonconst
  static_assert(std::is_same<value_type,const_value_type>::value ||
                std::is_same<typename src_field_type::value_type,non_const_value_type>::value,
                "Error! Cannot create a nonconst field from a const field.\n");
  if (&src!=*this) {
    m_header    = src.get_header_ptr();
    m_view      = src.get_view();
    m_allocated = src.is_allocated();
  }
}

template<typename ScalarType, typename MemSpace, typename MemManagement>
void Field<ScalarType,MemSpace,MemManagement>::allocate_view ()
{
  // Short names
  const auto& id = m_header->get_identifier();
  auto& alloc_prop = m_header->get_alloc_properties();

  // Check the identifier has all the dimensions set
  error::runtime_check(id.are_dimensions_set(), "Error! Cannot create a field until all the field's dimensions are set.\n");

  // Commit the allocation properties
  alloc_prop.commit();

  // Create the view, by quering allocation properties for the allocation size
  const int view_dim = alloc_prop.get_alloc_size() / sizeof(value_type);

  m_view = view_type(id.name(),view_dim);

  m_allocated = true;
}

} // namespace scream

#endif // SCREAM_FIELD_HPP
