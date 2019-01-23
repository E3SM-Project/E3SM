#ifndef SCREAM_FIELD_HPP
#define SCREAM_FIELD_HPP

#include "field_header.hpp"
#include <share/scream_types.hpp>
#include <share/scream_kokkos_meta.hpp>
#include <share/util/scream_std_meta.hpp>
#include <share/util/scream_kokkos_utils.hpp>

#include <memory>   // For std::shared_ptr

namespace scream
{

// ======================== FIELD ======================== //

// A field should be composed of metadata info (the header) and a pointer to the view

template<typename ScalarType, typename D=DefaultDevice>
class Field {
private:
public:

  using device_type          = D;
  using view_type            = typename KokkosTypes<device_type>::template view<ScalarType*>;
  using value_type           = typename view_type::traits::value_type;
  using const_value_type     = typename view_type::traits::const_value_type;
  using non_const_value_type = typename view_type::traits::non_const_value_type;
  using const_field_type     = Field<const_value_type, device_type>;
  using header_type          = FieldHeader;
  using identifier_type      = header_type::identifier_type;

  // Statically check that ScalarType is not an array.
  static_assert(view_type::Rank==1, "Error! ScalarType should not be an array type.\n");

  // Constructor(s)
  Field () = delete;
  explicit Field (const identifier_type& id);
  explicit Field (const header_type& header);

  // This constructor allows const->const, nonconst->nonconst, and nonconst->const copies
  template<typename SrcDT>
  Field (const Field<SrcDT,device_type>& src);

  // Assignment: allows const->const, nonconst->nonconst, and nonconst->const copies
  template<typename SrcDT>
  Field& operator= (const Field<SrcDT,device_type>& src);

  // ---- Getters ---- //
  const header_type& get_header () const { return *m_header; }
        header_type& get_header ()       { return *m_header; }
  const std::shared_ptr<header_type>& get_header_ptr () const { return m_header; }

  const view_type&   get_view   () const { return  m_view;   }

  template<typename DT>
  ko::Unmanaged<typename KokkosTypes<device_type>::template view<DT> >
  get_reshaped_view () const;

  bool is_allocated () const { return m_allocated; }

  // ---- Setters ---- //

  // Allocate the actual view
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

template<typename ScalarType, typename D>
Field<ScalarType,D>::
Field (const identifier_type& id)
 : m_header    (new header_type(id))
 , m_allocated (false)
{
  // At the very least, the allocation properties need to accommodate this field's value_type.
  m_header->get_alloc_properties().request_value_type_allocation<value_type>();
}

template<typename ScalarType, typename D>
template<typename SrcScalarType>
Field<ScalarType,D>::
Field (const Field<SrcScalarType,D>& src)
 : m_header    (src.get_header_ptr())
 , m_view      (src.get_view())
 , m_allocated (src.is_allocated())
{
  using src_field_type = Field<SrcScalarType,D>;

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

template<typename ScalarType, typename D>
template<typename SrcScalarType>
Field<ScalarType,D>&
Field<ScalarType,D>::
operator= (const Field<SrcScalarType,D>& src) {

  using src_field_type = decltype(src);
#ifndef CUDA_BUILD // TODO Figure out why nvcc isn't like this bit of code.
  // Check that underlying value type
  static_assert(std::is_same<non_const_value_type,
                             typename src_field_type::non_const_value_type
                            >::value,
                "Error! Cannot use copy constructor if the underlying value_type is different.\n");
  // Check that destination is const or source is nonconst
  static_assert(std::is_same<value_type,const_value_type>::value ||
                std::is_same<typename src_field_type::value_type,non_const_value_type>::value,
                "Error! Cannot create a nonconst field from a const field.\n");
#endif
  if (&src!=*this) {
    m_header    = src.get_header_ptr();
    m_view      = src.get_view();
    m_allocated = src.is_allocated();
  }

  return *this;
}

template<typename ScalarType, typename D>
template<typename DT>
ko::Unmanaged<typename KokkosTypes<D>::template view<DT> >
Field<ScalarType,D>::get_reshaped_view () const {
  // The dst value types
  using DstValueType = typename util::ValueType<DT>::type;

  // Get src details
  const auto& alloc_prop = m_header->get_alloc_properties();
  const auto& id = m_header->get_identifier();

  // Make sure input field is allocated
  error::runtime_check(m_allocated, "Error! Cannot reshape a field that has not been allocated yet.\n");

  // Make sure DstDT has an eligible rank: can only reinterpret if the data type rank does not change or if either src or dst have rank 1.
  constexpr int DstRank = util::GetRanks<DT>::rank==1;

  // Check the reinterpret cast makes sense for the two value types (need integer sizes ratio)
  error::runtime_check(alloc_prop.template is_allocation_compatible_with_value_type<DstValueType>(),
                       "Error! Source field allocation is not compatible with the destination field's value type.\n");

  // The destination view type
  using DstView = ko::Unmanaged<typename KokkosTypes<D>::template view<DT> >;
  typename DstView::traits::array_layout layout;

  const int num_values = alloc_prop.get_alloc_size() / sizeof(DstValueType);
  if (DstRank==1) {
    // We are staying 1d, possibly changing the data type
    layout.dimension[0] = num_values;
  } else {
    int num_last_dim_values = num_values;
    // The destination data type is a multi-dimensional array.
    for (int i=0; i<id.rank()-1; ++i) {
      layout.dimension[i] = id.dim(i);

      // Safety check: id.dim(0)*...*id.dim(id.rank()-2) should divide num_values, so we check
      error::runtime_check(num_last_dim_values % id.dim(i) == 0, "Error! Something is wrong with the allocation properties.\n");
      num_last_dim_values /= id.dim(i);
    }
    layout.dimension[id.rank()-1] = num_last_dim_values;
  }

  return DstView (reinterpret_cast<DstValueType*>(m_view.data()),layout);
}

template<typename ScalarType, typename D>
void Field<ScalarType,D>::allocate_view ()
{
  // Not sure if simply returning would be safe enough. Re-allocating
  // would definitely be error prone (someone may have already gotten
  // a subview of the field). However, it *seems* suspicious to call
  // this method twice, and I think it's more likely than not that
  // such a scenario would indicate a bug. Therefore, I am prohibiting it.
  error::runtime_check(!m_allocated, "Error! View was already allocated.\n");

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
