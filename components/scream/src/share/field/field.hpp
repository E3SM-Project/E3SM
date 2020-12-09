#ifndef SCREAM_FIELD_HPP
#define SCREAM_FIELD_HPP

#include "share/field/field_header.hpp"
#include "share/field/field_property_check.hpp"
#include "share/util/pointer_list.hpp"
#include "share/scream_types.hpp"

#include "ekat/ekat_type_traits.hpp"
#include "ekat/kokkos/ekat_kokkos_meta.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"

#include <memory>   // For std::shared_ptr

namespace scream
{

template<typename FieldType>
struct is_scream_field : public std::false_type {};

// ======================== FIELD ======================== //

// A field is composed of metadata info (the header) and a pointer to a view.
// Fields are always stored as 1D arrays of real-valued data. The associated
// view can be reshaped as needed to match a desired layout for a given client.
template<typename RealType>
class Field {
public:

  using value_type           = RealType;
  using device_type          = DefaultDevice;
  using header_type          = FieldHeader;
  using identifier_type      = header_type::identifier_type;
  using view_type            = typename KokkosTypes<device_type>::template view<value_type*>;
  using host_view_type       = typename KokkosTypes<HostDevice>::template view<value_type*>;
  using const_value_type     = typename std::add_const<value_type>::type;
  using non_const_value_type = typename std::remove_const<value_type>::type;

  using field_type           = Field<value_type>;
  using const_field_type     = Field<const_value_type>;
  using nonconst_field_type  = Field<non_const_value_type>;

  // Statically check that RealType is not an array.
  static_assert(view_type::Rank==1, "Error! RealType should not be an array type.\n");

  // A Field maintains a list of shared_ptrs to FieldPropertyChecks that can
  // determine whether it satisfies certain properties. We use the PointerList
  // class to provide simple access to the property checks.
  using property_check_type = FieldPropertyCheck<RealType>;
  using property_check_list = PointerList<std::shared_ptr<property_check_type>,
                                           property_check_type>;
  using property_check_iterator = typename property_check_list::iterator;

  // Constructor(s)
  Field () = default;
  explicit Field (const identifier_type& id);

  // This constructor allows const->const, nonconst->nonconst, and nonconst->const copies
  template<typename SrcDT>
  Field (const Field<SrcDT>& src);

  // Assignment: allows const->const, nonconst->nonconst, and nonconst->const copies
  template<typename SrcDT>
  Field& operator= (const Field<SrcDT>& src);

  // ---- Getters and const methods---- //
  const header_type& get_header () const { return *m_header; }
        header_type& get_header ()       { return *m_header; }
  const std::shared_ptr<header_type>& get_header_ptr () const { return m_header; }

  const view_type& get_view   () const { return  m_view;   }

  // Returns a const_field_type copy of this field
  const_field_type get_const () const { return const_field_type(*this); }

  // Adds a propery check to this field.
  void add_property_check(std::shared_ptr<property_check_type> property_check) {
    m_prop_checks->append(property_check);
  }

  // These (forward) iterators allow access to the set of property checks on the
  // field.
  property_check_iterator property_check_begin() const {
    return m_prop_checks->begin();
  }
  property_check_iterator property_check_end() const {
    return m_prop_checks->end();
  }

  // Allows to get the underlying view, reshaped for a different data type.
  // The class will check that the requested data type is compatible with the
  // allocation. This allows each field to be stored as a 1d array, but then
  // be reshaped to the desired layout before being used.
  template<typename DT>
  ekat::Unmanaged<typename KokkosTypes<device_type>::template view<DT> >
  get_reshaped_view () const;

  // Checks whether the underlying view has been already allocated.
  bool is_allocated () const { return m_allocated; }

  // ---- Setters and non-const methods ---- //

  // Allocate the actual view
  void allocate_view ();

protected:

  // Metadata (name, rank, dims, customere/providers, time stamp, ...)
  std::shared_ptr<header_type>    m_header;

  // Actual data.
  view_type                       m_view;

  // Keep track of whether the field has been allocated
  bool                            m_allocated;

  // List of property checks for this field.
  std::shared_ptr<property_check_list> m_prop_checks;
};

template<typename RealType>
struct is_scream_field<Field<RealType> > : public std::true_type {};

template<typename RealType>
bool operator< (const Field<RealType>& f1, const Field<RealType>& f2) {
  return f1.get_header().get_identifier() < f2.get_header().get_identifier();
}

// ================================= IMPLEMENTATION ================================== //

template<typename RealType>
Field<RealType>::
Field (const identifier_type& id)
 : m_header    (new header_type(id))
 , m_allocated (false)
 , m_prop_checks(new property_check_list)
{
  // At the very least, the allocation properties need to accommodate this field's value_type.
  m_header->get_alloc_properties().request_value_type_allocation<value_type>();
}

template<typename RealType>
template<typename SrcRealType>
Field<RealType>::
Field (const Field<SrcRealType>& src)
 : m_header    (src.get_header_ptr())
 , m_view      (src.get_view())
 , m_allocated (src.is_allocated())
{
  using src_field_type = Field<SrcRealType>;

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

template<typename RealType>
template<typename SrcRealType>
Field<RealType>&
Field<RealType>::
operator= (const Field<SrcRealType>& src) {

  using src_field_type = typename std::remove_reference<decltype(src)>::type;
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

  // If the field has a valid header, we only allow assignment of fields with
  // the *same* identifier.
  EKAT_REQUIRE_MSG(m_header->get_identifier().get_id_string()=="" ||
                     m_header->get_identifier()==src.get_header().get_identifier(),
                     "Error! Assignment of fields with different (and non-null) identifiers is prohibited.\n");

  // Since the type of *this and src may be different, we cannot do the usual
  // `if (this!=&src)`, cause the compiler cannot compare those pointers.
  // Therefore, we compare the stored pointers. Note that we don't compare
  // the 'm_allocated' member, cause its superfluous.
  // If either header or view are different, we copy everything
  if (m_header!=src.get_header_ptr() ||
      m_view!=src.get_view()) {
    m_header    = src.get_header_ptr();
    m_view      = src.get_view();
    m_allocated = src.is_allocated();
  }

  return *this;
}

template<typename RealType>
template<typename DT>
ekat::Unmanaged<typename KokkosTypes<DefaultDevice>::template view<DT> >
Field<RealType>::get_reshaped_view () const {
  // The dst value types
  using DstValueType = typename ekat::ValueType<DT>::type;

  // Get src details
  const auto& alloc_prop = m_header->get_alloc_properties();
  const auto& field_layout = m_header->get_identifier().get_layout();

  // Make sure input field is allocated
  EKAT_REQUIRE_MSG(m_allocated, "Error! Cannot reshape a field that has not been allocated yet.\n");

  // Make sure DstDT has an eligible rank: can only reinterpret if the data type rank does not change or if either src or dst have rank 1.
  constexpr int DstRank = ekat::GetRanks<DT>::rank;

  // Check the reinterpret cast makes sense for the two value types (need integer sizes ratio)
  EKAT_REQUIRE_MSG(alloc_prop.template is_allocation_compatible_with_value_type<DstValueType>(),
                       "Error! Source field allocation is not compatible with the destination field's value type.\n");

  // The destination view type
  using DstView = ekat::Unmanaged<typename KokkosTypes<DefaultDevice>::template view<DT> >;
  typename DstView::traits::array_layout kokkos_layout;

  const int num_values = alloc_prop.get_alloc_size() / sizeof(DstValueType);
  if (DstRank==1) {
    // We are staying 1d, possibly changing the data type
    kokkos_layout.dimension[0] = num_values;
  } else {
    int num_last_dim_values = num_values;
    // The destination data type is a multi-dimensional array.
    for (int i=0; i<field_layout.rank()-1; ++i) {
      kokkos_layout.dimension[i] = field_layout.dim(i);

      // Safety check: field_layout.dim(0)*...*field_layout.dim(field_layout.rank()-2) should divide num_values, so we check
      EKAT_REQUIRE_MSG(num_last_dim_values % field_layout.dim(i) == 0, "Error! Something is wrong with the allocation properties.\n");
      num_last_dim_values /= field_layout.dim(i);
    }
    kokkos_layout.dimension[field_layout.rank()-1] = num_last_dim_values;
  }

  return DstView (reinterpret_cast<DstValueType*>(m_view.data()),kokkos_layout);
}

template<typename RealType>
void Field<RealType>::allocate_view ()
{
  // Not sure if simply returning would be safe enough. Re-allocating
  // would definitely be error prone (someone may have already gotten
  // a subview of the field). However, it *seems* suspicious to call
  // this method twice, and I think it's more likely than not that
  // such a scenario would indicate a bug. Therefore, I am prohibiting it.
  EKAT_REQUIRE_MSG(!m_allocated, "Error! View was already allocated.\n");

  // Short names
  const auto& id     = m_header->get_identifier();
  const auto& layout = id.get_layout();
  auto& alloc_prop   = m_header->get_alloc_properties();

  // Check the identifier has all the dimensions set
  EKAT_REQUIRE_MSG(layout.are_dimensions_set(), "Error! Cannot allocate the view until all the field's dimensions are set.\n");

  // Commit the allocation properties
  alloc_prop.commit(layout);

  // Create the view, by quering allocation properties for the allocation size
  const int view_dim = alloc_prop.get_alloc_size() / sizeof(value_type);

  m_view = view_type(id.name(),view_dim);

  m_allocated = true;
}

} // namespace scream

#endif // SCREAM_FIELD_HPP
