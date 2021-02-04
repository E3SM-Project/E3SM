#ifndef SCREAM_FIELD_HPP
#define SCREAM_FIELD_HPP

#include "share/field/field_header.hpp"
#include "share/field/field_property_check.hpp"
#include "share/util/pointer_list.hpp"
#include "share/scream_types.hpp"

#include "ekat/ekat_type_traits.hpp"
#include "ekat/kokkos/ekat_kokkos_meta.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/kokkos/ekat_subview_utils.hpp"

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

  // The syntax of std::enable_if is way too long...
  template<bool c, typename T, typename F>
  using cond_t = typename std::conditional<c,T,F>::type;
  template<bool c, typename T>
  using if_t = typename std::enable_if<c,T>::type;

  // Various kokkos-related types
  using device_type = DefaultDevice;
  using kt          = KokkosTypes<device_type>;

  template<typename ViewT>
  using HM = typename ViewT::HostMirror;

  // Given a type Space, check if it is the device mem/exec space,
  // or even the device itself.
  template<typename Space>
  struct OnDevice {
    static constexpr bool value =
          std::is_same<Space,device_type::memory_space>::value ||
          std::is_same<Space,device_type::execution_space>::value ||
          std::is_same<Space,device_type>::value;
  };

  // Given a device view type V and a mem/exec space, returns V if the space
  // is on device, otherwise returns the host mirror of V
  template<typename DevViewT, typename Space>
  using get_view_type = cond_t<OnDevice<Space>::value,DevViewT,HM<DevViewT>>;

  template<typename DT>
  using view_type = typename kt::template view<DT>;

  template<typename DT>
  using uview_type = ekat::Unmanaged<view_type<DT>>;

  template<typename T, int N>
  using view_ND_type = uview_type<typename ekat::DataND<T,N>::type>;

  // These are used a few times in ctor and assignment, so use a shortcut
  template<typename T>
  using Const = typename std::add_const<T>::type;
  template<typename T>
  using NonConst = typename std::remove_const<T>::type;

  // Field stack classes types
  using RT                   = RealType;
  using header_type          = FieldHeader;
  using identifier_type      = FieldIdentifier;
  using field_type           = Field<RT>;
  using const_field_type     = Field<Const<RT>>;
  using nonconst_field_type  = Field<NonConst<RT>>;

  static constexpr int MaxRank = 6;

  // To make assignment impl simpler
  template<typename SrcRT>
  friend class Field;

  // Statically check that RealType is a legit numeric type.
  static_assert(std::is_floating_point<RT>::value,
                "Error! RealType should be a floating point type.\n");

  // A Field maintains a list of shared_ptrs to FieldPropertyChecks that can
  // determine whether it satisfies certain properties. We use the PointerList
  // class to provide simple access to the property checks.
  using property_check_type = FieldPropertyCheck<RealType>;
  using property_check_list = PointerList<std::shared_ptr<property_check_type>,
                                           property_check_type>;
  using property_check_iterator = typename property_check_list::iterator;

  // Constructor(s)
  Field () = default;
  Field (const std::shared_ptr<header_type>& h, const view_type<RT*>& v);
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

  template<typename Space = device_type::memory_space>
  const get_view_type<view_type<RT*>,Space>&
  get_view () const { return  get_view_impl<Space>();   }

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
  template<typename DT, typename Space = device_type::memory_space>
  get_view_type<uview_type<DT>,Space>
  get_reshaped_view () const;

  // If someone needs the host view, some sync routines might be needed.
  // Note: this class takes no responsibility in keeping track of whether
  //       a sync is required in either direction. Mainly because we expect
  //       host views to be seldom used, and even less frequently modified.
  void sync_to_host () const;
  void sync_to_dev () const;

  // Returns a subview of this field, slicing at entry k along dimension idim
  // NOTES:
  //   - the output field stores *the same* 1d view as this field. In order
  //     to get the N-1 dimensional view, call get_reshaped_view<DT>(), using
  //     the correct N-1 dimensional data type DT.
  //   - when calling get_reshaped_view<DT>() on the N-1 dimensional subfield,
  //     we first get an N-dimensional view, then subview it at index k along
  //     dimension idim.
  //   - idim must be either 0 or 1. This is b/c we cannot subview an N-dim
  //     view along idim=2+ while keeping LayoutRight. Kokkos would force the
  //     resulting view to have layout stride, which would conflict with the
  //     return type of get_reshaped_view<DT>().
  //   - If the field rank is 2, then idim cannot be 1. This is b/c Kokkos
  //     specializes view's traits for LayoutRight of rank 1, not allowing
  //     to store a stride for the slowest dimension.
  field_type subfield (const std::string& sf_name, const int idim, const int k) const;
  field_type subfield (const int idim, const int k) const;

  // Checks whether the underlying view has been already allocated.
  bool is_allocated () const { return m_view_d.data()!=nullptr; }

  // ---- Setters and non-const methods ---- //

  // Allocate the actual view
  void allocate_view ();

protected:

  template<typename Space>
  const if_t<OnDevice<Space>::value,view_type<RT*>>&
  get_view_impl   () const {
    EKAT_REQUIRE_MSG (is_allocated (),
        "Error! View was not yet allocated.\n");
    return  m_view_d;
  }

  template<typename Space>
  const if_t<!OnDevice<Space>::value,HM<view_type<RT*>>>&
  get_view_impl   () const {
    ensure_host_view ();
    return  m_view_h;
  }

  void ensure_host_view () const {
    EKAT_REQUIRE_MSG (is_allocated (),
        "Error! View was not yet allocated.\n");
    if (m_view_h.data()==nullptr) {
      m_view_h = Kokkos::create_mirror_view(m_view_d);
    }
  }

  // These SFINAE impl of get_subview are needed since subview_1 does not
  // exist for rank2 (or less) views.
  template<typename Space, typename T, int N>
  if_t<(N>2),
       get_view_type<view_ND_type<T,N-1>,Space>>
  get_subview_1 (const get_view_type<view_ND_type<T,N>,Space>& v, const int k) const {
    return ekat::subview_1(v,k);
  }


  template<typename Space, typename T, int N>
  if_t<(N<=2),
       get_view_type<view_ND_type<T,N-1>,Space>>
  get_subview_1 (const get_view_type<view_ND_type<T,N>,Space>&, const int) const {
    get_view_type<view_ND_type<T,N-1>,Space> sv;
    EKAT_ERROR_MSG ("Error! Cannot subview a rank2 view along the second dimension without losing LayoutRight.\n");
    return sv;
  }

  template<typename Space,typename T,int N>
  if_t<N==MaxRank,
       get_view_type<view_ND_type<T,N>,Space>>
  get_ND_view () const;

  template<typename Space,typename T,int N>
  if_t<(N<MaxRank),
       get_view_type<view_ND_type<T,N>,Space>>
  get_ND_view () const;

  // Metadata (name, rank, dims, customere/providers, time stamp, ...)
  std::shared_ptr<header_type>            m_header;

  // Actual data.
  view_type<NonConst<RT*>>                m_view_d;

  // Host mirror of the data. Mutable, to allow lazy creation
  mutable HM<view_type<NonConst<RT*>>>    m_view_h;

  // Keep track of whether the field has been allocated
  bool                                    m_allocated;

  // List of property checks for this field.
  std::shared_ptr<property_check_list>    m_prop_checks;
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
 : m_header     (create_header(id))
 , m_prop_checks(new property_check_list)
{
  // At the very least, the allocation properties need to accommodate this field's real type.
  m_header->get_alloc_properties().request_allocation<RT>();
}

template<typename RealType>
Field<RealType>::
Field (const std::shared_ptr<header_type>& h, const view_type<RT*>& v)
  : m_header(h)
  , m_view_d(v)
  , m_prop_checks(new property_check_list)
{
  // Nothing to do here
}

template<typename RealType>
template<typename SrcRealType>
Field<RealType>::
Field (const Field<SrcRealType>& src)
 : m_header (src.get_header_ptr())
 , m_view_d (src.get_view())
{
  using src_field_type = Field<SrcRealType>;

  // NOTE: the following checks might be redundant, since Kokkos::View copy
  //       constructor might already perform analogue checks in order to
  //       assign pointers.

  // Check that underlying value type
  static_assert(std::is_same<NonConst<RT>,NonConst<typename src_field_type::RT>>::value,
                "Error! Cannot use copy constructor if the underlying real type is different.\n");
  // Check that destination is const or source is nonconst
  static_assert(std::is_same<RT,Const<RT>>::value ||
                std::is_same<typename src_field_type::RT,NonConst<RT>>::value,
                "Error! Cannot create a nonconst field from a const field.\n");
}

template<typename RealType>
template<typename SrcRealType>
Field<RealType>&
Field<RealType>::
operator= (const Field<SrcRealType>& src) {

  using src_field_type = typename std::remove_reference<decltype(src)>::type;
#ifndef CUDA_BUILD // TODO Figure out why nvcc doesn't like this bit of code.
  // Check that underlying value type
  static_assert(
      std::is_same<NonConst<RT>,NonConst<typename src_field_type::RT>>::value,
      "Error! Cannot use copy constructor if the underlying value_type is different.\n");
  // Check that destination is const or source is nonconst
  static_assert(
      std::is_same<RT,Const<RT>>::value ||
      std::is_same<typename src_field_type::RT,NonConst<RT>>::value,
      "Error! Cannot create a nonconst field from a const field.\n");
#endif

  // If the field has a valid header (i.e., name!=""), then
  // we only allow assignment of fields with the *same* identifier.
  EKAT_REQUIRE_MSG(
      m_header->get_identifier().get_id_string()=="" ||
      m_header->get_identifier()==src.get_header().get_identifier(),
      "Error! Assignment of fields with different (and non-null) identifiers is prohibited.\n");

  // Since the type of *this and src may be different, we cannot do the usual
  // `if (this!=&src)`, cause the compiler cannot compare those pointers.
  // Therefore, we compare the stored pointers.
  // If either header or view are different, we copy everything
  if (m_header!=src.m_header ||
      m_view_d.data()!=src.m_view_d.data()) {
    m_header = src.m_header;
    m_view_d = src.m_view_d;
    m_view_h = src.m_view_h;
  }

  return *this;
}

template<typename RealType>
template<typename DT, typename Space>
auto Field<RealType>::get_reshaped_view () const
 -> get_view_type<uview_type<DT>,Space>
{
  // The destination view type on correct mem space
  using DstView = get_view_type<uview_type<DT>,Space>;
  // The dst value types
  using DstValueType = typename DstView::traits::value_type;
  // The ViewDimension object from the Dst View (used to check validity of possible compile-time extents)
  using dims_type = typename DstView::traits::dimension;

  // Get src details
  const auto& alloc_prop = m_header->get_alloc_properties();
  const auto& field_layout = m_header->get_identifier().get_layout();

  // We only allow to reshape to another 1d view (possibly to change the value type
  // to something like a pack), or to a view of the correct rank
  constexpr int DstRank = DstView::rank;

  EKAT_REQUIRE_MSG(DstRank==field_layout.rank(),
      "Error! You can only reshape to a view of the correct rank (equal to the FieldLayout's one).\n");

  // Check the reinterpret cast makes sense for the Dst value types (need integer sizes ratio)
  EKAT_REQUIRE_MSG(alloc_prop.template is_compatible<DstValueType>(),
      "Error! Source field allocation is not compatible with the requested value type.\n");

  // Make sure input field is allocated
  EKAT_REQUIRE_MSG(is_allocated(),
      "Error! Cannot reshape a field that has not been allocated yet.\n");

  // Start by reshaping into a view with all dyn extents
  const auto view_ND = get_ND_view<Space,DstValueType,DstRank>();

  constexpr int DstRankDynamic= DstView::rank_dynamic;

  using dyn_DT = typename decltype(view_ND)::traits::data_type;
  if (!std::is_same<dyn_DT,DT>::value) {
    // The user requested some compile-time dimensions.
    // Let's check that they are correct
    for (int i=DstRankDynamic; i<DstRank; ++i) {
      EKAT_REQUIRE_MSG(view_ND.extent(i)==dims_type::static_extent(i),
          "Error! The template DataType contains some invalid compile-time dimensions.\n");
    }
  }

  // Before building the DstView from view_ND, we have one more check:
  // if DstRankDynamic==0, kokkos specializes the view offset struct,
  // assuming *no* stride. That's fine, as long as this field alloc
  // props ensure that there is no stride
  EKAT_REQUIRE_MSG (DstRankDynamic>0 || alloc_prop.contiguous(),
      "Error! Cannot use all compile-time dimensions for strided views.\n");

  return DstView(view_ND);
}

template<typename RealType>
void Field<RealType>::
sync_to_host () const {
  // Sanity check
  EKAT_REQUIRE_MSG (is_allocated(),
      "Error! Input field must be allocated in order to sync host and device views.\n");

  // Ensure host view was created (lazy construction)
  ensure_host_view ();
  Kokkos::deep_copy(m_view_h,m_view_d);
}

template<typename RealType>
void Field<RealType>::
sync_to_dev () const {
  // Sanity check
  EKAT_REQUIRE_MSG (is_allocated(),
      "Error! Input field must be allocated in order to sync host and device views.\n");

  // Ensure host view was created (lazy construction)
  ensure_host_view ();
  Kokkos::deep_copy(m_view_d,m_view_h);
}

template<typename RealType>
Field<RealType> Field<RealType>::
subfield (const std::string& sf_name, const int idim, const int k) const {
  const auto& id = m_header->get_identifier();
  const auto& lt = id.get_layout();

  // Sanity checks
  EKAT_REQUIRE_MSG (is_allocated(),
      "Error! Input field must be allocated in order to subview it.\n");
  EKAT_REQUIRE_MSG (idim==0 || idim==1,
        "Error! Subview dimension index must be either 0 or 1.\n");

  // Create identifier for subfield
  std::vector<FieldTag> tags = lt.tags();
  std::vector<int> dims = lt.dims();
  tags.erase(tags.begin()+idim);
  dims.erase(dims.begin()+idim);
  FieldLayout sf_lt(tags,dims);
  FieldIdentifier sf_id(sf_name,sf_lt,id.get_units(),id.get_grid_name());

  // Create header
  auto sv_h = create_header(sf_id,m_header,idim,k);

  // Return subfield
  return field_type(sv_h,m_view_d);
}

template<typename RealType>
Field<RealType> Field<RealType>::
subfield (const int idim, const int k) const {
  return subfield(m_header->get_identifier().name(),idim,k);
}

template<typename RealType>
void Field<RealType>::allocate_view ()
{
  // Not sure if simply returning would be safe enough. Re-allocating
  // would definitely be error prone (someone may have already gotten
  // a subview of the field). However, it *seems* suspicious to call
  // this method twice, and I think it's more likely than not that
  // such a scenario would indicate a bug. Therefore, I am prohibiting it.
  EKAT_REQUIRE_MSG(!is_allocated(), "Error! View was already allocated.\n");

  // Short names
  const auto& id     = m_header->get_identifier();
  const auto& layout = id.get_layout_ptr();
  auto& alloc_prop   = m_header->get_alloc_properties();

  // Check the identifier has all the dimensions set
  EKAT_REQUIRE_MSG(layout->are_dimensions_set(),
      "Error! Cannot allocate the view until all the field's dimensions are set.\n");

  // Commit the allocation properties
  alloc_prop.commit(layout);

  // Create the view, by quering allocation properties for the allocation size
  const int view_dim = alloc_prop.get_alloc_size() / sizeof(RT);

  m_view_d = decltype(m_view_d)(id.name(),view_dim);
  Kokkos::deep_copy(m_view_d,0.0);
}

template<typename RealType>
template<typename Space,typename T,int N>
auto Field<RealType>::get_ND_view () const ->
  if_t<(N<MaxRank),get_view_type<view_ND_type<T,N>,Space>>
{
  const auto& fl = m_header->get_identifier().get_layout();
  EKAT_REQUIRE_MSG (N==1 || N==fl.rank(),
      "Error! Input Rank must either be 1 (flat array) or the actual field rank.\n");

  // Check if this field is a subview of another field
  const auto parent = m_header->get_parent().lock();
  if (parent!=nullptr) {
    // Parent field has correct layout to reinterpret the view into N+1-dim view
    Field<RealType> f(parent,m_view_d);
    auto v_np1 = f.get_ND_view<Space,T,N+1>();

    // Now we can subview v_np1 at the correct slice
    const auto& idx = m_header->get_alloc_properties().get_subview_idx();
    const int idim = idx.first;
    const int k    = idx.second;

    // So far we can only subview at first or second dimension.
    EKAT_REQUIRE_MSG (idim==0 || idim==1,
        "Error! Subview dimension index is out of bounds.\n");

    EKAT_REQUIRE_MSG (idim==0 || N>1,
        "Error! Cannot subview a rank-2 (or less) view along 2nd dimension without losing LayoutRight.\n");

    // Use SFINAE-ed get_subview helper function to pick correct
    // subview impl. If N+1<=2 and idim!=0, the code craps out in the check above.
    if (idim==0) {
      return ekat::subview(v_np1,k);
    } else {
      return get_subview_1<Space,T,N+1>(v_np1,k);
    }
  }

  // Compute extents from FieldLayout
  const auto& alloc_prop = m_header->get_alloc_properties();
  size_t num_values = alloc_prop.get_alloc_size() / sizeof(T);
  Kokkos::LayoutRight kl;
  for (int i=0; i<N-1; ++i) {
    kl.dimension[i] = fl.dim(i);
    num_values /= fl.dim(i);
  }
  kl.dimension[N-1] = num_values;
  auto ptr = reinterpret_cast<T*>(get_view<Space>().data());

  using ret_type = get_view_type<view_ND_type<T,N>,Space>;
  return ret_type (ptr,kl);
}

template<typename RealType>
template<typename Space,typename T,int N>
auto Field<RealType>::get_ND_view () const ->
  if_t<N==MaxRank,get_view_type<view_ND_type<T,N>,Space>>
{
  const auto& fl = m_header->get_identifier().get_layout();
  EKAT_REQUIRE_MSG (N==1 || N==fl.rank(),
      "Error! Input Rank must either be 1 (flat array) or the actual field rank.\n");

  // Given that N==MaxRank, this field cannot be a subview of another field
  EKAT_REQUIRE_MSG (m_header->get_parent().expired(),
      "Error! A view of rank 8 should not be the subview of another field.\n");

  // Compute extents from FieldLayout
  const auto& alloc_prop = m_header->get_alloc_properties();
  size_t num_values = alloc_prop.get_alloc_size() / sizeof(T);
  Kokkos::LayoutRight kl;
  for (int i=0; i<N-1; ++i) {
    kl.dimension[i] = fl.dim(i);
    num_values /= fl.dim(i);
  }
  kl.dimension[N-1] = num_values;
  auto ptr = reinterpret_cast<T*>(get_view<Space>().data());

  using ret_type = get_view_type<view_ND_type<T,N>,Space>;
  return ret_type (ptr,kl);
}

} // namespace scream

#endif // SCREAM_FIELD_HPP
