#ifndef SCREAM_FIELD_UTILS_HPP
#define SCREAM_FIELD_UTILS_HPP

#include <share/field/field.hpp>
#include <share/scream_pack.hpp>
#include <typeinfo>

namespace scream {

// ----------------- Free function for field reshaping ----------------- //

// This free functions allows to get and reshape the view of a field,
// changing the data type (with some restrictions).
// There are three use cases:
//  1) change the layout of the view
//  2) change the value_type of the view
//  3) a combination of 1) and 2)
//
// For 1), this is needed cause all Field's store a 1d view internally,
// so if you want a N-dimensional view, you are forced to reshape the
// view stored by the field.
// For 2), this is needed if you need to use the Pack structure to
// enhance vectorization. Fields are stored in an 'unpacked' way
// in the FieldRepository structure, so in order to access Pack's
// features, you first need to reshape the field.
//
// IMPORTANT NOTE
// 
// When changing layout while keeping the same value_type,
// we only need to make sure the layout has the correct
// size (i.e., product of all dimensions should equal the
// result of the size() method of the FieldIdentifier.
// However, when we change the value_type we need to make
// some additional checks, related to padding, which may
// have occurred due to some calls to 'request_value_type_allocation'
// on the FieldAllocProp object stored in the field.
// For instance, let's say the field identifier has rank 2,
// with dimensions (3,10). If one wants support for Pack<Real,4>
// on the field, he/she will call
//   f.get_header().get_alloc_properties().request_value_type_allocation<Pack<Real,4>>()
// This call will make sure that the allocation size is such that
// the last dimension of the field is a multiple of sizeof(Pack<Real,4>).
// In other words, we get the same allocation that we would get
// if the field had dimensions (3,12). And so far, so good.
// The problem arises when we call
//  v = get_reshaped_view<Real[3][10]>(f);
// This call will fail. That's because such view would have
// stride 10 between different rows, while in memory the stride
// is 12 (because an allocation accommodating Pack<Real,4> was
// requested). To overcome this problem, one can keep the
// dimensions non-static (i.e., runtime), and use
//  v = get_reshaped_view<Real**>(f);
// This will 

// 1) get an N-dimensional view from the 1d view in the field.
//    Here, N *must* match the rank of the field, as stated in
//    the FieldIdentifier stored in the field.
// 2) changed the value type of the view, going from a packed to an
//    unpacked view (or viceversa). This changes the 'extent' of the
//    view
//    
//    (and viceversa). The main use of this routine is for the FieldRepository,
//    which always stores 1d fields (otherwise it would need
//    one internal container for each layout). When the fields
//    are extracted from the repo, for efficiency, they might
//    need to be converted to a data type fully expanded.
//    E.g., the repo may store a field as double*, even if
//    its layout is (3,4,6). Upon extraction from the repo,
//    one can call reinterpret_field<double[3][4][6]>(f),
//    to have access to a view that is fully compile-time shaped.
// 2) go from a packed to an unpacked scalar type (or viceversa),
//    or to a packed scalar type with different pack length.
//    E.g., go from
//        Field<Pack<double,4>*[4]>
//    to
//        Field<double*[16],MS>
//    or even
//        Field<Pack<double,2>*[8]>
//
// 3) a combination of 1) and 2).

template<typename DstDT, typename ScalarType, typename MemSpace, typename MemManagement>
ViewType<DstDT,MemSpace,MemoryUnmanaged>
get_reshape_view (const Field<ScalarType,MemSpace,MemManagement>& src) {
  // The dst value types
  using DstValueType = typename util::ValueType<DstDT>::type;

  // Get src details
  const auto& src_alloc_prop = src.get_header().get_alloc_properties();
  const auto& src_id = src.get_header().get_identifier();

  // Make sure input field is allocated
  error::runtime_check(src.is_allocated(), "Error! Cannot reshape a field that has not been allocated yet.\n");

  // Make sure DstDT has an eligible rank: can only reinterpret if the data type rank does not change or if either src or dst have rank 1.
  constexpr int DstRank = util::GetRanks<DstDT>::rank==1;
  constexpr int SrcRank = util::GetRanks<SrcDT>::rank==1;
  static_assert(DstRank==SrcRank || DstRank==1 || SrcRank==1,
                "Error! You can only reinterpret_field if the src and dst ranks match or if one of the two is 1.\n");

  // Check the reinterpret cast makes sense for the two value types (need integer sizes ratio)
  error::runtime_check(src_alloc_prop.template is_allocation_compatible_with_value_type<DstValueType>(),
                       "Error! Source field allocation is not compatible with the destination field's value type.\n");

  // The destination field and view types
  using DstField = Field<DstDT,MemSpace,MemoryUnmanaged>;
  using DstView  = typename DstField::view_type;

  // The value types may differ, so we'll need to adjust the view's layout dimensions
  typename DstView::traits::array_layout layout;

  if (DstRank==1) {
    // We are going to 1d, possibly changing the data type
    layout.dimension[0] = src_alloc_prop.get_alloc_size() / sizeof(DstValueType);
  } else {
    // The destination data type is a multi-dimensional array.
    for (int i=0; i<src_id.rank()-1; ++i) {
      layout.dimension[i] = src_id.dim(i);
    }
    layout.dimension[src_id.rank()-1] = src_alloc_prop.get_fast_index_alloc_size() / sizeof(DstValueType);
  }

  DstView dst_view (reinterpret_cast<DstValueType*>(src.get_view().data()),layout);

  return DstField(src.get_header_ptr(),dst_view);
}

} // namespace scream

#endif // SCREAM_FIELD_UTILS_HPP
