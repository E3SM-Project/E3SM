#ifndef SCREAM_FIELD_UTILS_HPP
#define SCREAM_FIELD_UTILS_HPP

#include <share/field/field.hpp>
#include <share/scream_pack.hpp>
#include <typeinfo>

namespace scream {

// ----------------- Free functions for field reshaping ----------------- //

// This free functions allows to create an unmanaged view of a Field
// changing the layout of the field or its data type.
// There are two use cases:
// 1) go from N-dimensional data type to a 1-dimensional data type
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
// NOTE: if you need to go from, say, Field<Pack<double,4>*[4]> to a 1d
//       Field<double*> field, do it in two steps

template<typename DstDT, typename SrcDT, typename MemSpace, typename MemManagement>
Field<DstDT,MemSpace,MemoryUnmanaged>
reinterpret_field (const Field<SrcDT,MemSpace,MemManagement>& src) {
  // The src and dst value types
  using SrcValueType = typename util::ValueType<SrcDT>::type;
  using DstValueType = typename util::ValueType<DstDT>::type;

  // Check if they are packs or not
  using IsSrcPack = util::is_pack<SrcValueType>;
  using IsDstPack = util::is_pack<DstValueType>;

  // Get the scalar type (if not packs, then it's the same as value type)
  using SrcScalarType = typename IsSrcPack::scalar_type;
  using DstScalarType = typename IsDstPack::scalar_type;

  // Make sure input field is allocated
  error::runtime_check(src.get_alloc_properties().allocated, "Error! Cannot reshape a field that has not been allocated yet.\n");

  // Check the reinterpret cast makes sense for the two value types (need integer sizes ratio)
  constexpr int srcValueTypeSize = sizeof(SrcValueType);
  constexpr int dstValueTypeSize = sizeof(DstValueType);
  constexpr bool sizesRatioIsInteger = srcValueTypeSize>dstValueTypeSize ?
                                       (srcValueTypeSize % dstValueTypeSize == 0) :
                                       (dstValueTypeSize % srcValueTypeSize == 0);
  static_assert (sizesRatioIsInteger, "Error! Cannot reinterpret a field to a data type whose value_type's size is not a multiple or a divisor of the original one.\n");

  constexpr int sizesRatio = (srcValueTypeSize>dstValueTypeSize ?
                              srcValueTypeSize/dstValueTypeSize :
                              dstValueTypeSize/srcValueTypeSize);

  // Detect if one of the fields just stores a 1d pointer
  constexpr bool isSrc1D = std::is_same<SrcDT,SrcValueType*>::value;
  constexpr bool isDst1D = std::is_same<DstDT,DstValueType*>::value;

  // We only allow to call this routine for a few scenarios:
  //  1) DstDT = SrcDT
  //  2) DstDT = Pack<T,N>, SrcDT = T, and the other way around
  //  3) SrcDT = Pack<T,N>, DstDT = Pack<T,M>, with N%M=0 or M%N=0
  //  4) SrcDT = DstValueType* or viceversa
  // Case 1) is trivial (simply return src). Case 2 corresponds
  // to the case where a packed field is flattened (or viceversa).
  // Case 3) is more niche (and may not even ever happen), and
  // corresponds to the case where we change the pack size.
  // Case 4) is needed to be able to use a single container to store
  // fields with different layout (by reshaping all of them to 1d fields)

  constexpr bool check = std::is_same<SrcDT,DstDT>::value ||                                        // Either the same DataType,...
                        (IsSrcPack::value && std::is_same<SrcScalarType,DstValueType>::value) ||   // ...or pack<scalar,N> -> scalar
                        (IsDstPack::value && std::is_same<DstScalarType,SrcValueType>::value) ||   // ...or scalar -> pack<scalar,N>,...
                        (IsDstPack::value && IsSrcPack::value &&                                   // ...or pack<scalar,M> -> pack<scalar,N>
                         std::is_same<DstScalarType,SrcScalarType>::value) ||                      //    (provided common scalar type),...
                        (isSrc1D && std::is_same<SrcDT,DstValueType*>::value) ||            // ...or 1d -> Nd array,...
                        (isDst1D && std::is_same<DstDT,SrcValueType*>::value);              // ...or Nd -> 1d array,...

  static_assert (check, "Error! Unsupported use of reinterpret_field.\n");

  // The value types differ, so we'll need to adjust the view's layout dimensions
  using DstField = Field<DstDT,MemSpace,MemoryUnmanaged>;
  using DstView  = typename DstField::view_type;
  typename DstView::traits::array_layout layout;
  const auto& srcId = src.get_header().get_identifier();

  if (isDst1D) {
    // We are going to 1d, possibly changing the data type
    if (srcValueTypeSize>dstValueTypeSize) {
      // Dst value type is smaller => more entries
      layout.dimension[0] = srcId.size() * sizesRatio; 
    } else {
      // Dst value type is larger (or equal)  => less (or equal) entries
      layout.dimension[0] = srcId.size() / sizesRatio; 
    }
  } else {
    // We are going to an Nd array, possibly changing the data type

    // The first rank-1 dimensions will be the ones from the identifier.
    // The last ones may need to be adjusted, if one of the two data types
    // is a pack (or they both are, but with different pack sizes).
    for (int idim=0; idim<srcId.rank()-1; ++idim) {
      layout.dimension[idim] = srcId.dim(idim);
    }
    if (dstValueTypeSize>srcValueTypeSize) {
      // And to a larger valye type
      layout.dimension[srcId.rank()-1] = srcId.dims().back() / sizesRatio;
    } else {
      // And to an either smaller or equal value_type
      layout.dimension[srcId.rank()-1] = srcId.dims().back() * sizesRatio;
    }
  }

  DstView dst_view (reinterpret_cast<DstValueType*>(src.get_view().data()),layout);

  return DstField(src.get_header_ptr(),src.get_alloc_properties(),dst_view);
}

} // namespace scream

#endif // SCREAM_FIELD_UTILS_HPP
