#ifndef SCREAM_FIELD_IMPL_HPP
#define SCREAM_FIELD_IMPL_HPP

#include "share/field/field.hpp"

#include <array>

namespace scream
{

template<typename DT, HostOrDevice HD>
auto Field::get_view () const
 -> get_view_type<DT,HD>
{
  // The destination view type on correct mem space
  using DstView = get_view_type<DT,HD>;
  // The dst value types
  using DstValueType = typename DstView::traits::value_type;
  // The ViewDimension object from the Dst View (used to check validity of
  // possible compile-time extents)
  using dims_type = typename DstView::traits::dimension;
  // We only allow to reshape to a view of the correct rank
  constexpr int DstRank = DstView::rank;
  constexpr int DstRankDynamic= DstView::rank_dynamic;

  // Make sure input field is allocated
  EKAT_REQUIRE_MSG(is_allocated(),
      "Error! Cannot extract a field's view before allocation happens.\n");

  EKAT_REQUIRE_MSG (not m_header->is_broadcast(),
      "Error! Broadcast fields can only be accessed via get_strided_view.\n"
      " - field name: " + name() + "\n");

  EKAT_REQUIRE_MSG (not m_is_read_only || std::is_const<DstValueType>::value,
      "Error! Cannot get a view to non-const data if the field is read-only.\n");

  // Get src details
  const auto& alloc_prop = m_header->get_alloc_properties();
  const auto& field_layout = m_header->get_identifier().get_layout();

  EKAT_REQUIRE_MSG(DstRank==field_layout.rank(),
      "Error! You can only reshape to a view of the correct rank (equal to the FieldLayout's one).\n");

  // Check the reinterpret cast makes sense for the Dst value types (need integer sizes ratio)
  EKAT_REQUIRE_MSG(alloc_prop.template is_compatible<DstValueType>(),
      "Error! Source field allocation is not compatible with the requested value type.\n"
      " - field name: " + name() + "\n"
      " - field data type: " + e2str(data_type()) + "\n"
      " - alloc pack size: " + std::to_string(alloc_prop.get_largest_pack_size()) + "\n"
      " - DstValueType: " + std::string(typeid(DstValueType).name()) + "\n");

  // Start by reshaping into a ND view with all dynamic extents
  const auto view_ND = get_ND_view<HD,DstValueType,DstRank>();

  using dyn_DT = typename decltype(view_ND)::traits::data_type;
  if (!std::is_same<dyn_DT,DT>::value) {
    // The user requested some compile-time dimensions.
    // Let's check that they are correct
    for (int i=DstRankDynamic; i<DstRank; ++i) {
      EKAT_REQUIRE_MSG(view_ND.extent(i)==dims_type::static_extent(i),
          "Error! The template DataType contains an invalid compile-time dimensions:\n"
          "    - field name: " + m_header->get_identifier().name() + "\n"
          "    - dim index: " + std::to_string(i) + "\n"
          "    - input compile time dimension: " + std::to_string(dims_type::static_extent(i)) + "\n"
          "    - field internal dimension: " + std::to_string(view_ND.extent(i)) + "\n");
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

template <typename DT, HostOrDevice HD>
auto Field::get_strided_view() const ->
get_strided_view_type<DT, HD> {
  // The destination view type on correct mem space
  using DstView = get_strided_view_type<DT, HD>;
  // The dst value types
  using DstValueType = typename DstView::traits::value_type;
  // We only allow to reshape to a view of the correct rank
  constexpr int DstRank = DstView::rank;

  if constexpr (DstRank > 0) {
    // Get src details
    const auto& alloc_prop = m_header->get_alloc_properties();
    const auto& fl = m_header->get_identifier().get_layout();

    // Checks
    EKAT_REQUIRE_MSG(
        is_allocated(),
        "Error! Cannot extract a field's view before allocation happens.\n");
    EKAT_REQUIRE_MSG(not m_is_read_only || std::is_const<DstValueType>::value,
                    "Error! Cannot get a view to non-const data if the field is "
                    "read-only.\n");
    EKAT_REQUIRE_MSG(alloc_prop.template is_compatible<DstValueType>(),
                    "Error! Source field allocation is not compatible with the "
                    "requested value type.\n");

    // Check if this field is a broadcast of another field
    if (m_header->is_broadcast()) {
      // Get the strided view of the source field (whose rank is smaller than DstRank),
      // and build a rank-DstRank view, with the same pointer, and stride 0 along
      // the broadcasted dims. NOTE: we MUST build the src field, since only the Field
      // structure can build the original view (from which we grab the strides)
      const auto& src = m_header->template get_extra_data<Field>("bcast_src");
      const auto& src_layout = src.get_header().get_identifier().get_layout();

      DstValueType* ptr = nullptr;
      std::array<size_t,MaxRank> src_strides = {0};
      auto get_src_strides = [&](const auto& v) {
        ptr = v.data();
        for (int i=0; i<src_layout.rank(); ++i) {
          src_strides[i] = v.stride(i);
        }
      };
      switch (src_layout.rank()) {
        case 0: get_src_strides(src.get_strided_view<data_nd_t<DstValueType,0>,HD>()); break;
        case 1: get_src_strides(src.get_strided_view<data_nd_t<DstValueType,1>,HD>()); break;
        case 2: get_src_strides(src.get_strided_view<data_nd_t<DstValueType,2>,HD>()); break;
        case 3: get_src_strides(src.get_strided_view<data_nd_t<DstValueType,3>,HD>()); break;
        case 4: get_src_strides(src.get_strided_view<data_nd_t<DstValueType,4>,HD>()); break;
        case 5: get_src_strides(src.get_strided_view<data_nd_t<DstValueType,5>,HD>()); break;
        default:
          EKAT_ERROR_MSG ("Error! Unexpected source rank in broadcast field.\n"
                          " - field name: " + name() + "\n");
      }

      auto tgt_layout = m_header->get_identifier().get_layout();
      // We already check broadcast works when we created this field
      auto dim_map = src_layout.broadcast_map(tgt_layout).value();
      Kokkos::LayoutStride kl;
      for (int i=0; i<DstRank; ++i) {
        kl.dimension[i] = fl.dim(i);
        kl.stride[i]    = dim_map[i]>=0 ? src_strides[dim_map[i]] : 0;
      }
      return DstView(ptr,kl);
    }

    // Check if this field is a subview of another field
    const auto parent = m_header->get_parent();
    if (parent != nullptr) {
      // Parent field has correct layout to reinterpret the view into N+1-dim view,
      // for single-slice subfield, and N-dim view for multi-slice subfield.
      // So create the parent field on the fly, use it to get the N+{1,0}-dim view,
      // then subview it. NOTE: we can set protected members, since f is the same
      // type of this class.
      Field f;
      f.m_header = parent;
      f.m_data = m_data;

      // get subview info to determine whether we are single- or multi-slicing
      const auto& sv_alloc_prop = m_header->get_alloc_properties();
      const auto& info = sv_alloc_prop.get_subview_info();
      const int idim = info.dim_idx;
      const int k = info.slice_idx;
      const int k_end = info.slice_idx_end;

      // k_end has not been set by a multi-slice subfield function
      if (k_end == -1) {
        // Take an (n + 1)-dimensional == DstRank (== 2D, in practice) view
        // with normal LayoutRight
        auto v_np1 = f.get_ND_view<HD, DstValueType, DstRank + 1>();

        // As of now, we can only single-slice subview at first or second dimension.
        EKAT_REQUIRE_MSG(idim == 0 || idim == 1,
                        "Error! Subview dimension index is out of bounds.\n");

        // Use correct subview utility
        if (idim == 0) {
          return DstView(ekat::subview(v_np1, k));
        } else {
          return DstView(ekat::subview_1(v_np1, k));
        }
      // k_end has been set, so we're multi-slicing
      } else if (k_end > 0) {
        // rank doesn't change for multi-slice
        EKAT_REQUIRE_MSG(DstRank == fl.rank(),
                        "Error! Destination view rank must be equal to parent "
                        "field's rank for multi-sliced subview.\n");
        auto v_fullsize = f.get_ND_view<HD, DstValueType, DstRank>();

        return DstView(ekat::subview(
            v_fullsize, Kokkos::make_pair<int, int>(k, k_end), idim));
      }
    }
  }
  // Either not a subfield or requesting a zero-D view from a
  // 0D or 1D subfield, so stride == 1, and we can create the
  // strided view from the LayoutRight 1d view.
  return DstView(get_ND_view<HD, DstValueType, DstRank>());
}

template<HostOrDevice HD, typename T, int N>
auto Field::get_ND_view () const
  -> std::enable_if_t<(N < MaxRank), get_view_type<data_nd_t<T, N>, HD>>
{
  const auto& fl = m_header->get_identifier().get_layout();
  EKAT_REQUIRE_MSG (N==1 || N==fl.rank(),
      "Error! Input Rank must either be 1 (flat array) or the actual field rank.\n");

  // Check if this field is a subview of another field
  const auto parent = m_header->get_parent();
  if (parent!=nullptr) {
    // Parent field has correct layout to reinterpret the view into N+1-dim view
    // So create the parent field on the fly, use it to get the N+1-dim view, then subview it.
    // NOTE: we can set protected members, since f is the same type of this class.
    Field f;
    f.m_header = parent;
    f.m_data   = m_data;

    auto v_np1 = f.get_ND_view<HD,T,N+1>();

    // Now we can subview v_np1 at the correct slice
    const auto& info = m_header->get_alloc_properties().get_subview_info();
    const int idim = info.dim_idx;
    const int k    = info.slice_idx;

    // So far we can only subview at first or second dimension.
    EKAT_REQUIRE_MSG (idim==0 || idim==1,
        "Error! Subview dimension index is out of bounds.\n");

    EKAT_REQUIRE_MSG (idim==0 || N>1,
        "Error! Cannot subview a rank-2 (or less) view along 2nd dimension "
        "without losing LayoutRight.\n");

    // Use SFINAE-ed get_subview helper function to pick correct
    // subview impl. If N+1<=2 and idim!=0, the code craps out in the check above.
    if (idim==0) {
      return ekat::subview(v_np1,k);
    } else {
      return get_subview_1<HD,T,N+1>(v_np1,k);
    }
  }

  // Compute extents from FieldLayout
  const auto& alloc_prop = m_header->get_alloc_properties();
  auto num_values = alloc_prop.get_alloc_size() / sizeof(T);
  Kokkos::LayoutRight kl;
  for (int i=0; i<N; ++i) {
    if (i==N-1) {
      kl.dimension[i] = num_values;
    } else {
      kl.dimension[i] = fl.dim(i);
      num_values = fl.dim(i)==0 ? 0 : num_values/fl.dim(i);
    }
  }

  auto ptr = get_internal_view_data<T,HD>();
  using ret_type = get_view_type<data_nd_t<T,N>,HD>;
  return ret_type (ptr,kl);
}

template<HostOrDevice HD,typename T,int N>
auto Field::get_ND_view () const
  -> std::enable_if_t<N==MaxRank,get_view_type<data_nd_t<T,N>,HD>>
{
  static_assert(HD==Host or HD==Device,
      "Invalid value for non-type template argument HD.\n");

  const auto& fl = m_header->get_identifier().get_layout();
  EKAT_REQUIRE_MSG (N==1 || N==fl.rank(),
      "Error! Input Rank must either be 1 (flat array) or the actual field rank.\n");

  // Given that N==MaxRank, this field cannot be a subview of another field
  EKAT_REQUIRE_MSG (m_header->get_parent()==nullptr,
      "Error! A view of rank " + std::to_string(MaxRank) + " should not be the subview of another field.\n");

  // Compute extents from FieldLayout
  const auto& alloc_prop = m_header->get_alloc_properties();
  auto num_values = alloc_prop.get_alloc_size() / sizeof(T);
  Kokkos::LayoutRight kl;
  for (int i=0; i<N-1; ++i) {
    kl.dimension[i] = fl.dim(i);
    num_values /= fl.dim(i);
  }
  kl.dimension[N-1] = num_values;

  auto ptr = get_internal_view_data<T,HD>();
  using ret_type = get_view_type<data_nd_t<T,N>,HD>;
  return ret_type (ptr,kl);
}

// NOTE: DO NOT USE--this circumvents compile-time issues with
// subview slicing in get_strided_view()
template<HostOrDevice HD,typename T,int N>
auto Field::get_ND_view () const
  -> std::enable_if_t<(N >= MaxRank + 1),get_view_type<data_nd_t<T,N>,HD>>
{
  EKAT_ERROR_MSG("Error! Cannot call get_ND_view for rank greater than "
                 "MaxRank = 6.\n"
                 "This should never be called at run time.\n"
                 "Please contact developer if this functionality is required\n");
  return get_view_type<data_nd_t<T,N>,HD>();
}

} // namespace scream

#endif // SCREAM_FIELD_IMPL_HPP
