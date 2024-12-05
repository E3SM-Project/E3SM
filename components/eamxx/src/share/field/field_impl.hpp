#ifndef SCREAM_FIELD_IMPL_HPP
#define SCREAM_FIELD_IMPL_HPP

#include "share/field/field.hpp"
#include "share/util/scream_array_utils.hpp"
#include "share/util/scream_universal_constants.hpp"

#include <ekat/ekat_type_traits.hpp>

namespace scream
{

template<typename ViewT, typename>
Field::
Field (const identifier_type& id,
       const ViewT& view_d)
 : Field(id)
{
  constexpr int N = ViewT::rank;
  using ScalarT  = typename ViewT::traits::value_type;
  using ExeSpace = typename ViewT::traits::execution_space;

  EKAT_REQUIRE_MSG ( (std::is_same<ExeSpace,typename device_t::execution_space>::value),
      "Error! This constructor of Field requires a view from device.\n");

  EKAT_REQUIRE_MSG (id.data_type()==get_data_type<ScalarT>(),
      "Error! Input view data type does not match what is stored in the field identifier.\n"
      " - field name: " + id.name() + "\n"
      " - field data type: " + e2str(id.data_type()) + "\n");

  const auto& fl = id.get_layout();
  EKAT_REQUIRE_MSG (N==fl.rank(),
      "Error! This constructor of Field requires a device view of the correct rank.\n"
      " - field name: " + id.name() + "\n"
      " - field rank: " + std::to_string(fl.rank()) + "\n"
      " - view rank : " + std::to_string(N) + "\n");
  for (int i=0; i<(N-1); ++i) {
    EKAT_REQUIRE_MSG (view_d.extent_int(i)==fl.dims()[i],
        "Error! Input view has the wrong i-th extent.\n"
        " - field name: " + id.name() + "\n"
        " - idim: " + std::to_string(i) + "\n"
        " - layout i-th dim: " + std::to_string(fl.dims()[i]) + "\n"
        " - view i-th dim: " + std::to_string(view_d.extent(i)) + "\n");
  }

  auto& alloc_prop = m_header->get_alloc_properties();
  if (N>0 and view_d.extent_int(N-1)!=fl.dims().back()) {
    EKAT_REQUIRE_MSG (view_d.extent_int(N-1)>=fl.dims()[N-1],
        "Error! Input view has the wrong last extent.\n"
        " - field name: " + id.name() + "\n"
        " - layout last dim: " + std::to_string(fl.dims()[N-1]) + "\n"
        " - view last dim: " + std::to_string(view_d.extent(N-1)) + "\n");

    // We have a padded view. We don't know what the pack size was, so we pick the largest
    // power of 2 that divides the last extent
    auto last_view_dim = view_d.extent_int(N-1);
    int last_fl_dim    = fl.dims().back();

    // This should get the smallest pow of 2 that gives npacks*pack_size==view_last_dim
    int ps = 1;
    int packed_length = 0;
    do {
      ps *= 2;
      auto npacks = (last_fl_dim + ps - 1) / ps;
      packed_length = ps*npacks;
    }
    while (packed_length!=last_view_dim);

    alloc_prop.request_allocation(ps);
  }
  alloc_prop.commit(fl);

  // Create an unmanaged dev view, and its host mirror
  const auto view_dim = alloc_prop.get_alloc_size();
  char* data = reinterpret_cast<char*>(view_d.data());
  m_data.d_view = decltype(m_data.d_view)(data,view_dim);
  m_data.h_view = Kokkos::create_mirror_view(m_data.d_view);

  // Since we created m_data.d_view from a raw pointer, we don't get any
  // ref counting from the kokkos view. Hence, to ensure that the input view
  // survives as long as this Field, we store it as extra data in the header
  m_header->set_extra_data("orig_view",view_d);
}

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

  EKAT_REQUIRE_MSG (not m_is_read_only || std::is_const<DstValueType>::value,
      "Error! Cannot get a view to non-const data if the field is read-only.\n");

  // Get src details
  const auto& alloc_prop = m_header->get_alloc_properties();
  const auto& field_layout = m_header->get_identifier().get_layout();

  EKAT_REQUIRE_MSG(DstRank==field_layout.rank(),
      "Error! You can only reshape to a view of the correct rank (equal to the FieldLayout's one).\n");

  // Check the reinterpret cast makes sense for the Dst value types (need integer sizes ratio)
  EKAT_REQUIRE_MSG(alloc_prop.template is_compatible<DstValueType>(),
      "Error! Source field allocation is not compatible with the requested value type.\n");

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

    // Check if this field is a subview of another field
    const auto parent = m_header->get_parent().lock();
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

template<typename ST, HostOrDevice From, HostOrDevice To>
void Field::sync_views_impl () const {
  // For all Kokkos::deep_copy() calls we will pass in an instance of the
  // device execution space so that we are asynchronous w.r.t. host.
  using DeviceExecSpace = typename Field::get_device<Device>::execution_space;

  // Rank 0 will always be contiguous. Copy and return early.
  if (rank() == 0) {
    Kokkos::deep_copy(DeviceExecSpace(), get_view<ST, To>(), get_view<const ST, From>());
    return;
  }

  const bool is_contiguous = get_header().get_alloc_properties().contiguous();
  if (is_contiguous) {
    // For contiguous fields, simply use Kokkos::deep_copy().
    switch (rank()) {
      case 1:
        Kokkos::deep_copy(DeviceExecSpace(), get_view<ST*, To>(), get_view<const ST*, From>());
        break;
      case 2:
        Kokkos::deep_copy(DeviceExecSpace(), get_view<ST**, To>(), get_view<const ST**, From>());
        break;
      case 3:
        Kokkos::deep_copy(DeviceExecSpace(), get_view<ST***, To>(), get_view<const ST***, From>());
        break;
      case 4:
        Kokkos::deep_copy(DeviceExecSpace(), get_view<ST****, To>(), get_view<const ST****, From>());
        break;
      case 5:
        Kokkos::deep_copy(DeviceExecSpace(), get_view<ST*****, To>(), get_view<const ST*****, From>());
        break;
      case 6:
        Kokkos::deep_copy(DeviceExecSpace(), get_view<ST******, To>(), get_view<const ST******, From>());
        break;
      default:
        EKAT_ERROR_MSG ("Error! Unsupported field rank in Field::sync_to_host.\n");
    }
  } else {
    auto sync_helper = [this] () {
      if constexpr (To==Host) m_contiguous_field->sync_to_host();
      else                    m_contiguous_field->sync_to_dev();
    };
    switch (rank()) {
      case 1:
        Kokkos::deep_copy(DeviceExecSpace(),
                          m_contiguous_field->get_view<ST*, From>(),
                          get_strided_view<const ST*, From>());
        sync_helper();
        Kokkos::deep_copy(DeviceExecSpace(),
                          get_strided_view<ST*, To>(),
                          m_contiguous_field->get_view<const ST*, To>());
        break;
      case 2:
        Kokkos::deep_copy(DeviceExecSpace(),
                          m_contiguous_field->get_view<ST**, From>(),
                          get_strided_view<const ST**, From>());
        sync_helper();
        Kokkos::deep_copy(DeviceExecSpace(),
                          get_strided_view<ST**, To>(),
                          m_contiguous_field->get_view<const ST**, To>());
        break;
      case 3:
        Kokkos::deep_copy(DeviceExecSpace(),
                          m_contiguous_field->get_view<ST***, From>(),
                          get_strided_view<const ST***, From>());
        sync_helper();
        Kokkos::deep_copy(DeviceExecSpace(),
                          get_strided_view<ST***, To>(),
                          m_contiguous_field->get_view<const ST***, To>());
        break;
      case 4:
        Kokkos::deep_copy(DeviceExecSpace(),
                          m_contiguous_field->get_view<ST****, From>(),
                          get_strided_view<const ST****, From>());
        sync_helper();
        Kokkos::deep_copy(DeviceExecSpace(),
                          get_strided_view<ST****, To>(),
                          m_contiguous_field->get_view<const ST****, To>());
        break;
      case 5:
        Kokkos::deep_copy(DeviceExecSpace(),
                          m_contiguous_field->get_view<ST*****, From>(),
                          get_strided_view<const ST*****, From>());
        sync_helper();
        Kokkos::deep_copy(DeviceExecSpace(),
                          get_strided_view<ST*****, To>(),
                          m_contiguous_field->get_view<const ST*****, To>());
        break;
      case 6:
        Kokkos::deep_copy(DeviceExecSpace(),
                          m_contiguous_field->get_view<ST******, From>(),
                          get_strided_view<const ST******, From>());
        sync_helper();
        Kokkos::deep_copy(DeviceExecSpace(),
                          get_strided_view<ST******, To>(),
                          m_contiguous_field->get_view<const ST******, To>());
        break;
      default:
        EKAT_ERROR_MSG ("Error! Unsupported field rank in Field::sync_to_host.\n");
    }
  }
}

template<HostOrDevice HD>
void Field::
deep_copy (const Field& src) const {
  EKAT_REQUIRE_MSG (not m_is_read_only,
      "Error! Cannot call deep_copy on read-only fields.\n");

  EKAT_REQUIRE_MSG (data_type()==src.data_type(),
      "Error! Cannot copy fields with different data type.\n");

  switch (data_type()) {
    case DataType::IntType:
      deep_copy_impl<HD,int>(src);
      break;
    case DataType::FloatType:
      deep_copy_impl<HD,float>(src);
      break;
    case DataType::DoubleType:
      deep_copy_impl<HD,double>(src);
      break;
    default:
      EKAT_ERROR_MSG ("Error! Unrecognized field data type in Field::deep_copy.\n");
  }
}

template<typename ST, HostOrDevice HD>
void Field::
deep_copy (const ST value) const {
  EKAT_REQUIRE_MSG (not m_is_read_only,
      "Error! Cannot call deep_copy on read-only fields.\n");

  const auto my_data_type = data_type();
  switch (my_data_type) {
    case DataType::IntType:
      EKAT_REQUIRE_MSG( (std::is_convertible<ST,int>::value),
          "Error! Input value type is not convertible to field data type.\n"
          "   - Input value type: " + ekat::ScalarTraits<ST>::name() + "\n"
          "   - Field data type : " + e2str(my_data_type) + "\n");
      deep_copy_impl<HD,int>(value);
      break;
    case DataType::FloatType:
      EKAT_REQUIRE_MSG( (std::is_convertible<ST,float>::value),
          "Error! Input value type is not convertible to field data type.\n"
          "   - Input value type: " + ekat::ScalarTraits<ST>::name() + "\n"
          "   - Field data type : " + e2str(my_data_type) + "\n");
      deep_copy_impl<HD,float>(value);
      break;
    case DataType::DoubleType:
      EKAT_REQUIRE_MSG( (std::is_convertible<ST,double>::value),
          "Error! Input value type is not convertible to field data type.\n"
          "   - Input value type: " + ekat::ScalarTraits<ST>::name() + "\n"
          "   - Field data type : " + e2str(my_data_type) + "\n");
      deep_copy_impl<HD,double>(value);
      break;
    default:
      EKAT_ERROR_MSG ("Error! Unrecognized field data type in Field::deep_copy.\n");
  }
}

template<HostOrDevice HD, typename ST>
void Field::
deep_copy_impl (const Field& src) const {

  const auto& layout     =     get_header().get_identifier().get_layout();
  const auto& layout_src = src.get_header().get_identifier().get_layout();
  EKAT_REQUIRE_MSG(layout==layout_src,
       "ERROR: Unable to copy field " + src.get_header().get_identifier().name() +
          " to field " + get_header().get_identifier().name() + ".  Layouts don't match.\n");
  const auto  rank = layout.rank();

  // For rank 0 view, we only need to copy a single value and return
  if (rank == 0) {
    auto v     =     get_view<      ST,HD>();
    auto v_src = src.get_view<const ST,HD>();
    Kokkos::deep_copy(v,v_src);
    return;
  }

  // Note: we can't just do a deep copy on get_view_impl<HD>(), since this
  //       field might be a subfield of another. We need the reshaped view.
  //       Also, don't call Kokkos::deep_copy if this field and src have
  //       different pack sizes.
  auto src_alloc_props = src.get_header().get_alloc_properties();
  auto tgt_alloc_props =     get_header().get_alloc_properties();

  using device_t = typename Field::get_device<HD>;
  using exec_space = typename device_t::execution_space;
  using RangePolicy = Kokkos::RangePolicy<exec_space>;

  auto policy = RangePolicy(0,layout.size());

  using extents_type = typename ekat::KokkosTypes<device_t>::template view_1d<int>;
  extents_type ext;
  if constexpr (HD==Device) {
    ext = layout.extents();
  } else {
    ext = layout.extents_h();
  }
  switch (rank) {
    case 1:
      {
        if (src_alloc_props.contiguous() and tgt_alloc_props.contiguous()) {
          auto v     =     get_view<      ST*,HD>();
          auto v_src = src.get_view<const ST*,HD>();
          Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const int idx) {
            v(idx) = v_src(idx);
          });
        } else {
          auto v     =     get_strided_view<      ST*,HD>();
          auto v_src = src.get_strided_view<const ST*,HD>();
          Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const int idx) {
            v(idx) = v_src(idx);
          });
        }
      }
      break;
    case 2:
      {
        if (src_alloc_props.contiguous() and tgt_alloc_props.contiguous()) {
          auto v     =     get_view<      ST**,HD>();
          auto v_src = src.get_view<const ST**,HD>();
          Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const int idx) {
            int i,j;
            unflatten_idx(idx,ext,i,j);
            v(i,j) = v_src(i,j);
          });
        }
        else {
          auto v     =     get_strided_view<      ST**,HD>();
          auto v_src = src.get_strided_view<const ST**,HD>();
          Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const int idx) {
            int i,j;
            unflatten_idx(idx,ext,i,j);
            v(i,j) = v_src(i,j);
          });
        }
      }
      break;
    case 3:
      {
        if (src_alloc_props.contiguous() and tgt_alloc_props.contiguous()) {
          auto v     =     get_view<      ST***,HD>();
          auto v_src = src.get_view<const ST***,HD>();
          Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const int idx) {
            int i,j,k;
            unflatten_idx(idx,ext,i,j,k);
            v(i,j,k) = v_src(i,j,k);
          });
        } else {
          auto v     =     get_strided_view<      ST***,HD>();
          auto v_src = src.get_strided_view<const ST***,HD>();
          Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const int idx) {
            int i,j,k;
            unflatten_idx(idx,ext,i,j,k);
            v(i,j,k) = v_src(i,j,k);
          });
        }
      }
      break;
    case 4:
      {
        if (src_alloc_props.contiguous() and tgt_alloc_props.contiguous()) {
          auto v     =     get_view<      ST****,HD>();
          auto v_src = src.get_view<const ST****,HD>();
          Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const int idx) {
            int i,j,k,l;
            unflatten_idx(idx,ext,i,j,k,l);
            v(i,j,k,l) = v_src(i,j,k,l);
          });
        } else {
          auto v     =     get_strided_view<      ST****,HD>();
          auto v_src = src.get_strided_view<const ST****,HD>();
          Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const int idx) {
            int i,j,k,l;
            unflatten_idx(idx,ext,i,j,k,l);
            v(i,j,k,l) = v_src(i,j,k,l);
          });
        }
      }
      break;
    case 5:
      {
        if (src_alloc_props.contiguous() and tgt_alloc_props.contiguous()) {
          auto v     =     get_view<      ST*****,HD>();
          auto v_src = src.get_view<const ST*****,HD>();
          Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const int idx) {
            int i,j,k,l,m;
            unflatten_idx(idx,ext,i,j,k,l,m);
            v(i,j,k,l,m) = v_src(i,j,k,l,m);
          });
        } else {
          auto v     =     get_view<      ST*****,HD>();
          auto v_src = src.get_view<const ST*****,HD>();
          Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const int idx) {
            int i,j,k,l,m;
            unflatten_idx(idx,ext,i,j,k,l,m);
            v(i,j,k,l,m) = v_src(i,j,k,l,m);
          });
        }
      }
      break;
    default:
      EKAT_ERROR_MSG ("Error! Unsupported field rank in 'deep_copy'.\n");
  }
}

template<HostOrDevice HD, typename ST>
void Field::deep_copy_impl (const ST value) const {

  // Note: we can't just do a deep copy on get_view_impl<HD>(), since this
  //       field might be a subfield of another. Instead, get the
  //       reshaped view first, based on the field rank.

  const auto& layout = get_header().get_identifier().get_layout();
  const auto  rank   = layout.rank();
  switch (rank) {
    case 0:
      {
        auto v = get_view<ST,HD>();
        Kokkos::deep_copy(v,value);
      }
      break;
    case 1:
      {
        if (m_header->get_alloc_properties().contiguous()) {
          auto v = get_view<ST*,HD>();
          Kokkos::deep_copy(v,value);
        } else {
          auto v = get_strided_view<ST*,HD>();
          Kokkos::deep_copy(v,value);
        }
      }
      break;
    case 2:
      {
        if (m_header->get_alloc_properties().contiguous()) {
          auto v = get_view<ST**,HD>();
          Kokkos::deep_copy(v,value);
        } else {
          auto v = get_strided_view<ST**,HD>();
          Kokkos::deep_copy(v,value);
        }
      }
      break;
    case 3:
      {
        if (m_header->get_alloc_properties().contiguous()) {
          auto v = get_view<ST***,HD>();
          Kokkos::deep_copy(v,value);
        } else {
          auto v = get_strided_view<ST***,HD>();
          Kokkos::deep_copy(v,value);
        }
      }
      break;
    case 4:
      {
        if (m_header->get_alloc_properties().contiguous()) {
          auto v = get_view<ST****,HD>();
          Kokkos::deep_copy(v,value);
        } else {
          auto v = get_strided_view<ST****,HD>();
          Kokkos::deep_copy(v,value);
        }
      }
      break;
    case 5:
      {
        if (m_header->get_alloc_properties().contiguous()) {
          auto v = get_view<ST*****,HD>();
          Kokkos::deep_copy(v,value);
        } else {
          auto v = get_strided_view<ST*****,HD>();
          Kokkos::deep_copy(v,value);
        }
      }
      break;
    case 6:
      {
        if (m_header->get_alloc_properties().contiguous()) {
          auto v = get_view<ST******,HD>();
          Kokkos::deep_copy(v,value);
        } else {
          auto v = get_strided_view<ST******,HD>();
          Kokkos::deep_copy(v,value);
        }
      }
      break;
    default:
      EKAT_ERROR_MSG ("Error! Unsupported field rank in 'deep_copy'.\n");
  }
}

template<HostOrDevice HD, typename ST>
void Field::
update (const Field& x, const ST alpha, const ST beta)
{
  const auto& dt = data_type();

  // Determine if there is a FillValue that requires extra treatment.
  ST fill_val = constants::DefaultFillValue<ST>().value;

  if (x.get_header().has_extra_data("mask_value")) {

    if (dt==DataType::IntType) {
      fill_val = x.get_header().get_extra_data<int>("mask_value");
    } else if (dt==DataType::FloatType) {
      fill_val = x.get_header().get_extra_data<float>("mask_value");
    } else if (dt==DataType::DoubleType) {
      fill_val = x.get_header().get_extra_data<double>("mask_value");
    } else {
      EKAT_ERROR_MSG ("Error! Unrecognized/unsupported field data type in Field::update.\n");
    }
  }

  // If user passes, say, double alpha/beta for an int field, we should error out, warning about
  // a potential narrowing rounding. The other way around, otoh, is allowed (even though
  // there's an upper limit to the int values that a double can store, it is unlikely the user
  // will use such large factors).
  const auto dt_st = get_data_type<ST>();
  EKAT_REQUIRE_MSG (not is_narrowing_conversion(dt_st,dt),
      "Error! Coefficients alpha/beta may be narrowed when converted to x/y data type.\n"
      " - x/y data type  : " + e2str(dt) + "\n"
      " - coeff data type: " + e2str(dt_st) + "\n");

  if (dt==DataType::IntType) {
    return update_impl<CombineMode::ScaleUpdate,HD,int>(x,alpha,beta,fill_val);
  } else if (dt==DataType::FloatType) {
    return update_impl<CombineMode::ScaleUpdate,HD,float>(x,alpha,beta,fill_val);
  } else if (dt==DataType::DoubleType) {
    return update_impl<CombineMode::ScaleUpdate,HD,double>(x,alpha,beta,fill_val);
  } else {
    EKAT_ERROR_MSG ("Error! Unrecognized/unsupported field data type in Field::update.\n");
  }
}

template<HostOrDevice HD, typename ST>
void Field::
scale (const ST beta)
{
  const auto& dt = data_type();

  // Determine if there is a FillValue that requires extra treatment.
  ST fill_val = constants::DefaultFillValue<ST>().value;
  if (get_header().has_extra_data("mask_value")) {
    fill_val = get_header().get_extra_data<ST>("mask_value");
  }

  // If user passes, say, double beta for an int field, we should error out, warning about
  // a potential narrowing rounding. The other way around, otoh, is allowed (even though
  // there's an upper limit to the int values that a double can store, it is unlikely the user
  // will use such large factors).
  const auto dt_st = get_data_type<ST>();
  EKAT_REQUIRE_MSG (not is_narrowing_conversion(dt_st,dt),
      "Error! Coefficients alpha/beta may be narrowed when converted to x/y data type.\n"
      " - x/y data type  : " + e2str(dt) + "\n"
      " - coeff data type: " + e2str(dt_st) + "\n");

  if (dt==DataType::IntType) {
    return update_impl<CombineMode::Rescale,HD,int>(*this,ST(0),beta,fill_val);
  } else if (dt==DataType::FloatType) {
    return update_impl<CombineMode::Rescale,HD,float>(*this,ST(0),beta,fill_val);
  } else if (dt==DataType::DoubleType) {
    return update_impl<CombineMode::Rescale,HD,double>(*this,ST(0),beta,fill_val);
  } else {
    EKAT_ERROR_MSG ("Error! Unrecognized/unsupported field data type in Field::scale.\n");
  }
}

template<HostOrDevice HD>
void Field::
scale (const Field& x)
{
  const auto& dt = data_type();
  if (dt==DataType::IntType) {
    int fill_val = constants::DefaultFillValue<int>().value;
    if (get_header().has_extra_data("mask_value")) {
      fill_val = get_header().get_extra_data<int>("mask_value");
    }
    return update_impl<CombineMode::Multiply,HD,int>(x,0,0,fill_val);
  } else if (dt==DataType::FloatType) {
    float fill_val = constants::DefaultFillValue<float>().value;
    if (get_header().has_extra_data("mask_value")) {
      fill_val = get_header().get_extra_data<float>("mask_value");
    }
    return update_impl<CombineMode::Multiply,HD,float>(x,0,0,fill_val);
  } else if (dt==DataType::DoubleType) {
    double fill_val = constants::DefaultFillValue<double>().value;
    if (get_header().has_extra_data("mask_value")) {
      fill_val = get_header().get_extra_data<double>("mask_value");
    }
    return update_impl<CombineMode::Multiply,HD,double>(x,0,0,fill_val);
  } else {
    EKAT_ERROR_MSG ("Error! Unrecognized/unsupported field data type in Field::scale.\n");
  }
}



template<CombineMode CM, HostOrDevice HD,typename ST>
void Field::
update_impl (const Field& x, const ST alpha, const ST beta, const ST fill_val)
{
  // Check x/y are allocated
  EKAT_REQUIRE_MSG (is_allocated(),
      "Error! Cannot update field, since it is not allocated.\n"
      " - field name: " + name() + "\n");
  EKAT_REQUIRE_MSG (x.is_allocated(),
      "Error! Cannot update field, since source field is not allocated.\n"
      " - field name: " + x.name() + "\n");

  // Check y is writable
  EKAT_REQUIRE_MSG (not is_read_only(),
      "Error! Cannot update field, as it is read-only.\n"
      " - field name: " + name() + "\n");

  // Check compatibility between fields data type
  const auto dt_x = x.data_type();
  const auto dt_y = data_type();

  // We could be more general, and do
  //   EKAT_REQUIRE_MSG (not is_narrowing_conversion(dt_x,dt_y),"");
  // but then we'd need N*N template instantiations of the impl function
  EKAT_REQUIRE_MSG (dt_x==dt_y,
      "Error! X and Y field must have the same data type.\n"
      " - x name: " + x.name() + "\n"
      " - y name: " + name() + "\n"
      " - x data type: " + e2str(dt_x) + "\n"
      " - y data type: " + e2str(dt_y) + "\n");

  const auto& y_l = get_header().get_identifier().get_layout();
  const auto& x_l = x.get_header().get_identifier().get_layout();
  EKAT_REQUIRE_MSG (y_l==x_l,
      "Error! Incompatible layouts for update_field.\n"
      " - x name: " + x.name() + "\n"
      " - y name: " + name() + "\n"
      " - x layout: " + x_l.to_string() + "\n"
      " - y layout: " + y_l.to_string() + "\n");

  using device_t = typename Field::get_device<HD>;
  using exec_space = typename device_t::execution_space;
  using RangePolicy = Kokkos::RangePolicy<exec_space>;
  auto policy = RangePolicy(0,x_l.size());

  using extents_type = typename ekat::KokkosTypes<device_t>::template view_1d<int>;
  extents_type ext;
  if constexpr (HD==Device) {
    ext = x_l.extents();
  } else {
    ext = x_l.extents_h();
  }

  switch (x_l.rank()) {
    case 0:
      {
        auto xv = x.get_view<const ST,HD>();
        auto yv =   get_view<      ST,HD>();
        Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const int /*idx*/) {
          combine_and_fill<CM>(xv(),yv(),fill_val,alpha,beta);
        });
      }
      break;
    case 1:
      {
        // Must handle the case where one of the two views is strided
        if (x.get_header().get_alloc_properties().contiguous() and
              get_header().get_alloc_properties().contiguous()) {
          auto xv = x.get_view<const ST*,HD>();
          auto yv =   get_view<      ST*,HD>();
          Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const int idx) {
            combine_and_fill<CM>(xv(idx),yv(idx),fill_val,alpha,beta);
          });
        } else {
          auto xv = x.get_strided_view<const ST*,HD>();
          auto yv =   get_strided_view<      ST*,HD>();
          Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const int idx) {
            combine_and_fill<CM>(xv(idx),yv(idx),fill_val,alpha,beta);
          });
        }
      }
      break;
    case 2:
      {
        auto xv = x.get_view<const ST**,HD>();
        auto yv =   get_view<      ST**,HD>();
        Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const int idx) {
          int i,j;
          unflatten_idx(idx,ext,i,j);
          combine_and_fill<CM>(xv(i,j),yv(i,j),fill_val,alpha,beta);
        });
      }
      break;
    case 3:
      {
        auto xv = x.get_view<const ST***,HD>();
        auto yv =   get_view<      ST***,HD>();
        Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const int idx) {
          int i,j,k;
          unflatten_idx(idx,ext,i,j,k);
          combine_and_fill<CM>(xv(i,j,k),yv(i,j,k),fill_val,alpha,beta);
        });
      }
      break;
    case 4:
      {
        auto xv = x.get_view<const ST****,HD>();
        auto yv =   get_view<      ST****,HD>();
        Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const int idx) {
          int i,j,k,l;
          unflatten_idx(idx,ext,i,j,k,l);
          combine_and_fill<CM>(xv(i,j,k,l),yv(i,j,k,l),fill_val,alpha,beta);
        });
      }
      break;
    case 5:
      {
        auto xv = x.get_view<const ST*****,HD>();
        auto yv =   get_view<      ST*****,HD>();
        Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const int idx) {
          int i,j,k,l,m;
          unflatten_idx(idx,ext,i,j,k,l,m);
          combine_and_fill<CM>(xv(i,j,k,l,m),yv(i,j,k,l,m),fill_val,alpha,beta);
        });
      }
      break;
    case 6:
      {
        auto xv = x.get_view<const ST******,HD>();
        auto yv =   get_view<      ST******,HD>();
        Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const int idx) {
          int i,j,k,l,m,n;
          unflatten_idx(idx,ext,i,j,k,l,m,n);
          combine_and_fill<CM>(xv(i,j,k,l,m,n),yv(i,j,k,l,m,n),fill_val,alpha,beta);
        });
      }
      break;
    default:
      EKAT_ERROR_MSG ("Error! Rank not supported in update_field.\n"
          " - x name: " + x.name() + "\n"
          " - y name: " + name() + "\n");
  }
  Kokkos::fence();
}

template<HostOrDevice HD, typename T, int N>
auto Field::get_ND_view () const
  -> if_t<(N < MaxRank), get_view_type<data_nd_t<T, N>, HD>>
{
  const auto& fl = m_header->get_identifier().get_layout();
  EKAT_REQUIRE_MSG (N==1 || N==fl.rank(),
      "Error! Input Rank must either be 1 (flat array) or the actual field rank.\n");

  // Check if this field is a subview of another field
  const auto parent = m_header->get_parent().lock();
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
  auto ptr = reinterpret_cast<T*>(get_view_impl<HD>().data());

  using ret_type = get_view_type<data_nd_t<T,N>,HD>;

  return ret_type (ptr,kl);
}

template<HostOrDevice HD,typename T,int N>
auto Field::get_ND_view () const
  -> if_t<N==MaxRank,get_view_type<data_nd_t<T,N>,HD>>
{
  const auto& fl = m_header->get_identifier().get_layout();
  EKAT_REQUIRE_MSG (N==1 || N==fl.rank(),
      "Error! Input Rank must either be 1 (flat array) or the actual field rank.\n");

  // Given that N==MaxRank, this field cannot be a subview of another field
  EKAT_REQUIRE_MSG (m_header->get_parent().expired(),
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
  auto ptr = reinterpret_cast<T*>(get_view_impl<HD>().data());

  using ret_type = get_view_type<data_nd_t<T,N>,HD>;
  return ret_type (ptr,kl);
}

// NOTE: DO NOT USE--this circumvents compile-time issues with
// subview slicing in get_strided_view()
template<HostOrDevice HD,typename T,int N>
auto Field::get_ND_view () const
  -> if_t<(N >= MaxRank + 1),get_view_type<data_nd_t<T,N>,HD>>
{
  EKAT_ERROR_MSG("Error! Cannot call get_ND_view for rank greater than "
                 "MaxRank = 6.\n"
                 "This should never be called at run time.\n"
                 "Please contact developer if this functionality is required\n");
}

} // namespace scream

#endif // SCREAM_FIELD_IMPL_HPP
