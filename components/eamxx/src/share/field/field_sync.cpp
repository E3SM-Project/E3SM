#include "share/field/field.hpp"

namespace scream
{

void Field::
sync_to_host (const bool fence) const {
  // Sanity check
  EKAT_REQUIRE_MSG (is_allocated(),
      "Error! Input field must be allocated in order to sync host and device views.\n");

  // Check for early return if Host and Device are the same memory space
  if (host_and_device_share_memory_space()) return;

  // We allow sync_to_host for constant fields. Temporarily disable read only flag.
  const bool original_read_only = m_is_read_only;
  m_is_read_only = false;

  switch (data_type()) {
    case DataType::IntType:
      sync_views_impl<int, Device, Host>();
      break;
    case DataType::FloatType:
      sync_views_impl<float, Device, Host>();
      break;
    case DataType::DoubleType:
      sync_views_impl<double, Device, Host>();
      break;
    default:
      EKAT_ERROR_MSG("Error! Unrecognized field data type in Field::sync_to_host.\n");
  }

  if (fence) Kokkos::fence();

  // Return field to read-only state
  m_is_read_only = original_read_only;
}

void Field::
sync_to_dev (const bool fence) const {
  // Sanity check
  EKAT_REQUIRE_MSG (is_allocated(),
      "Error! Input field must be allocated in order to sync host and device views.\n");

  // Check for early return if Host and Device are the same memory space
  if (host_and_device_share_memory_space()) return;

  switch (data_type()) {
    case DataType::IntType:
      sync_views_impl<int, Host, Device>();
      break;
    case DataType::FloatType:
      sync_views_impl<float, Host, Device>();
      break;
    case DataType::DoubleType:
      sync_views_impl<double, Host, Device>();
      break;
    default:
      EKAT_ERROR_MSG("Error! Unrecognized field data type in Field::sync_to_dev.\n");
  }

  if (fence) Kokkos::fence();
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

} // namespace scream
