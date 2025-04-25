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

} // namespace scream
