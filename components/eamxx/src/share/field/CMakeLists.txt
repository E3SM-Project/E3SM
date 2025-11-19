add_library(eamxx_field
  field_alloc_prop.cpp
  field_identifier.cpp
  field_header.cpp
  field_layout.cpp
  field_tracking.cpp
  field.cpp
  field_group.cpp
  field_sync.cpp
  eti/field_deep_copy_double.cpp
  eti/field_deep_copy_float.cpp
  eti/field_deep_copy_int.cpp
  eti/field_get_strided_view_double_device.cpp
  eti/field_get_strided_view_float_device.cpp
  eti/field_get_strided_view_int_device.cpp
  eti/field_get_view_double_device.cpp
  eti/field_get_view_float_device.cpp
  eti/field_get_view_int_device.cpp
  eti/field_update_double.cpp
  eti/field_update_float.cpp
  eti/field_update_impl_double.cpp
  eti/field_update_impl_fill_double.cpp
  eti/field_update_impl_fill_float.cpp
  eti/field_update_impl_fill_int.cpp
  eti/field_update_impl_float.cpp
  eti/field_update_impl_int.cpp
  eti/field_update_int.cpp
  eti/field_update_max_min_double.cpp
  eti/field_update_max_min_float.cpp
  eti/field_update_max_min_impl_double.cpp
  eti/field_update_max_min_impl_fill_double.cpp
  eti/field_update_max_min_impl_fill_float.cpp
  eti/field_update_max_min_impl_fill_int.cpp
  eti/field_update_max_min_impl_float.cpp
  eti/field_update_max_min_impl_int.cpp
  eti/field_update_max_min_int.cpp
  utils/compute_mask.cpp
  utils/horiz_contraction.cpp
  utils/perturb.cpp
  utils/print_hyperslab.cpp
  utils/randomize.cpp
  utils/reductions.cpp
  utils/transpose.cpp
  utils/vert_contraction.cpp
  utils/views_are_equal.cpp
)

if (EAMXX_ENABLE_GPU)
  # The HostOrDevice::Device is different from HostOrDevice::Host,
  # so we need to ETI also these fcns for HostOrDevice::Host
  # NOTE: we DON'T do ETI for the update/max/min utils, as they are
  #       only used on device (at least for now). OTOH, the methods
  #       to extract views are used in a few places
  target_sources (eamxx_field PRIVATE
    eti/field_get_strided_view_double_host.cpp
    eti/field_get_strided_view_float_host.cpp
    eti/field_get_strided_view_int_host.cpp
    eti/field_get_view_double_host.cpp
    eti/field_get_view_float_host.cpp
    eti/field_get_view_int_host.cpp
  )
endif()

target_link_libraries (eamxx_field PUBLIC
  eamxx_core
  eamxx_utils)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
