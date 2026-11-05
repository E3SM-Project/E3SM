add_library(eamxx_field
  field_alloc_prop.cpp
  field_identifier.cpp
  field_header.cpp
  field_layout.cpp
  field_tracking.cpp
  field.cpp
  field_group.cpp
  field_sync.cpp
  field_reader.cpp
  # Unity build files combine multiple ETI instantiations into single translation units.
  # This reduces header parsing overhead (headers parsed once per unity file instead of
  # once per individual ETI file) while keeping template instantiation distributed.
  # See eti/field_*_device_unity.cpp files which each #include their respective ETI sources.
  eti/field_get_strided_view_device_unity.cpp
  eti/field_get_view_device_unity.cpp
  eti/field_deep_copy_impl_device_unity.cpp
  eti/field_deep_copy_masked_device_unity.cpp
  eti/field_update_impl_device_unity.cpp
  eti/field_update_fill_aware_device_unity.cpp
  eti/field_update_masked_device_unity.cpp
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
    eti/field_get_strided_view_host_unity.cpp
    eti/field_get_view_host_unity.cpp
  )
endif()

target_link_libraries (eamxx_field PUBLIC
  eamxx_core
  eamxx_utils
  eamxx_scorpio_interface)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
