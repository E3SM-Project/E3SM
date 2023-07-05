# Add external packages

# Add the YAKL library
add_subdirectory(${E3SM_EXTERNALS_ROOT}/YAKL ${CMAKE_CURRENT_BINARY_DIR}/YAKL)
if (CMAKE_VERSION VERSION_GREATER "3.18.0")
  set_property(TARGET yakl PROPERTY CUDA_ARCHITECTURES OFF)
endif()
