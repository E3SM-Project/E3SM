set(SCREAM_CONTROL_SOURCES
  atmosphere_driver.cpp
  surface_coupling.cpp
  surface_coupling_export.cpp
  fvphyshack.cpp
)

set(SCREAM_CONTROL_HEADERS
  atmosphere_driver.hpp
  surface_coupling.hpp
)

add_library(scream_control ${SCREAM_CONTROL_SOURCES})
target_link_libraries(scream_control PUBLIC scream_share scream_io)

if (Kokkos_ENABLE_CUDA)
  # This is to silence some nvcc warning that is a CUDA compiler bug
  # See https://github.com/kokkos/kokkos/issues/1473 for more details
  target_compile_options (scream_control PUBLIC
      $<$<COMPILE_LANGUAGE:CXX>:-Xcudafe --diag_suppress=esa_on_defaulted_function_ignored>)
endif()

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
