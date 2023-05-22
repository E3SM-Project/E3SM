set (dynLibName)
if ("${SCREAM_DYNAMICS_DYCORE}" STREQUAL "HOMME")
  # Need these two to be linked, due to homme deps

  # Recall that src/dynamics/homme/CMakeLists.txt does not build a dyn lib.
  # It only defines a macro to build it, so that different tests can build
  # a dyn lib if necessary.
  # Here, we ask to create a dyn lib depending on the config options that
  # were set during e3sm configuration

  # This is a list of cmake vars whose values are then used when calling
  # CreateDynamicsLib, to correctly build the dynamics library within Homme.
  # We set them to a default, but each compset should set its values anyways.
  set (SCREAM_NP 4 CACHE STRING "The number of Gauss points per element.")
  set (SCREAM_NUM_TRACERS 10 CACHE STRING "The max number of tracers.")
  set (SCREAM_USE_PIO FALSE CACHE STRING "Whether Homme can use PIO.")

  CreateDynamicsLib ("theta-l_kokkos" ${SCREAM_NP} ${SCREAM_NUM_VERTICAL_LEV} ${SCREAM_NUM_TRACERS})
endif()

set(SCREAM_CONTROL_SOURCES
  atmosphere_driver.cpp
  fvphyshack.cpp
  atmosphere_surface_coupling_importer.cpp
  atmosphere_surface_coupling_exporter.cpp
  surface_coupling_utils.cpp
)

set(SCREAM_CONTROL_HEADERS
  atmosphere_driver.hpp
  atmosphere_surface_coupling.hpp
  surface_coupling_utils.hpp
)

add_library(scream_control ${SCREAM_CONTROL_SOURCES})
target_link_libraries(scream_control PUBLIC scream_share scream_io ${dynLibName})

if (Kokkos_ENABLE_CUDA)
  # This is to silence some nvcc warning that is a CUDA compiler bug
  # See https://github.com/kokkos/kokkos/issues/1473 for more details
  target_compile_options (scream_control PUBLIC
      $<$<COMPILE_LANGUAGE:CXX>:-Xcudafe --diag_suppress=esa_on_defaulted_function_ignored>)
endif()

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
