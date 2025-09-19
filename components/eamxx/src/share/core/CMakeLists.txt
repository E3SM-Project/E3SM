# Generate eamxx_config.h and eamxx_config.f
include (EkatUtils)
EkatConfigFile(${CMAKE_CURRENT_SOURCE_DIR}/eamxx_config.h.in
               ${CMAKE_CURRENT_BINARY_DIR}/eamxx_config.h
               F90_FILE ${CMAKE_CURRENT_BINARY_DIR}/eamxx_config.f)

# Generate eamxx_version.h
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/eamxx_version.h.in
               ${CMAKE_CURRENT_BINARY_DIR}/eamxx_version.h)

add_library(eamxx_core
  eamxx_config.cpp
  eamxx_session.cpp
)

target_include_directories(eamxx_core PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR} # for config/version headers
  ${SCREAM_SRC_DIR}           # for all eamxx headers
)

# This used to be in `eamxx_config.h/eamxx_config.f`, but we did accidentally remove
# the includes at least twice. To avoid this from happening again, we add this crucial
# macro as a compile definition of the eamxx_core target, which is linked against by
# all other eamxx libraries.
if (SCREAM_CIME_BUILD)
  target_compile_definitions (eamxx_core PUBLIC SCREAM_CIME_BUILD)
endif()

# Attach custom f90 flags, so they are picked up by any EAMxx target
target_compile_options(eamxx_core PUBLIC
  $<$<COMPILE_LANGUAGE:Fortran>:${SCREAM_Fortran_FLAGS}>
)

# When dealing with subtle bugs in read/write of unit tests baselines,
# using ascii baselines can help inspect that the file indeed looks as expected
option (EAMXX_ASCII_BASELINES "Whether physics units tests should assume ascii (instead of binary) format for baselines" OFF)
if (EAMXX_ASCII_BASELINES)
  target_compile_definitions (eamxx_core PUBLIC EAMXX_ASCII_BASELINES)
endif()

# We have some issues with RDC and repeated libraries in the link line.
# It's a known issue, and nvcc_wrapper has a flag for handling this.
if (Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE)
  if(CMAKE_VERSION VERSION_LESS "3.13.0")
    target_link_libraries(eamxx_core PUBLIC "--remove-duplicate-link-files")
  else()
    target_link_options(eamxx_core PUBLIC "--remove-duplicate-link-files")
  endif()
endif()

# Only link what this library need. Other eamxx libs will link other tpls
target_link_libraries(eamxx_core PUBLIC ekat::Core ekat::KokkosUtils)

# Handle the case where we have python support
if (EAMXX_ENABLE_PYTHON)
  if (NOT Python_EXECUTABLE)
    message (FATAL_ERROR "You must set Python_EXECUTABLE to point to a valid python3 interpreter")
  endif()

  execute_process(
    COMMAND "${Python_EXECUTABLE}" -c "import sys; assert sys.version_info >= (3,9), 'Python version must be >= 3.9'"
    RESULT_VARIABLE python_version_check_result
    ERROR_VARIABLE python_version_check_error
  )
  if (python_version_check_result)
      message(FATAL_ERROR "Python version check failed: ${python_version_check_error}")
  endif()

  # Need pybind11
  execute_process(
    COMMAND "${Python_EXECUTABLE}" -c "import pybind11; print(pybind11.get_cmake_dir())"
    OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE pybind11_ROOT)

  find_package(pybind11 REQUIRED HINTS ${pybind11_ROOT})

  target_link_libraries(eamxx_core PUBLIC pybind11::embed)
  target_compile_definitions (eamxx_core PUBLIC EAMXX_HAS_PYTHON)

  # Used in eamxx_pysession.hpp to get current path
  target_link_libraries(eamxx_core PUBLIC stdc++fs)
endif()

if (NOT SCREAM_LIB_ONLY)
  # Nano-library, providing ekat_initialize/finalize_test_session implementations
  # This lib is linked to ALL tests created by CreateUnitTest from cmake/ScreamUtils.cmake
  add_library(eamxx_test_session
    eamxx_test_session.cpp
  )
  target_include_directories(eamxx_test_session PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
  )
  target_link_libraries(eamxx_test_session PUBLIC eamxx_core)

  add_subdirectory(tests)
endif()
