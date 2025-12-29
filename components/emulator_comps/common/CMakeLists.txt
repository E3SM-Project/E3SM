cmake_minimum_required(VERSION 3.16)

#==============================================================================
# Emulator Common Library
#
# Contains shared framework code for all emulator components (atm, ocn, ice).
#==============================================================================

set(EMULATOR_COMMON_SOURCES
  src/emulator_comp.cpp
  src/emulator_config.cpp
  src/coupling_fields.cpp
  src/emulator_io.cpp
  src/emulator_logger.cpp
  src/emulator_diagnostics.cpp
  src/emulator_output_stream.cpp
  src/emulator_output_manager.cpp
  # Diagnostics subfolder (modular derived diagnostics)
  src/diagnostics/horiz_avg_diagnostic.cpp
  src/diagnostics/vert_slice_diagnostic.cpp
  src/diagnostics/diagnostic_factory.cpp
  # Inference backends
  src/inference/inference_factory.cpp
  src/inference/stub_backend.cpp
)

set(EMULATOR_COMMON_HEADERS
  src/emulator_context.hpp
  src/emulator_comp.hpp
  src/emulator_config.hpp
  src/coupling_fields.hpp
  src/emulator_io.hpp
  src/emulator_logger.hpp
  src/inference/inference_backend.hpp
  src/inference/stub_backend.hpp
)

# Create shared library
add_library(emulator_common 
  ${EMULATOR_COMMON_SOURCES}
)

target_include_directories(emulator_common PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${CMAKE_CURRENT_SOURCE_DIR}/src/diagnostics
  ${CMAKE_CURRENT_SOURCE_DIR}/src/inference
)

# Link MPI
find_package(MPI REQUIRED)
target_link_libraries(emulator_common PUBLIC MPI::MPI_CXX MPI::MPI_Fortran)

# C++17 for std::any
target_compile_features(emulator_common PUBLIC cxx_std_17)

# Suppress warnings from yaml-cpp headers (std::iterator deprecated in C++17)
# This affects our code that includes yaml-cpp headers
target_compile_options(emulator_common PRIVATE -Wno-deprecated-declarations)

# PIO integration for I/O - required for emulator_io.hpp
# For CIME builds:
#   - Source headers (pio.h etc) are in externals/scorpio/src/clib
#   - Generated headers (pio_config.h) are in INSTALL_SHAREDPATH/include
# For Standalone builds:
#   - Source headers are in externals/scorpio/src/clib
#   - Generated headers (pio_config.h) are in ${CMAKE_BINARY_DIR}/externals/scorpio/src/clib
set(SCORPIO_CLIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../../externals/scorpio/src/clib)

if(EMULATOR_STANDALONE_BUILD)
  # Standalone build: use both source and build directories
  # pio.h is in source dir, pio_config.h is generated in build dir
  set(SCORPIO_BUILD_CLIB_DIR ${CMAKE_BINARY_DIR}/externals/scorpio/src/clib)
  message(STATUS "  PIO: Using scorpio (standalone build)")
  message(STATUS "    Source: ${SCORPIO_CLIB_DIR}")
  message(STATUS "    Build:  ${SCORPIO_BUILD_CLIB_DIR}")
  target_include_directories(emulator_common PUBLIC 
    ${SCORPIO_CLIB_DIR}
    ${SCORPIO_BUILD_CLIB_DIR}
  )
  target_compile_definitions(emulator_common PUBLIC EMULATOR_HAS_PIO)
elseif(DEFINED INSTALL_SHAREDPATH AND EXISTS ${SCORPIO_CLIB_DIR}/pio.h)
  # CIME build: need both source and install paths
  message(STATUS "  PIO: Using scorpio (CIME build)")
  message(STATUS "    Source: ${SCORPIO_CLIB_DIR}")
  message(STATUS "    Config: ${INSTALL_SHAREDPATH}/include")
  target_include_directories(emulator_common PUBLIC 
    ${SCORPIO_CLIB_DIR}
    ${INSTALL_SHAREDPATH}/include
  )
  target_compile_definitions(emulator_common PUBLIC EMULATOR_HAS_PIO)
elseif(EXISTS ${SCORPIO_CLIB_DIR}/pio.h)
  # Fallback - just source path (may fail if pio_config.h needed)
  message(STATUS "  PIO: Using scorpio from: ${SCORPIO_CLIB_DIR}")
  target_include_directories(emulator_common PUBLIC ${SCORPIO_CLIB_DIR})
  target_compile_definitions(emulator_common PUBLIC EMULATOR_HAS_PIO)
else()
  message(WARNING "  PIO headers not found - I/O functionality will be limited")
endif()

# Link against csm_share if available (provides PIO/MCT/ESMF and proper link order)
find_package(CsmShare QUIET)
if(TARGET csm_share)
  message(STATUS "  Linking against csm_share (provides PIO/MCT/ESMF)")
  target_link_libraries(emulator_common PUBLIC csm_share)
elseif(TARGET pioc)
  # Standalone build - link to pioc target from SCORPIO
  message(STATUS "  Linking against pioc target (standalone)")
  target_link_libraries(emulator_common PUBLIC pioc)
elseif(TARGET spio)
  # Fallback to spio
  message(STATUS "  Linking against spio target")
  target_link_libraries(emulator_common PUBLIC spio)
endif()

# Optional LibTorch integration
option(EMULATOR_HAS_LIBTORCH "Enable LibTorch integration" OFF)
if(EMULATOR_HAS_LIBTORCH)
  find_package(Torch REQUIRED)
  
  # Add LibTorch backend sources
  target_sources(emulator_common PRIVATE
    src/inference/libtorch_backend.cpp
  )
  
  target_link_libraries(emulator_common PUBLIC ${TORCH_LIBRARIES})
  target_include_directories(emulator_common PUBLIC ${TORCH_INCLUDE_DIRS})
  target_compile_definitions(emulator_common PUBLIC EMULATOR_HAS_LIBTORCH EMULATOR_HAS_MPI)
  message(STATUS "  LibTorch found: ${TORCH_LIBRARIES}")
  message(STATUS "  LibTorch includes: ${TORCH_INCLUDE_DIRS}")
endif()

#------------------------------------------------------------------------------
# yaml-cpp (REQUIRED for YAML configuration)
# build_emulator_comps.cmake may have already found/built yaml-cpp
# If not, we build it from EKAT submodule
#------------------------------------------------------------------------------
if(EMULATOR_YAML_CPP_TARGET)
  # yaml-cpp was already found by build_emulator_comps.cmake
  # Use find_package to get the include dirs, then link by library name
  find_package(yaml-cpp QUIET)
  if(yaml-cpp_FOUND AND TARGET yaml-cpp)
    message(STATUS "  yaml-cpp: Using system yaml-cpp target")
    target_link_libraries(emulator_common PUBLIC yaml-cpp)
  elseif(yaml-cpp_FOUND)
    # Old-style find, use variables directly
    message(STATUS "  yaml-cpp: Using system yaml-cpp (legacy)")
    target_include_directories(emulator_common PUBLIC ${YAML_CPP_INCLUDE_DIR})
    target_link_libraries(emulator_common PUBLIC ${YAML_CPP_LIBRARIES})
  else()
    # Just try the library name
    message(STATUS "  yaml-cpp: Linking against yaml-cpp")
    target_link_libraries(emulator_common PUBLIC yaml-cpp)
  endif()
  target_compile_definitions(emulator_common PUBLIC EMULATOR_HAS_YAML_CPP)
else()
  # Build yaml-cpp from EKAT submodule (standalone build or fallback)
  set(YAML_CPP_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../../externals/ekat/extern/yaml-cpp)
  
  if(EXISTS ${YAML_CPP_SOURCE_DIR}/CMakeLists.txt)
    message(STATUS "  yaml-cpp: Building from EKAT submodule (fallback)")
    message(STATUS "    Source: ${YAML_CPP_SOURCE_DIR}")
    
    # Collect yaml-cpp source files directly
    file(GLOB YAML_CPP_SOURCES ${YAML_CPP_SOURCE_DIR}/src/*.cpp)
    
    # Add yaml-cpp sources to emulator_common
    target_sources(emulator_common PRIVATE ${YAML_CPP_SOURCES})
    
    # Add yaml-cpp include directories
    target_include_directories(emulator_common PUBLIC ${YAML_CPP_SOURCE_DIR}/include)
    
    target_compile_definitions(emulator_common PUBLIC EMULATOR_HAS_YAML_CPP)
  else()
    message(FATAL_ERROR "yaml-cpp not found! Checked:\n"
            "  - EMULATOR_YAML_CPP_TARGET (from build_emulator_comps)\n"
            "  - EKAT submodule at ${YAML_CPP_SOURCE_DIR}")
  endif()
endif()
