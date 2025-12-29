cmake_minimum_required(VERSION 3.16)

#==============================================================================
# Emulator Components
#
# Generalized framework for emulated E3SM components (atmosphere, ocean, ice)
#==============================================================================

project(EmulatorComps LANGUAGES CXX C Fortran)

#------------------------------------------------------------------------------
# Build Options
#------------------------------------------------------------------------------
option(EMULATOR_STANDALONE_BUILD "Build emulator_comps standalone (without CIME)" OFF)
option(EMULATOR_BUILD_TESTS "Build emulator component tests" OFF)

# Inference backend: stub or libtorch
set(EATM_INFERENCE_BACKEND "stub" CACHE STRING "Inference backend for EATM (stub|libtorch)")
set_property(CACHE EATM_INFERENCE_BACKEND PROPERTY STRINGS "stub" "libtorch")

# Create C++ definitions based on inference backend
string(TOUPPER "${EATM_INFERENCE_BACKEND}" EATM_INFERENCE_BACKEND_UPPER)
add_compile_definitions(EATM_INFERENCE_BACKEND_${EATM_INFERENCE_BACKEND_UPPER})

# Enable backend dependencies automatically
if(EATM_INFERENCE_BACKEND STREQUAL "libtorch")
  set(EMULATOR_HAS_LIBTORCH ON CACHE BOOL "Enable LibTorch integration" FORCE)
endif()

#------------------------------------------------------------------------------
# Standalone Build Configuration
#------------------------------------------------------------------------------
if(EMULATOR_STANDALONE_BUILD)
  message(STATUS "Configuring standalone build...")
  include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/EmulatorStandalone.cmake)
endif()

#------------------------------------------------------------------------------
# Build Libraries
#------------------------------------------------------------------------------

# Build common framework first
add_subdirectory(common)

# Build component-specific emulators
add_subdirectory(eatm)

# Future components:
# add_subdirectory(eocn)
# add_subdirectory(eice)

#------------------------------------------------------------------------------
# Testing
#------------------------------------------------------------------------------
if(EMULATOR_BUILD_TESTS)
  enable_testing()
  
  # Include test utilities
  include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/EmulatorTestUtils.cmake)
  
  message(STATUS "")
  message(STATUS "=== Configuring Tests ===")
  
  # Add test directories for each component
  add_subdirectory(common/tests)
  add_subdirectory(eatm/tests)
  
  message(STATUS "=== Tests Configured ===")
  message(STATUS "")
endif()
