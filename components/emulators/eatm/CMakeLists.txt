cmake_minimum_required(VERSION 3.16)

set(EATM_CXX_SOURCES
  src/atm.cpp
  src/atm_f2c.cpp
)

set(EATM_CXX_HEADERS
  src/atm.hpp
)

set(EATM_F90_SOURCES
  src/atm_comp_mct.F90
  src/atm_f2c.F90
)

if(STANDALONE_MODE)
  # Build only C++ sources in standalone mode
  add_library(eatm ${EATM_CXX_SOURCES})
else()
  # Build library with both C++ and Fortran sources
  add_library(eatm ${EATM_CXX_SOURCES} ${EATM_F90_SOURCES})

  # Find and link CIME shared library
  find_package(CsmShare REQUIRED)
  target_link_libraries(eatm PUBLIC csm_share)

  # Set up Fortran module directory
  set_target_properties(eatm PROPERTIES 
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)
  target_include_directories(eatm PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/modules)
endif()
# Include paths
target_include_directories(eatm PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${CMAKE_CURRENT_SOURCE_DIR}/../common/src
)

# Link to common emulator framework
target_link_libraries(eatm PUBLIC emulator_common)

# C++17 required
target_compile_features(eatm PUBLIC cxx_std_17)

#------------------------------------------------------------------------------
# Status summary
#------------------------------------------------------------------------------
message(STATUS "EATM configuration summary:")
message(STATUS "  Core sources: ${EATM_CXX_SOURCES}")
