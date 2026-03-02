cmake_minimum_required(VERSION 3.16)

set(EMULATORATM_CXX_SOURCES
  src/atm.cpp
  src/atm_f2c.cpp
)

set(EMULATORATM_CXX_HEADERS
  src/atm.hpp
)

set(EMULATORATM_F90_SOURCES
  src/atm_comp_mct.F90
)

if(STANDALONE_MODE)
  # Build only C++ sources in standalone mode
  add_library(emulatoratm ${EMULATORATM_CXX_SOURCES})
else()
  # Build library with both C++ and Fortran sources
  add_library(emulatoratm ${EMULATORATM_CXX_SOURCES} ${EMULATORATM_F90_SOURCES})

  # Find and link CIME shared library
  find_package(CsmShare REQUIRED)
  target_link_libraries(emulatoratm PUBLIC csm_share)

  # Set up Fortran module directory
  set_target_properties(emulatoratm PROPERTIES 
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)
  target_include_directories(emulatoratm PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/modules)
endif()
# Include paths
target_include_directories(emulatoratm PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${CMAKE_CURRENT_SOURCE_DIR}/../common/src
)

# Link to common emulator framework
target_link_libraries(emulatoratm PUBLIC emulator_common)

# MPI is required for both standalone and full builds
find_package(MPI REQUIRED)
target_link_libraries(emulatoratm PUBLIC MPI::MPI_CXX)

# C++17 required
target_compile_features(emulatoratm PUBLIC cxx_std_17)
