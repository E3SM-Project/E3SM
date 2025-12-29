cmake_minimum_required(VERSION 3.16)

#==============================================================================
# Emulated Atmosphere Component (EATM)
#==============================================================================

# Options for inference backends are now handled in emulator_common
# See components/emulator_comps/common/CMakeLists.txt

#------------------------------------------------------------------------------
# Core EATM sources (always built)
#------------------------------------------------------------------------------
set(EATM_CXX_SOURCES
  src/emulator_atm.cpp
  src/emulator_atm_interface.cpp
  src/impl/atm_field_manager.cpp
  src/impl/atm_coupling.cpp
  src/impl/atm_io.cpp
  src/impl/atm_field_data_provider.cpp
)

set(EATM_CXX_HEADERS
  src/emulator_atm.hpp
  src/impl/atm_field_manager.hpp
  src/impl/atm_coupling.hpp
  src/impl/atm_io.hpp
  src/impl/atm_field_data_provider.hpp
)

#------------------------------------------------------------------------------
# Fortran MCT wrapper sources (only for CIME builds, not standalone)
#------------------------------------------------------------------------------
set(EATM_F90_SOURCES
  src/atm_comp_mct.F90
  src/emulator_atm_f2c.F90
)

#------------------------------------------------------------------------------
# Create EATM library
#------------------------------------------------------------------------------

if(EMULATOR_STANDALONE_BUILD)
  # Standalone mode: C++ only, no Fortran MCT wrappers
  add_library(eatm ${EATM_CXX_SOURCES})
  message(STATUS "  EATM: Standalone mode (C++ only, no MCT wrappers)")
else()
  # CIME build: Include Fortran MCT wrappers
  add_library(eatm
    ${EATM_CXX_SOURCES}
    ${EATM_F90_SOURCES}
  )
  
  # Link to E3SM shared infrastructure (csm_share)
  # This provides ESMF (or esmf.mod stub), MCT, and derived types required by F90 wrapper
  find_package(CsmShare REQUIRED)
  target_link_libraries(eatm PUBLIC csm_share)
  
  # Set Fortran Module Directory
  set_target_properties(eatm PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)
  
  # Expose modules to consumers (cpl/driver)
  target_include_directories(eatm PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/modules)
  
  message(STATUS "  EATM: CIME build (with MCT wrappers)")
endif()

# Include paths
target_include_directories(eatm PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${CMAKE_CURRENT_SOURCE_DIR}/src/impl
  ${CMAKE_CURRENT_SOURCE_DIR}/../common/src
)

# Link to common emulator framework
target_link_libraries(eatm PUBLIC emulator_common)

# C++17 required
target_compile_features(eatm PUBLIC cxx_std_17)

# Create alias 'atm' so build_model can link against it (only for CIME builds)
if(NOT EMULATOR_STANDALONE_BUILD)
  if(NOT TARGET atm)
    add_library(atm ALIAS eatm)
  endif()
endif()

#------------------------------------------------------------------------------
# Status summary
#------------------------------------------------------------------------------
message(STATUS "EATM configuration summary:")
message(STATUS "  Core sources: ${EATM_CXX_SOURCES}")

