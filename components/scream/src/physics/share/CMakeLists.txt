set(CIMEROOT ${SCREAM_BASE_DIR}/../../cime)
list(APPEND CMAKE_MODULE_PATH ${CIMEROOT}/src/CMake)

set(GENF90 ${CIMEROOT}/src/externals/genf90/genf90.pl)
set(ENABLE_GENF90 True)
include(genf90_utils)
include(Sourcelist_utils)

set(PHYSICS_SHARE_SRCS
  physics_share_f2c.F90
  physics_share.cpp
  physics_only_grids_manager.cpp
  physics_test_data.cpp
  ${SCREAM_BASE_DIR}/../cam/src/physics/cam/physics_utils.F90
  ${SCREAM_BASE_DIR}/../cam/src/physics/cam/scream_abortutils.F90
)

# Add ETI source files if not on CUDA
if (NOT CUDA_BUILD)
  list(APPEND PHYSICS_SHARE_SRCS
  physics_saturation.cpp)
endif()

add_library(physics_share ${PHYSICS_SHARE_SRCS})
target_include_directories(physics_share PUBLIC ${SCREAM_INCLUDE_DIRS})
target_include_directories(physics_share SYSTEM PUBLIC ${SCREAM_TPL_INCLUDE_DIRS} ${CIMEROOT}/src/share/include)
set_target_properties(physics_share PROPERTIES
  Fortran_MODULE_DIRECTORY ${SCREAM_F90_MODULES})
target_link_libraries(physics_share scream_share ${SCREAM_TPL_LIBRARIES})

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()

