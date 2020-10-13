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
  ${SCREAM_BASE_DIR}/../eam/src/physics/cam/physics_utils.F90
  ${SCREAM_BASE_DIR}/../eam/src/physics/cam/scream_abortutils.F90
  ${SCREAM_BASE_DIR}/../eam/src/physics/cam/wv_sat_scream.F90
  ${SCREAM_BASE_DIR}/../eam/src/physics/cam/micro_p3_utils.F90
  ${SCREAM_BASE_DIR}/../eam/src/physics/cam/debug_info.F90
)

# Add ETI source files if not on CUDA
if (NOT CUDA_BUILD)
  list(APPEND PHYSICS_SHARE_SRCS
  physics_saturation.cpp)
endif()

add_library(physics_share ${PHYSICS_SHARE_SRCS})
target_compile_options(physics_share PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${SCREAM_Fortran_FLAGS}>)
set_target_properties(physics_share PROPERTIES Fortran_MODULE_DIRECTORY ${SCREAM_F90_MODULES})
target_link_libraries(physics_share scream_share)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
