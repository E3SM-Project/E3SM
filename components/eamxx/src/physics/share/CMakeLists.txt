set(PHYSICS_SHARE_SRCS
  physics_share_f2c.F90
  physics_share.cpp
  physics_test_data.cpp
  scream_trcmix.cpp
  ${SCREAM_BASE_DIR}/../eam/src/physics/cam/physics_utils.F90
  ${SCREAM_BASE_DIR}/../eam/src/physics/cam/scream_abortutils.F90
  ${SCREAM_BASE_DIR}/../eam/src/physics/cam/wv_sat_scream.F90
  ${SCREAM_BASE_DIR}/../eam/src/physics/cam/micro_p3_utils.F90
  ${SCREAM_BASE_DIR}/../eam/src/physics/cam/debug_info.F90
)

# Add ETI source files if not on CUDA/HIP
if (NOT EAMXX_ENABLE_GPU)
  list(APPEND PHYSICS_SHARE_SRCS
  physics_saturation.cpp)
endif()

add_library(physics_share ${PHYSICS_SHARE_SRCS})
set_target_properties(physics_share PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules
)
target_include_directories(physics_share PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
)
target_link_libraries(physics_share scream_share)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
