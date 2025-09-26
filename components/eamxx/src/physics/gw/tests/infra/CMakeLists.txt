set(PATH_TO_LEGACY_GW ${SCREAM_BASE_DIR}/../eam/src/physics/cam/gw)
set(INFRA_SRCS
  gw_test_data.cpp
  gw_iso_c.f90
  ${PATH_TO_LEGACY_GW}/gw_utils.F90
  ${PATH_TO_LEGACY_GW}/gw_common.F90
  ${PATH_TO_LEGACY_GW}/gw_convect.F90
  ${PATH_TO_LEGACY_GW}/gw_diffusion.F90
  ${PATH_TO_LEGACY_GW}/gw_front.F90
  ${PATH_TO_LEGACY_GW}/gw_oro.F90
  ${PATH_TO_LEGACY_GW}/../vdiff_lu_solver.F90
)

set_target_properties(gw PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules
)

add_library(gw_test_infra ${INFRA_SRCS})
target_link_libraries(gw_test_infra PUBLIC gw scream_test_support)
target_include_directories(gw_test_infra PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  ${PATH_TO_LEGACY_GW}
  ${PATH_TO_LEGACY_GW}/..
)
