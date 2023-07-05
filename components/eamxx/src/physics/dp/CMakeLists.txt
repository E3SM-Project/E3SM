set(DP_SRCS
  dp_f90.cpp
  dp_iso_c.f90
  #${SCREAM_BASE_DIR}/../eam/src/control/apply_iop_forcing.F90
  #${SCREAM_BASE_DIR}/../eam/src/dynamics/se/se_iop_intr_mod.F90",
  #${SCREAM_BASE_DIR}/../eam/src/control/iop_data_mod.F90",
  #${SCREAM_BASE_DIR}/../eam/src/control/history_iop.F90"
)

if (NOT SCREAM_LIB_ONLY)
  list(APPEND DP_SRCS
    dp_functions_f90.cpp
    ) # Add f90 bridges needed for testing
endif()

# Add ETI source files if not on CUDA/HIP
if (NOT EAMXX_ENABLE_GPU OR Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE OR Kokkos_ENABLE_HIP_RELOCATABLE_DEVICE_CODE)
  list(APPEND DP_SRCS
    eti/dp_advance_iop_forcing.cpp
    eti/dp_advance_iop_nudging.cpp
    eti/dp_advance_iop_subsidence.cpp
    eti/dp_iop_setinitial.cpp
    eti/dp_iop_broadcast.cpp
    eti/dp_apply_iop_forcing.cpp
    eti/dp_iop_domain_relaxation.cpp
    ) # DP ETI SRCS
endif()

add_library(dp ${DP_SRCS})
set_target_properties(dp PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules
)
target_include_directories(dp PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/../share
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/impl
  ${SCREAM_BASE_DIR}/../eam/src/physics/cam
)
target_link_libraries(dp physics_share scream_share)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
