# Add RRTMGP source files. Right now, this builds the Fortran source for RRTMGP.
# Once the C++ port of RRTMGP is complete we should switch over to that.
set(RRTMGP_SRCS
  atmosphere_radiation.cpp
  rrtmgp.cpp
  scream_rrtmgp_interface.cpp
  scream_rrtmgp_interface.F90
  ${SCREAM_BASE_DIR}/../cam/src/physics/rrtmgp/external/rrtmgp/kernels/mo_gas_optics_kernels.F90
  ${SCREAM_BASE_DIR}/../cam/src/physics/rrtmgp/external/rrtmgp/kernels/mo_rrtmgp_util_reorder_kernels.F90
  ${SCREAM_BASE_DIR}/../cam/src/physics/rrtmgp/external/rrtmgp/mo_gas_concentrations.F90
  ${SCREAM_BASE_DIR}/../cam/src/physics/rrtmgp/external/rrtmgp/mo_gas_optics.F90
  ${SCREAM_BASE_DIR}/../cam/src/physics/rrtmgp/external/rrtmgp/mo_gas_optics_rrtmgp.F90
  ${SCREAM_BASE_DIR}/../cam/src/physics/rrtmgp/external/rrtmgp/mo_rrtmgp_constants.F90
  ${SCREAM_BASE_DIR}/../cam/src/physics/rrtmgp/external/rrtmgp/mo_rrtmgp_util_reorder.F90
  ${SCREAM_BASE_DIR}/../cam/src/physics/rrtmgp/external/rrtmgp/mo_rrtmgp_util_string.F90
  ${SCREAM_BASE_DIR}/../cam/src/physics/rrtmgp/external/rte/mo_fluxes.F90
  ${SCREAM_BASE_DIR}/../cam/src/physics/rrtmgp/external/rte/mo_optical_props.F90
  ${SCREAM_BASE_DIR}/../cam/src/physics/rrtmgp/external/rte/mo_rte_kind.F90
  ${SCREAM_BASE_DIR}/../cam/src/physics/rrtmgp/external/rte/mo_rte_lw.F90
  ${SCREAM_BASE_DIR}/../cam/src/physics/rrtmgp/external/rte/mo_rte_sw.F90
  ${SCREAM_BASE_DIR}/../cam/src/physics/rrtmgp/external/rte/mo_rte_util_array.F90
  ${SCREAM_BASE_DIR}/../cam/src/physics/rrtmgp/external/rte/mo_source_functions.F90
  ${SCREAM_BASE_DIR}/../cam/src/physics/rrtmgp/external/rte/kernels/mo_fluxes_broadband_kernels.F90
  ${SCREAM_BASE_DIR}/../cam/src/physics/rrtmgp/external/rte/kernels/mo_optical_props_kernels.F90
  ${SCREAM_BASE_DIR}/../cam/src/physics/rrtmgp/external/rte/kernels/mo_rte_solver_kernels.F90
  ${SCREAM_BASE_DIR}/../cam/src/physics/rrtmgp/external/extensions/mo_fluxes_byband.F90
  ${SCREAM_BASE_DIR}/../cam/src/physics/rrtmgp/external/extensions/mo_fluxes_byband_kernels.F90
  ${SCREAM_BASE_DIR}/../cam/src/physics/rrtmgp/external/extensions/mo_fluxes_bygpoint.F90
  ${SCREAM_BASE_DIR}/../cam/src/physics/rrtmgp/external/extensions/mo_heating_rates.F90
  ${SCREAM_BASE_DIR}/../cam/src/physics/rrtmgp/external/extensions/mo_rrtmgp_clr_all_sky.F90
)

set(RRTMGP_HEADERS
  rrtmgp.hpp
  scream_rrtmgp_interface.hpp
  atmosphere_radiation.hpp
)

add_library(rrtmgp ${RRTMGP_SRCS})
target_include_directories(rrtmgp PUBLIC ${SCREAM_INCLUDE_DIRS})
target_include_directories(rrtmgp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../share)
target_include_directories(rrtmgp SYSTEM PUBLIC ${SCREAM_TPL_INCLUDE_DIRS} ${CIMEROOT}/src/share/include)
set_target_properties(rrtmgp PROPERTIES
      Fortran_MODULE_DIRECTORY ${SCREAM_F90_MODULES})
target_link_libraries(rrtmgp physics_share scream_share ${SCREAM_TPL_LIBRARIES})

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()

