# Add RRTMGP source files. Right now, this builds the Fortran source for RRTMGP.
# Once the C++ port of RRTMGP is complete we should switch over to that.
set(RRTMGP_SRCS
  rrtmgp.cpp
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
)

add_library(rrtmgp ${RRTMGP_SRCS})

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
