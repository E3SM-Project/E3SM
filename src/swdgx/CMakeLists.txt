# Relevant Directories within the HOMME repository
# holding source files needed for building the library
# and executable sweqx.

SET(SRC_BASE    ${CMAKE_CURRENT_SOURCE_DIR}/..)
SET(SRC_SHARE   ${SRC_BASE}/share)
SET(SRC_UTILS   ${SRC_BASE}/../utils/csm_share)

# Make INCLUDE_DIRS global so the tests can access it
SET (SWDGX_INCLUDE_DIRS ${PIO_INCLUDE_DIRS}
                        ${TIMING_INCLUDE_DIRS} CACHE INTERNAL "include directories for swdgx")

SET(SWDGX_SRCS
    ${SRC_BASE}/jrio.c
    ${SRC_SHARE}/kinds.F90
    ${SRC_SHARE}/dimensions_mod.F90
    ${SRC_SHARE}/control_mod.F90
    ${SRC_SHARE}/physical_constants.F90
    ${SRC_SHARE}/quadrature_mod.F90
    ${SRC_SHARE}/parallel_mod.F90
    ${SRC_SHARE}/hybrid_mod.F90
    ${SRC_SHARE}/reduction_mod.F90
    ${SRC_SHARE}/params_mod.F90
    ${SRC_SHARE}/gridgraph_mod.F90
    ${SRC_SHARE}/coordinate_systems_mod.F90
    ${SRC_SHARE}/edge_mod.F90
    ${SRC_SHARE}/element_mod.F90
    ${SRC_SHARE}/global_norms_mod.F90
    ${SRC_SHARE}/thread_mod.F90
    ${SRC_SHARE}/ll_mod.F90
    ${SRC_SHARE}/metagraph_mod.F90
    ${SRC_SHARE}/schedule_mod.F90
    ${SRC_SHARE}/schedtype_mod.F90
    ${SRC_SHARE}/bndry_mod.F90
    ${SRC_SHARE}/derivative_mod.F90
    ${SRC_SHARE}/viscosity_mod.F90
    ${SRC_SHARE}/time_mod.F90
    ${SRC_BASE}/common_io_mod.F90
    ${SRC_BASE}/pio_io_mod.F90
    ${SRC_BASE}/netcdf_io_mod.F90
    ${SRC_SHARE}/spacecurve_mod.F90
    ${SRC_SHARE}/cube_mod.F90
    ${SRC_SHARE}/interpolate_mod.F90
    ${SRC_SHARE}/spelt_mod.F90
    ${SRC_SHARE}/fvm_mod.F90
    ${SRC_SHARE}/fvm_control_volume_mod.F90
    ${SRC_SHARE}/fvm_line_integrals_mod.F90 
    ${SRC_SHARE}/fvm_line_integrals_flux_mod.F90 
    ${SRC_SHARE}/fvm_analytic_mod.F90 
    ${SRC_SHARE}/fvm_filter_mod.F90 
    ${SRC_SHARE}/fvm_reconstruction_mod.F90 
    ${SRC_BASE}/interp_movie_mod.F90
    ${SRC_UTILS}/shr_kind_mod.F90
    ${SRC_UTILS}/shr_mpi_mod.F90
    ${SRC_UTILS}/shr_sys_mod.F90
    ${SRC_UTILS}/shr_file_mod.F90
    ${SRC_UTILS}/shr_vmath_mod.F90
    ${SRC_UTILS}/shr_vmath_fwrap.c
    ${SRC_SHARE}/domain_mod.F90
    ${SRC_BASE}/restart_io_mod.F90
    ${SRC_SHARE}/cg_mod.F90
    ${SRC_SHARE}/linear_algebra_mod.F90
    ${SRC_SHARE}/solver_mod.F90
    ${SRC_SHARE}/filter_mod.F90
    ${SRC_BASE}/dg_flux_mod.F90
    ${SRC_BASE}/dg_core_mod.F90
    ${SRC_BASE}/dg_tests_mod.F90
    ${SRC_BASE}/ref_state_mod.F90
    ${SRC_BASE}/dg_tvdrk_mod.F90
    ${SRC_BASE}/state_mod.F90
    ${SRC_BASE}/dg_sweq_mod.F90
    ${SRC_BASE}/common_movie_mod.F90
    ${SRC_SHARE}/namelist_mod.F90
    ${SRC_SHARE}/mass_matrix_mod.F90
    ${SRC_SHARE}/dof_mod.F90
    ${SRC_SHARE}/perfmodel_mod.F90
    ${SRC_SHARE}/metis_mod.F90
    ${SRC_BASE}/checksum_mod.F90
    ${SRC_BASE}/repro_sum_mod.F90
    ${SRC_BASE}/restart_mod.F90
    ${SRC_BASE}/init_mod.F90
    ${SRC_BASE}/dg_movie_mod.F90
    ${SRC_SHARE}/mesh_mod.F90
)

SET(SWDGX_SRCS ${SWDGX_SRCS} ${SRC_BASE}/dg_main.F90 CACHE INTERNAL "Sources for swdgx")

# Set up defaults
IF (NOT SWDGX_NP)
  SET (SWDGX_NP 6)
ENDIF ()
IF (NOT SWDGX_PLEV)
  SET (SWDGX_PLEV 1)
ENDIF ()
IF (NOT SWDGX_USE_PIO)
  SET (SWDGX_USE_PIO FALSE) 
ENDIF ()
IF (NOT SWDGX_USE_ENERGY)
  SET (SWDGX_USE_ENERGY FALSE) 
ENDIF ()
IF (NOT SWDGX_NC)
  SET (SWDGX_NC 4)
ENDIF ()

MESSAGE(STATUS "Building swdgx with:")
MESSAGE(STATUS "  NP = ${SWDGX_NP}")
MESSAGE(STATUS "  NC = ${SWDGX_NC}")
MESSAGE(STATUS "  PLEV = ${SWDGX_PLEV}")
MESSAGE(STATUS "  PIO = ${SWDGX_USE_PIO}")
MESSAGE(STATUS "  ENERGY = ${SWDGX_USE_ENERGY}")

# The following macro handles the compilation with the appropriate variables
############################################################################
# createTestExec(exec_name exec_type NP NC PLEV USE_PIO USE_ENERGY)
############################################################################
createTestExec(swdgx swdgx ${SWDGX_NP} ${SWDGX_NC} ${SWDGX_PLEV}
               ${SWDGX_USE_PIO} ${SWDGX_USE_ENERGY})

