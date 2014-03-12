# Relevant Directories within the HOMME repository
# holding source files needed for building the library
# and executable sweqx.

SET(SRC_BASE    ${CMAKE_CURRENT_SOURCE_DIR}/..)
SET(SRC_SHARE   ${SRC_BASE}/share)
SET(SRC_UTILS   ${SRC_BASE}/../utils/csm_share)
SET(PHYSICS_DIR ${HOMME_SOURCE_DIR}/physics)
SET(UTILS_TIMING_DIR ${HOMME_SOURCE_DIR}/utils/timing)

# Make INCLUDE_DIRS global so the tests can access it
SET (PRIMDGX_INCLUDE_DIRS ${PIO_INCLUDE_DIRS}
                        ${TIMING_INCLUDE_DIRS} CACHE INTERNAL "include directories for primdgx")

SET(PRIMDGX_SRCS
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
    ${SRC_SHARE}/spacecurve_mod.F90
    ${SRC_SHARE}/cube_mod.F90
    ${SRC_SHARE}/interpolate_mod.F90
    ${SRC_SHARE}/spelt_mod.F90
    ${SRC_SHARE}/fvm_control_volume_mod.F90
    ${SRC_SHARE}/fvm_reconstruction_mod.F90
    ${SRC_SHARE}/fvm_filter_mod.F90
    ${SRC_SHARE}/fvm_mod.F90
    ${SRC_SHARE}/fvm_line_integrals_mod.F90 
    ${SRC_SHARE}/fvm_line_integrals_flux_mod.F90 
    ${SRC_SHARE}/fvm_analytic_mod.F90 
    ${PHYSICS_DIR}/physics_types_mod.F90
    ${SRC_SHARE}/dof_mod.F90
    ${SRC_BASE}/common_io_mod.F90
    ${SRC_BASE}/pio_io_mod.F90
    ${PHYSICS_DIR}/physics_io_mod.F90
    ${SRC_BASE}/aquaplanet.F90
    ${SRC_BASE}/aquaplanet_io_mod.F90
    ${SRC_SHARE}/hybvcoord_mod.F90
    ${SRC_SHARE}/time_mod.F90
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
    ${SRC_BASE}/common_movie_mod.F90
    ${SRC_SHARE}/namelist_mod.F90
    ${SRC_SHARE}/mass_matrix_mod.F90
    ${SRC_SHARE}/perfmodel_mod.F90
    ${SRC_SHARE}/metis_mod.F90
    ${UTILS_TIMING_DIR}/perf_mod.F90
    ${UTILS_TIMING_DIR}/perf_utils.F90
    ${SRC_BASE}/checksum_mod.F90
    ${SRC_BASE}/dg3d_vertical_mod.F90
    ${SRC_BASE}/ref_state_mod.F90
    ${SRC_BASE}/dg3d_errors_mod.F90
    ${SRC_BASE}/dg_core_mod.F90
    ${SRC_BASE}/dg3d_tests_mod.F90
    ${SRC_BASE}/dg3d_remap_mod.F90
    ${SRC_BASE}/dg3d_core_mod.F90
    ${SRC_BASE}/dg3d_dynamics_mod.F90
    ${SRC_BASE}/dg3d_prim_mod.F90
    ${SRC_BASE}/state_mod.F90
    ${SRC_BASE}/repro_sum_mod.F90
    ${SRC_BASE}/dg3d_primeq_mod.F90
    ${SRC_BASE}/restart_mod.F90
    ${SRC_BASE}/init_mod.F90
    ${SRC_BASE}/netcdf_io_mod.F90
    ${SRC_SHARE}/mesh_mod.F90 
)

SET(PRIMDGX_SRCS ${PRIMDGX_SRCS} ${SRC_BASE}/dg3d_main.F90 CACHE INTERNAL "Sources for primdgx")

# Set up defaults
IF (NOT PRIMDGX_NP)
  SET (PRIMDGX_NP 6)
ENDIF ()
IF (NOT PRIMDGX_PLEV)
  SET (PRIMDGX_PLEV 26)
ENDIF ()
IF (NOT PRIMDGX_USE_PIO)
  SET (PRIMDGX_USE_PIO FALSE) 
ENDIF ()
IF (NOT PRIMDGX_USE_ENERGY)
  SET (PRIMDGX_USE_ENERGY FALSE) 
ENDIF ()
IF (NOT PRIMDGX_NC)
  SET (PRIMDGX_NC 4)
ENDIF ()

MESSAGE(STATUS "Building primdgx with:")
MESSAGE(STATUS "  NP = ${PRIMDGX_NP}")
MESSAGE(STATUS "  NC = ${PRIMDGX_NC}")
MESSAGE(STATUS "  PLEV = ${PRIMDGX_PLEV}")
MESSAGE(STATUS "  PIO = ${PRIMDGX_USE_PIO}")
MESSAGE(STATUS "  ENERGY = ${PRIMDGX_USE_ENERGY}")

# The following macro handles the compilation with the appropriate variables
############################################################################
# createTestExec(exec_name exec_type NP NC PLEV USE_PIO USE_ENERGY)
############################################################################
createTestExec(primdgx primdgx ${PRIMDGX_NP} ${PRIMDGX_NC} ${PRIMDGX_PLEV}
               ${PRIMDGX_USE_PIO} ${PRIMDGX_USE_ENERGY})
