# Relevant Directories within the HOMME repository
# holding source files needed for building the library
# and executable sweqx.

SET(SRC_BASE    ${CMAKE_CURRENT_SOURCE_DIR}/..)
SET(SRC_SHARE   ${SRC_BASE}/share)
SET(SRC_UTILS   ${SRC_BASE}/../utils/csm_share)


# Make INCLUDE_DIRS global so the tests can access it
SET (SWEQX_INCLUDE_DIRS ${PIO_INCLUDE_DIRS}
                        ${TIMING_INCLUDE_DIRS} CACHE INTERNAL "include directories for sweqx")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}
                    ${SWEQX_INCLUDE_DIRS})

CONFIGURE_FILE(config.h.cmake.in config.h)
ADD_DEFINITIONS(-DHAVE_CONFIG_H)

SET(SWEQX_SRCS
               ${SRC_BASE}/advance_mod.F90
               ${SRC_BASE}/checksum_mod.F90
               ${SRC_BASE}/common_io_mod.F90
               ${SRC_BASE}/common_movie_mod.F90
               ${SRC_BASE}/fvm_bsp_mod.F90
               ${SRC_BASE}/fvm_transformation_mod.F90
               ${SRC_BASE}/init_mod.F90
               ${SRC_BASE}/interp_movie_mod.F90
               ${SRC_BASE}/netcdf_io_mod.F90
               ${SRC_BASE}/pio_io_mod.F90
               ${SRC_BASE}/repro_sum_mod.F90
               ${SRC_BASE}/restart_io_mod.F90
               ${SRC_BASE}/restart_mod.F90
               ${SRC_BASE}/rk_mod.F90
               ${SRC_BASE}/shallow_water_mod.F90
               ${SRC_BASE}/shal_movie_mod.F90
               ${SRC_BASE}/state_mod.F90
               ${SRC_BASE}/sweq_mod.F90
               ${SRC_BASE}/types_mod.F90
               ${SRC_SHARE}/bndry_mod.F90
               ${SRC_SHARE}/cg_mod.F90
               ${SRC_SHARE}/control_mod.F90
               ${SRC_SHARE}/coordinate_systems_mod.F90
               ${SRC_SHARE}/cube_mod.F90
               ${SRC_SHARE}/derivative_mod.F90
               ${SRC_SHARE}/dimensions_mod.F90
               ${SRC_SHARE}/dof_mod.F90
               ${SRC_SHARE}/domain_mod.F90
               ${SRC_SHARE}/edge_mod.F90
               ${SRC_SHARE}/element_mod.F90
               ${SRC_SHARE}/filter_mod.F90
               ${SRC_SHARE}/fvm_analytic_mod.F90
               ${SRC_SHARE}/fvm_control_volume_mod.F90
               ${SRC_SHARE}/fvm_filter_mod.F90
               ${SRC_SHARE}/fvm_line_integrals_mod.F90
               ${SRC_SHARE}/fvm_line_integrals_flux_mod.F90
               ${SRC_SHARE}/fvm_mod.F90
               ${SRC_SHARE}/fvm_reconstruction_mod.F90
               ${SRC_SHARE}/global_norms_mod.F90
               ${SRC_SHARE}/gridgraph_mod.F90
               ${SRC_SHARE}/hybrid_mod.F90
               ${SRC_SHARE}/interpolate_mod.F90
               ${SRC_SHARE}/spelt_mod.F90
               ${SRC_SHARE}/kinds.F90
               ${SRC_SHARE}/linear_algebra_mod.F90
               ${SRC_SHARE}/ll_mod.F90
               ${SRC_SHARE}/mass_matrix_mod.F90
               ${SRC_SHARE}/mesh_mod.F90
               ${SRC_SHARE}/metagraph_mod.F90
               ${SRC_SHARE}/metis_mod.F90
               ${SRC_SHARE}/namelist_mod.F90
               ${SRC_SHARE}/parallel_mod.F90
               ${SRC_SHARE}/params_mod.F90
               ${SRC_SHARE}/physical_constants.F90
               ${SRC_SHARE}/quadrature_mod.F90
               ${SRC_SHARE}/reduction_mod.F90
               ${SRC_SHARE}/schedule_mod.F90
               ${SRC_SHARE}/solver_mod.F90
               ${SRC_SHARE}/spacecurve_mod.F90
               ${SRC_SHARE}/thread_mod.F90
               ${SRC_SHARE}/time_mod.F90
               ${SRC_SHARE}/viscosity_mod.F90
               ${SRC_UTILS}/shr_file_mod.F90
               ${SRC_UTILS}/shr_kind_mod.F90
               ${SRC_UTILS}/shr_mpi_mod.F90
               ${SRC_UTILS}/shr_sys_mod.F90
               ${SRC_UTILS}/shr_vmath_mod.F90
               ${SRC_BASE}/ref_state_mod.F90
               ${SRC_BASE}/jrio.c

)

SET(SWEQX_SRCS ${SWEQX_SRCS} ${SRC_BASE}/main.F90 CACHE INTERNAL "Sources for sweqx")

# Set up defaults
IF (NOT SWEQX_NP)
  SET (SWEQX_NP 8) 
ENDIF ()
IF (NOT SWEQX_PLEV)
  SET (SWEQX_PLEV 4) 
ENDIF ()
IF (NOT SWEQX_USE_PIO)
  SET (SWEQX_USE_PIO FALSE) 
ENDIF ()
IF (NOT SWEQX_USE_ENERGY)
  SET (SWEQX_USE_ENERGY FALSE) 
ENDIF ()
IF (NOT SWEQX_NC)
  SET (SWEQX_NC 4)
ENDIF ()

MESSAGE(STATUS "Building sweqx with:")
MESSAGE(STATUS "  NP = ${SWEQX_NP}")
MESSAGE(STATUS "  NC = ${SWEQX_NC}")
MESSAGE(STATUS "  PLEV = ${SWEQX_PLEV}")
MESSAGE(STATUS "  PIO = ${SWEQX_USE_PIO}")
MESSAGE(STATUS "  ENERGY = ${SWEQX_USE_ENERGY}")

# The following macro handles the compilation with the appropriate variables
############################################################################
# createTestExec(exec_name exec_type NP PLEV USE_PIO USE_ENERGY)
############################################################################
createTestExec(sweqx sweqx ${SWEQX_NP} ${SWEQX_NC} ${SWEQX_PLEV} 
               ${SWEQX_USE_PIO} ${SWEQX_USE_ENERGY})

