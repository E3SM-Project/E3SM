# Relevant Directories within the HOMME repository
# holding source files needed for building the library
# and executable sweqx.

macro(sweqx_setup)

SET(TARGET_DIR        ${HOMME_SOURCE_DIR}/src/sweqx)
SET(SRC_BASE          ${HOMME_SOURCE_DIR}/src)
SET(SRC_SHARE         ${HOMME_SOURCE_DIR}/src/share)
SET(SRC_UTILS         ${HOMME_SOURCE_DIR}/utils/csm_share)
SET(UTILS_TIMING_DIR  ${CMAKE_BINARY_DIR}/utils/timing)

# Find F90 files in target directory
FILE(GLOB TARGET_F90  ${TARGET_DIR}/*.F90)

# Make INCLUDE_DIRS global so the tests can access it
SET (EXEC_INCLUDE_DIRS ${PIO_INCLUDE_DIRS} ${UTILS_TIMING_DIR} )

SET(SWEQX_SRCS
  ${TARGET_F90}
  ${SRC_BASE}/checksum_mod.F90
  ${SRC_BASE}/common_io_mod.F90
  ${SRC_BASE}/common_movie_mod.F90
  ${SRC_BASE}/interp_movie_mod.F90
  ${SRC_BASE}/netcdf_io_mod.F90
  ${SRC_BASE}/pio_io_mod.F90
  ${SRC_BASE}/repro_sum_mod.F90
  ${SRC_BASE}/restart_io_mod.F90
  ${SRC_SHARE}/bndry_mod_base.F90
  ${SRC_SHARE}/cg_mod.F90
  ${SRC_SHARE}/control_mod.F90
  ${SRC_SHARE}/coordinate_systems_mod.F90
  ${SRC_SHARE}/cube_mod.F90
  ${SRC_SHARE}/derivative_mod_base.F90
  ${SRC_SHARE}/dimensions_mod.F90
  ${SRC_SHARE}/dof_mod.F90
  ${SRC_SHARE}/domain_mod.F90
  ${SRC_SHARE}/edge_mod_base.F90
  ${SRC_SHARE}/element_mod.F90
  ${SRC_SHARE}/fvm_bsp_mod.F90
  ${SRC_SHARE}/fvm_transformation_mod.F90
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
  ${SRC_SHARE}/schedtype_mod.F90
  ${SRC_SHARE}/edgetype_mod.F90
  ${SRC_SHARE}/solver_init_mod_base.F90
  ${SRC_SHARE}/solver_mod.F90
  ${SRC_SHARE}/spacecurve_mod.F90
  ${SRC_SHARE}/thread_mod.F90
  ${SRC_SHARE}/time_mod.F90
  ${SRC_SHARE}/unit_tests_mod.F90
  ${SRC_SHARE}/viscosity_mod_base.F90
  ${SRC_UTILS}/shr_file_mod.F90
  ${SRC_UTILS}/shr_kind_mod.F90
  ${SRC_UTILS}/shr_mpi_mod.F90
  ${SRC_UTILS}/shr_sys_mod.F90
  ${SRC_UTILS}/shr_vmath_mod.F90
  ${SRC_UTILS}/shr_const_mod.F90
  ${SRC_UTILS}/shr_spfn_mod.F90
  ${SRC_BASE}/ref_state_mod.F90
  ${SRC_BASE}/jrio.c
)

SET(EXEC_SOURCES ${SWEQX_SRCS} )

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
IF (NOT QSIZE_D)
  SET (QSIZE_D 4)
ENDIF ()

endmacro(sweqx_setup)

############################################################################
# createTestExec(exec_name exec_type NP PLEV USE_PIO USE_ENERGY)
############################################################################

sweqx_setup()
createTestExec(sweqx sweqx ${SWEQX_NP} ${SWEQX_NC} ${SWEQX_PLEV} ${SWEQX_USE_PIO} ${SWEQX_USE_ENERGY} ${QSIZE_D})

