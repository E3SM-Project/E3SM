# Relevant Directories within the HOMME repository
# holding source files needed for building the library
# and executable spelt.

SET(UTILS_SHARE_DIR ${HOMME_SOURCE_DIR}/utils/csm_share)
SET(SRC_DIR         ${HOMME_SOURCE_DIR}/src)
SET(SRC_SHARE_DIR   ${HOMME_SOURCE_DIR}/src/share)
SET(PHYSICS_DIR     ${HOMME_SOURCE_DIR}/physics)
SET(UTILS_TIMING_DIR ${HOMME_SOURCE_DIR}/utils/timing)

# Make INCLUDE_DIRS global so the tests can access it
SET (SPELT_INCLUDE_DIRS ${PIO_INCLUDE_DIRS}
                        ${TIMING_INCLUDE_DIRS}
                        ${SRC_SHARE_DIR} CACHE INTERNAL "include directories for spelt")

SET(SPELT_SRCS ${SRC_SHARE_DIR}/bndry_mod.F90
               ${SRC_SHARE_DIR}/spelt_mod.F90
               ${SRC_DIR}/spelt_bench_mod.F90
               ${SRC_SHARE_DIR}/fvm_analytic_mod.F90
               ${SRC_DIR}/fvm_bsp_mod.F90
               ${SRC_SHARE_DIR}/fvm_control_volume_mod.F90
               ${SRC_SHARE_DIR}/fvm_reconstruction_mod.F90
               ${SRC_SHARE_DIR}/fvm_filter_mod.F90
               ${SRC_DIR}/fvm_init_mod.F90
               ${SRC_DIR}/fvm_transformation_mod.F90
               ${SRC_SHARE_DIR}/cg_mod.F90
               ${SRC_DIR}/checksum_mod.F90
               ${SRC_DIR}/netcdf_io_mod.F90
               ${SRC_DIR}/common_movie_mod.F90
               ${SRC_DIR}/common_io_mod.F90
               ${SRC_SHARE_DIR}/control_mod.F90
               ${SRC_SHARE_DIR}/coordinate_systems_mod.F90
               ${SRC_SHARE_DIR}/cube_mod.F90
               ${SRC_SHARE_DIR}/derivative_mod.F90
               ${SRC_SHARE_DIR}/dimensions_mod.F90
               ${SRC_SHARE_DIR}/dof_mod.F90
               ${SRC_SHARE_DIR}/domain_mod.F90
               ${SRC_SHARE_DIR}/edge_mod.F90
               ${SRC_SHARE_DIR}/element_mod.F90
               ${SRC_SHARE_DIR}/gridgraph_mod.F90
               ${SRC_SHARE_DIR}/hybrid_mod.F90
               ${SRC_DIR}/interp_movie_mod.F90
               ${SRC_DIR}/netcdf_io_mod.F90
               ${SRC_DIR}/shal_movie_mod.F90
               ${SRC_SHARE_DIR}/interpolate_mod.F90
               ${SRC_SHARE_DIR}/kinds.F90
               ${SRC_SHARE_DIR}/ll_mod.F90
               ${SRC_SHARE_DIR}/mass_matrix_mod.F90
               ${SRC_SHARE_DIR}/metagraph_mod.F90
               ${SRC_SHARE_DIR}/metis_mod.F90
               ${SRC_SHARE_DIR}/mesh_mod.F90
               ${SRC_SHARE_DIR}/namelist_mod.F90
               ${SRC_DIR}/netcdf_io_mod.F90
               ${SRC_SHARE_DIR}/parallel_mod.F90
               ${SRC_SHARE_DIR}/params_mod.F90
               ${SRC_SHARE_DIR}/physical_constants.F90
               ${SRC_SHARE_DIR}/quadrature_mod.F90
               ${SRC_SHARE_DIR}/reduction_mod.F90
               ${SRC_DIR}/repro_sum_mod.F90
               ${SRC_DIR}/restart_mod.F90
               ${SRC_DIR}/restart_io_mod.F90
               ${SRC_SHARE_DIR}/schedule_mod.F90
               ${SRC_SHARE_DIR}/spacecurve_mod.F90
               ${SRC_DIR}/surfaces_mod.F90
               ${SRC_DIR}/state_mod.F90
               ${SRC_DIR}/repro_sum_mod.F90
               ${SRC_SHARE_DIR}/thread_mod.F90
               ${SRC_SHARE_DIR}/time_mod.F90
               ${SRC_SHARE_DIR}/viscosity_mod.F90
               ${SRC_SHARE_DIR}/global_norms_mod.F90
               ${SRC_DIR}/pio_io_mod.F90
               ${UTILS_SHARE_DIR}/shr_kind_mod.F90
               ${UTILS_SHARE_DIR}/shr_mpi_mod.F90
               ${UTILS_SHARE_DIR}/shr_sys_mod.F90
               ${UTILS_SHARE_DIR}/shr_file_mod.F90
               ${UTILS_SHARE_DIR}/shr_vmath_mod.F90 
)

# Make SRCS global so the tests can access it
SET(SPELT_SRCS ${SPELT_SRCS} ${SRC_DIR}/spelt_main.F90 CACHE INTERNAL "Sources for spelt")

# Set up defaults
IF (NOT SPELT_NP)
  SET (SPELT_NP 8) 
ENDIF ()
IF (NOT SPELT_PLEV)
  SET (SPELT_PLEV 20) 
ENDIF ()
IF (NOT SPELT_USE_PIO)
  SET (SPELT_USE_PIO FALSE) 
ENDIF ()
IF (NOT SPELT_USE_ENERGY)
  SET (SPELT_USE_ENERGY FALSE) 
ENDIF ()
IF (NOT SPELT_NC)
  SET (SPELT_NC 4)
ENDIF ()

MESSAGE(STATUS "Building spelt with:")
MESSAGE(STATUS "  NP = ${SPELT_NP}")
MESSAGE(STATUS "  NC = ${SPELT_NC}")
MESSAGE(STATUS "  PLEV = ${SPELT_PLEV}")
MESSAGE(STATUS "  PIO = ${SPELT_USE_PIO}")
MESSAGE(STATUS "  ENERGY = ${SPELT_USE_ENERGY}")

# The following macro handles the compilation with the appropriate variables
############################################################################
# createTestExec(exec_name exec_type SRCS NP PLEV USE_PIO USE_ENERGY)
############################################################################
createTestExec(spelt spelt ${SPELT_NP} ${SPELT_NC} ${SPELT_PLEV} 
               ${SPELT_USE_PIO}  ${SPELT_USE_ENERGY})
