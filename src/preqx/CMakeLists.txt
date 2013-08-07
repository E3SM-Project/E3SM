# Relevant Directories within the HOMME repository
# holding source files needed for building the library
# and executable preqx.

SET(UTILS_SHARE_DIR ${HOMME_SOURCE_DIR}/utils/csm_share)
SET(SRC_DIR         ${HOMME_SOURCE_DIR}/src)
SET(SRC_SHARE_DIR   ${HOMME_SOURCE_DIR}/src/share)
SET(PHYSICS_DIR     ${HOMME_SOURCE_DIR}/physics)
SET(UTILS_TIMING_DIR ${HOMME_SOURCE_DIR}/utils/timing)

# Make INCLUDE_DIRS global so the tests can access it
SET (PREQX_INCLUDE_DIRS ${PIO_INCLUDE_DIRS}
                        ${TIMING_INCLUDE_DIRS}
                        ${SRC_SHARE_DIR} CACHE INTERNAL "include directories for preqx")

SET(PREQX_SRCS ${SRC_SHARE_DIR}/kinds.F90
               ${SRC_SHARE_DIR}/dimensions_mod.F90
               ${SRC_SHARE_DIR}/control_mod.F90
               ${SRC_SHARE_DIR}/physical_constants.F90
               ${SRC_SHARE_DIR}/coordinate_systems_mod.F90
               ${SRC_SHARE_DIR}/params_mod.F90
               ${SRC_SHARE_DIR}/gridgraph_mod.F90
               ${SRC_SHARE_DIR}/edge_mod.F90
               ${SRC_SHARE_DIR}/cuda_mod.F90
               ${SRC_SHARE_DIR}/element_mod.F90
               ${SRC_SHARE_DIR}/parallel_mod.F90
               ${SRC_SHARE_DIR}/hybrid_mod.F90
               ${SRC_SHARE_DIR}/reduction_mod.F90
               ${SRC_SHARE_DIR}/quadrature_mod.F90
               ${SRC_SHARE_DIR}/global_norms_mod.F90
               ${SRC_SHARE_DIR}/ll_mod.F90
               ${SRC_SHARE_DIR}/metagraph_mod.F90
               ${SRC_SHARE_DIR}/schedule_mod.F90
               ${SRC_SHARE_DIR}/thread_mod.F90
               ${SRC_SHARE_DIR}/bndry_mod.F90
               ${SRC_SHARE_DIR}/derivative_mod.F90
               ${SRC_SHARE_DIR}/viscosity_mod.F90
               ${SRC_SHARE_DIR}/spacecurve_mod.F90
               ${SRC_SHARE_DIR}/cube_mod.F90
               ${SRC_SHARE_DIR}/interpolate_mod.F90
               ${PHYSICS_DIR}/physics_types_mod.F90
               ${SRC_SHARE_DIR}/dof_mod.F90
               ${SRC_DIR}/common_io_mod.F90
               ${SRC_SHARE_DIR}/spelt_mod.F90
               ${SRC_SHARE_DIR}/fvm_mod.F90
               ${SRC_SHARE_DIR}/fvm_control_volume_mod.F90
               ${SRC_SHARE_DIR}/fvm_line_integrals_mod.F90 
               ${SRC_SHARE_DIR}/fvm_line_integrals_flux_mod.F90 
               ${SRC_SHARE_DIR}/fvm_analytic_mod.F90 
               ${SRC_SHARE_DIR}/fvm_filter_mod.F90 
               ${SRC_SHARE_DIR}/fvm_reconstruction_mod.F90 
               ${SRC_DIR}/pio_io_mod.F90 
               ${SRC_DIR}/repro_sum_mod.F90 
               ${PHYSICS_DIR}/physics_io_mod.F90 
               ${SRC_SHARE_DIR}/time_mod.F90 
               ${PHYSICS_DIR}/gravity_wave_drag_mod.F90 
               ${SRC_SHARE_DIR}/physics_mod.F90 
               ${SRC_SHARE_DIR}/hybvcoord_mod.F90 
               ${SRC_DIR}/aquaplanet.F90 
               ${SRC_DIR}/aquaplanet_io_mod.F90 
               ${SRC_DIR}/interp_movie_mod.F90 
               ${UTILS_SHARE_DIR}/shr_kind_mod.F90 
               ${UTILS_SHARE_DIR}/shr_mpi_mod.F90 
               ${UTILS_SHARE_DIR}/shr_sys_mod.F90 
               ${UTILS_SHARE_DIR}/shr_file_mod.F90 
               ${UTILS_SHARE_DIR}/shr_vmath_mod.F90 
               ${UTILS_SHARE_DIR}/shr_vmath_fwrap.c 
               ${UTILS_TIMING_DIR}/perf_mod.F90 
               ${UTILS_TIMING_DIR}/perf_utils.F90 
               ${SRC_SHARE_DIR}/domain_mod.F90 
               ${SRC_DIR}/restart_io_mod.F90 
               ${SRC_SHARE_DIR}/cg_mod.F90 
               ${SRC_SHARE_DIR}/linear_algebra_mod.F90 
               ${SRC_SHARE_DIR}/solver_mod.F90 
               ${SRC_SHARE_DIR}/filter_mod.F90 
               ${SRC_DIR}/common_movie_mod.F90 
               ${SRC_SHARE_DIR}/namelist_mod.F90 
               ${SRC_SHARE_DIR}/mass_matrix_mod.F90 
               ${SRC_DIR}/checksum_mod.F90 
               ${SRC_DIR}/baroclinic_inst_mod.F90 
               ${SRC_SHARE_DIR}/diffusion_mod.F90 
               ${SRC_SHARE_DIR}/prim_si_mod.F90 
               ${SRC_DIR}/asp_tests.F90 
               ${SRC_SHARE_DIR}/prim_state_mod.F90 
               ${SRC_DIR}/prim_restart_mod.F90 
               ${SRC_SHARE_DIR}/prim_si_ref_mod.F90 
               ${SRC_SHARE_DIR}/prim_advance_mod.F90 
               ${SRC_SHARE_DIR}/prim_advection_mod.F90 
               ${SRC_DIR}/forcing_mod.F90 
               ${PHYSICS_DIR}/emanuel/convect43c.F90 
               ${PHYSICS_DIR}/heldsuarez/held_suarez_mod.F90 
               ${PHYSICS_DIR}/column_model_mod.F90 
               ${PHYSICS_DIR}/column_types_mod.F90 
               ${SRC_SHARE_DIR}/metis_mod.F90 
               ${SRC_SHARE_DIR}/prim_driver_mod.F90 
               ${SRC_DIR}/vertical_mod.F90 
               ${PHYSICS_DIR}/emanuel/mjo_mod.F90 
               ${PHYSICS_DIR}/multicloud_mod.F90 
               ${PHYSICS_DIR}/eigenmodes_mod.F90 
               ${SRC_DIR}/prim_movie_mod.F90 
               ${SRC_DIR}/surfaces_mod.F90 
               ${SRC_DIR}/repro_sum_mod.F90 
               ${SRC_DIR}/interpolate_driver_mod.F90 
               ${SRC_DIR}/netcdf_io_mod.F90 
               ${SRC_SHARE_DIR}/mesh_mod.F90 
               ${SRC_SHARE_DIR}/Manager.F90
)

# Make SRCS global so the tests can access it
SET(PREQX_SRCS ${PREQX_SRCS} ${SRC_DIR}/prim_main.F90 CACHE INTERNAL "Sources for preqx")

# Set up defaults
IF (NOT PREQX_NP)
  SET (PREQX_NP 4) 
ENDIF ()
IF (NOT PREQX_PLEV)
  SET (PREQX_PLEV 20) 
ENDIF ()
IF (NOT PREQX_USE_PIO)
  SET (PREQX_USE_PIO FALSE) 
ENDIF ()
IF (NOT PREQX_USE_ENERGY)
  SET (PREQX_USE_ENERGY FALSE) 
ENDIF ()
IF (NOT PREQX_NC)
  SET (PREQX_NC 4)
ENDIF ()

MESSAGE(STATUS "Building preqx with:")
MESSAGE(STATUS "  NP = ${PREQX_NP}")
MESSAGE(STATUS "  NC = ${PREQX_NC}")
MESSAGE(STATUS "  PLEV = ${PREQX_PLEV}")
MESSAGE(STATUS "  PIO = ${PREQX_USE_PIO}")
MESSAGE(STATUS "  ENERGY = ${PREQX_USE_ENERGY}")

IF (${PREQX_USE_CUDA_FORTRAN})
  SET(USE_CUDA_FORTRAN TRUE)
  MESSAGE(STATUS "  Using CUDA Fortran")
ENDIF ()

IF (${ENABLE_HORIZ_OPENMP})
  SET(HORIZ_OPENMP TRUE)
  MESSAGE(STATUS "  Using HORIZ_OPENMP")
ENDIF ()

IF (${ENABLE_VERT_OPENMP})
  SET(VERT_OPENMP TRUE)
  MESSAGE(STATUS "  Using VERT_OPENMP")
ENDIF ()

IF (${ENABLE_HORIZ_OPENMP} AND ${ENABLE_VERT_OPENMP})
  SET(NESTED_OPENMP TRUE)
ENDIF ()

# The following macro handles the compilation with the appropriate variables
############################################################################
# createTestExec(exec_name exec_type NP PLEV USE_PIO USE_ENERGY)
############################################################################
createTestExec(preqx preqx ${PREQX_NP} ${PREQX_NC} ${PREQX_PLEV} 
               ${PREQX_USE_PIO}  ${PREQX_USE_ENERGY})

IF (${BUILD_PREQX_FVM})
  MESSAGE(STATUS "  Also building preqx_fvm")
  createTestExec(preqx_fvm preqx ${PREQX_NP} ${PREQX_NC} ${PREQX_PLEV} 
                 ${PREQX_USE_PIO}  ${PREQX_USE_ENERGY})
  SET_TARGET_PROPERTIES(preqx_fvm PROPERTIES COMPILE_DEFINITIONS "_FVM")
ENDIF ()

IF (${BUILD_PREQX_SPELT})
  MESSAGE(STATUS "  Also building preqx_spelt")
  createTestExec(preqx_spelt preqx ${PREQX_NP} ${PREQX_NC} ${PREQX_PLEV} 
                 ${PREQX_USE_PIO}  ${PREQX_USE_ENERGY})
  SET_TARGET_PROPERTIES(preqx_spelt PROPERTIES COMPILE_DEFINITIONS "_SPELT")
ENDIF ()

