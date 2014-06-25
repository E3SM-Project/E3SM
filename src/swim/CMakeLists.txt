# Relevant Directories within the HOMME repository
# holding source files needed for building the library
# and executable swim.

INCLUDE_DIRECTORIES ( ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS} )
LINK_DIRECTORIES    ( ${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS} )

SET(UTILS_SHARE_DIR  ${HOMME_SOURCE_DIR}/utils/csm_share)
SET(SRC_DIR          ${HOMME_SOURCE_DIR}/src)
SET(SRC_SHARE_DIR    ${HOMME_SOURCE_DIR}/src/share)
SET(UTILS_TIMING_DIR ${HOMME_SOURCE_DIR}/utils/timing)
SET(UTILS_TRILINOS_DIR ${HOMME_SOURCE_DIR}/utils/trilinos)

# Make INCLUDE_DIRS global so the tests can access it
SET (SWIM_INCLUDE_DIRS ${PIO_INCLUDE_DIRS}
                       ${TIMING_INCLUDE_DIRS}
                       ${SRC_SHARE_DIR} CACHE INTERNAL "include directories for swim") 

# Make SRCS global so the tests can access it
SET(SWIM_SRCS ${SRC_SHARE_DIR}/kinds.F90
               ${SRC_DIR}/advance_mod.F90
               ${SRC_DIR}/checksum_mod.F90
               ${SRC_DIR}/common_io_mod.F90
               ${SRC_DIR}/common_movie_mod.F90
               ${SRC_DIR}/init_mod.F90
               ${SRC_DIR}/interp_movie_mod.F90
               ${SRC_DIR}/netcdf_io_mod.F90
               ${SRC_DIR}/pio_io_mod.F90
               ${SRC_DIR}/repro_sum_mod.F90
               ${SRC_DIR}/restart_io_mod.F90
               ${SRC_DIR}/restart_mod.F90
               ${SRC_DIR}/rk_mod.F90
               ${SRC_DIR}/shallow_water_mod.F90
               ${SRC_DIR}/shal_movie_mod.F90
               ${SRC_DIR}/state_mod.F90
               ${SRC_DIR}/sweq_mod.F90
               ${SRC_DIR}/types_mod.F90
               ${SRC_SHARE_DIR}/fvm_bsp_mod.F90
               ${SRC_SHARE_DIR}/fvm_transformation_mod.F90
               ${SRC_SHARE_DIR}/bndry_mod.F90
               ${SRC_SHARE_DIR}/cg_mod.F90
               ${SRC_SHARE_DIR}/control_mod.F90
               ${SRC_SHARE_DIR}/coordinate_systems_mod.F90
               ${SRC_SHARE_DIR}/cube_mod.F90
               ${SRC_SHARE_DIR}/derivative_mod.F90
               ${SRC_SHARE_DIR}/dimensions_mod.F90
               ${SRC_SHARE_DIR}/dof_mod.F90
               ${SRC_SHARE_DIR}/domain_mod.F90
               ${SRC_SHARE_DIR}/edge_mod.F90
               ${SRC_SHARE_DIR}/element_mod.F90
               ${SRC_SHARE_DIR}/filter_mod.F90
               ${SRC_SHARE_DIR}/fvm_analytic_mod.F90
               ${SRC_SHARE_DIR}/fvm_control_volume_mod.F90
               ${SRC_SHARE_DIR}/fvm_filter_mod.F90
               ${SRC_SHARE_DIR}/fvm_line_integrals_mod.F90
               ${SRC_SHARE_DIR}/fvm_line_integrals_flux_mod.F90
               ${SRC_SHARE_DIR}/fvm_mod.F90
               ${SRC_SHARE_DIR}/fvm_reconstruction_mod.F90
               ${SRC_SHARE_DIR}/global_norms_mod.F90
               ${SRC_SHARE_DIR}/gridgraph_mod.F90
               ${SRC_SHARE_DIR}/hybrid_mod.F90
               ${SRC_SHARE_DIR}/interpolate_mod.F90
               ${SRC_SHARE_DIR}/spelt_mod.F90
               ${SRC_SHARE_DIR}/linear_algebra_mod.F90
               ${SRC_SHARE_DIR}/ll_mod.F90
               ${SRC_SHARE_DIR}/mass_matrix_mod.F90
               ${SRC_SHARE_DIR}/mesh_mod.F90
               ${SRC_SHARE_DIR}/metagraph_mod.F90
               ${SRC_SHARE_DIR}/metis_mod.F90
               ${SRC_SHARE_DIR}/namelist_mod.F90
               ${SRC_SHARE_DIR}/parallel_mod.F90
               ${SRC_SHARE_DIR}/params_mod.F90
               ${SRC_SHARE_DIR}/physical_constants.F90
               ${SRC_SHARE_DIR}/quadrature_mod.F90
               ${SRC_SHARE_DIR}/reduction_mod.F90
               ${SRC_SHARE_DIR}/schedule_mod.F90
               ${SRC_SHARE_DIR}/schedtype_mod.F90
               ${SRC_SHARE_DIR}/solver_mod.F90
               ${SRC_SHARE_DIR}/spacecurve_mod.F90
               ${SRC_SHARE_DIR}/thread_mod.F90
               ${SRC_SHARE_DIR}/time_mod.F90
               ${SRC_SHARE_DIR}/viscosity_mod.F90
               ${UTILS_SHARE_DIR}/shr_file_mod.F90
               ${UTILS_SHARE_DIR}/shr_kind_mod.F90
               ${UTILS_SHARE_DIR}/shr_mpi_mod.F90
               ${UTILS_SHARE_DIR}/shr_sys_mod.F90
               ${UTILS_SHARE_DIR}/shr_vmath_mod.F90
               ${UTILS_SHARE_DIR}/shr_vmath_fwrap.c
               ${UTILS_SHARE_DIR}/shr_const_mod.F90
               ${UTILS_SHARE_DIR}/shr_spfn_mod.F90
               ${SRC_DIR}/ref_state_mod.F90
               ${SRC_SHARE_DIR}/perfmodel_mod.F90
               ${SRC_DIR}/jrio.c
               ${UTILS_TRILINOS_DIR}/trilinosNoxSolver.cpp
               ${UTILS_TRILINOS_DIR}/trilinosModelEvaluator.cpp
               ${UTILS_TRILINOS_DIR}/block_precon_interface.cpp
               ${UTILS_TRILINOS_DIR}/precon_interface.cpp
               ${SRC_DIR}/derived_type_mod.F90
               ${SRC_DIR}/precon_type_mod.F90
               ${SRC_DIR}/implicit_mod.F90
               ${SRC_DIR}/precon_mod.F90
               ${UTILS_TRILINOS_DIR}/belos_helm.cpp
)


SET(SWIM_SRCS ${SWIM_SRCS} ${SRC_DIR}/main.F90 CACHE INTERNAL "Sources for swim")

#Set up defaults
IF (NOT SWIM_NP)
  SET (SWIM_NP 4)
ENDIF ()
IF (NOT SWIM_NC)
  SET (SWIM_NC 4)
ENDIF ()
IF (NOT SWIM_PLEV)
  SET (SWIM_PLEV 1)
ENDIF ()
IF (NOT SWIM_USE_PIO)
  SET (SWIM_USE_PIO FALSE)
ENDIF ()
IF (NOT SWIM_USE_ENERGY)
  SET (SWIM_USE_ENERGY FALSE)
ENDIF ()

MESSAGE(STATUS "Building swim with:")
MESSAGE(STATUS "  NP = ${SWIM_NP}")
MESSAGE(STATUS "  NC = ${SWIM_NC}")
MESSAGE(STATUS "  PLEV = ${SWIM_PLEV}")
MESSAGE(STATUS "  PIO = ${SWIM_USE_PIO}")
MESSAGE(STATUS "  ENERGY = ${SWIM_USE_ENERGY}")

# The following macro handles the compilation with the appropriate variables
############################################################################
# createTestExec(exec_name exec_type NP PLEV USE_PIO USE_ENERGY)
############################################################################
createTestExec(swim swim ${SWIM_NP} ${SWIM_NC} ${SWIM_PLEV} ${SWIM_USE_PIO} ${SWIM_USE_ENERGY})

SET_PROPERTY(TARGET swim PROPERTY LINKER_LANGUAGE Fortran)

TARGET_LINK_LIBRARIES(swim ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})

