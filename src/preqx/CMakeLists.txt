# Relevant Directories within the HOMME repository
# holding source files needed for building the library
# and executable preqx.

INCLUDE_DIRECTORIES ( ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS} )
LINK_DIRECTORIES    ( ${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS} )

SET(UTILS_SHARE_DIR ${HOMME_SOURCE_DIR}/utils/csm_share)
SET(SRC_DIR         ${HOMME_SOURCE_DIR}/src)
SET(SRC_SHARE_DIR   ${HOMME_SOURCE_DIR}/src/share)
SET(PHYSICS_DIR     ${HOMME_SOURCE_DIR}/physics)
SET(UTILS_TIMING_DIR ${HOMME_SOURCE_DIR}/utils/timing)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}
                    ${PIO_INCLUDE_DIRS}
                    ${TIMING_INCLUDE_DIRS}
                    ${SRC_SHARE_DIR})

CONFIGURE_FILE(config.h.cmake.in config.h)
ADD_DEFINITIONS(-DHAVE_CONFIG_H)

SET(SRCS_PREQX ${SRC_SHARE_DIR}/kinds.F90
					     ${SRC_SHARE_DIR}/dimensions_mod.F90
					     ${SRC_SHARE_DIR}/control_mod.F90
					     ${SRC_SHARE_DIR}/physical_constants.F90
					     ${SRC_SHARE_DIR}/coordinate_systems_mod.F90
					     ${SRC_SHARE_DIR}/params_mod.F90
					     ${SRC_SHARE_DIR}/gridgraph_mod.F90
					     ${SRC_SHARE_DIR}/edge_mod.F90
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

ADD_LIBRARY(homme_preqx ${SRCS_PREQX})
TARGET_LINK_LIBRARIES(homme_preqx pio timing ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

ADD_EXECUTABLE(preqx ${SRC_DIR}/prim_main.F90)

ADD_DEPENDENCIES(preqx blas lapack)
TARGET_LINK_LIBRARIES(preqx homme_preqx)

IF (NOT HOMME_FIND_BLASLAPACK)
  TARGET_LINK_LIBRARIES(preqx lapack blas)
  ADD_DEPENDENCIES(preqx blas lapack)
ENDIF()

INSTALL(TARGETS preqx
        RUNTIME DESTINATION bin)

# Macro to create the individual tests
macro(createPreqxTest execName macroNP macroPLEV macroUSE_PIO
                               macroWITH_ENERGY)

  # Backup the cmake variables
  SET(tempNP ${NUM_POINTS})
  SET(tempPLEV ${NUM_PLEV})
  SET(tempUSE_PIO ${PIO})
  SET(tempWITH_ENERGY ${ENERGY_DIAGNOSTICS})

  # Set the variable to the macro variables
  SET(NUM_POINTS ${macroNP})
  SET(NUM_PLEV ${macroPLEV})
  SET(PIO ${macroUSE_PIO})
  SET(ENERGY_DIAGNOSTICS ${macroWITH_ENERGY})

  # This is needed to compile the test executables with the correct options
  CONFIGURE_FILE(config.h.cmake.in
                ${CMAKE_CURRENT_BINARY_DIR}/${execName}_config/config.h)

  # Reset the include_directories
  SET_DIRECTORY_PROPERTIES(PROPERTY INCLUDE_DIRECTORIES "")

  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/${execName}_config
                      ${PIO_INCLUDE_DIRS}
                      ${TIMING_INCLUDE_DIRS}
                      ${SRC_SHARE_DIR}
                      )

  ADD_EXECUTABLE(${execName} ${SRCS_PREQX} ${SRC_DIR}/prim_main.F90)

  # More thought needs to go into these options
  #IF (macroUSE_PIO)
  #ELSE ()
  #ENDIF ()

  #IF (macroWITH_ENERGY)
  #ENDIF ()

  TARGET_LINK_LIBRARIES(${execName} pio timing ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

  # Move the module files out of the way so the parallel build 
  # doesn't have a race condition
  SET_TARGET_PROPERTIES(${execName} PROPERTIES
     Fortran_MODULE_DIRECTORY ${execName}_modules)

  IF (NOT HOMME_FIND_BLASLAPACK)
    TARGET_LINK_LIBRARIES(${execName} lapack blas)
    ADD_DEPENDENCIES(${execName} blas lapack)
  ENDIF()

  INSTALL(TARGETS ${execName}
          RUNTIME DESTINATION tests)

  # Restore the original the cmake variables
  SET(NUM_POINTS ${tempNP})
  SET(NUM_PLEV ${tempPLEV})
  SET(PIO ${tempUSE_PIO})
  SET(ENERGY_DIAGNOSTICS ${tempWITH_ENERGY})

endmacro(createPreqxTest)

# There are only three separate executables that need to be compiled
#   separately
#                     NP  PLEV  USE_PIO  WITH_ENERGY
createPreqxTest(baroA  8    20    FALSE        TRUE)
createPreqxTest(baroB  8    20     TRUE       FALSE)
createPreqxTest(baroC  4    26    FALSE        TRUE)

