# Relevant Directories within the HOMME repository
# holding source files needed for building the library
# and executable sweqx.

SET(SRC_BASE    ${CMAKE_CURRENT_SOURCE_DIR}/..)
SET(SRC_SHARE   ${SRC_BASE}/share)
SET(SRC_UTILS   ${SRC_BASE}/../utils/csm_share)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}
                    ${PIO_INCLUDE_DIRS}
                    ${TIMING_INCLUDE_DIRS}
                    ${CMAKE_CURRENT_SOURCE_DIR})

CONFIGURE_FILE(config.h.cmake.in config.h)
ADD_DEFINITIONS(-DHAVE_CONFIG_H)

SET(SRCS_F90
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

ADD_LIBRARY(homme_sweqx ${SRCS_F90})
TARGET_LINK_LIBRARIES(homme_sweqx pio timing ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

ADD_EXECUTABLE(sweqx ${SRC_BASE}/main.F90)

TARGET_LINK_LIBRARIES(sweqx homme_sweqx)

IF (NOT HOMME_FIND_BLASLAPACK)
  TARGET_LINK_LIBRARIES(sweqx lapack blas)
  ADD_DEPENDENCIES(sweqx blas lapack)
ENDIF()

INSTALL(TARGETS sweqx
        RUNTIME DESTINATION bin)

# Macro to create the individual tests
macro(createSweqxTest execName macroNP macroPLEV macroUSE_PIO
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
  CONFIGURE_FILE(config.h.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/${execName}_config/config.h)

  # Reset the include_directories
  SET_DIRECTORY_PROPERTIES(PROPERTY INCLUDE_DIRECTORIES "")

  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/${execName}_config
                      ${PIO_INCLUDE_DIRS}
                      ${TIMING_INCLUDE_DIRS}
                      ${SRC_SHARE_DIR}
                      )

  ADD_EXECUTABLE(${execName} ${SRCS_F90} ${SRC_BASE}/main.F90)

  # Some more thought needs to go into these entries
  #IF (macroUSE_PIO)
  #ELSE ()
  #ENDIF ()

  #IF (macroWITH_ENERGY)
  #ENDIF ()

  TARGET_LINK_LIBRARIES(${execName} pio timing ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

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

endmacro(createSweqxTest)

# There are only three separate executables that need to be compiled
#   separately
#                     NP  PLEV  USE_PIO  WITH_ENERGY
createSweqxTest(swtcA  8     4    FALSE       FALSE)
createSweqxTest(swtcB  4     1    FALSE       FALSE)
createSweqxTest(swtcC  8     1    FALSE       FALSE)



