# Relevant Directories within the HOMME repository
# holding source files needed for building the library
# and executable swim.

INCLUDE_DIRECTORIES ( ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS} )
LINK_DIRECTORIES    ( ${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS} )

SET(SRC_BASE    ${CMAKE_CURRENT_SOURCE_DIR}/..)
SET(SRC_SHARE   ${SRC_BASE}/share)
SET(SRC_UTILS   ${SRC_BASE}/../utils/csm_share)
SET(SRC_TRIL    ${SRC_BASE}/../utils/trilinos)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}
                    ${CMAKE_CURRENT_BINARY_DIR}/share
                    ${PIO_INCLUDE_DIRS}
                    ${TIMING_INCLUDE_DIRS})

# Set up swim specfic defaults 
IF (NOT SWIM_NP)
  SET (SWIM_NP 4 CACHE STRING "Default number of points (NP) in the element (default value is 4)")
ENDIF ()
IF (NOT SWIM_NC)
  SET (SWIM_NC 4 CACHE STRING "Default number of cells (NC*NC) in each element (default value is 4)")
ENDIF ()
IF (NOT SWIM_PLEV)
  SET (SWIM_PLEV 1 CACHE STRING "Default number of levels (PLEV) in each element (default value is 1)")
ENDIF ()
IF (NOT SWIM_NTRAC)
  SET (SWIM_NTRAC 1 CACHE STRING "Default number of tracers (NTRAC) (default value is 1)")
ENDIF ()

MESSAGE(STATUS "Building swim with:")
MESSAGE(STATUS "  NP = ${SWIM_NP}")
MESSAGE(STATUS "  NC = ${SWIM_NC}")
MESSAGE(STATUS "  PLEV = ${SWIM_PLEV}")
MESSAGE(STATUS "  NTRAC = ${SWIM_NTRAC}")

# For the configuration
SET(NUM_POINTS ${SWIM_NP})
SET(NUM_CELLS ${SWIM_NC})
SET(NUM_PLEV ${SWIM_PLEV})
SET(NUM_TRACERS ${SWIM_NTRAC})

IF (${macroUSE_PIO})
  SET(PIO TRUE)
  SET(PIO_INTERP)
ELSE ()
  SET(PIO)
  SET(PIO_INTERP TRUE)
ENDIF ()

SET(ENERGY_DIAGNOSTICS ${macroWITH_ENERGY})

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
               ${SRC_UTILS}/shr_vmath_fwrap.c
               ${SRC_BASE}/ref_state_mod.F90
               ${SRC_SHARE}/perfmodel_mod.F90
               ${SRC_BASE}/jrio.c
)

SET(SRCS_SWIM     ${SRC_TRIL}/donoxloca.cpp
                  ${SRC_TRIL}/block_precon_interface.cpp
                  ${SRC_TRIL}/noxlocainterface.cpp
                  ${SRC_TRIL}/precon_interface.cpp
                  ${SRC_BASE}/derived_type_mod.F90
                  ${SRC_BASE}/implicit_mod.F90
                  ${SRC_BASE}/precon_mod.F90
)

ADD_LIBRARY(homme_swim ${SRCS_F90} ${SRCS_SWIM})
ADD_DEPENDENCIES(homme_swim pio timing)

ADD_EXECUTABLE(swim ${SRC_BASE}/main.F90)

ADD_DEPENDENCIES(swim homme_swim pio timing)

SET_PROPERTY(TARGET swim PROPERTY LINKER_LANGUAGE Fortran)

TARGET_LINK_LIBRARIES(swim homme_swim pio timing ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})

INSTALL(TARGETS swim
        RUNTIME DESTINATION bin)
