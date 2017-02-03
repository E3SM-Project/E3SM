# Relevant Directories within the HOMME repository
# holding source files needed for building the library
# and executable prim.
macro(prim_setup)

INCLUDE_DIRECTORIES ( ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS} )
LINK_DIRECTORIES    ( ${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS} )

SET(UTILS_SHARE_DIR ${HOMME_SOURCE_DIR}/utils/csm_share)
SET(SRC_DIR         ${HOMME_SOURCE_DIR}/src)
SET(SRC_SHARE_DIR   ${HOMME_SOURCE_DIR}/src/share)
SET(PHYSICS_DIR     ${HOMME_SOURCE_DIR}/physics)
SET(UTILS_TIMING_DIR ${CMAKE_BINARY_DIR}/utils/timing)
SET(UTILS_TRILINOS_DIR ${HOMME_SOURCE_DIR}/utils/trilinos)

# Make INCLUDE_DIRS global so the tests can access it
SET (EXEC_INCLUDE_DIRS ${PIO_INCLUDE_DIRS}
                       ${UTILS_TIMING_DIR} )

FILE(GLOB SB81_SRC_F90  ${SRC_SHARE_DIR}/sb81/*.F90)

# Make SRCS global so the tests can access it
SET(PRIM_SRCS 
  ${SB81_SRC_F90}
  ${SRC_SHARE_DIR}/kinds.F90
  ${SRC_DIR}/test_src/asp_tests.F90
  ${SRC_DIR}/test_src/baroclinic_inst_mod.F90
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
  ${SRC_SHARE_DIR}/schedtype_mod.F90
  ${SRC_SHARE_DIR}/edgetype_mod.F90
  ${SRC_SHARE_DIR}/thread_mod.F90
  ${SRC_SHARE_DIR}/bndry_mod.F90
  ${SRC_SHARE_DIR}/derivative_mod.F90
  ${SRC_SHARE_DIR}/viscosity_mod.F90
  ${SRC_SHARE_DIR}/unit_tests_mod.F90
  ${SRC_SHARE_DIR}/spacecurve_mod.F90
  ${SRC_SHARE_DIR}/cube_mod.F90
  ${SRC_SHARE_DIR}/interpolate_mod.F90
  ${SRC_SHARE_DIR}/dof_mod.F90
  ${SRC_DIR}/common_io_mod.F90
  ${SRC_DIR}/pio_io_mod.F90 
  ${SRC_DIR}/repro_sum_mod.F90 
  ${SRC_SHARE_DIR}/time_mod.F90 
  ${SRC_SHARE_DIR}/physics_mod.F90 
  ${SRC_SHARE_DIR}/hybvcoord_mod.F90 
  ${SRC_DIR}/interp_movie_mod.F90 
  ${UTILS_SHARE_DIR}/shr_kind_mod.F90 
  ${UTILS_SHARE_DIR}/shr_mpi_mod.F90 
  ${UTILS_SHARE_DIR}/shr_sys_mod.F90 
  ${UTILS_SHARE_DIR}/shr_file_mod.F90 
  ${UTILS_SHARE_DIR}/shr_vmath_mod.F90 
  ${UTILS_SHARE_DIR}/shr_vmath_fwrap.c 
  ${UTILS_SHARE_DIR}/shr_spfn_mod.F90 
  ${UTILS_SHARE_DIR}/shr_const_mod.F90 
  ${SRC_SHARE_DIR}/domain_mod.F90 
  ${SRC_DIR}/restart_io_mod.F90 
  ${SRC_SHARE_DIR}/cg_mod.F90 
  ${SRC_SHARE_DIR}/linear_algebra_mod.F90 
  ${SRC_SHARE_DIR}/solver_mod.F90 
  ${SRC_DIR}/common_movie_mod.F90 
  ${SRC_SHARE_DIR}/namelist_mod.F90 
  ${SRC_SHARE_DIR}/mass_matrix_mod.F90 
  ${SRC_DIR}/checksum_mod.F90 
  ${SRC_SHARE_DIR}/prim_si_mod.F90 
  ${SRC_SHARE_DIR}/prim_state_mod.F90
  ${SRC_DIR}/prim_restart_mod.F90 
  ${SRC_SHARE_DIR}/prim_advance_mod.F90 
  ${SRC_SHARE_DIR}/prim_advection_mod.F90
  ${SRC_SHARE_DIR}/prim_implicit_mod.F90
  ${SRC_SHARE_DIR}/metis_mod.F90 
  ${SRC_SHARE_DIR}/prim_driver_base.F90 
  ${SRC_DIR}/prim_movie_mod.F90 
  ${SRC_DIR}/surfaces_mod.F90 
  ${SRC_DIR}/repro_sum_mod.F90 
  ${SRC_DIR}/interpolate_driver_mod.F90 
  ${SRC_DIR}/netcdf_io_mod.F90 
  ${SRC_SHARE_DIR}/mesh_mod.F90 
  ${SRC_SHARE_DIR}/prim_derived_type_mod.F90 
  ${UTILS_TRILINOS_DIR}/trilinosNoxSolver.cpp 
  ${UTILS_TRILINOS_DIR}/trilinosModelEvaluator.cpp 
  ${UTILS_TRILINOS_DIR}/block_precon_interface.cpp 
  ${UTILS_TRILINOS_DIR}/precon_interface.cpp 
  ${UTILS_TRILINOS_DIR}/belos_helm.cpp 
  ${SRC_DIR}/derived_type_mod.F90 
  ${PHYSICS_DIR}/physics_types_mod.F90
  ${PHYSICS_DIR}/physics_io_mod.F90 
  ${PHYSICS_DIR}/column_model_mod.F90 
  ${PHYSICS_DIR}/column_types_mod.F90 
  ${PHYSICS_DIR}/heldsuarez/held_suarez_mod.F90 
  ${SRC_DIR}/test_mod.F90
)


SET(EXEC_SOURCES ${PRIM_SRCS} ${SRC_DIR}/prim_main.F90 )

# Set up defaults
IF (NOT PRIM_NP)
  SET (PRIM_NP 4) 
ENDIF ()
IF (NOT PRIM_PLEV)
  SET (PRIM_PLEV 26) 
ENDIF ()
IF (NOT PRIM_USE_PIO)
  SET (PRIM_USE_PIO FALSE) 
ENDIF ()
IF (NOT PRIM_USE_ENERGY)
  SET (PRIM_USE_ENERGY FALSE) 
ENDIF ()
IF (NOT PRIM_NC)
  SET (PRIM_NC 4)
ENDIF ()
IF (NOT QSIZE_D)
  SET (QSIZE_D 4)
ENDIF ()

endmacro(prim_setup)

prim_setup()
# createTestExec(exec_name exec_type NP PLEV USE_PIO USE_ENERGY)
createTestExec(prim prim ${PRIM_NP} ${PRIM_NC} ${PRIM_PLEV} 
               ${PRIM_USE_PIO}  ${PRIM_USE_ENERGY} ${QSIZE_D})

SET_PROPERTY(TARGET prim PROPERTY LINKER_LANGUAGE Fortran)

TARGET_LINK_LIBRARIES(prim ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
