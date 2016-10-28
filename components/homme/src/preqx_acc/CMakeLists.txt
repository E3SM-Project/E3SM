# Relevant Directories within the HOMME repository
# holding source files needed for building the library
# and executable preqx.

macro(preqx_acc_setup)

SET(UTILS_SHARE_DIR ${HOMME_SOURCE_DIR}/utils/csm_share)
SET(SRC_DIR         ${HOMME_SOURCE_DIR}/src)
SET(SRC_SHARE_DIR   ${HOMME_SOURCE_DIR}/src/share)
SET(PHYSICS_DIR     ${HOMME_SOURCE_DIR}/physics)
SET(UTILS_TIMING_DIR ${CMAKE_BINARY_DIR}/utils/timing)

# Make INCLUDE_DIRS global so the tests can access it
SET (EXEC_INCLUDE_DIRS ${PIO_INCLUDE_DIRS}
                            ${UTILS_TIMING_DIR} )

SET(PREQX_SRCS_F90
  ${SRC_DIR}/test_src/asp_tests.F90
  ${SRC_DIR}/test_src/baroclinic_inst_mod.F90
  ${SRC_DIR}/test_src/dcmip2012_test1_2_3.F90
  ${SRC_DIR}/test_src/dcmip2012_test4.F90
  ${SRC_DIR}/test_src/dcmip2012_test5.F90
  ${SRC_SHARE_DIR}/bndry_mod_base.F90
  ${SRC_SHARE_DIR}/derivative_mod_base.F90
  ${SRC_SHARE_DIR}/edge_mod_base.F90
  ${SRC_SHARE_DIR}/prim_advection_mod_base.F90 
  ${SRC_SHARE_DIR}/solver_init_mod_base.F90
  ${SRC_SHARE_DIR}/vertremap_mod_base.F90
  ${SRC_SHARE_DIR}/viscosity_mod_base.F90
  ${SRC_DIR}/preqx/dcmip_tests.F90
  ${SRC_DIR}/preqx_acc/bndry_mod.F90
  ${SRC_DIR}/preqx_acc/derivative_mod.F90
  ${SRC_DIR}/preqx_acc/edge_mod.F90
  ${SRC_DIR}/preqx_acc/openacc_utils_mod.F90
  ${SRC_DIR}/preqx_acc/prim_advection_mod.F90 
  ${SRC_DIR}/preqx_acc/solver_init_mod.F90
  ${SRC_DIR}/preqx_acc/vertremap_mod.F90
  ${SRC_DIR}/preqx_acc/viscosity_mod.F90
  ${SRC_SHARE_DIR}/kinds.F90
  ${SRC_SHARE_DIR}/dimensions_mod.F90
  ${SRC_SHARE_DIR}/control_mod.F90
  ${SRC_SHARE_DIR}/physical_constants.F90
  ${SRC_SHARE_DIR}/coordinate_systems_mod.F90
  ${SRC_SHARE_DIR}/params_mod.F90
  ${SRC_SHARE_DIR}/gridgraph_mod.F90
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
  ${SRC_SHARE_DIR}/unit_tests_mod.F90
  ${SRC_SHARE_DIR}/spacecurve_mod.F90
  ${SRC_SHARE_DIR}/cube_mod.F90
  ${SRC_SHARE_DIR}/interpolate_mod.F90
  ${SRC_SHARE_DIR}/fvm_bsp_mod.F90
  ${SRC_SHARE_DIR}/fvm_transformation_mod.F90
  ${SRC_SHARE_DIR}/dof_mod.F90
  ${SRC_DIR}/common_io_mod.F90
  ${SRC_SHARE_DIR}/fvm_mod.F90
  ${SRC_SHARE_DIR}/fvm_control_volume_mod.F90
  ${SRC_SHARE_DIR}/fvm_line_integrals_mod.F90 
  ${SRC_SHARE_DIR}/fvm_line_integrals_flux_mod.F90 
  ${SRC_SHARE_DIR}/fvm_analytic_mod.F90 
  ${SRC_SHARE_DIR}/fvm_filter_mod.F90 
  ${SRC_SHARE_DIR}/fvm_reconstruction_mod.F90 
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
  ${UTILS_SHARE_DIR}/shr_const_mod.F90 
  ${UTILS_SHARE_DIR}/shr_spfn_mod.F90 
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
  ${SRC_SHARE_DIR}/prim_si_ref_mod.F90 
  ${SRC_SHARE_DIR}/prim_advance_mod.F90 
  ${SRC_DIR}/forcing_mod.F90 
  ${SRC_SHARE_DIR}/metis_mod.F90 
  ${SRC_SHARE_DIR}/prim_driver_mod.F90 
  ${SRC_DIR}/vertical_mod.F90 
  ${SRC_DIR}/prim_movie_mod.F90 
  ${SRC_DIR}/surfaces_mod.F90 
  ${SRC_DIR}/interpolate_driver_mod.F90 
  ${SRC_DIR}/netcdf_io_mod.F90 
  ${SRC_SHARE_DIR}/mesh_mod.F90 
  ${SRC_SHARE_DIR}/Manager.F90
  ${SRC_DIR}/prim_main.F90 
  ${SRC_DIR}/ref_state_mod.F90 
  ${PHYSICS_DIR}/physics_types_mod.F90
  ${PHYSICS_DIR}/physics_io_mod.F90 
  ${PHYSICS_DIR}/heldsuarez/held_suarez_mod.F90 
  ${PHYSICS_DIR}/column_model_mod.F90 
  ${PHYSICS_DIR}/column_types_mod.F90 
  ${SRC_DIR}/test_mod.F90
  )

SET(PREQX_SRCS_C
  ${UTILS_SHARE_DIR}/shr_vmath_fwrap.c
  ${SRC_DIR}/jrio.c
)

# If the user specified a file for custom compiler options use those
IF (DEFINED PREQX_CUSTOM_FLAGS_FILE)
  setCustomCompilerFlags(PREQX_CUSTOM_FLAGS_FILE PREQX_SRCS_F90)
ENDIF ()

# Make SRCS global so the tests can access it
SET(EXEC_SOURCES ${PREQX_ACC_SRCS} ${PREQX_SRCS_C} ${PREQX_SRCS_F90})

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
IF (NOT QSIZE_D)
  SET (QSIZE_D 4)
ENDIF ()


SET(USE_OPENACC TRUE)
SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OPENACC_Fortran_FLAGS}")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OPENACC_Linker_FLAGS}")

# nested openMP would never be used with openACC
SET(NESTED_OPENMP FALSE )
SET(COLUMN_OPENMP FALSE )

endmacro(preqx_acc_setup)




############################################################################
# createTestExec(exec_name exec_type NP PLEV USE_PIO USE_ENERGY)
############################################################################

preqx_acc_setup()
createTestExec(preqx_acc preqx_acc ${PREQX_NP} ${PREQX_NC} ${PREQX_PLEV} 
               ${PREQX_USE_PIO}  ${PREQX_USE_ENERGY} ${QSIZE_D})


