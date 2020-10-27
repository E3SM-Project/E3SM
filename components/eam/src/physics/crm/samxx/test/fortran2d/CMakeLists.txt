
add_executable(fortran2d ../shr_const_mod.F90 ../perf_mod.F90 ../fortran_driver.F90
               ../../../crmdims.F90
               ../../../params_kind.F90
               ../../../RNG_MT.F90
               ../../../crm_input_module.F90
               ../../../crm_output_module.F90
               ../../../crm_rad_module.F90
               ../../../crm_state_module.F90
               ../../../openacc_utils.F90
               ${FORTRAN_SRC})
target_link_libraries(fortran2d ${NCFLAGS})
set_property(TARGET fortran2d APPEND PROPERTY COMPILE_FLAGS ${DEFS2D} )

