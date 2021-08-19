
add_executable(cpp3d ../dmdf.F90 ../cpp_driver.F90
               ../../../crmdims.F90
               ../../../params_kind.F90
               ../../../crm_input_module.F90
               ../../../crm_output_module.F90
               ../../../crm_rad_module.F90
               ../../../crm_state_module.F90
               ../../../crm_ecpp_output_module.F90
               ../../../ecppvars.F90
               ../../../openacc_utils.F90
               ${CPP_SRC})
target_link_libraries(cpp3d yakl ${NCFLAGS})
set_property(TARGET cpp3d APPEND PROPERTY COMPILE_FLAGS ${DEFS3D} )

include(${YAKL_HOME}/process_cxx_source_files.cmake)
process_cxx_source_files("${CUDA_SRC}")
include_directories(${YAKL_BIN})

