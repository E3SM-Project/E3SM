
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
set_property(TARGET cpp3d PROPERTY LINK_FLAGS "-Wl,--defsym,main=MAIN__  -lifcore")
set_property(TARGET cpp3d PROPERTY LINKER_LANGUAGE CXX)

include(${YAKL_HOME}/yakl_utils.cmake)
yakl_process_target(cpp3d)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/../yakl)

