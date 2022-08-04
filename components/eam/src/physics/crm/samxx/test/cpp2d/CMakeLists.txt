
add_executable(cpp2d ../dmdf.F90 ../cpp_driver.F90
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
target_link_libraries(cpp2d yakl ${NCFLAGS})
set_property(TARGET cpp2d APPEND PROPERTY COMPILE_FLAGS ${DEFS2D} )
set_property(TARGET cpp2d PROPERTY LINK_FLAGS "-Wl,--defsym,main=MAIN__  -lifcore")
set_property(TARGET cpp2d PROPERTY LINKER_LANGUAGE CXX)

include(${YAKL_HOME}/yakl_utils.cmake)
yakl_process_target(cpp2d)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/../yakl)

