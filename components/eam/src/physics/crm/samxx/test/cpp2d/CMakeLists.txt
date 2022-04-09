
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
target_link_libraries(cpp2d yakl ekat p3 shoc physics_share scream_share ${NCFLAGS})
set_property(TARGET cpp2d APPEND PROPERTY COMPILE_FLAGS ${DEFS2D} )
set_property(TARGET cpp2d PROPERTY LINKER_LANGUAGE Fortran)

include(${YAKL_HOME}/yakl_utils.cmake)
yakl_process_cxx_source_files("${CUDA_SRC}")

include_directories(${CMAKE_CURRENT_BINARY_DIR}/../yakl)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/../scream/src/physics/p3)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/../scream/src/physics/shoc)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/../scream/src/physics/share)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/../scream/src/share)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/../scream)

include_directories(${YAKL_HOME})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../scream/src/physics/p3)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../scream/src/physics/shoc)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../scream/src/physics/share)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../scream/src/share)




