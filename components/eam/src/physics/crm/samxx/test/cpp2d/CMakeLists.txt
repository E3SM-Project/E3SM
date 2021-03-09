
add_executable(cpp2d ../dmdf.F90 ../cpp_driver.F90
               ../../../crmdims.F90
               ../../../params_kind.F90
               ../../../crm_input_module.F90
               ../../../crm_output_module.F90
               ../../../crm_rad_module.F90
               ../../../crm_state_module.F90
               ../../../openacc_utils.F90
               ${CPP_SRC})
target_link_libraries(cpp2d yakl ${NCFLAGS})
set_property(TARGET cpp2d APPEND PROPERTY COMPILE_FLAGS ${DEFS2D} )

set_source_files_properties(${CUDA_SRC} PROPERTIES COMPILE_FLAGS "${YAKL_CXX_FLAGS}")
if ("${ARCH}" STREQUAL "CUDA")
  set_source_files_properties(${CUDA_SRC} PROPERTIES LANGUAGE CUDA)
  set_source_files_properties(${CUDA_SRC} PROPERTIES COMPILE_FLAGS "${CUDA_FLAGS}")
  include_directories(${YAKL_CUB_HOME})
endif()
include_directories(${CMAKE_CURRENT_BINARY_DIR}/../yakl)

include_directories(${YAKL_HOME})

