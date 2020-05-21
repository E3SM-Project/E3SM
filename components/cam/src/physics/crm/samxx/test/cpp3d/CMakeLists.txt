
add_executable(cpp3d ../dmdf.F90 ../cpp_driver.F90 ${CPP_SRC})
target_link_libraries(cpp3d yakl ${NCFLAGS})
set_property(TARGET cpp3d APPEND PROPERTY COMPILE_FLAGS ${DEFS3D} )

set_source_files_properties(${CUDA_SRC} PROPERTIES COMPILE_FLAGS "${YAKL_CXX_FLAGS}")
if ("${ARCH}" STREQUAL "CUDA")
  enable_language(CUDA)
  set_source_files_properties(${CUDA_SRC} PROPERTIES LANGUAGE CUDA)
  include_directories(${YAKL_CUB_HOME})
endif()

include_directories(${YAKL_HOME})

