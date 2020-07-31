###############################################################################
## This assumes you have already run add_subdirectory(${YAKL_HOME} ./yakl)
## It also assumes you've enabled CXX and C as languages
###############################################################################

file(GLOB SAMXX_SRC *.F90 *.cpp)
file(GLOB CUDA_SRC  *.cpp)

add_subdirectory(${YAKL_HOME} ./yakl  )
include_directories(${YAKL_HOME})

add_executable(cpp2d ../dmdf.F90 ../cpp_driver.F90 ../../fft.F90 ${CPP_SRC})
target_link_libraries(cpp2d yakl ${NCFLAGS})
set_property(TARGET cpp2d APPEND PROPERTY COMPILE_FLAGS ${DEFS2D} )

set_source_files_properties(${CUDA_SRC} PROPERTIES COMPILE_FLAGS "${YAKL_CXX_FLAGS}")
if ("${ARCH}" STREQUAL "CUDA")
  enable_language(CUDA)
  set_source_files_properties(${CUDA_SRC} PROPERTIES LANGUAGE CUDA)
set_source_files_properties(${CUDA_SRC} PROPERTIES COMPILE_FLAGS "${CUDA_FLAGS}")
  include_directories(${YAKL_CUB_HOME})
endif()
include_directories(${CMAKE_CURRENT_BINARY_DIR}/../yakl)

include_directories(${YAKL_HOME})

