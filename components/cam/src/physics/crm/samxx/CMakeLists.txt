###############################################################################
## This assumes you have already run add_subdirectory(${YAKL_HOME} ${YAKL_BIN})
## It also assumes you've enabled CXX and C as languages
###############################################################################

file(GLOB F90_SRC   *.F90)
file(GLOB CXX_SRC   *.cpp)
file(GLOB SAMXX_SRC *.cpp *.F90)

add_library(samxx STATIC ${SAMXX_SRC})
set_target_properties(samxx PROPERTIES LINKER_LANGUAGE CUDA CUDA_SEPARABLE_COMPILATION OFF CUDA_RESOLVE_DEVICE_SYMBOLS ON)
target_link_libraries(samxx yakl)

set_source_files_properties(${F90_SRC} PROPERTIES COMPILE_FLAGS "${CPPDEFS} ${FFLAGS}")
if ("${ARCH}" STREQUAL "CUDA")
  set_source_files_properties(${CXX_SRC} PROPERTIES LANGUAGE CUDA)
  set_source_files_properties(${CXX_SRC} PROPERTIES COMPILE_FLAGS "${CPPDEFS} ${YAKL_CXX_FLAGS}")
  include_directories(${YAKL_CUB_HOME})
endif()
include_directories(${YAKL_BIN})
include_directories(${YAKL_HOME})

