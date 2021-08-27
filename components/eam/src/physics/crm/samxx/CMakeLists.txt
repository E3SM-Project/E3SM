###############################################################################
## This assumes you have already run add_subdirectory(${YAKL_HOME} ${YAKL_BIN})
## It also assumes you've enabled CXX and C as languages (and CUDA if used)
###############################################################################

file(GLOB F90_SRC   cpp_interface_mod.F90 fft.F90 params.F90 fftpack5.F90 fftpack5_1d.F90)
file(GLOB CXX_SRC   *.cpp)
file(GLOB SAMXX_SRC ${CXX_SRC} ${F90_SRC})

# samxx will be a static library
add_library(samxx STATIC ${SAMXX_SRC})
if (${USE_CUDA})
  # samxx will be CUDA-linked with device symbols resolved at library creation
  set_target_properties(samxx PROPERTIES LINKER_LANGUAGE CUDA CUDA_SEPARABLE_COMPILATION OFF CUDA_RESOLVE_DEVICE_SYMBOLS ON)
endif()
# samxx needs to link with the yakl library
target_link_libraries(samxx yakl)
target_compile_features(samxx PUBLIC cxx_std_14)

# Set fortran compiler flags
set_source_files_properties(${F90_SRC} PROPERTIES COMPILE_FLAGS "${CPPDEFS} ${FFLAGS}")

# Set YAKL compiler flags
include(${YAKL_HOME}/yakl_utils.cmake)
yakl_process_cxx_source_files("${CXX_SRC}")
message(STATUS "YAKL Flags: ${YAKL_COMPILER_FLAGS}")

# Include YAKL source and library directories
include_directories(${YAKL_BIN})

