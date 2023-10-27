cmake_minimum_required(VERSION 3.18)
project(GPTL C Fortran)

include(FortranCInterface)
FortranCInterface_HEADER(cmake_fortran_c_interface.h
  MACRO_NAMESPACE  "FCI_")

include(${CASEROOT}/Macros.cmake)

set(SRCS_C
  GPTLget_memusage.c
  GPTLprint_memusage.c
  GPTLutil.c
  f_wrappers.c
  gptl.c
  gptl_papi.c)

set(SRCS_F90
  perf_mod.F90
  perf_utils.F90)

set(TIMING_INCLUDE_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${GPTL_DIR}
  ${INSTALL_SHAREDPATH}/include)

set(TIMING_CPPDEFS "-DINCLUDE_CMAKE_FCI")
if (NOT MPILIB STREQUAL "mpi-serial")
  string(APPEND TIMING_CPPDEFS " -DHAVE_MPI")
endif()
if (compile_threaded)
  string(APPEND TIMING_CPPDEFS " -DTHREADED_OMP")
endif()
string(TOUPPER COMP_INTERFACE UPVAR)
string(APPEND TIMING_CPPDEFS " -D${UPVAR}_INTERFACE")
if (DEBUG)
  string(APPEND TIMING_CPPDEFS " -DDEBUG")
endif()

if (CPRE)
  # Not sure what to do here or if this is really needed
endif()

if (FREEFLAGS)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FREEFLAGS}")
endif()

add_library(timing ${SRCS_F90} ${SRCS_C})
target_include_directories(timing PRIVATE ${TIMING_INCLUDE_DIRS})
target_compile_definitions(timing PRIVATE ${TIMING_CPPDEFS})
set_target_properties(timing PROPERTIES PUBLIC_HEADER "gptl.h;${CMAKE_BINARY_DIR}/perf_mod.mod;${CMAKE_BINARY_DIR}/perf_utils.mod;")

install(TARGETS timing
  PUBLIC_HEADER DESTINATION include)
