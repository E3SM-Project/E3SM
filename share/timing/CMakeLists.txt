cmake_minimum_required(VERSION 3.18)

# Needs to happen before project call
function(set_compilers_e3sm)
  # Grab CXX compiler from CIME
  include(${CASEROOT}/Macros.cmake)

  set(CMAKE_CXX_COMPILER     ${CMAKE_CXX_COMPILER}     CACHE STRING "The CXX compiler")
  set(CMAKE_C_COMPILER       ${CMAKE_C_COMPILER}       CACHE STRING "The C   compiler")
  set(CMAKE_Fortran_COMPILER ${CMAKE_Fortran_COMPILER} CACHE STRING "The Fortran compiler")

endfunction()
set_compilers_e3sm()

include(${CASEROOT}/Macros.cmake)

project(GPTL C Fortran)

set(CMAKE_VERBOSE_MAKEFILE TRUE)

include(${CASEROOT}/Macros.cmake)

include(FortranCInterface)
FortranCInterface_HEADER(cmake_fortran_c_interface.h
  MACRO_NAMESPACE  "FCI_")

# Just assume all files are free
set(CMAKE_Fortran_FORMAT FREE)

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

set(CPPDEFS "${CPPDEFS} ${TIMING_CPPDEFS}")
separate_arguments(CPPDEFS_LIST UNIX_COMMAND "${CPPDEFS}")

add_library(gptl ${SRCS_F90} ${SRCS_C})
target_include_directories(gptl PRIVATE ${TIMING_INCLUDE_DIRS})
target_compile_definitions(gptl PRIVATE ${CPPDEFS_LIST})
set_target_properties(gptl PROPERTIES PUBLIC_HEADER "gptl.h;${CMAKE_BINARY_DIR}/perf_mod.mod;${CMAKE_BINARY_DIR}/perf_utils.mod;")

install(TARGETS gptl
  PUBLIC_HEADER DESTINATION include)
