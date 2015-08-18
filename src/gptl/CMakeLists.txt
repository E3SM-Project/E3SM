cmake_minimum_required (VERSION 2.8.10)
project (GPTL C Fortran)
include (FortranCInterface)
FortranCInterface_HEADER(cmake_fortran_c_interface.h
                         MACRO_NAMESPACE "FCI_")

#==============================================================================
#  DEFINE THE TARGET
#==============================================================================

set (GPTL_C_SRCS GPTLget_memusage.c
                 GPTLprint_memusage.c
                 GPTLutil.c
                 f_wrappers.c
                 gptl.c
                 gptl_papi.c
                 threadutil.c)

set (GPTL_Fortran_SRCS perf_mod.F90
                       perf_utils.F90)

add_library (gptl EXCLUDE_FROM_ALL ${GPTL_Fortran_SRCS} ${GPTL_C_SRCS})

set (GPTL_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} 
                       ${CMAKE_CURRENT_BINARY_DIR} 
     CACHE STRING "GPTL include directories")
set (GPTL_DEFINITIONS -DINCLUDE_CMAKE_FCI
     CACHE STRING "GPTL compiler definitions")

#==============================================================================
#  DEFINE THE DEPENDENCIES
#==============================================================================

#===== MPI =====
find_package (MPI REQUIRED)
if (MPI_C_FOUND AND MPI_Fortran_FOUND)
    list (APPEND GPTL_DEFINITIONS -DHAVE_MPI)
else ()
    message (FATAL_ERROR "GPTL needs the C and Fortran MPI libraries")
endif ()

#==============================================================================
#  ADD INCLUDES AND DEFINITIONS
#==============================================================================

include_directories(${GPTL_INCLUDE_DIRS})
add_definitions (${GPTL_DEFINITIONS})
