cmake_minimum_required (VERSION 2.8.12)
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

set (GPTL_Fortran_MODS ${CMAKE_CURRENT_BINARY_DIR}/perf_mod.mod
                       ${CMAKE_CURRENT_BINARY_DIR}/perf_utils.mod)

add_library (gptl ${GPTL_Fortran_SRCS} ${GPTL_C_SRCS})

target_include_directories (gptl
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
    PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

target_compile_definitions (gptl
    PUBLIC INCLUDE_CMAKE_FCI)
    
#==============================================================================
#  DEFINE THE INSTALL
#==============================================================================

# Library
install (TARGETS gptl DESTINATION lib)

# Header/Include File
install (FILES ${CMAKE_CURRENT_SOURCE_DIR}/gptl.h DESTINATION include)

# Fortran Modules
install (FILES ${GPTL_Fortran_MODS} DESTINATION include)

#==============================================================================
#  DEFINE THE DEPENDENCIES
#==============================================================================

#===== MPI =====
find_package (MPI REQUIRED)
if (MPI_C_FOUND AND MPI_Fortran_FOUND)
    target_compile_definitions (gptl
        PUBLIC HAVE_MPI)
else ()
    message (FATAL_ERROR "GPTL needs the C and Fortran MPI libraries")
endif ()
