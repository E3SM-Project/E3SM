cmake_minimum_required (VERSION 3.0.0)
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

target_include_directories (gptl
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
    PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

target_compile_definitions (gptl
    PRIVATE INCLUDE_CMAKE_FCI)

#==============================================================================
#  DEFINE THE DEPENDENCIES
#==============================================================================

#===== MPI =====
find_package (MPI)
if (MPI_C_FOUND AND MPI_Fortran_FOUND)
    target_compile_definitions (gptl 
        PRIVATE HAVE_MPI)
    target_include_directories (gptl
        PRIVATE ${MPI_C_INCLUDE_PATH}
        PRIVATE ${MPI_Fortran_INCLUDE_PATH})
    target_link_libraries (gptl
        ${MPI_C_LIBRARIES}
        ${MPI_Fortran_LIBRARIES})
else ()
    message (FATAL_ERROR "GPTL needs the C and Fortran MPI libraries")
endif ()
