cmake_minimum_required (VERSION 2.8.12)
project (GPTL C Fortran)
include (CheckFunctionExists)
include (FortranCInterface)
FortranCInterface_HEADER(cmake_fortran_c_interface.h
                         MACRO_NAMESPACE "FCI_")
include (LibCheck)

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

#===== System Checks =====
if (CMAKE_SYSTEM_NAME MATCHES "AIX")
    target_compile_definitions (gptl
        PUBLIC _AIX)
else ()
    target_compile_definitions (gptl
        PUBLIC ${CMAKE_SYSTEM_NAME_CAPS})
endif ()

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

#===== PAPI =====
if (ENABLE_PAPI)
    find_package (PAPI)
    if (PAPI_FOUND)
        target_compile_definitions (gptl
            PUBLIC HAVE_PAPI)
        target_include_directories (gptl
            PUBLIC ${PAPI_INCLUDE_DIRECTORIES})
        target_link_libraries (gptl
            PUBLIC ${PAPI_LIBRARIES})
    endif ()
endif ()

#===== LIBRT =====
if (ENABLE_LIBRT)
    find_package (LIBRT)
    if (LIBRT_FOUND)
        target_compile_definitions (gptl
            PUBLIC HAVE_LIBRT)
        target_include_directories (gptl
            PUBLIC ${LIBRT_INCLUDE_DIRECTORIES})
        target_link_libraries (gptl
            PUBLIC ${LIBRT_LIBRARIES})
    endif () 
endif ()

#===== MPI =====
find_package (MPI REQUIRED)
if (MPI_C_FOUND AND MPI_Fortran_FOUND)
    target_compile_definitions (gptl
        PUBLIC HAVE_MPI)
    target_include_directories (gptl
        PUBLIC ${MPI_C_INCLUDE_PATH}
        PUBLIC ${MPI_Fortran_INCLUDE_PATH})
    target_link_libraries (gptl
        PUBLIC ${MPI_C_LIBRARIES}
        PUBLIC ${MPI_Fortran_LIBRARIES})
        
    # Check library for varn functions
    set (CMAKE_REQUIRED_LIBRARIES ${MPI_C_LIBRARIES})
    check_function_exists (MPI_Comm_f2c MPI_COMM_F2C)
    if (NOT MPI_COMM_F2C)
        target_compile_definitions(gptl
            PUBLIC NO_COMM_F2C)
    endif()

    # Check for MPI Fortran module
    find_path(MPIMOD_PATH
              NAMES mpi.mod
              HINTS ${MPI_Fortran_INCLUDE_PATH})
    if (MPIMOD_PATH)
        check_macro (MPIMOD_WORKS
                     NAME TryMPIMod.f90
                     HINTS ${CMAKE_MODULE_PATH}
                     DEFINITIONS -I${MPIMOD_PATH}
                     COMMENT "whether MPI Fortran module is supported")
        if (${MPIMOD_WORKS})
            message (STATUS "Working MPI Fortran module detected and enabled.")
        else ()
            message (STATUS "Working MPI Fortran module not detected and therefore disabled.")
            target_compile_definitions (gptl
                PUBLIC NO_MPIMOD)
        endif ()
    else ()
        message (STATUS "MPI Fortran module not detected and therefore disabled.")
        target_compile_definitions (gptl
            PUBLIC NO_MPIMOD)
    endif ()

else ()
    message (FATAL_ERROR "GPTL needs the C and Fortran MPI libraries")
endif ()

#===== GetTimeOfDay =====
get_target_property (GPTL_LINK_LIBRARIES gptl LINK_LIBRARIES)
set (CMAKE_REQUIRED_LIBRARIES ${GPTL_LINK_LIBRARIES})
check_function_exists (gettimeofday GETTIMEOFDAY)
if (NOT GETTIMEOFDAY)
    target_compile_definitions(gptl
        PUBLIC NO_GETTIMEOFDAY)
endif ()
