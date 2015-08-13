#
# - Functions to check/validate the PnetCDF package(s) that were found
#

#
# - Check PnetCDF_C
#
function (check_PnetCDF_C)

    if (NOT PnetCDF_C_CHECKED)

        # Get version string
        try_run (PnetCDF_C_VERSION_RUNVAR PnetCDF_C_VERSION_COMPVAR
                 ${CMAKE_CURRENT_BINARY_DIR}/tryPnetCDF_C_VERSION
                 ${CMAKE_SOURCE_DIR}/cmake/TryPnetCDF_VERSION.c
                 COMPILE_DEFINITIONS -I${PnetCDF_C_INCLUDE_DIR}
                 COMPILE_OUTPUT_VARIABLE TryPnetCDF_OUT
                 RUN_OUTPUT_VARIABLE PnetCDF_C_VERSION)
        if (PnetCDF_C_VERSION)
            if (PnetCDF_C_VERSION VERSION_LESS PnetCDF_FIND_VERSION)
                message (FATAL_ERROR "PnetCDF_C version insufficient")
            else ()
                message (STATUS "Found PnetCDF_C version ${PnetCDF_C_VERSION}")
            endif ()
        else ()
            message (STATUS "Could not find PnetCDF_C version")
        endif ()
        set (PnetCDF_C_VERSION ${PnetCDF_C_VERSION} PARENT_SCOPE)

        set (PnetCDF_C_CHECKED TRUE CACHE BOOL "PnetCDF_C checked")
        
    endif ()
        
endfunction ()

#
# - Check NetCDF_Fortran
#
function (check_PnetCDF_Fortran)

    if (NOT PnetCDF_Fortran_CHECKED)

        # Get version string
        set (COMP_DEFS)
        foreach (incdir IN LISTS PnetCDF_Fortran_INCLUDE_DIRS)
            list (APPEND COMP_DEFS "-I${incdir}")
        endforeach ()
        try_run (PnetCDF_Fortran_VERSION_RUNVAR PnetCDF_Fortran_VERSION_COMPVAR
                 ${CMAKE_CURRENT_BINARY_DIR}/tryPnetCDF_Fortran_VERSION
                 ${CMAKE_SOURCE_DIR}/cmake/TryPnetCDF_VERSION.f90
                 COMPILE_DEFINITIONS ${COMP_DEFS}
                 LINK_LIBRARIES ${PnetCDF_Fortran_LIBRARIES}
                 COMPILE_OUTPUT_VARIABLE TryPnetCDF_OUT
                 RUN_OUTPUT_VARIABLE PnetCDF_Fortran_VERSION)
        if (PnetCDF_Fortran_VERSION)
            string (STRIP ${PnetCDF_Fortran_VERSION} PnetCDF_Fortran_VERSION)
            if (PnetCDF_Fortran_VERSION VERSION_LESS PnetCDF_FIND_VERSION)
                message (FATAL_ERROR "PnetCDF_Fortran version insufficient")
            else ()
                message (STATUS "Found PnetCDF_Fortran version ${PnetCDF_Fortran_VERSION}")
            endif ()
        else ()
            message (STATUS "Could not find PnetCDF_Fortran version")
        endif ()
        set (PnetCDF_Fortran_VERSION ${PnetCDF_Fortran_VERSION} PARENT_SCOPE)

        set (PnetCDF_Fortran_CHECKED TRUE CACHE BOOL "PnetCDF_Fortran checked")
    
    endif ()
    
endfunction ()