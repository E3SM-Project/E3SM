#
# - Functions to check/validate the PnetCDF package(s) that were found
#   These functions should be run from the FindPnetCDF module

#
# - Check PnetCDF Version
#
function (check_PnetCDF_VERSION INC_DIR)

    if (NOT DEFINED PnetCDF_VERSION)
    
        # Get version string
        try_run (RUN_RESULT COMPILE_RESULT
                 ${CMAKE_CURRENT_BINARY_DIR}/tryPnetCDF_VERSION
                 ${CMAKE_SOURCE_DIR}/cmake/TryPnetCDF_VERSION.c
                 COMPILE_DEFINITIONS -I${INC_DIR}
                 COMPILE_OUTPUT_VARIABLE TryPnetCDF_OUT
                 RUN_OUTPUT_VARIABLE PnetCDF_VERSION)
        if (COMPILE_RESULT AND NOT RUN_RESULT STREQUAL FAILED_TO_RUN)
            if (PnetCDF_VERSION VERSION_LESS PnetCDF_FIND_VERSION)
                message (FATAL_ERROR "PnetCDF version insufficient")
            else ()
                message (STATUS "Found PnetCDF version ${PnetCDF_VERSION}")
            endif ()
            set (PnetCDF_VERSION ${PnetCDF_VERSION}
                 CACHE STRING "PnetCDF version string")
        else ()
            message (STATUS "PnetCDF version could not be found")
        endif ()
        
    endif ()
        
endfunction ()
