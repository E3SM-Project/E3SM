#
# - Functions to check/validate the NetCDF package(s) that were found
#   These functions should be used from withing the FindNetCDF module.
#

#
# - Check NetCDF Version
#
function (check_NetCDF_VERSION META_HINTS)

    if (NOT DEFINED NetCDF_VERSION)
    
        find_path (NetCDF_META_DIR
                   NAMES netcdf_meta.h
                   HINTS ${META_HINTS})
        if (NetCDF_META_DIR)
        
            # Get version string
            try_run (RUN_RESULT COMPILE_RESULT
                     ${CMAKE_CURRENT_BINARY_DIR}/tryNetCDF_VERSION
                     ${CMAKE_SOURCE_DIR}/cmake/TryNetCDF_VERSION.c
                     COMPILE_DEFINITIONS -I${NetCDF_META_DIR}
                     COMPILE_OUTPUT_VARIABLE TryNetCDF_OUT
                     RUN_OUTPUT_VARIABLE NetCDF_VERSION)
            if (COMPILE_RESULT AND NOT RUN_RESULT STREQUAL FAILED_TO_RUN)
                if (NetCDF_VERSION VERSION_LESS NetCDF_FIND_VERSION)
                    message (FATAL_ERROR "NetCDF version insufficient")
                else ()
                    message (STATUS "Found NetCDF version ${NetCDF_VERSION}")
                endif ()
                set (NetCDF_VERSION ${NetCDF_VERSION}
                     CACHE STRING "NetCDF version string")
            else ()
                message (STATUS "NetCDF version could not be found")
            endif ()
                        
        else ()
            unset (NetCDF_META_DIR CACHE)
        endif ()
        
    endif ()
        
endfunction ()

#
# - Check if NetCDF has DAP support
#
function (check_NetCDF_HAS_DAP META_HINTS)

    if (NOT DEFINED NetCDF_HAS_DAP)
    
        find_path (NetCDF_META_DIR
                   NAMES netcdf_meta.h
                   HINTS ${META_HINTS})
        if (NetCDF_META_DIR)

            # Test for DAP support (requires CURL)
            try_compile(COMPILE_RESULT 
                        ${CMAKE_CURRENT_BINARY_DIR}/tryNetCDF_DAP
                        SOURCES ${CMAKE_SOURCE_DIR}/cmake/TryNetCDF_DAP.c
                        COMPILE_DEFINITIONS -I${NetCDF_META_DIR}
                        OUTPUT_VARIABLE TryNetCDF_OUT)
            if (COMPILE_RESULT)
                message (STATUS "NetCDF has DAP support")
            else ()
                message (STATUS "NetCDF does not have DAP support")
            endif ()            
            set (NetCDF_HAS_DAP ${COMPILE_RESULT}
                 CACHE BOOL "Whether NetCDF has DAP support enabled")

        else ()
            unset (NetCDF_META_DIR CACHE)
        endif ()
    
    endif ()
    
endfunction ()

#
# - Check if NetCDF has parallel support
#
function (check_NetCDF_HAS_PARALLEL META_HINTS)

    if (NOT DEFINED NetCDF_HAS_PARALLEL)
    
        find_path (NetCDF_META_DIR
                   NAMES netcdf_meta.h
                   HINTS ${META_HINTS})
        if (NetCDF_META_DIR)

            # Test for PARALLEL support
            try_compile(COMPILE_RESULT 
                        ${CMAKE_CURRENT_BINARY_DIR}/tryNetCDF_PARALLEL
                        SOURCES ${CMAKE_SOURCE_DIR}/cmake/TryNetCDF_PARALLEL.c
                        COMPILE_DEFINITIONS -I${NetCDF_META_DIR}
                        OUTPUT_VARIABLE TryNetCDF_OUT)
            if (COMPILE_RESULT)
                message (STATUS "NetCDF has parallel support enabled")
            else ()
                message (STATUS "NetCDF has parallel support disabled")
            endif ()
            set (NetCDF_HAS_PARALLEL ${COMPILE_RESULT}
                 CACHE BOOL "Whether NetCDF has parallel support enabled")

        else ()
            unset (NetCDF_META_DIR CACHE)
        endif ()
    
    endif ()
    
endfunction ()

#
# - Check if NetCDF has PnetCDF support
#
function (check_NetCDF_HAS_PNETCDF META_HINTS)

    if (NOT DEFINED NetCDF_HAS_PNETCDF)
    
        find_path (NetCDF_META_DIR
                   NAMES netcdf_meta.h
                   HINTS ${META_HINTS})
        if (NetCDF_META_DIR)

            # Test for PNETCDF support
            try_compile(COMPILE_RESULT
                        ${CMAKE_CURRENT_BINARY_DIR}/tryNetCDF_PNETCDF
                        SOURCES ${CMAKE_SOURCE_DIR}/cmake/TryNetCDF_PNETCDF.c
                        COMPILE_DEFINITIONS -I${NetCDF_C_META_DIR}
                        OUTPUT_VARIABLE TryNetCDF_OUT)
            if (COMPILE_RESULT)
                message (STATUS "NetCDF requires PnetCDF")
            else ()
                message (STATUS "NetCDF does not require PnetCDF")
            endif ()
            
            set (NetCDF_HAS_PNETCDF ${COMPILE_RESULT}
                 CACHE BOOL "Whether NetCDF was build with PnetCDF support")

        else ()
            unset (NetCDF_META_DIR CACHE)
        endif ()
    
    endif ()
    
endfunction ()
