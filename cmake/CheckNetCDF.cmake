#
# - Functions to check/validate the NetCDF package(s) that were found
#

#
# - Check NetCDF_C
#
function (check_NetCDF_C)

    if (NOT NetCDF_C_CHECKED)
    
        find_path (NetCDF_C_META_DIR
                   NAMES netcdf_meta.h
                   HINTS ${NetCDF_C_INCLUDE_DIRS})
        if (NetCDF_C_META_DIR)
        
            # Get version string
            try_run (NetCDF_C_VERSION_RUNVAR NetCDF_C_VERSION_COMPVAR
                     ${CMAKE_CURRENT_BINARY_DIR}/tryNetCDF_C_VERSION
                     ${CMAKE_SOURCE_DIR}/cmake/TryNetCDF_VERSION.c
                     COMPILE_DEFINITIONS -I${NetCDF_C_META_DIR}
                     COMPILE_OUTPUT_VARIABLE TryNetCDF_OUT
                     RUN_OUTPUT_VARIABLE NetCDF_C_VERSION)
            if (NetCDF_C_VERSION)
                if (NetCDF_C_VERSION VERSION_LESS NetCDF_FIND_VERSION)
                    message (FATAL_ERROR "NetCDF_C version insufficient")
                else ()
                    message (STATUS "Found NetCDF_C version ${NetCDF_C_VERSION}")
                endif ()
            endif ()
            set (NetCDF_C_VERSION ${NetCDF_C_VERSION} PARENT_SCOPE)
            
            # Test for DAP support (requires CURL)
            try_compile(NetCDF_C_HAS_DAP 
                        ${CMAKE_CURRENT_BINARY_DIR}/tryNetCDF_DAP
                        SOURCES ${CMAKE_SOURCE_DIR}/cmake/TryNetCDF_DAP.c
                        COMPILE_DEFINITIONS -I${NetCDF_C_META_DIR}
                        OUTPUT_VARIABLE TryNetCDF_OUT)
            if (NetCDF_C_HAS_DAP)
                message (STATUS "NetCDF_C has DAP support")
            else ()
                message (STATUS "NetCDF_C does not have DAP support")
            endif ()
            set (NetCDF_C_HAS_DAP ${NetCDF_C_HAS_DAP} PARENT_SCOPE)
    
            # Test for PARALLEL support
            try_compile(NetCDF_C_HAS_PARALLEL 
                        ${CMAKE_CURRENT_BINARY_DIR}/tryNetCDF_PARALLEL
                        SOURCES ${CMAKE_SOURCE_DIR}/cmake/TryNetCDF_PARALLEL.c
                        COMPILE_DEFINITIONS -I${NetCDF_C_META_DIR}
                        OUTPUT_VARIABLE TryNetCDF_OUT)
            if (NetCDF_C_HAS_PARALLEL)
                message (STATUS "NetCDF_C has parallel support")
            else ()
                message (STATUS "NetCDF_C does not have parallel support")
            endif ()
            set (NetCDF_C_HAS_PARALLEL ${NetCDF_C_HAS_PARALLEL} PARENT_SCOPE)
                
            # Test for PNETCDF support
            try_compile(NetCDF_C_HAS_PNETCDF
                        ${CMAKE_CURRENT_BINARY_DIR}/tryNetCDF_PNETCDF
                        SOURCES ${CMAKE_SOURCE_DIR}/cmake/TryNetCDF_PNETCDF.c
                        COMPILE_DEFINITIONS -I${NetCDF_C_META_DIR}
                        OUTPUT_VARIABLE TryNetCDF_OUT)
            if (NetCDF_C_HAS_PNETCDF)
                message (STATUS "NetCDF_C requires PnetCDF")
                find_package(PnetCDF COMPONENTS C)
                if (PnetCDF_C_FOUND)
                    set (temp ${NetCDF_C_INCLUDE_DIRS} ${PnetCDF_C_INCLUDE_DIRS})
                    set (NetCDF_C_INCLUDE_DIRS ${temp} PARENT_SCOPE)
                    set (temp ${NetCDF_C_LIBRARIES} ${PnetCDF_C_LIBRARIES})
                    set (NetCDF_C_LIBRARIES ${temp} PARENT_SCOPE)
                endif ()
            else ()
                message (STATUS "NetCDF_C does not require PnetCDF")
            endif ()
            set (NetCDF_C_HAS_PNETCDF ${NetCDF_C_HAS_PNETCDF} PARENT_SCOPE)
            
            
        else ()
            message (WARNING "Could not find netcdf_meta.h")
        endif ()
        
        set (NetCDF_C_CHECKED TRUE CACHE BOOL "NetCDF_C checked")
        
    endif ()
        
endfunction ()

#
# - Check NetCDF_Fortran
#
function (check_NetCDF_Fortran)

    if (NOT NetCDF_Fortran_CHECKED)
    
        # Check NetCDF_C if it needs PnetCDF
        if (NetCDF_C_HAS_PNETCDF)
            find_package(PnetCDF COMPONENTS C)
            if (PnetCDF_C_FOUND)
                set (temp ${NetCDF_Fortran_INCLUDE_DIRS} ${PnetCDF_C_INCLUDE_DIRS})
                set (NetCDF_Fortran_INCLUDE_DIRS ${temp} PARENT_SCOPE)
                set (temp ${NetCDF_Fortran_LIBRARIES} ${PnetCDF_C_LIBRARIES})
                set (NetCDF_Fortran_LIBRARIES ${temp} PARENT_SCOPE)
            endif ()
        endif ()
    
        # Get version string
        set (COMP_DEFS)
        foreach (incdir IN LISTS NetCDF_Fortran_INCLUDE_DIRS)
            list (APPEND COMP_DEFS "-I${incdir}")
        endforeach ()
        message ("COMP_DEFS = ${COMP_DEFS}")
        message ("NetCDF_Fortran_LIBRARIES = ${NetCDF_Fortran_LIBRARIES}")
        try_run (NetCDF_Fortran_VERSION_RUNVAR
                 NetCDF_Fortran_VERSION_COMPVAR
                 ${CMAKE_CURRENT_BINARY_DIR}/tryNetCDF_Fortran_VERSION
                 ${CMAKE_SOURCE_DIR}/cmake/TryNetCDF_VERSION.f90
                 COMPILE_DEFINITIONS ${COMP_DEFS}
                 LINK_LIBRARIES ${NetCDF_Fortran_LIBRARIES}
                 COMPILE_OUTPUT_VARIABLE TryNetCDF_OUT
                 RUN_OUTPUT_VARIABLE NetCDF_Fortran_VERSION)
        if (NetCDF_Fortran_VERSION)
            string (STRIP ${NetCDF_Fortran_VERSION} NetCDF_Fortran_VERSION)
            if (NetCDF_Fortran_VERSION VERSION_LESS NetCDF_FIND_VERSION)
                message (FATAL_ERROR "NetCDF_Fortan version insufficient")
            else ()
                message (STATUS "Found NetCDF_Fortran version ${NetCDF_Fortran_VERSION}")
            endif ()
        else ()
            message (STATUS "Could not find NetCDF_Fortran version")
        endif ()
        set (NetCDF_Fortran_VERSION ${NetCDF_Fortran_VERSION} PARENT_SCOPE)

        set (NetCDF_Fortran_CHECKED TRUE CACHE BOOL "NetCDF_Fortran checked")
    
    endif ()
    
endfunction ()