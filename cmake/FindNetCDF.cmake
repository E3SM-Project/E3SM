# - Try to find NetCDF
#
# This can be controlled by setting the NetCDF_DIR (or, equivalently, the 
# NETCDF environment variable), or NetCDF_<lang>_DIR CMake variables, where
# <lang> is the COMPONENT language one needs.
#
# Once done, this will define:
#
#   NetCDF_<lang>_FOUND        (BOOL) - system has NetCDF
#   NetCDF_<lang>_IS_SHARED    (BOOL) - whether library is shared/dynamic
#   NetCDF_<lang>_INCLUDE_DIR  (PATH) - Location of the C header file
#   NetCDF_<lang>_INCLUDE_DIRS (LIST) - the NetCDF include directories
#   NetCDF_<lang>_LIBRARY      (FILE) - Path to the C library file
#   NetCDF_<lang>_LIBRARIES    (LIST) - link these to use NetCDF
#   NetCDF_<lang>_DEFINITIONS  (LIST) - preprocessor macros to use with NetCDF
#   NetCDF_<lang>_OPTIONS      (LIST) - compiler options to use NetCDF
#
# The available COMPONENTS are: C Fortran
# If no components are specified, it assumes only C
include (LibFindLibraryMacros)

# Define NetCDF C Component
define_package_component (NetCDF DEFAULT
                          COMPONENT C
                          INCLUDE_NAMES netcdf.h
                          LIBRARY_NAMES netcdf)

# Define NetCDF Fortran Component
define_package_component (NetCDF
                          COMPONENT Fortran
                          INCLUDE_NAMES netcdf.mod netcdf.inc
                          LIBRARY_NAMES netcdff)
                       
# Search for list of valid components requested
find_valid_components (NetCDF)

# SEARCH FOR VALIDATED COMPONENTS
foreach (comp IN LISTS NetCDF_FIND_VALID_COMPONENTS)

    # If not found already, search...
    if (NOT NetCDF_${comp}_FOUND)

        # Manually add the MPI include and library dirs to search paths
        if (MPI_${comp}_FOUND)
            set (NetCDF_${comp}_INCLUDE_HINTS ${MPI_${comp}_INCLUDE_PATH})
            set (NetCDF_${comp}_LIBRARY_HINTS)
            foreach (lib IN LISTS ${MPI_${comp}_LIBRARIES})
                get_filename_component (libdir ${lib} PATH)
                list (APPEND NetCDF_${comp}_LIBRARY_HINTS ${libdir})
                unset (libdir)
            endforeach ()
        endif ()
        
        message ("NetCDF_${comp}_LIBRARY_HINTS = ${NetCDF_${comp}_LIBRARY_HINTS}")

        # Search for the package component    
        find_package_component(NetCDF COMPONENT ${comp}
                               INCLUDE_NAMES ${NetCDF_${comp}_INCLUDE_NAMES}
                               INCLUDE_HINTS ${NetCDF_${comp}_INCLUDE_HINTS}
                               LIBRARY_NAMES ${NetCDF_${comp}_LIBRARY_NAMES}
                               LIBRARY_HINTS ${NetCDF_${comp}_LIBRARY_HINTS})
        
    endif ()
    
endforeach ()

# SEARCH FOR DEPENDENCIES (only if SHARED libraries were found)
foreach (comp IN LISTS NetCDF_FIND_VALID_COMPONENTS)

    # If the component was found, and it is a static library...
    if (NetCDF_${comp}_FOUND AND NOT NetCDF_${comp}_IS_SHARED)
        
        # Search only if dependencies for this component were not already found
        if (NOT NetCDF_${comp}_DEPENDENCIES_SEARCHED)

            # COMPONENT: C
            if (comp STREQUAL C)
        
                # Look in netcdf_meta.h include file
                if (NOT NetCDF_C_META_DIR)
                
                    find_path (NetCDF_C_META_DIR
                               NAMES netcdf_meta.h
                               HINTS ${NetCDF_C_INCLUDE_DIRS})
                    if (NetCDF_C_META_DIR)
                    
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
                        
                    else ()
                    
                        message (WARNING "Could not find netcdf_meta.h")
                         
                    endif ()
                    
                endif ()
        
                # DEPENDENCY: HDF5
                find_package (HDF5 COMPONENTS HL C)
                if (HDF5_C_FOUND)
                    list (APPEND NetCDF_C_INCLUDE_DIRS ${HDF5_C_INCLUDE_DIRS}
                                                       ${HDF5_HL_INCLUDE_DIRS})
                    list (APPEND NetCDF_C_LIBRARIES ${HDF5_C_LIBRARIES}
                                                    ${HDF5_HL_LIBRARIES})
                endif ()

                # DEPENDENCY: CURL (If DAP enabled)
                if (NetCDF_C_HAS_DAP)
                    find_package (CURL)
                    if (CURL_FOUND)
                        list (APPEND NetCDF_C_INCLUDE_DIRS ${CURL_INCLUDE_DIRS})
                        list (APPEND NetCDF_C_LIBRARIES ${CURL_LIBRARIES}
                                                        ${CURL_LIBRARIES})
                    endif ()
                endif ()
                                
                # DEPENDENCY: LIBDL Math
                list (APPEND NetCDF_C_LIBRARIES -ldl -lm)
        
            # COMPONENT: Fortran
            elseif (comp STREQUAL Fortran)
                            
                # DEPENDENCY: NetCDF
                find_package (NetCDF COMPONENTS C)
                if (NetCDF_C_FOUND)
                    list (APPEND NetCDF_Fortran_INCLUDE_DIRS ${NetCDF_C_INCLUDE_DIRS})
                    list (APPEND NetCDF_Fortran_LIBRARIES ${NetCDF_C_LIBRARIES})
                endif ()
        
            endif ()
            
            set (NetCDF_${comp}_DEPENDENCIES_SEARCHED TRUE)
            
        endif ()
        
    endif ()

endforeach ()
