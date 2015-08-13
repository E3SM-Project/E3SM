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
include (CheckNetCDF)

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

#==============================================================================
# SEARCH FOR VALIDATED COMPONENTS
foreach (NetCDF_comp IN LISTS NetCDF_FIND_VALID_COMPONENTS)

    # If not found already, search...
    if (NOT NetCDF_${NetCDF_comp}_FOUND)

        # Manually add the MPI include and library dirs to search paths
        if (MPI_${NetCDF_comp}_FOUND)
            set (NetCDF_${NetCDF_comp}_INCLUDE_HINTS ${MPI_${NetCDF_comp}_INCLUDE_PATH})
            set (NetCDF_${NetCDF_comp}_LIBRARY_HINTS)
            foreach (lib IN LISTS MPI_${NetCDF_comp}_LIBRARIES)
                get_filename_component (libdir ${lib} PATH)
                list (APPEND NetCDF_${NetCDF_comp}_LIBRARY_HINTS ${libdir})
                unset (libdir)
            endforeach ()
        endif ()
        
        # Search for the package component    
        find_package_component(NetCDF COMPONENT ${NetCDF_comp}
                               INCLUDE_HINTS ${NetCDF_${NetCDF_comp}_INCLUDE_HINTS}
                               LIBRARY_HINTS ${NetCDF_${NetCDF_comp}_LIBRARY_HINTS})

        # Dependencies
        if (NetCDF_comp STREQUAL C AND NOT NetCDF_C_IS_SHARED)
        
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
                    list (APPEND NetCDF_C_LIBRARIES ${CURL_LIBRARIES})
                endif ()
            endif ()
            
            # DEPENDENCY: PnetCDF (if PnetCDF enabled)
            if (NetCDF_C_HAS_PNETCDF)
                find_package (PnetCDF COMPONENTS C)
                if (CURL_FOUND)
                    list (APPEND NetCDF_C_INCLUDE_DIRS ${PnetCDF_C_INCLUDE_DIRS})
                    list (APPEND NetCDF_C_LIBRARIES ${PnetCDF_C_LIBRARIES})
                endif ()
            endif ()
                            
            # DEPENDENCY: LIBDL Math
            list (APPEND NetCDF_C_LIBRARIES -ldl -lm)

        elseif (NetCDF_comp STREQUAL Fortran AND NOT NetCDF_Fortran_IS_SHARED)
        
            # DEPENDENCY: NetCDF
            set (orig_comps ${NetCDF_FIND_VALID_COMPONENTS})
            find_package (NetCDF COMPONENTS C)
            set (NetCDF_FIND_VALID_COMPONENTS ${orig_comps})
            if (NetCDF_C_FOUND)
                list (APPEND NetCDF_Fortran_INCLUDE_DIRS ${NetCDF_C_INCLUDE_DIRS})
                list (APPEND NetCDF_Fortran_LIBRARIES ${NetCDF_C_LIBRARIES})
            endif ()
            
        endif ()

    endif ()
    
endforeach ()

#==============================================================================
# CHECKS AND DEPENDENCIES
foreach (NetCDF_comp IN LISTS NetCDF_FIND_VALID_COMPONENTS)
    if (NetCDF_comp STREQUAL C AND NetCDF_C_FOUND)
        check_NetCDF_C ()
    elseif (NetCDF_comp STREQUAL Fortran AND NetCDF_Fortran_FOUND)
        check_NetCDF_Fortran ()
    endif ()
endforeach ()
