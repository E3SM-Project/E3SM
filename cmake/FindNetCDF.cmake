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
# The available COMPONENTS are: C, CXX, CXX4, Fortran
# If no components are specified, it assumes only C
include (LibFindLibraryMacros)

set (NetCDF_VALID_COMPONENTS C CXX Fortran)

if (NOT NetCDF_FIND_COMPONENTS)
    set (NetCDF_FIND_COMPONENTS C)
endif ()

set (NetCDF_FIND_VALID_COMPONENTS)
foreach (comp IN LISTS NetCDF_FIND_COMPONENTS)
    if (";${NetCDF_VALID_COMPONENTS};" MATCHES ";${comp};")
        list (APPEND NetCDF_FIND_VALID_COMPONENTS ${comp})
    endif ()
endforeach ()

set (NetCDF_C_INCLUDE_NAMES netcdf.h)
set (NetCDF_CXX_INCLUDE_NAMES netcdf)
set (NetCDF_CXX4_INCLUDE_NAMES netcdf)
set (NetCDF_Fortran_INCLUDE_NAMES netcdf.mod netcdf.inc)

set (NetCDF_C_LIBRARY_NAMES netcdf)
set (NetCDF_CXX_LIBRARY_NAMES netcdf_c++4 netcdf_c++)
set (NetCDF_Fortran_LIBRARY_NAMES netcdff)

foreach (comp IN LISTS NetCDF_FIND_VALID_COMPONENTS)

    find_package_component(NetCDF COMPONENT ${comp}
                           INCLUDE_NAMES ${NetCDF_${comp}_INCLUDE_NAMES}
                           LIBRARY_NAMES ${NetCDF_${comp}_LIBRARY_NAMES}
                           PRE_SEARCH_HINTS ${MPI_${comp}_INCLUDE_PATH})
    
    # Handle Dependencies, if static
    if (NOT NetCDF_${comp}_IS_SHARED)
    
        if (comp STREQUAL C)
        
            # DEPENDENCY: HDF5
            find_package (HDF5 REQUIRED COMPONENTS C HL)
            if (HDF5_FOUND)
                list (APPEND NetCDF_C_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS})
                list (APPEND NetCDF_C_LIBRARIES ${HDF5_C_LIBRARIES}
                                                ${HDF5_HL_LIBRARIES})
            endif ()
        
            # DEPENDENCY: CURL
            find_package (CURL REQUIRED)
            if (CURL_FOUND)
                list (APPEND NetCDF_C_INCLUDE_DIRS ${CURL_INCLUDE_DIRS})
                list (APPEND NetCDF_C_LIBRARIES ${CURL_LIBRARIES})
            endif ()
            
        elseif (comp STREQUAL CXX OR comp STREQUAL Fortran)
        
            # DEPENDENCY: NetCDF_C
            find_package (NetCDF REQUIRED COMPONENTS C)
            if (NetCDF_C_FOUND)
                list (APPEND NetCDF_CXX_INCLUDE_DIRS ${NetCDF_C_INCLUDE_DIRS})
                list (APPEND NetCDF_CXX_LIBRARIES ${NetCDF_C_LIBRARIES})
            endif ()

        endif ()

    endif ()
    
endforeach ()
