#
# - Wrapper for finding static libraries ONLY
#
macro (find_static_library)
    set (_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})
    set (CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX})
    find_library(${ARGN})
    set (CMAKE_FIND_LIBRARY_SUFFIXES ${_CMAKE_FIND_LIBRARY_SUFFIXES})
endmacro ()

#
# - Wrapper for finding shared/dynamic libraries ONLY
#
macro (find_shared_library)
    set (_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})
    set (CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_SHARED_LIBRARY_SUFFIX})
    find_library(${ARGN})
    set (CMAKE_FIND_LIBRARY_SUFFIXES ${_CMAKE_FIND_LIBRARY_SUFFIXES})
endmacro ()

#
# - Basic find package macro for a specific component
#
include (CMakeParseArguments)
include(FindPackageHandleStandardArgs)
macro (find_package_component PKG)

    # Parse the input arguments
    set(options)
    set(oneValueArgs COMPONENT)
    set(multiValueArgs INCLUDE_NAMES LIBRARY_NAMES)
    cmake_parse_arguments(${PKG} "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})    

    # Determine include dir search order
    set (${PKG}_${${PKG}_COMPONENT}_INCLUDE_HINTS)
    if (MPI_${${PKG}_COMPONENT}_FOUND)
        list (APPEND ${PKG}_${${PKG}_COMPONENT}_INCLUDE_HINTS ${MPI_${${PKG}_COMPONENT}_INCLUDE_PATH})
    endif ()
    if (${PKG}_${${PKG}_COMPONENT}_DIR)
        list (APPEND ${PKG}_${${PKG}_COMPONENT}_INCLUDE_HINTS ${${PKG}_${${PKG}_COMPONENT}_DIR}/include)
    endif ()
    if (${PKG}_DIR)
        list (APPEND ${PKG}_${${PKG}_COMPONENT}_INCLUDE_HINTS ${${PKG}_DIR}/include)
    endif ()
    if (DEFINED ENV{NETCDF})
        list (APPEND ${PKG}_${${PKG}_COMPONENT}_INCLUDE_HINTS $ENV{NETCDF}/include)
    endif ()
    
    # Search for include file
    find_path (${PKG}_${${PKG}_COMPONENT}_INCLUDE_DIR
               NAMES ${${PKG}_${${PKG}_COMPONENT}_INCLUDE_NAMES}
               HINTS ${${PKG}_${${PKG}_COMPONENT}_INCLUDE_HINTS})
               
    # Unset include search variables
    unset (${PKG}_${${PKG}_COMPONENT}_INCLUDE_HINTS)
    
    # Determine library dir search order
    set (${PKG}_${${PKG}_COMPONENT}_LIBRARY_HINTS)
    if (MPI_${${PKG}_COMPONENT}_FOUND)
        list (APPEND ${PKG}_${${PKG}_COMPONENT}_LIBRARY_HINTS ${MPI_${${PKG}_COMPONENT}_LIBRARIES})
    endif ()
    if (${PKG}_${${PKG}_COMPONENT}_DIR)
        list (APPEND ${PKG}_${${PKG}_COMPONENT}_LIBRARY_HINTS ${${PKG}_${${PKG}_COMPONENT}_DIR}/lib)
    endif ()
    if (${PKG}_DIR)
        list (APPEND ${PKG}_${${PKG}_COMPONENT}_LIBRARY_HINTS ${${PKG}_DIR}/lib)
    endif ()
    if (DEFINED ENV{NETCDF})
        list (APPEND ${PKG}_${${PKG}_COMPONENT}_LIBRARY_HINTS $ENV{NETCDF}/lib)
    endif ()
    
    # Search for library file
    include (LibFindLibraryMacros)
    set (${PKG}_${${PKG}_COMPONENT}_IS_SHARED TRUE)
    if (PREFER_SHARED)
        find_shared_library (${PKG}_${${PKG}_COMPONENT}_LIBRARY
                             NAMES ${${PKG}_${${PKG}_COMPONENT}_LIBRARY_NAMES}
                             HINTS ${${PKG}_${${PKG}_COMPONENT}_LIBRARY_HINTS})
        if (NOT ${PKG}_${${PKG}_COMPONENT}_LIBRARY)
            find_static_library (${PKG}_${${PKG}_COMPONENT}_LIBRARY
                                 NAMES ${${PKG}_${${PKG}_COMPONENT}_LIBRARY_NAMES}
                                 HINTS ${${PKG}_${${PKG}_COMPONENT}_LIBRARY_HINTS})
            if (${PKG}_${${PKG}_COMPONENT}_LIBRARY)
                set (${PKG}_${${PKG}_COMPONENT}_IS_SHARED FALSE)
            endif ()
        endif ()
    else ()
        find_static_library (${PKG}_${${PKG}_COMPONENT}_LIBRARY
                             NAMES ${${PKG}_${${PKG}_COMPONENT}_LIBRARY_NAMES}
                             HINTS ${${PKG}_${${PKG}_COMPONENT}_LIBRARY_HINTS})
        if (${PKG}_${${PKG}_COMPONENT}_LIBRARY)
            set (${PKG}_${${PKG}_COMPONENT}_IS_SHARED FALSE)
        else ()
            find_shared_library (${PKG}_${${PKG}_COMPONENT}_LIBRARY
                                 NAMES ${${PKG}_${${PKG}_COMPONENT}_LIBRARY_NAMES}
                                 HINTS ${${PKG}_${${PKG}_COMPONENT}_LIBRARY_HINTS})
        endif ()
    endif ()
    
    # Unset include search variables
    unset (${PKG}_${${PKG}_COMPONENT}_LIBRARY_HINTS)

    # handle the QUIETLY and REQUIRED arguments and 
    # set NetCDF_C_FOUND to TRUE if all listed variables are TRUE
    find_package_handle_standard_args (${PKG}_${${PKG}_COMPONENT} DEFAULT_MSG
                                       ${PKG}_${${PKG}_COMPONENT}_LIBRARY 
                                       ${PKG}_${${PKG}_COMPONENT}_INCLUDE_DIR)
    mark_as_advanced (${PKG}_${${PKG}_COMPONENT}_INCLUDE_DIR ${PKG}_${${PKG}_COMPONENT}_LIBRARY)
    
    # HACK For bug in CMake v3.0:
    set (NetCDF_C_FOUND ${NETCDF_C_FOUND})

    # Set return variables
    if (${PKG}_${${PKG}_COMPONENT}_FOUND)
        set (${PKG}_${${PKG}_COMPONENT}_INCLUDE_DIRS ${${PKG}_${${PKG}_COMPONENT}_INCLUDE_DIR})
        set (${PKG}_${${PKG}_COMPONENT}_LIBRARIES ${${PKG}_${${PKG}_COMPONENT}_LIBRARY})
        set (${PKG}_${${PKG}_COMPONENT}_DEFINITIONS)
        set (${PKG}_${${PKG}_COMPONENT}_OPTIONS)
    endif ()
endmacro ()