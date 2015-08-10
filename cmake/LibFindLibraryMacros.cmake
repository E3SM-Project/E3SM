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
macro (find_package_component PKG COMP)

    # Parse the input arguments
    

    # Determine include dir search order
    set (${PKG}_${COMP}_INCLUDE_HINTS)
    if (MPI_${COMP}_FOUND)
        list (APPEND ${PKG}_${COMP}_INCLUDE_HINTS ${MPI_${COMP}_INCLUDE_PATH})
    endif ()
    if (${PKG}_${COMP}_DIR)
        list (APPEND ${PKG}_${COMP}_INCLUDE_HINTS ${${PKG}_${COMP}_DIR}/include)
    endif ()
    if (${PKG}_DIR)
        list (APPEND ${PKG}_${COMP}_INCLUDE_HINTS ${${PKG}_DIR}/include)
    endif ()
    if (DEFINED ENV{NETCDF})
        list (APPEND ${PKG}_${COMP}_INCLUDE_HINTS $ENV{NETCDF}/include)
    endif ()
    
    # Search for include file
    find_path (${PKG}_${COMP}_INCLUDE_DIR
               NAMES netcdf.h
               HINTS ${${PKG}_${COMP}_INCLUDE_HINTS})
               
    # Unset include search variables
    unset (${PKG}_${COMP}_INCLUDE_HINTS)
    
    # Determine library dir search order
    set (${PKG}_${COMP}_LIBRARY_HINTS)
    if (MPI_${COMP}_FOUND)
        list (APPEND ${PKG}_${COMP}_LIBRARY_HINTS ${MPI_${COMP}_LIBRARIES})
    endif ()
    if (${PKG}_${COMP}_DIR)
        list (APPEND ${PKG}_${COMP}_LIBRARY_HINTS ${${PKG}_${COMP}_DIR}/lib)
    endif ()
    if (${PKG}_DIR)
        list (APPEND ${PKG}_${COMP}_LIBRARY_HINTS ${${PKG}_DIR}/lib)
    endif ()
    if (DEFINED ENV{NETCDF})
        list (APPEND ${PKG}_${COMP}_LIBRARY_HINTS $ENV{NETCDF}/lib)
    endif ()
    
    # Search for library file
    include (LibFindLibraryMacros)
    set (${PKG}_${COMP}_IS_SHARED TRUE)
    if (PREFER_SHARED)
        find_shared_library (${PKG}_${COMP}_LIBRARY
                             NAMES netcdf
                             HINTS ${${PKG}_${COMP}_LIBRARY_HINTS})
        if (NOT ${PKG}_${COMP}_LIBRARY)
            find_static_library (${PKG}_${COMP}_LIBRARY
                                 NAMES netcdf
                                 HINTS ${${PKG}_${COMP}_LIBRARY_HINTS})
            if (${PKG}_${COMP}_LIBRARY)
                set (${PKG}_${COMP}_IS_SHARED FALSE)
            endif ()
        endif ()
    else ()
        find_static_library (${PKG}_${COMP}_LIBRARY
                             NAMES netcdf
                             HINTS ${${PKG}_${COMP}_LIBRARY_HINTS})
        if (${PKG}_${COMP}_LIBRARY)
            set (${PKG}_${COMP}_IS_SHARED FALSE)
        else ()
            find_shared_library (${PKG}_${COMP}_LIBRARY
                                 NAMES netcdf
                                 HINTS ${${PKG}_${COMP}_LIBRARY_HINTS})
        endif ()
    endif ()
    
    # Unset include search variables
    unset (${PKG}_${COMP}_LIBRARY_HINTS)
    
    # Set return variables
    set (${PKG}_${COMP}_INCLUDE_DIRS ${${PKG}_${COMP}_INCLUDE_DIR})
    set (${PKG}_${COMP}_LIBRARIES ${${PKG}_${COMP}_LIBRARY})
    set (${PKG}_${COMP}_DEFINITIONS)
    set (${PKG}_${COMP}_OPTIONS)
endmacro ()