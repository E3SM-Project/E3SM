include (CMakeParseArguments)
include(FindPackageHandleStandardArgs)

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
# - Basic find package macro for a specific COMPonent
#
macro (find_package_component PKG)

    # Parse the input arguments
    set (options)
    set (oneValueArgs COMPONENT)
    set (multiValueArgs INCLUDE_NAMES INCLUDE_HINTS LIBRARY_NAMES LIBRARY_HINTS)
    cmake_parse_arguments (${PKG} "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})    
    set (COMP ${${PKG}_COMPONENT})
    if (COMP)
        set (PKGCOMP ${PKG}_${COMP})
    else ()
        set (PKGCOMP ${PKG})
    endif ()
    string (TOUPPER ${PKG} PKGUP)
    string (TOUPPER ${PKGCOMP} PKGCOMPUP)
    
    # Determine include dir search order
    set (INCLUDE_HINTS)
    if (${PKG}_INCLUDE_HINTS)
        list (APPEND INCLUDE_HINTS ${${PKG}_INCLUDE_HINTS})
    endif ()
    if (${PKGCOMP}_DIR)
        list (APPEND INCLUDE_HINTS ${${PKGCOMP}_DIR}/include)
    endif ()
    if (${PKG}_DIR)
        list (APPEND INCLUDE_HINTS ${${PKG}_DIR}/include)
    endif ()
    if (DEFINED ENV{${PKGUP}})
        list (APPEND INCLUDE_HINTS $ENV{${PKGUP}}/include)
    endif ()
    
    # Search for include file
    find_path (${PKGCOMP}_INCLUDE_DIR
               NAMES ${${PKG}_INCLUDE_NAMES}
               HINTS ${INCLUDE_HINTS})
               
    # Unset include search variables
    unset (INCLUDE_HINTS)
    
    # Determine library dir search order
    set (LIBRARY_HINTS)
    if (${PKG}_LIBRARY_HINTS)
        list (APPEND LIBRARY_HINTS ${${PKG}_LIBRARY_HINTS})
    endif ()
    if (${PKGCOMP}_DIR)
        list (APPEND LIBRARY_HINTS ${${PKGCOMP}_DIR}/lib)
    endif ()
    if (${PKG}_DIR)
        list (APPEND LIBRARY_HINTS ${${PKG}_DIR}/lib)
    endif ()
    if (DEFINED ENV{${PKGUP}})
        list (APPEND LIBRARY_HINTS $ENV{${PKGUP}}/lib)
    endif ()
    
    # Search for library file
    set (${PKGCOMP}_IS_SHARED TRUE)
    if (PREFER_SHARED)
        find_shared_library (${PKGCOMP}_LIBRARY
                             NAMES ${${PKG}_LIBRARY_NAMES}
                             HINTS ${LIBRARY_HINTS})
        if (NOT ${PKGCOMP}_LIBRARY)
            find_static_library (${PKGCOMP}_LIBRARY
                                 NAMES ${${PKG}_LIBRARY_NAMES}
                                 HINTS ${LIBRARY_HINTS})
            if (${PKGCOMP}_LIBRARY)
                set (${PKGCOMP}_IS_SHARED FALSE)
            endif ()
        endif ()
    else ()
        find_static_library (${PKGCOMP}_LIBRARY
                             NAMES ${${PKG}_LIBRARY_NAMES}
                             HINTS ${LIBRARY_HINTS})
        if (${PKGCOMP}_LIBRARY)
            set (${PKGCOMP}_IS_SHARED FALSE)
        else ()
            find_shared_library (${PKGCOMP}_LIBRARY
                                 NAMES ${${PKG}_LIBRARY_NAMES}
                                 HINTS ${LIBRARY_HINTS})
        endif ()
    endif ()
    
    # Unset include search variables
    unset (LIBRARY_HINTS)

    # handle the QUIETLY and REQUIRED arguments and 
    # set NetCDF_C_FOUND to TRUE if all listed variables are TRUE
    find_package_handle_standard_args (${PKGCOMP} DEFAULT_MSG
                                       ${PKGCOMP}_LIBRARY 
                                       ${PKGCOMP}_INCLUDE_DIR)
    mark_as_advanced (${PKGCOMP}_INCLUDE_DIR ${PKGCOMP}_LIBRARY)
    
    # HACK For bug in CMake v3.0:
    set (${PKGCOMP}_FOUND ${${PKGCOMPUP}_FOUND})

    # Set return variables
    if (${PKGCOMP}_FOUND)
        set (${PKGCOMP}_INCLUDE_DIRS ${${PKGCOMP}_INCLUDE_DIR})
        set (${PKGCOMP}_LIBRARIES ${${PKGCOMP}_LIBRARY})
        set (${PKGCOMP}_DEFINITIONS)
        set (${PKGCOMP}_OPTIONS)
    endif ()
    
    unset (COMP)
    unset (PKGCOMPUP)
    unset (PKGUP)

endmacro ()



