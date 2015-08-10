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
include (CMakeParseArguments)
include(FindPackageHandleStandardArgs)
macro (find_package_component PKG)

    # Parse the input arguments
    set (options)
    set (oneValueArgs COMPONENT)
    set (multiValueArgs INCLUDE_NAMES PRE_SEARCH_HINTS POST_SEARCH_HINTS LIBRARY_NAMES)
    cmake_parse_arguments (${PKG} "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})    
    set (COMP ${${PKG}_COMPONENT})
    string (TOUPPER ${PKG} PKGUP)
    string (TOUPPER ${COMP} COMPUP)
    
    
    # Determine include dir search order
    set (INCLUDE_HINTS)
    if (${PKG}_PRE_SEARCH_HINTS)
        foreach (hint IN LISTS ${PKG}_PRE_SEARCH_HINTS)
            list (APPEND INCLUDE_HINTS ${hint}/include)
        endforeach ()
    endif ()
    if (${PKG}_${COMP}_DIR)
        list (APPEND INCLUDE_HINTS ${${PKG}_${COMP}_DIR}/include)
    endif ()
    if (${PKG}_DIR)
        list (APPEND INCLUDE_HINTS ${${PKG}_DIR}/include)
    endif ()
    if (DEFINED ENV{${PKGUP}})
        list (APPEND INCLUDE_HINTS $ENV{${PKGUP}}/include)
    endif ()
    if (${PKG}_POST_SEARCH_HINTS)
        foreach (hint IN LISTS ${PKG}_POST_SEARCH_HINTS)
            list (APPEND INCLUDE_HINTS ${hint}/include)
        endforeach ()
    endif ()
    
    # Search for include file
    find_path (${PKG}_${COMP}_INCLUDE_DIR
               NAMES ${${PKG}_INCLUDE_NAMES}
               HINTS ${INCLUDE_HINTS})
               
    # Unset include search variables
    unset (INCLUDE_HINTS)
    
    # Determine library dir search order
    set (LIBRARY_HINTS)
    if (${PKG}_PRE_SEARCH_HINTS)
        foreach (hint IN LISTS ${PKG}_PRE_SEARCH_HINTS)
            list (APPEND INCLUDE_HINTS ${hint}/lib)
        endforeach ()
    endif ()
    if (${PKG}_${COMP}_DIR)
        list (APPEND LIBRARY_HINTS ${${PKG}_${COMP}_DIR}/lib)
    endif ()
    if (${PKG}_DIR)
        list (APPEND LIBRARY_HINTS ${${PKG}_DIR}/lib)
    endif ()
    if (DEFINED ENV{${PKGUP}})
        list (APPEND LIBRARY_HINTS $ENV{${PKGUP}}/lib)
    endif ()
    if (${PKG}_POST_SEARCH_HINTS)
        foreach (hint IN LISTS ${PKG}_POST_SEARCH_HINTS)
            list (APPEND INCLUDE_HINTS ${hint}/lib)
        endforeach ()
    endif ()
    
    # Search for library file
    include (LibFindLibraryMacros)
    set (${PKG}_${COMP}_IS_SHARED TRUE)
    if (PREFER_SHARED)
        find_shared_library (${PKG}_${COMP}_LIBRARY
                             NAMES ${${PKG}_LIBRARY_NAMES}
                             HINTS ${LIBRARY_HINTS})
        if (NOT ${PKG}_${COMP}_LIBRARY)
            find_static_library (${PKG}_${COMP}_LIBRARY
                                 NAMES ${${PKG}_LIBRARY_NAMES}
                                 HINTS ${LIBRARY_HINTS})
            if (${PKG}_${COMP}_LIBRARY)
                set (${PKG}_${COMP}_IS_SHARED FALSE)
            endif ()
        endif ()
    else ()
        find_static_library (${PKG}_${COMP}_LIBRARY
                             NAMES ${${PKG}_LIBRARY_NAMES}
                             HINTS ${LIBRARY_HINTS})
        if (${PKG}_${COMP}_LIBRARY)
            set (${PKG}_${COMP}_IS_SHARED FALSE)
        else ()
            find_shared_library (${PKG}_${COMP}_LIBRARY
                                 NAMES ${${PKG}_LIBRARY_NAMES}
                                 HINTS ${LIBRARY_HINTS})
        endif ()
    endif ()
    
    # Unset include search variables
    unset (LIBRARY_HINTS)

    # handle the QUIETLY and REQUIRED arguments and 
    # set NetCDF_C_FOUND to TRUE if all listed variables are TRUE
    find_package_handle_standard_args (${PKG}_${COMP} DEFAULT_MSG
                                       ${PKG}_${COMP}_LIBRARY 
                                       ${PKG}_${COMP}_INCLUDE_DIR)
    mark_as_advanced (${PKG}_${COMP}_INCLUDE_DIR ${PKG}_${COMP}_LIBRARY)
    
    # HACK For bug in CMake v3.0:
    set (${PKG}_${COMP}_FOUND ${${PKGUP}_${COMPUP}_FOUND})

    # Set return variables
    if (${PKG}_${COMP}_FOUND)
        set (${PKG}_${COMP}_INCLUDE_DIRS ${${PKG}_${COMP}_INCLUDE_DIR})
        set (${PKG}_${COMP}_LIBRARIES ${${PKG}_${COMP}_LIBRARY})
        set (${PKG}_${COMP}_DEFINITIONS)
        set (${PKG}_${COMP}_OPTIONS)
    endif ()
    
    unset (COMP)
    unset (COMPUP)
    unset (PKGUP)

endmacro ()