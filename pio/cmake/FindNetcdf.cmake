# library for reading and writing self describing array data.
#
# This module invokes the NETCDF wrapper compiler that should be installed
# alongside NETCDF.  Depending upon the NETCDF Configuration, the wrapper compiler
# is called either h5cc or h5pcc.  If this succeeds, the module will then call
# the compiler with the -show argument to see what flags are used when compiling
# an NETCDF client application.
#
# The module will optionally accept the COMPONENTS argument.  If no COMPONENTS
# are specified, then the find module will default to finding only the NETCDF C
# library.  If one or more COMPONENTS are specified, the module will attempt to
# find the language bindings for the specified components.  Currently, the only
# valid components are C, CXX, FORTRAN and F90.
#
# On UNIX systems, this module will read the variable NETCDF_USE_STATIC_LIBRARIES
# to determine whether or not to prefer a static link to a dynamic link for NETCDF
# and all of it's dependencies.  To use this feature, make sure that the
# NETCDF_USE_STATIC_LIBRARIES variable is set before the call to find_package.
#
# To provide the module with a hint about where to find your NETCDF installation,
# you can set the environment variable NETCDF_ROOT.  The Find module will then
# look in this path when searching for NETCDF executables, paths, and libraries.
#
# In addition to finding the includes and libraries required to compile an NETCDF
# client application, this module also makes an effort to find tools that come
# with the NETCDF distribution that may be useful for regression testing.
#
# This module will define the following variables:
#  NETCDF_INCLUDE_DIRS - Location of the NETCDF includes
#  NETCDF_INCLUDE_DIR - Location of the NETCDF includes (deprecated)
#  NETCDF_DEFINITIONS - Required compiler definitions for NETCDF
#  NETCDF_C_LIBRARIES - Required libraries for the NETCDF C bindings.
#  NETCDF_CXX_LIBRARIES - Required libraries for the NETCDF C++ bindings
#  NETCDF_FORTRAN_LIBRARIES - Required libraries for the NETCDF FORTRAN bindings
#  NETCDF_F90_LIBRARIES - Required libraries for the NETCDF FORTRAN 90 bindings
#  NETCDF_LIBRARIES - Required libraries for all requested bindings
#  NETCDF_FOUND - true if NETCDF was found on the system
#  NETCDF_LIBRARY_DIRS - the full set of library directories
#  NETCDF_IS_PARALLEL - Whether or not NETCDF was found with parallel IO support
#  NETCDF_CONFIG_EXECUTABLE - the path to the NC-CONFIG tool

#=============================================================================
# Copyright 2009 Kitware, Inc.
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)

# This module is maintained by Thijs Heus <thijs.heus@zmaw.de>.

include(SelectLibraryConfigurations)
include(FindPackageHandleStandardArgs)

# List of the valid NETCDF components
set( NETCDF_VALID_COMPONENTS
    FORTRAN
    F90
    CXX
    C
)

# Invoke the NETCDF wrapper compiler.  The compiler return value is stored to the
# return_value argument, the text output is stored to the output variable.
macro( _NETCDF_CONFIG flag output return_value )
    if( NETCDF_CONFIG_EXECUTABLE )
        exec_program( ${NETCDF_CONFIG_EXECUTABLE}
            ARGS ${flag}
            OUTPUT_VARIABLE ${output}
            RETURN_VALUE ${return_value}
        )
        if( ${${return_value}} EQUAL 0 )
            # do nothing
        else()
            message( STATUS
              "Unable to determine ${flag} from NC-CONFIG." )
        endif()
    endif()
endmacro()
#
# try to find the NETCDF wrapper compilers
find_program( NETCDF_CONFIG_EXECUTABLE
    NAMES nc-config
    HINTS ENV NETCDF_ROOT
    PATH_SUFFIXES bin Bin
    DOC "NETCDF CONFIG PROGRAM.  Used only to detect NETCDF compile flags." )
mark_as_advanced( NETCDF_CONFIG_EXECUTABLE )
set(output "no")
_NETCDF_CONFIG (--has-hdf5 output return)
set(HAS_HDF5 FALSE)
set(NETCDF_IS_PARALLEL FALSE)

if(${output} STREQUAL yes)
  set(HAS_HDF5 TRUE)
  if(NETCDF_FIND_QUIETLY OR NOT NETCDF_FIND_REQUIRED)
    find_package(HDF5)
  else()
    find_package(HDF5 REQUIRED)
  endif()
#        list( APPEND NETCDF_LIBRARIES_DEBUG
#            ${HDF5_LIBRARIES_DEBUG} )
#        list( APPEND NETCDF_LIBRARIES_RELEASE
#            ${HDF5_LIBRARIES_RELEASE}

  TRY_COMPILE(NC4_SUPPORT ${CMAKE_CURRENT_BINARY_DIR}/tryNC4
                          ${CMAKE_CURRENT_SOURCE_DIR}/cmake/TryNC4.f90)
  IF( ${NC4_SUPPORT})
    set (NETCDF_IS_PARALLEL ${HDF5_IS_PARALLEL})
  else()
    set (NETCDF_IS_PARALLEL FALSE)
  endif()
endif()

#_NETCDF_CONFIG (--has-pnetcdf output return)
#if(${output} STREQUAL yes)
#  set (NETCDF_IS_PARALLEL TRUE)
#else()
#   set(NETCDF_IS_PARALLEL FALSE)
#endif()
#set( NETCDF_IS_PARALLEL TRUE CACHE BOOL
#    "NETCDF library compiled with parallel IO support" )


if( NETCDF_INCLUDE_DIRS AND NETCDF_LIBRARIES )
    # Do nothing: we already have NETCDF_INCLUDE_PATH and NETCDF_LIBRARIES in the
    # cache, it would be a shame to override them
else()
    if( NOT NETCDF_FIND_COMPONENTS )
        set( NETCDF_LANGUAGE_BINDINGS "C" )
    else()
        # add the extra specified components, ensuring that they are valid.
        foreach( component ${NETCDF_FIND_COMPONENTS} )
            list( FIND NETCDF_VALID_COMPONENTS ${component} component_location )
            if( ${component_location} EQUAL -1 )
                message( FATAL_ERROR
                    "\"${component}\" is not a valid NETCDF component." )
            else()
                list( APPEND NETCDF_LANGUAGE_BINDINGS ${component} )
            endif()
        endforeach()
    endif()

    # seed the initial lists of libraries to find with items we know we need
    set( NETCDF_C_INCLUDE_NAMES netcdf.h )
    set( NETCDF_CXX_INCLUDE_NAMES netcdfcpp.h ${NETCDF_C_INCLUDE_NAMES} )
    set( NETCDF_FORTRAN_INCLUDE_NAMES ${NETCDF_C_INCLUDE_NAMES} )
    set( NETCDF_F90_INCLUDE_NAMES netcdf.mod typesizes.mod ${NETCDF_C_INCLUDE_NAMES} )

    set( NETCDF_C_LIBRARY_NAMES netcdf)
    set( NETCDF_CXX_LIBRARY_NAMES netcdf_c++ ${NETCDF_C_LIBRARY_NAMES} )
    set( NETCDF_FORTRAN_LIBRARY_NAMES netcdff ${NETCDF_C_LIBRARY_NAMES})
    set( NETCDF_F90_LIBRARY_NAMES ${NETCDF_FORTRAN_LIBRARY_NAMES} )

    set( NETCDF_REQUIRED netcdf.h netcdfcpp.h netcdf.mod typesizes.mod netcdf netcdf_c++)

    foreach( LANGUAGE ${NETCDF_LANGUAGE_BINDINGS} )
        foreach( INC ${NETCDF_${LANGUAGE}_INCLUDE_NAMES} )
          find_path( NETCDF_${INC}_INCLUDE_DIR ${INC}
              HINTS
                  ${NETCDF_${LANGUAGE}_INCLUDE_FLAGS}
                      ENV NETCDF_ROOT
              PATHS
              PATH_SUFFIXES
                  include
                  Include
          )
          mark_as_advanced( NETCDF_${INC}_INCLUDE_DIR )
          if (NETCDF_${INC}_INCLUDE_DIR)
            list( APPEND NETCDF_INCLUDE_DIRS ${NETCDF_${INC}_INCLUDE_DIR} )
          else()
            list( FIND NETCDF_REQUIRED ${INC} location )
            if( ${location} EQUAL -1 )
            else()
              if(NETCDF_FIND_REQUIRED)
                message( SEND_ERROR "\"${INC}\" is not found." )
              endif()
            else()
            endif()
          endif()
        endforeach()
        # find the NETCDF libraries
        foreach( LIB ${NETCDF_${LANGUAGE}_LIBRARY_NAMES} )
            if( UNIX AND NETCDF_USE_STATIC_LIBRARIES )
                # According to bug 1643 on the CMake bug tracker, this is the
                # preferred method for searching for a static library.
                # See http://www.cmake.org/Bug/view.php?id=1643.  We search
                # first for the full static library name, but fall back to a
                # generic search on the name if the static search fails.
                set( THIS_LIBRARY_SEARCH_DEBUG lib${LIB}d.a ${LIB}d )
                set( THIS_LIBRARY_SEARCH_RELEASE lib${LIB}.a ${LIB} )
            else()
                set( THIS_LIBRARY_SEARCH_DEBUG ${LIB}d )
                set( THIS_LIBRARY_SEARCH_RELEASE ${LIB} )
            endif()
            find_library( NETCDF_${LIB}_LIBRARY_DEBUG
                NAMES ${THIS_LIBRARY_SEARCH_DEBUG}
                HINTS ${NETCDF_${LANGUAGE}_LIBRARY_DIRS}
                ENV NETCDF_ROOT
                PATH_SUFFIXES lib64 Lib64 lib Lib)
            find_library( NETCDF_${LIB}_LIBRARY_RELEASE
                NAMES ${THIS_LIBRARY_SEARCH_RELEASE}
                HINTS ${NETCDF_${LANGUAGE}_LIBRARY_DIRS}
                ENV NETCDF_ROOT
                PATH_SUFFIXES lib64 Lib64 lib Lib )
            select_library_configurations( NETCDF_${LIB} )
            # even though we adjusted the individual library names in
            # select_library_configurations, we still need to distinguish
            # between debug and release variants because NETCDF_LIBRARIES will
            # need to specify different lists for debug and optimized builds.
            # We can't just use the NETCDF_${LIB}_LIBRARY variable (which was set
            # up by the selection macro above) because it may specify debug and
            # optimized variants for a particular library, but a list of
            # libraries is allowed to specify debug and optimized only once.
          if (NETCDF_${LIB}_LIBRARY_RELEASE)
            list( APPEND NETCDF_LIBRARIES_RELEASE ${NETCDF_${LIB}_LIBRARY_RELEASE} )
          endif()
          if (NETCDF_${LIB}_LIBRARY_DEBUG)
            list( APPEND NETCDF_LIBRARIES_DEBUG ${NETCDF_${LIB}_LIBRARY_DEBUG} )
          endif()
          if (NETCDF_${LIB}_LIBRARY_RELEASE OR NETCDF_${LIB}_LIBRARY_DEBUG )
          else()
#             message( STATUS "\"${LIB}\" is not found." )
            list( FIND NETCDF_REQUIRED ${LIB} location )
            if( ${location} EQUAL -1 )
            else()
              if(NETCDF_FIND_REQUIRED)
                message( SEND_ERROR "\"${LIB}\" is not found." )
              endif()
           endif()
          endif()
        endforeach()
        list( APPEND NETCDF_LIBRARY_DIRS ${NETCDF_${LANGUAGE}_LIBRARY_DIRS} )

        # Append the libraries for this language binding to the list of all
        # required libraries.
        list( APPEND NETCDF_LIBRARIES_DEBUG
            ${NETCDF_${LANGUAGE}_LIBRARIES_DEBUG} )
        list( APPEND NETCDF_LIBRARIES_RELEASE
            ${NETCDF_${LANGUAGE}_LIBRARIES_RELEASE} )
    endforeach()

    # We may have picked up some duplicates in various lists during the above
    # process for the language bindings (both the C and C++ bindings depend on
    # libz for example).  Remove the duplicates.
    if( NETCDF_INCLUDE_DIRS )
        list( REMOVE_DUPLICATES NETCDF_INCLUDE_DIRS )
    endif()
    if( NETCDF_LIBRARIES_DEBUG )
        list( REMOVE_DUPLICATES NETCDF_LIBRARIES_DEBUG )
    endif()
    if( NETCDF_LIBRARIES_RELEASE )
        list( REMOVE_DUPLICATES NETCDF_LIBRARIES_RELEASE )
    endif()
    if( NETCDF_LIBRARY_DIRS )
        list( REMOVE_DUPLICATES NETCDF_LIBRARY_DIRS )
    endif()

    # Construct the complete list of NETCDF libraries with debug and optimized
    # variants when the generator supports them.
    if( CMAKE_CONFIGURATION_TYPES OR CMAKE_BUILD_TYPE )
        set( NETCDF_LIBRARIES
            debug ${NETCDF_LIBRARIES_DEBUG}
            optimized ${NETCDF_LIBRARIES_RELEASE} )
    else()
        set( NETCDF_LIBRARIES ${NETCDF_LIBRARIES_RELEASE} )
    endif()
endif()
message ("NETCDF_IS_PARALLEL ${NETCDF_IS_PARALLEL} here")
find_package_handle_standard_args( NETCDF DEFAULT_MSG
    NETCDF_LIBRARIES
    NETCDF_INCLUDE_DIRS
)

mark_as_advanced(
    NETCDF_INCLUDE_DIRS
    NETCDF_LIBRARIES
    NETCDF_LIBRARY_DIRS
)

# For backwards compatibility we set NETCDF_INCLUDE_DIR to the value of
# NETCDF_INCLUDE_DIRS
set( NETCDF_INCLUDE_DIR "${NETCDF_INCLUDE_DIRS}" )

