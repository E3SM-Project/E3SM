# - Try to find PnetCDF
#
# This can be controlled by setting the PnetCDF_DIR (or, equivalently, the 
# PNETCDF environment variable), or PnetCDF_<lang>_DIR CMake variables, where
# <lang> is the COMPONENT language one needs.
#
# Once done, this will define:
#
#   PnetCDF_<lang>_FOUND        (BOOL) - system has PnetCDF
#   PnetCDF_<lang>_IS_SHARED    (BOOL) - whether library is shared/dynamic
#   PnetCDF_<lang>_INCLUDE_DIR  (PATH) - Location of the C header file
#   PnetCDF_<lang>_INCLUDE_DIRS (LIST) - the PnetCDF include directories
#   PnetCDF_<lang>_LIBRARY      (FILE) - Path to the C library file
#   PnetCDF_<lang>_LIBRARIES    (LIST) - link these to use PnetCDF
#   PnetCDF_<lang>_DEFINITIONS  (LIST) - preprocessor macros to use with PnetCDF
#   PnetCDF_<lang>_OPTIONS      (LIST) - compiler options to use PnetCDF
#
# The available COMPONENTS are: C, CXX, Fortran

# Default to C component if no components specified
if (NOT PnetCDF_FIND_COMPONENTS)
    set (PnetCDF_FIND_COMPONENTS C)
endif ()

# COMPONENT: C
if (";${PnetCDF_FIND_COMPONENTS};" MATCHES ";C;")

    # Determine include dir search order
    set (PnetCDF_C_INCLUDE_HINTS)
    if (PnetCDF_C_DIR)
        list (APPEND PnetCDF_C_INCLUDE_HINTS ${PnetCDF_C_DIR}/include)
    endif ()
    if (PnetCDF_DIR)
        list (APPEND PnetCDF_C_INCLUDE_HINTS ${PnetCDF_DIR}/include)
    endif ()
    if (DEFINED ENV{PNETCDF})
        list (APPEND PnetCDF_C_INCLUDE_HINTS $ENV{PNETCDF}/include)
    endif ()
    
    # Search for include file
    find_path (PnetCDF_C_INCLUDE_DIR
               NAMES pnetcdf.h
               HINTS ${PnetCDF_C_INCLUDE_HINTS})
               
    # Unset include search variables
    unset (PnetCDF_C_INCLUDE_HINTS)
    
    # Determine library dir search order
    set (PnetCDF_C_LIBRARY_HINTS)
    if (PnetCDF_C_DIR)
        list (APPEND PnetCDF_C_LIBRARY_HINTS ${PnetCDF_C_DIR}/lib)
    endif ()
    if (PnetCDF_DIR)
        list (APPEND PnetCDF_C_LIBRARY_HINTS ${PnetCDF_DIR}/lib)
    endif ()
    if (DEFINED ENV{PNETCDF})
        list (APPEND PnetCDF_C_LIBRARY_HINTS $ENV{PNETCDF}/lib)
    endif ()
    
    # Search for library file
    include (LibFindLibraryMacros)
    set (PnetCDF_C_IS_SHARED TRUE)
    if (PREFER_SHARED)
        find_shared_library (PnetCDF_C_LIBRARY
                             NAMES pnetcdf
                             HINTS ${PnetCDF_C_LIBRARY_HINTS})
        if (NOT PnetCDF_C_LIBRARY)
            find_static_library (PnetCDF_C_LIBRARY
                                 NAMES pnetcdf
                                 HINTS ${PnetCDF_C_LIBRARY_HINTS})
            if (PnetCDF_C_LIBRARY)
                set (PnetCDF_C_IS_SHARED FALSE)
            endif ()
        endif ()
    else ()
        find_static_library (PnetCDF_C_LIBRARY
                             NAMES pnetcdf
                             HINTS ${PnetCDF_C_LIBRARY_HINTS})
        if (PnetCDF_C_LIBRARY)
            set (PnetCDF_C_IS_SHARED FALSE)
        else ()
            find_shared_library (PnetCDF_C_LIBRARY
                                 NAMES pnetcdf
                                 HINTS ${PnetCDF_C_LIBRARY_HINTS})
        endif ()
    endif ()
    
    # Unset include search variables
    unset (PnetCDF_C_LIBRARY_HINTS)
    
    # Set return variables
    set (PnetCDF_C_INCLUDE_DIRS ${PnetCDF_C_INCLUDE_DIR})
    set (PnetCDF_C_LIBRARIES ${PnetCDF_C_LIBRARY})
    set (PnetCDF_C_DEFINITIONS)
    set (PnetCDF_C_OPTIONS)
    
    include(FindPackageHandleStandardArgs)
    # handle the QUIETLY and REQUIRED arguments and 
    # set PnetCDF_C_FOUND to TRUE if all listed variables are TRUE
    find_package_handle_standard_args (PnetCDF_C DEFAULT_MSG
                                       PnetCDF_C_LIBRARY PnetCDF_C_INCLUDE_DIR)
    mark_as_advanced (PnetCDF_C_INCLUDE_DIR PnetCDF_C_LIBRARY)
    
    # HACK For bug in CMake v3.0:
    set (PnetCDF_C_FOUND ${PNETCDF_C_FOUND})

endif ()

# COMPONENT: CXX
if (";${PnetCDF_FIND_COMPONENTS};" MATCHES ";CXX;")

    # Determine include dir search order
    set (PnetCDF_CXX_INCLUDE_HINTS)
    if (PnetCDF_CXX_DIR)
        list (APPEND PnetCDF_CXX_INCLUDE_HINTS ${PnetCDF_CXX_DIR}/include)
    endif ()
    if (PnetCDF_DIR)
        list (APPEND PnetCDF_CXX_INCLUDE_HINTS ${PnetCDF_DIR}/include)
    endif ()
    if (DEFINED ENV{PNETCDF})
        list (APPEND PnetCDF_CXX_INCLUDE_HINTS $ENV{PNETCDF}/include)
    endif ()
    
    # Search for include file
    find_path (PnetCDF_CXX_INCLUDE_DIR
               NAMES pnetcdf
               HINTS ${PnetCDF_CXX_INCLUDE_HINTS})
               
    # Unset include search variables
    unset (PnetCDF_CXX_INCLUDE_HINTS)
    
    # Determine library dir search order
    set (PnetCDF_CXX_LIBRARY_HINTS)
    if (PnetCDF_CXX_DIR)
        list (APPEND PnetCDF_CXX_LIBRARY_HINTS ${PnetCDF_CXX_DIR}/lib)
    endif ()
    if (PnetCDF_DIR)
        list (APPEND PnetCDF_CXX_LIBRARY_HINTS ${PnetCDF_DIR}/lib)
    endif ()
    if (DEFINED ENV{PNETCDF})
        list (APPEND PnetCDF_CXX_LIBRARY_HINTS $ENV{PNETCDF}/lib)
    endif ()
    
    # Search for library file
    include (LibFindLibraryMacros)
    set (PnetCDF_CXX_IS_SHARED TRUE)
    if (PREFER_SHARED)
        find_shared_library (PnetCDF_CXX_LIBRARY
                             NAMES pnetcdf
                             HINTS ${PnetCDF_CXX_LIBRARY_HINTS})
        if (NOT PnetCDF_CXX_LIBRARY)
            find_static_library (PnetCDF_CXX_LIBRARY
                                 NAMES pnetcdf
                                 HINTS ${PnetCDF_CXX_LIBRARY_HINTS})
            if (PnetCDF_CXX_LIBRARY)
                set (PnetCDF_CXX_IS_SHARED FALSE)
            endif ()
        endif ()
    else ()
        find_static_library (PnetCDF_CXX_LIBRARY
                             NAMES pnetcdf
                             HINTS ${PnetCDF_CXX_LIBRARY_HINTS})
        if (PnetCDF_CXX_LIBRARY)
            set (PnetCDF_CXX_IS_SHARED FALSE)
        else ()
            find_shared_library (PnetCDF_CXX_LIBRARY
                                 NAMES pnetcdf
                                 HINTS ${PnetCDF_CXX_LIBRARY_HINTS})
        endif ()
    endif ()
    
    # Unset include search variables
    unset (PnetCDF_CXX_LIBRARY_HINTS)
    
    # Set return variables
    set (PnetCDF_CXX_INCLUDE_DIRS ${PnetCDF_CXX_INCLUDE_DIR})
    set (PnetCDF_CXX_LIBRARIES ${PnetCDF_CXX_LIBRARY})
    set (PnetCDF_CXX_DEFINITIONS)
    set (PnetCDF_CXX_OPTIONS)
        
    include(FindPackageHandleStandardArgs)
    # handle the QUIETLY and REQUIRED arguments and 
    # set PnetCDF_CXX_FOUND to TRUE if all listed variables are TRUE
    find_package_handle_standard_args (PnetCDF_CXX DEFAULT_MSG
                                       PnetCDF_CXX_LIBRARY PnetCDF_CXX_INCLUDE_DIR)
    mark_as_advanced (PnetCDF_CXX_INCLUDE_DIR PnetCDF_CXX_LIBRARY)
    
    # HACK For bug in CMake v3.0:
    set (PnetCDF_CXX_FOUND ${PNETCDF_CXX_FOUND})

endif ()

# COMPONENT: Fortran
if (";${PnetCDF_FIND_COMPONENTS};" MATCHES ";Fortran;")

    # Determine include dir search order
    set (PnetCDF_Fortran_INCLUDE_HINTS)
    if (PnetCDF_Fortran_DIR)
        list (APPEND PnetCDF_Fortran_INCLUDE_HINTS ${PnetCDF_Fortran_DIR}/include)
    endif ()
    if (PnetCDF_DIR)
        list (APPEND PnetCDF_Fortran_INCLUDE_HINTS ${PnetCDF_DIR}/include)
    endif ()
    if (DEFINED ENV{PNETCDF})
        list (APPEND PnetCDF_Fortran_INCLUDE_HINTS $ENV{PNETCDF}/include)
    endif ()
    
    # Search for include file
    find_path (PnetCDF_Fortran_INCLUDE_DIR
               NAMES pnetcdf.mod pnetcdf.inc
               HINTS ${PnetCDF_Fortran_INCLUDE_HINTS})
               
    # Unset include search variables
    unset (PnetCDF_Fortran_INCLUDE_HINTS)
    
    # Determine library dir search order
    set (PnetCDF_Fortran_LIBRARY_HINTS)
    if (PnetCDF_Fortran_DIR)
        list (APPEND PnetCDF_Fortran_LIBRARY_HINTS ${PnetCDF_Fortran_DIR}/lib)
    endif ()
    if (PnetCDF_DIR)
        list (APPEND PnetCDF_Fortran_LIBRARY_HINTS ${PnetCDF_DIR}/lib)
    endif ()
    if (DEFINED ENV{PNETCDF})
        list (APPEND PnetCDF_Fortran_LIBRARY_HINTS $ENV{PNETCDF}/lib)
    endif ()
    
    # Search for library file
    include (LibFindLibraryMacros)
    set (PnetCDF_Fortran_IS_SHARED TRUE)
    if (PREFER_SHARED)
        find_shared_library (PnetCDF_Fortran_LIBRARY
                             NAMES pnetcdf
                             HINTS ${PnetCDF_Fortran_LIBRARY_HINTS})
        if (NOT PnetCDF_Fortran_LIBRARY)
            find_static_library (PnetCDF_Fortran_LIBRARY
                                 NAMES pnetcdf
                                 HINTS ${PnetCDF_Fortran_LIBRARY_HINTS})
            if (PnetCDF_Fortran_LIBRARY)
                set (PnetCDF_Fortran_IS_SHARED FALSE)
            endif ()
        endif ()
    else ()
        find_static_library (PnetCDF_Fortran_LIBRARY
                             NAMES pnetcdf
                             HINTS ${PnetCDF_Fortran_LIBRARY_HINTS})
        if (PnetCDF_Fortran_LIBRARY)
            set (PnetCDF_Fortran_IS_SHARED FALSE)
        else ()
            find_shared_library (PnetCDF_Fortran_LIBRARY
                                 NAMES pnetcdf
                                 HINTS ${PnetCDF_Fortran_LIBRARY_HINTS})
        endif ()
    endif ()
    
    # Unset include search variables
    unset (PnetCDF_Fortran_LIBRARY_HINTS)
    
    # Set return variables
    set (PnetCDF_Fortran_INCLUDE_DIRS ${PnetCDF_Fortran_INCLUDE_DIR})
    set (PnetCDF_Fortran_LIBRARIES ${PnetCDF_Fortran_LIBRARY})
    set (PnetCDF_Fortran_DEFINITIONS)
    set (PnetCDF_Fortran_OPTIONS)
    
    include(FindPackageHandleStandardArgs)
    # handle the QUIETLY and REQUIRED arguments and 
    # set PnetCDF_Fortran_FOUND to TRUE if all listed variables are TRUE
    find_package_handle_standard_args (PnetCDF_Fortran DEFAULT_MSG
                                       PnetCDF_Fortran_LIBRARY PnetCDF_Fortran_INCLUDE_DIR)
    mark_as_advanced (PnetCDF_Fortran_INCLUDE_DIR PnetCDF_Fortran_LIBRARY)
    
    # HACK For bug in CMake v3.0:
    set (PnetCDF_Fortran_FOUND ${NETCDF_FORTRAN_FOUND})

endif ()
