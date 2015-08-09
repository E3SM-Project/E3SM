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

# COMPONENT: C
if (";${NetCDF_FIND_COMPONENTS};" MATCHES ";C;")

    # Determine include dir search order
    set (NetCDF_C_INCLUDE_HINTS)
    if (NetCDF_C_DIR)
        list (APPEND NetCDF_C_INCLUDE_HINTS ${NetCDF_C_DIR}/include)
    endif ()
    if (NetCDF_DIR)
        list (APPEND NetCDF_C_INCLUDE_HINTS ${NetCDF_DIR}/include)
    endif ()
    if (DEFINED ENV{NETCDF})
        list (APPEND NetCDF_C_INCLUDE_HINTS $ENV{NETCDF}/include)
    endif ()
    
    # Search for include file
    find_path (NetCDF_C_INCLUDE_DIR
               NAMES netcdf.h
               HINTS ${NetCDF_C_INCLUDE_HINTS})
               
    # Unset include search variables
    unset (NetCDF_C_INCLUDE_HINTS)
    
    # Determine library dir search order
    set (NetCDF_C_LIBRARY_HINTS)
    if (NetCDF_C_DIR)
        list (APPEND NetCDF_C_LIBRARY_HINTS ${NetCDF_C_DIR}/lib)
    endif ()
    if (NetCDF_DIR)
        list (APPEND NetCDF_C_LIBRARY_HINTS ${NetCDF_DIR}/lib)
    endif ()
    if (DEFINED ENV{NETCDF})
        list (APPEND NetCDF_C_LIBRARY_HINTS $ENV{NETCDF}/lib)
    endif ()
    
    # Search for library file
    include (LibFindLibraryMacros)
    set (NetCDF_C_IS_SHARED TRUE)
    if (PREFER_SHARED)
        find_shared_library (NetCDF_C_LIBRARY
                             NAMES netcdf
                             HINTS ${NetCDF_C_LIBRARY_HINTS})
        if (NOT NetCDF_C_LIBRARY)
            find_static_library (NetCDF_C_LIBRARY
                                 NAMES netcdf
                                 HINTS ${NetCDF_C_LIBRARY_HINTS})
            if (NetCDF_C_LIBRARY)
                set (NetCDF_C_IS_SHARED FALSE)
            endif ()
        endif ()
    else ()
        find_static_library (NetCDF_C_LIBRARY
                             NAMES netcdf
                             HINTS ${NetCDF_C_LIBRARY_HINTS})
        if (NetCDF_C_LIBRARY)
            set (NetCDF_C_IS_SHARED FALSE)
        else ()
            find_shared_library (NetCDF_C_LIBRARY
                                 NAMES netcdf
                                 HINTS ${NetCDF_C_LIBRARY_HINTS})
        endif ()
    endif ()
    
    # Unset include search variables
    unset (NetCDF_C_LIBRARY_HINTS)
    
    # Set return variables
    set (NetCDF_C_INCLUDE_DIRS ${NetCDF_C_INCLUDE_DIR})
    set (NetCDF_C_LIBRARIES ${NetCDF_C_LIBRARY})
    set (NetCDF_C_DEFINITIONS)
    set (NetCDF_C_OPTIONS)
    
    # If static, look for dependencies (REQUIRED)
    if (NOT NetCDF_C_IS_SHARED)
    
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
         
    endif ()
    
    include(FindPackageHandleStandardArgs)
    # handle the QUIETLY and REQUIRED arguments and 
    # set NetCDF_C_FOUND to TRUE if all listed variables are TRUE
    find_package_handle_standard_args (NetCDF_C DEFAULT_MSG
                                       NetCDF_C_LIBRARY NetCDF_C_INCLUDE_DIR)
    mark_as_advanced (NetCDF_C_INCLUDE_DIR NetCDF_C_LIBRARY)
    
    # HACK For bug in CMake v3.0:
    set (NetCDF_C_FOUND ${NETCDF_C_FOUND})

endif ()

# COMPONENT: CXX
if (";${NetCDF_FIND_COMPONENTS};" MATCHES ";CXX;")

    # Determine include dir search order
    set (NetCDF_CXX_INCLUDE_HINTS)
    if (NetCDF_CXX_DIR)
        list (APPEND NetCDF_CXX_INCLUDE_HINTS ${NetCDF_CXX_DIR}/include)
    endif ()
    if (NetCDF_DIR)
        list (APPEND NetCDF_CXX_INCLUDE_HINTS ${NetCDF_DIR}/include)
    endif ()
    if (DEFINED ENV{NETCDF})
        list (APPEND NetCDF_CXX_INCLUDE_HINTS $ENV{NETCDF}/include)
    endif ()
    
    # Search for include file
    find_path (NetCDF_CXX_INCLUDE_DIR
               NAMES netcdfcpp.h
               HINTS ${NetCDF_CXX_INCLUDE_HINTS})
               
    # Unset include search variables
    unset (NetCDF_CXX_INCLUDE_HINTS)
    
    # Determine library dir search order
    set (NetCDF_CXX_LIBRARY_HINTS)
    if (NetCDF_CXX_DIR)
        list (APPEND NetCDF_CXX_LIBRARY_HINTS ${NetCDF_CXX_DIR}/lib)
    endif ()
    if (NetCDF_DIR)
        list (APPEND NetCDF_CXX_LIBRARY_HINTS ${NetCDF_DIR}/lib)
    endif ()
    if (DEFINED ENV{NETCDF})
        list (APPEND NetCDF_CXX_LIBRARY_HINTS $ENV{NETCDF}/lib)
    endif ()
    
    # Search for library file
    include (LibFindLibraryMacros)
    set (NetCDF_CXX_IS_SHARED TRUE)
    if (PREFER_SHARED)
        find_shared_library (NetCDF_CXX_LIBRARY
                             NAMES netcdf_c++
                             HINTS ${NetCDF_CXX_LIBRARY_HINTS})
        if (NOT NetCDF_CXX_LIBRARY)
            find_static_library (NetCDF_CXX_LIBRARY
                                 NAMES netcdf_c++
                                 HINTS ${NetCDF_CXX_LIBRARY_HINTS})
            if (NetCDF_CXX_LIBRARY)
                set (NetCDF_CXX_IS_SHARED FALSE)
            endif ()
        endif ()
    else ()
        find_static_library (NetCDF_CXX_LIBRARY
                             NAMES netcdf_c++
                             HINTS ${NetCDF_CXX_LIBRARY_HINTS})
        if (NetCDF_CXX_LIBRARY)
            set (NetCDF_CXX_IS_SHARED FALSE)
        else ()
            find_shared_library (NetCDF_CXX_LIBRARY
                                 NAMES netcdf_c++
                                 HINTS ${NetCDF_CXX_LIBRARY_HINTS})
        endif ()
    endif ()
    
    # Unset include search variables
    unset (NetCDF_CXX_LIBRARY_HINTS)
    
    # Set return variables
    set (NetCDF_CXX_INCLUDE_DIRS ${NetCDF_CXX_INCLUDE_DIR})
    set (NetCDF_CXX_LIBRARIES ${NetCDF_CXX_LIBRARY})
    set (NetCDF_CXX_DEFINITIONS)
    set (NetCDF_CXX_OPTIONS)
        
    include(FindPackageHandleStandardArgs)
    # handle the QUIETLY and REQUIRED arguments and 
    # set NetCDF_CXX_FOUND to TRUE if all listed variables are TRUE
    find_package_handle_standard_args (NetCDF_CXX DEFAULT_MSG
                                       NetCDF_CXX_LIBRARY NetCDF_CXX_INCLUDE_DIR)
    mark_as_advanced (NetCDF_CXX_INCLUDE_DIR NetCDF_CXX_LIBRARY)
    
    # HACK For bug in CMake v3.0:
    set (NetCDF_CXX_FOUND ${NETCDF_CXX_FOUND})

endif ()

# COMPONENT: CXX4
if (";${NetCDF_FIND_COMPONENTS};" MATCHES ";CXX4;")

    # Determine include dir search order
    set (NetCDF_CXX4_INCLUDE_HINTS)
    if (NetCDF_CXX4_DIR)
        list (APPEND NetCDF_CXX4_INCLUDE_HINTS ${NetCDF_CXX4_DIR}/include)
    endif ()
    if (NetCDF_DIR)
        list (APPEND NetCDF_CXX4_INCLUDE_HINTS ${NetCDF_DIR}/include)
    endif ()
    if (DEFINED ENV{NETCDF})
        list (APPEND NetCDF_CXX4_INCLUDE_HINTS $ENV{NETCDF}/include)
    endif ()
    
    # Search for include file
    find_path (NetCDF_CXX4_INCLUDE_DIR
               NAMES netcdfcpp.h
               HINTS ${NetCDF_CXX4_INCLUDE_HINTS})
               
    # Unset include search variables
    unset (NetCDF_CXX4_INCLUDE_HINTS)
    
    # Determine library dir search order
    set (NetCDF_CXX4_LIBRARY_HINTS)
    if (NetCDF_CXX4_DIR)
        list (APPEND NetCDF_CXX4_LIBRARY_HINTS ${NetCDF_CXX4_DIR}/lib)
    endif ()
    if (NetCDF_DIR)
        list (APPEND NetCDF_CXX4_LIBRARY_HINTS ${NetCDF_DIR}/lib)
    endif ()
    if (DEFINED ENV{NETCDF})
        list (APPEND NetCDF_CXX4_LIBRARY_HINTS $ENV{NETCDF}/lib)
    endif ()
    
    # Search for library file
    include (LibFindLibraryMacros)
    set (NetCDF_CXX4_IS_SHARED TRUE)
    if (PREFER_SHARED)
        find_shared_library (NetCDF_CXX4_LIBRARY
                             NAMES netcdf_c++4
                             HINTS ${NetCDF_CXX4_LIBRARY_HINTS})
        if (NOT NetCDF_CXX4_LIBRARY)
            find_static_library (NetCDF_CXX4_LIBRARY
                                 NAMES netcdf_c++4
                                 HINTS ${NetCDF_CXX4_LIBRARY_HINTS})
            if (NetCDF_CXX4_LIBRARY)
                set (NetCDF_CXX4_IS_SHARED FALSE)
            endif ()
        endif ()
    else ()
        find_static_library (NetCDF_CXX4_LIBRARY
                             NAMES netcdf_c++4
                             HINTS ${NetCDF_CXX4_LIBRARY_HINTS})
        if (NetCDF_CXX4_LIBRARY)
            set (NetCDF_CXX4_IS_SHARED FALSE)
        else ()
            find_shared_library (NetCDF_CXX4_LIBRARY
                                 NAMES netcdf_c++4
                                 HINTS ${NetCDF_CXX4_LIBRARY_HINTS})
        endif ()
    endif ()
    
    # Unset include search variables
    unset (NetCDF_CXX4_LIBRARY_HINTS)
    
    # Set return variables
    set (NetCDF_CXX4_INCLUDE_DIRS ${NetCDF_CXX4_INCLUDE_DIR})
    set (NetCDF_CXX4_LIBRARIES ${NetCDF_CXX4_LIBRARY})
    set (NetCDF_CXX4_DEFINITIONS)
    set (NetCDF_CXX4_OPTIONS)
        
    include(FindPackageHandleStandardArgs)
    # handle the QUIETLY and REQUIRED arguments and 
    # set NetCDF_CXX4_FOUND to TRUE if all listed variables are TRUE
    find_package_handle_standard_args (NetCDF_CXX4 DEFAULT_MSG
                                       NetCDF_CXX4_LIBRARY NetCDF_CXX4_INCLUDE_DIR)
    mark_as_advanced (NetCDF_CXX4_INCLUDE_DIR NetCDF_CXX4_LIBRARY)
    
    # HACK For bug in CMake v3.0:
    set (NetCDF_CXX4_FOUND ${NETCDF_CXX4_FOUND})

endif ()

# COMPONENT: Fortran
if (";${NetCDF_FIND_COMPONENTS};" MATCHES ";Fortran;")

    # Determine include dir search order
    set (NetCDF_Fortran_INCLUDE_HINTS)
    if (NetCDF_Fortran_DIR)
        list (APPEND NetCDF_Fortran_INCLUDE_HINTS ${NetCDF_Fortran_DIR}/include)
    endif ()
    if (NetCDF_DIR)
        list (APPEND NetCDF_Fortran_INCLUDE_HINTS ${NetCDF_DIR}/include)
    endif ()
    if (DEFINED ENV{NETCDF})
        list (APPEND NetCDF_Fortran_INCLUDE_HINTS $ENV{NETCDF}/include)
    endif ()
    
    # Search for include file
    find_path (NetCDF_Fortran_INCLUDE_DIR
               NAMES netcdf.h
               HINTS ${NetCDF_Fortran_INCLUDE_HINTS})
               
    # Unset include search variables
    unset (NetCDF_Fortran_INCLUDE_HINTS)
    
    # Determine library dir search order
    set (NetCDF_Fortran_LIBRARY_HINTS)
    if (NetCDF_Fortran_DIR)
        list (APPEND NetCDF_Fortran_LIBRARY_HINTS ${NetCDF_Fortran_DIR}/lib)
    endif ()
    if (NetCDF_DIR)
        list (APPEND NetCDF_Fortran_LIBRARY_HINTS ${NetCDF_DIR}/lib)
    endif ()
    if (DEFINED ENV{NETCDF})
        list (APPEND NetCDF_Fortran_LIBRARY_HINTS $ENV{NETCDF}/lib)
    endif ()
    
    # Search for library file
    include (LibFindLibraryMacros)
    set (NetCDF_Fortran_IS_SHARED TRUE)
    if (PREFER_SHARED)
        find_shared_library (NetCDF_Fortran_LIBRARY
                             NAMES netcdf
                             HINTS ${NetCDF_Fortran_LIBRARY_HINTS})
        if (NOT NetCDF_Fortran_LIBRARY)
            find_static_library (NetCDF_Fortran_LIBRARY
                                 NAMES netcdf
                                 HINTS ${NetCDF_Fortran_LIBRARY_HINTS})
            if (NetCDF_Fortran_LIBRARY)
                set (NetCDF_Fortran_IS_SHARED FALSE)
            endif ()
        endif ()
    else ()
        find_static_library (NetCDF_Fortran_LIBRARY
                             NAMES netcdf
                             HINTS ${NetCDF_Fortran_LIBRARY_HINTS})
        if (NetCDF_Fortran_LIBRARY)
            set (NetCDF_Fortran_IS_SHARED FALSE)
        else ()
            find_shared_library (NetCDF_Fortran_LIBRARY
                                 NAMES netcdf
                                 HINTS ${NetCDF_Fortran_LIBRARY_HINTS})
        endif ()
    endif ()
    
    # Unset include search variables
    unset (NetCDF_Fortran_LIBRARY_HINTS)
    
    # Set return variables
    set (NetCDF_Fortran_INCLUDE_DIRS ${NetCDF_Fortran_INCLUDE_DIR})
    set (NetCDF_Fortran_LIBRARIES ${NetCDF_Fortran_LIBRARY})
    set (NetCDF_Fortran_DEFINITIONS)
    set (NetCDF_Fortran_OPTIONS)
    
    # If static, look for dependencies (REQUIRED)
    if (NOT NetCDF_Fortran_IS_SHARED)
    
        # DEPENDENCY: HDF5
        find_package (HDF5 REQUIRED COMPONENTS C HL)
        if (HDF5_FOUND)
            list (APPEND NetCDF_Fortran_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS})
            list (APPEND NetCDF_Fortran_LIBRARIES ${HDF5_Fortran_LIBRARIES}
                                                  ${HDF5_Fortran_HL_LIBRARIES})
        endif ()
    
        # DEPENDENCY: CURL
        find_package (CURL REQUIRED)
        if (CURL_FOUND)
            list (APPEND NetCDF_Fortran_INCLUDE_DIRS ${CURL_INCLUDE_DIRS})
            list (APPEND NetCDF_Fortran_LIBRARIES ${CURL_LIBRARIES})
        endif ()
         
    endif ()
    
    include(FindPackageHandleStandardArgs)
    # handle the QUIETLY and REQUIRED arguments and 
    # set NetCDF_Fortran_FOUND to TRUE if all listed variables are TRUE
    find_package_handle_standard_args (NetCDF_Fortran DEFAULT_MSG
                                       NetCDF_Fortran_LIBRARY NetCDF_Fortran_INCLUDE_DIR)
    mark_as_advanced (NetCDF_Fortran_INCLUDE_DIR NetCDF_Fortran_LIBRARY)
    
    # HACK For bug in CMake v3.0:
    set (NetCDF_Fortran_FOUND ${NETCDF_FORTRAN_FOUND})

endif ()
