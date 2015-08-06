# - Try to find the math library
#
# Once done, this will define:
#
#   LIBM_FOUND        (BOOL) - system has NetCDF
#   LIBM_INCLUDE_DIRS (LIST) - the NetCDF include directories
#   LIBM_LIBRARIES    (LIST) - link these to use NetCDF

# Search for include file
find_path (LIBM_INCLUDE_DIR
           NAMES math.h)

# Search for library file
if (BUILD_SHARED_LIBS)
    find_library (LIBM_LIBRARY
                  NAMES m)
else ()
    find_library (LIBM_LIBRARY
                  NAMES libm.a)
endif ()

# Set return variables
set (LIBM_INCLUDE_DIRS ${LIBM_INCLUDE_DIR})
set (LIBM_LIBRARIES ${LIBM_LIBRARY})

# Handle standard find_package arguments
include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and 
# set NetCDF_C_FOUND to TRUE if all listed variables are TRUE
find_package_handle_standard_args (LIBM DEFAULT_MSG
                                   LIBM_LIBRARY LIBM_INCLUDE_DIR)
mark_as_advanced (LIBM_INCLUDE_DIR LIBM_LIBRARY)
