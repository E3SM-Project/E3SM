include(FindPackageHandleStandardArgs)

FIND_PATH(PNETCDF_INCLUDE_DIR 
          pnetcdf.h
          HINTS ${PNETCDF_DIR}/include)


IF (${PREFER_SHARED})
  FIND_LIBRARY(PNETCDF_LIBRARY 
               NAMES pnetcdf
               HINTS ${PNETCDF_DIR}/lib)

ELSE ()
  FIND_LIBRARY(PNETCDF_LIBRARY 
               NAMES libpnetcdf.a pnetcdf
               HINTS ${PNETCDF_DIR}/lib)
ENDIF ()

SET(PNETCDF_LIBRARIES ${PNETCDF_LIBRARY} )
SET(PNETCDF_INCLUDE_DIRS ${PNETCDF_INCLUDE_DIR} )

# Handle QUIETLY and REQUIRED.
find_package_handle_standard_args(pnetcdf DEFAULT_MSG
  PNETCDF_LIBRARY PNETCDF_INCLUDE_DIR  )
