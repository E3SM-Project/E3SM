

find_path(Netcdf_INCLUDE_DIR 
          NAMES netcdf.h
          PATHS ${NETCDF_DIR} ${Homme_NETCDF_DIR}
          PATH_SUFFIXES include
          NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)

find_library(Netcdf_LIBRARY 
             NAMES libnetcdf.a netcdf
             HINTS ${Netcdf_INCLUDE_DIR}/../lib
             NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)

find_library(NetcdfF_LIBRARY 
             NAMES libnetcdff.a netcdff
             HINTS ${Netcdf_INCLUDE_DIR}/../lib
             NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)

find_library(Pnetcdf_LIBRARY 
             NAMES libpnetcdf.a netcdff
             HINTS ${Netcdf_INCLUDE_DIR}/../lib
             NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)

find_path(Netcdf_NC_CONFIG_BIN
          NAMES nc-config
          HINTS ${Netcdf_INCLUDE_DIR}/../bin
          NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)

# Includes only need to be Netcdf (for now)
set(Netcdf_INCLUDE_DIRS ${Netcdf_INCLUDE_DIR})

# Store libraries in Netcdf_LIBRARIES
set(Netcdf_LIBRARIES ${NetcdfF_LIBRARY} ${Netcdf_LIBRARY} ${Pnetcdf_LIBRARY})

IF (NOT ${Netcdf_NC_CONFIG_BIN} STREQUAL Netcdf_NC_CONFIG_BIN-NOTFOUND)

  # Probe nc-config to determine dependencies of Netcdf
  MESSAGE(STATUS "nc-config found at ${Netcdf_NC_CONFIG_BIN}")

  # use nc-confg --has-nc4 to determine if Netcdf depends upon HDF5
  EXECUTE_PROCESS(COMMAND ${Netcdf_NC_CONFIG_BIN}/nc-config --has-nc4
    RESULT_VARIABLE Homme_result
    OUTPUT_VARIABLE Homme_output
    ERROR_VARIABLE Homme_error
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  IF (${Homme_error})
    MESSAGE(FATAL_ERROR "${Netcdf_NC_CONFIG_BIN}/nc-config --has-nc4 produced an error")
  ELSE ()
    IF (${Homme_output} STREQUAL yes)
      SET (NETCDF_REQUIRE_HDF5 TRUE)
      MESSAGE(STATUS "nc-config: Netcdf depends upon HDF5")
    ELSE ()
      SET (NETCDF_REQUIRE_HDF5 FALSE)
      MESSAGE(STATUS "nc-config: Netcdf does not depend upon HDF5")
    ENDIF ()
  ENDIF ()

  # use nc-confg --has-dap to determine if Netcdf depends upon CURL
  EXECUTE_PROCESS(COMMAND ${Netcdf_NC_CONFIG_BIN}/nc-config --has-dap
    RESULT_VARIABLE Homme_result
    OUTPUT_VARIABLE Homme_output
    ERROR_VARIABLE Homme_error
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  IF (${Homme_error})
    MESSAGE(FATAL_ERROR "${Netcdf_NC_CONFIG_BIN}/nc-config --has-dap produced an error")
  ELSE ()
    IF (${Homme_output} STREQUAL yes)
      SET (NETCDF_REQUIRE_CURL TRUE)
      MESSAGE(STATUS "nc-config: Netcdf depends upon CURL")
    ELSE ()
      SET (NETCDF_REQUIRE_CURL FALSE)
      MESSAGE(STATUS "nc-config: Netcdf does not depend upon CURL")
    ENDIF ()
  ENDIF ()

ELSE ()
  SET (NETCDF_REQUIRE_HDF5 TRUE)
  SET (NETCDF_REQUIRE_CURL TRUE)
  MESSAGE(STATUS "nc-config not found assuming hdf5 and curl dependencies")
ENDIF ()

IF (${NETCDF_REQUIRE_CURL}) 
  
  # For some reasone CURL uses CURL_ROOT rather than CURL_DIR
  #   - change the variable for consistency
  SET(CURL_ROOT ${CURL_DIR})
  find_package(CURL)

  IF (${CURL_FOUND})
    MESSAGE(STATUS "Found CURL: ${CURL_LIBRARY}")
    set(Netcdf_LIBRARIES ${Netcdf_LIBRARIES} ${CURL_LIBRARY})
  ELSE ()
    MESSAGE(FATAL_ERROR "CURL Not found")
  ENDIF ()
ENDIF ()

IF (${NETCDF_REQUIRE_HDF5}) 

  find_path(HDF5_INCLUDE_DIR
            hdf5.h
            PATHS ${HDF5_DIR} ${Homme_HDF5_DIR}
            PATH_SUFFIXES include
            NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)

  find_library(HDF5_LIBRARY
               NAMES libhdf5.a hdf5
               PATHS ${HDF5_DIR} ${Homme_HDF5_DIR}
               PATH_SUFFIXES lib lib64
               NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)

  find_library(HDF5hl_LIBRARY
               NAMES libhdf5_hl.a hdf5_hl
               PATHS ${HDF5_DIR} ${Homme_HDF5_DIR}
               PATH_SUFFIXES lib lib64
               NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)


  if(${HDF5_LIBRARY} STREQUAL "HDF5_LIBRARY-NOTFOUND" OR ${HDF5hl_LIBRARY} STREQUAL "HDF5hl_LIBRARY-NOTFOUND")
    set(HDF5_FOUND OFF)
    MESSAGE(FATAL_ERROR "HDF5 not found, set HDF5_DIR to appropriate installation")
  else()
    set(HDF5_FOUND ON)
    MESSAGE(STATUS "Found HDF5: ${HDF5hl_LIBRARY} ${HDF5_LIBRARY}")
    set(Netcdf_LIBRARIES ${Netcdf_LIBRARIES} ${HDF5hl_LIBRARY} ${HDF5_LIBRARY})
  endif()

  # Check to see which dependencies (ZLIB, SZIP) hdf5 has
  INCLUDE(CheckSymbolExists)
  CHECK_SYMBOL_EXISTS(H5_HAVE_SZLIB_H "${HDF5_INCLUDE_DIR}/H5pubconf.h" HDF5_REQUIRE_SZ)
  CHECK_SYMBOL_EXISTS(H5_HAVE_ZLIB_H "${HDF5_INCLUDE_DIR}/H5pubconf.h" HDF5_REQUIRE_ZLIB)

  IF (${HDF5_REQUIRE_ZLIB})

    MESSAGE(STATUS "This HDF5 library requires ZLIB")
    # Find package always finds the shared object
    #find_package(ZLIB REQUIRED)
    find_library(ZLIB_LIBRARY
                 NAMES libz.a z
                 PATHS ${ZLIB_DIR} ${Homme_ZLIB_DIR}
                 PATH_SUFFIXES lib lib64
                 NO_SYSTEM_ENVIRONMENT_PATH)

    IF(${ZLIB_LIBRARY} STREQUAL "ZLIB_LIBRARY-NOTFOUND")
      SET(ZLIB_FOUND OFF)
      MESSAGE(FATAL_ERROR "ZLIB Not found")
    ELSE()
      SET(ZLIB_FOUND ON)
    ENDIF ()
    MESSAGE(STATUS "Found ZLIB_LIBRARY: ${ZLIB_LIBRARY}")
    MESSAGE(STATUS "ZLIB_LIBRARIES: ${ZLIB_LIBRARIES}")
    IF (${ZLIB_FOUND})
      MESSAGE(STATUS "Found ZLIB: ${ZLIB_LIBRARY}")
      set(Netcdf_LIBRARIES ${Netcdf_LIBRARIES} ${ZLIB_LIBRARY})
    ELSE ()
      MESSAGE(FATAL_ERROR "ZLIB Not found, set ZLIB_DIR to appropriate installation")
    ENDIF ()
  ENDIF ()

  IF (${HDF5_REQUIRE_SZ})
    MESSAGE(STATUS "This HDF5 library requires SZIP")

    find_library(SZIP_LIBRARY
                 NAMES libsz.a sz
                 PATHS ${SZIP_DIR} ${Homme_SZIP_DIR}
                 PATH_SUFFIXES lib lib64
                 NO_SYSTEM_ENVIRONMENT_PATH)

    IF(${SZIP_LIBRARY} STREQUAL "SZIP_LIBRARY-NOTFOUND")
      SET(SZIP_FOUND OFF)
      MESSAGE(FATAL_ERROR "SZIP Not found, set SZIP_DIR to appropriate installation")
    ELSE()
      SET(SZIP_FOUND ON)
      MESSAGE(STATUS "Found SZIP: ${SZIP_LIBRARY}")
      SET(Netcdf_LIBRARIES ${Netcdf_LIBRARIES} ${SZIP_LIBRARY})
    ENDIF()
  ENDIF ()

ENDIF()

# If we haven't had a FATAL_ERROR then we have found Netcdf and all dependencies
set(Netcdf_FOUND ON)

MESSAGE(STATUS "Found Netcdf:")
MESSAGE(STATUS "  Libraries: ${Netcdf_LIBRARIES}")
MESSAGE(STATUS "  Includes:  ${Netcdf_INCLUDE_DIRS}")

IF(Netcdf_FIND_REQUIRED AND NOT Netcdf_FOUND)
  MESSAGE(FATAL_ERROR "Did not find required library Netcdf")
ELSEIF(Netcdf_FOUND)
  MESSAGE(STATUS "Found Netcdf and all dependencies")
ENDIF()
