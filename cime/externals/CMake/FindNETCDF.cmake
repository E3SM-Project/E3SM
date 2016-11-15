

find_path(Netcdf_INCLUDE_DIR
          NAMES netcdf.h
          PATHS ${NETCDF_DIR}
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

find_path(Netcdf_NC_CONFIG_BIN
          NAMES nc-config
          HINTS ${Netcdf_INCLUDE_DIR}/../bin
          NO_CMAKE_SYSTEM_PATH)

find_file(NETCDF4_PAR_H netcdf_par.h
          HINTS ${Netcdf_INCLUDE_DIR}
          NO_DEFAULT_PATH )

if(NOT NETCDF4_PAR_H)
  set(NETCDF4_PARALLEL "no")
  MESSAGE("NETCDF built without MPIIO")
else()
  set(NETCDF4_PARALLEL "yes")
  MESSAGE("NETCDF built with hdf5 MPIIO support")
endif()

# Store libraries in Netcdf_LIBRARIES
set(Netcdf_LIBRARIES ${NetcdfF_LIBRARY} ${Netcdf_LIBRARY})

IF (NOT ${Netcdf_NC_CONFIG_BIN} STREQUAL Netcdf_NC_CONFIG_BIN-NOTFOUND)

  # Probe nc-config to determine dependencies of Netcdf
  MESSAGE(STATUS "nc-config found at ${Netcdf_NC_CONFIG_BIN}")

  # use nc-config --has-pnetcdf to determine if Netcdf depends on pnetcdf
  EXECUTE_PROCESS(COMMAND ${Netcdf_NC_CONFIG_BIN}/nc-config --has-pnetcdf
    RESULT_VARIABLE NCCONFIG_RESULT
    OUTPUT_VARIABLE NCCONFIG_OUTPUT
    ERROR_VARIABLE NCCONFIG_ERROR
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  IF (${NCCONFIG_ERROR})
    MESSAGE(WARNING "${Netcdf_NC_CONFIG_BIN}/nc-config --has-pnetcdf produced an error. Assuming no pnetcdf")
    SET (NETCDF_REQUIRE_PNETCDF FALSE)
  ELSE ()
    IF (${NCCONFIG_OUTPUT} STREQUAL yes)
      SET (NETCDF_REQUIRE_PNETCDF TRUE)
      MESSAGE(STATUS "nc-config: Netcdf depends upon Pnetcdf")
    ELSE ()
      SET (NETCDF_REQUIRE_PNETCDF FALSE)
      MESSAGE(STATUS "nc-config: Netcdf does not depend upon Pnetcdf")
    ENDIF ()
  ENDIF ()

  # use nc-confg --has-nc4 to determine if Netcdf depends upon HDF5
  EXECUTE_PROCESS(COMMAND ${Netcdf_NC_CONFIG_BIN}/nc-config --has-nc4
    RESULT_VARIABLE NCCONFIG_RESULT
    OUTPUT_VARIABLE NCCONFIG_OUTPUT
    ERROR_VARIABLE NCCONFIG_ERROR
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  IF (${NCCONFIG_ERROR})
    MESSAGE(WARNING "${Netcdf_NC_CONFIG_BIN}/nc-config --has-nc4 produced an error. Assuming hdf5")
    SET (NETCDF_REQUIRE_HDF5 TRUE)
  ELSE ()
    IF (${NCCONFIG_OUTPUT} STREQUAL yes)
      SET (NETCDF_REQUIRE_HDF5 TRUE)
      MESSAGE(STATUS "nc-config: Netcdf depends upon HDF5")
    ELSE ()
      SET (NETCDF_REQUIRE_HDF5 FALSE)
      MESSAGE(STATUS "nc-config: Netcdf does not depend upon HDF5")
    ENDIF ()
  ENDIF ()

  # use nc-confg --has-dap to determine if Netcdf depends upon CURL
  EXECUTE_PROCESS(COMMAND ${Netcdf_NC_CONFIG_BIN}/nc-config --has-dap
    RESULT_VARIABLE NCCONFIG_RESULT
    OUTPUT_VARIABLE NCCONFIG_OUTPUT
    ERROR_VARIABLE NCCONFIG_ERROR
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  IF (${NCCONFIG_ERROR})
    MESSAGE(WARNING "${Netcdf_NC_CONFIG_BIN}/nc-config --has-dap produced an error. Assuming curl.")
    SET (NETCDF_REQUIRE_CURL TRUE)
  ELSE ()
    IF (${NCCONFIG_OUTPUT} STREQUAL yes)
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

IF (${NETCDF_REQUIRE_PNETCDF})

  find_library(Pnetcdf_LIBRARY
               NAMES libpnetcdf.a pnetcdf
               HINTS ${PNETCDF_DIR} ${NETCDF_DIR}
               PATH_SUFFIXES lib
               NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)

  find_path(Pnetcdf_INCLUDE_DIR
            pnetcdf.h
            PATHS ${PNETCDF_DIR} ${NETCDF_DIR}
            PATH_SUFFIXES include
            NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)

  IF (${Pnetcdf_LIBRARY} STREQUAL "Pnetcdf_LIBRARY-NOTFOUND" OR ${Pnetcdf_INCLUDE_DIR} STREQUAL "Pnetcdf_INCLUDE_DIR-NOTFOUND")
    MESSAGE(FATAL_ERROR "Pnetcdf not found, set PNETCDF_DIR to appropriate installation or install Pnetcdf alongside Netcdf")
  ELSE ()
    MESSAGE(STATUS "Found Pnetcdf: ${Pnetcdf_LIBRARY}")
    set(Netcdf_LIBRARIES ${Netcdf_LIBRARIES} ${Pnetcdf_LIBRARY})
    set(Netcdf_INCLUDE_DIR ${Netcdf_INCLUDE_DIR} ${Pnetcdf_INCLUDE_DIR})
  ENDIF()

ENDIF()

IF (${NETCDF_REQUIRE_HDF5})

  IF (HDF5_DIR)
    find_path(HDF5_INCLUDE_DIR
              hdf5.h
              PATHS ${HDF5_DIR}
              PATH_SUFFIXES include
              NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)

    find_library(HDF5_LIBRARY
                 NAMES libhdf5.a hdf5
                 PATHS ${HDF5_DIR}
                 PATH_SUFFIXES lib lib64
                 NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)

    find_library(HDF5hl_LIBRARY
                 NAMES libhdf5_hl.a hdf5_hl
                 PATHS ${HDF5_DIR}
                 PATH_SUFFIXES lib lib64
                 NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)


    if(${HDF5_LIBRARY} STREQUAL "HDF5_LIBRARY-NOTFOUND" OR ${HDF5hl_LIBRARY} STREQUAL "HDF5hl_LIBRARY-NOTFOUND")
      MESSAGE(FATAL_ERROR "HDF5 not found, set HDF5_DIR to appropriate installation")
    else()
      MESSAGE(STATUS "Found HDF5: ${HDF5hl_LIBRARY} ${HDF5_LIBRARY}")
      set(Netcdf_LIBRARIES ${Netcdf_LIBRARIES} ${HDF5hl_LIBRARY} ${HDF5_LIBRARY})
    endif()

  else()
    FIND_PACKAGE(HDF5 COMPONENTS C HL)

    if(${HDF5_FOUND})
      MESSAGE(STATUS "Adding hdf5 libraries ")
      set(NETCDF_C_LIBRARY ${NETCDF_C_LIBRARY} ${HDF5_LIBRARIES})
    endif()
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
                 PATHS ${ZLIB_DIR}
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
                 PATHS ${SZIP_DIR}
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

# Set all caps vars too
set(NETCDF_FOUND ON)
set(NETCDF_INCLUDE_DIR ${Netcdf_INCLUDE_DIR})
set(NETCDF_LIBRARIES ${Netcdf_LIBRARIES})

MESSAGE(STATUS "Found Netcdf:")
MESSAGE(STATUS "  Libraries: ${Netcdf_LIBRARIES}")
MESSAGE(STATUS "  Includes:  ${Netcdf_INCLUDE_DIR}")

IF(Netcdf_FIND_REQUIRED AND NOT Netcdf_FOUND)
  MESSAGE(FATAL_ERROR "Did not find required library Netcdf")
ELSEIF(Netcdf_FOUND)
  MESSAGE(STATUS "Found Netcdf and all dependencies")
ENDIF()
