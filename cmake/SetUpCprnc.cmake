
INCLUDE(ExternalProject)

# CPRNC locations
SET (CPRNC_INSTALL_DIR ${CMAKE_BINARY_DIR}/utils/cprnc)
SET (CPRNC_BINARY ${CMAKE_BINARY_DIR}/utils/cprnc/bin/cprnc)

######################################################################################
# If there is a queue then we cannot use the parallel fortran compiler to compile cprnc
#   since we then would have to submit to the queue to run the exectuable. Instead, we 
#   can just feed cprnc a separate serial fortran compiler. 
# Of course doing this automatically is not simple.
######################################################################################
IF (${HOMME_QUEUING}) # If we have a queue
  IF (NOT FORTRAN_SER) # If the user has not specified a serial fortran compiler
    ##################################################################################
    # Note that the following looks totally insane. However there is a rationale here.
    ##################################################################################
    # If we are using a parallel version of a compiler, then it makes sense that we 
    #   would also have the serial version of that compiler - so first, search for that. 
    #   If the serial version of that compiler is not found, then prefer other fortan 
    #   compilers in the following order (gnu,intel,pgi,pathscale)
    ##################################################################################
    IF (${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
      FIND_PROGRAM(FORTRAN_SER ifort)
      IF (NOT FORTRAN_SER)
        FIND_PROGRAM(FORTRAN_SER gfortran)
      ENDIF ()
      IF (NOT FORTRAN_SER)
        FIND_PROGRAM(FORTRAN_SER pgf90)
      ENDIF ()
      IF (NOT FORTRAN_SER)
        FIND_PROGRAM(FORTRAN_SER pathf90)
      ENDIF ()
    ELSEIF (${CMAKE_Fortran_COMPILER_ID} STREQUAL GNU)
      FIND_PROGRAM(FORTRAN_SER gfortran)
      IF (NOT FORTRAN_SER)
        FIND_PROGRAM(FORTRAN_SER ifort)
      ENDIF ()
      IF (NOT FORTRAN_SER)
        FIND_PROGRAM(FORTRAN_SER pgf90)
      ENDIF ()
      IF (NOT FORTRAN_SER)
        FIND_PROGRAM(FORTRAN_SER pathf90)
      ENDIF ()
    ELSEIF (${CMAKE_Fortran_COMPILER_ID} STREQUAL PGI)
      FIND_PROGRAM(FORTRAN_SER pgf90)
      IF (NOT FORTRAN_SER)
        FIND_PROGRAM(FORTRAN_SER gfortran)
      ENDIF ()
      IF (NOT FORTRAN_SER)
        FIND_PROGRAM(FORTRAN_SER ifort)
      ENDIF ()
      IF (NOT FORTRAN_SER)
        FIND_PROGRAM(FORTRAN_SER pathf90)
      ENDIF ()
    ELSEIF (${CMAKE_Fortran_COMPILER_ID} STREQUAL PathScale)
      FIND_PROGRAM(FORTRAN_SER pathf90)
      IF (NOT FORTRAN_SER)
        FIND_PROGRAM(FORTRAN_SER gfortran)
      ENDIF ()
      IF (NOT FORTRAN_SER)
        FIND_PROGRAM(FORTRAN_SER ifort)
      ENDIF ()
      IF (NOT FORTRAN_SER)
        FIND_PROGRAM(FORTRAN_SER pgf90)
      ENDIF ()
    ENDIF ()
  ENDIF ()
ELSE ()
  # We don't have a queue so we can (presumably use the primary fortran compiler)
  SET (FORTRAN_SER ${CMAKE_Fortran_COMPILER})
ENDIF ()

IF (FORTRAN_SER)
  MESSAGE(STATUS "Using serial fortran compiler ${FORTRAN_SER} to compile CPRNC")
ELSE ()
  MESSAGE(FATAL_ERROR "Could not find a serial fortran compiler to compile CPRNC")
ENDIF ()
  

# The following builds cprnc 
ExternalProject_Add(
  cprnc
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/utils/cprnc
  CMAKE_ARGS
    -DCMAKE_Fortran_COMPILER=${FORTRAN_SER}
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DNETCDF_DIR=${CPRNC_NETCDF_DIR}
    -DHDF5_DIR=${CPRNC_HDF5_DIR}
    -DSZIP_DIR=${CPRNC_SZIP_DIR}
    -DZLIB_DIR=${CPRNC_ZLIB_DIR}
    -DCURL_DIR=${CPRNC_CURL_DIR}
    -DCMAKE_INSTALL_PREFIX=${CPRNC_INSTALL_DIR}
  INSTALL_DIR ${CPRNC_INSTALL_DIR}
  LOG_CONFIGURE ON)



