##############################################################################
# Compiler specific options
##############################################################################

SET(CMAKE_Fortran_FLAGS "")

MESSAGE(STATUS "CMAKE_Fortran_COMPILER_ID = ${CMAKE_Fortran_COMPILER_ID}")
# Need this for a fix in repro_sum_mod
IF (${CMAKE_Fortran_COMPILER_ID} STREQUAL XL)
  ADD_DEFINITIONS(-DnoI8)
ENDIF ()

IF (DEFINED BASE_FFLAGS)
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${BASE_FFLAGS}")
ELSE ()
  IF (CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
    if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "10")
      SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
    ELSE ()
      SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -w -fallow-argument-mismatch -ffree-line-length-none")
    endif()
    ADD_DEFINITIONS(-DCPRGNU)
  ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL PGI)
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mextend -Mflushz")
    # Needed by csm_share
    ADD_DEFINITIONS(-DCPRPGI)
  ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL PathScale)
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -extend-source")
  ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -assume byterecl")
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fp-model fast -ftz")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fp-model fast -ftz")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fp-model fast -ftz")
    #SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fp-model fast -qopt-report=5 -ftz")
    #SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mP2OPT_hpo_matrix_opt_framework=0 -fp-model fast -qopt-report=5 -ftz")

    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -diag-disable 8291")

    # remark #8291: Recommended relationship between field width 'W' and the number of fractional digits 'D' in this edit descriptor is 'W>=D+7'.

    # Needed by csm_share
    ADD_DEFINITIONS(-DCPRINTEL)
  ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL XL)
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -WF,-C! -qstrict -qnosave")
    # Needed by csm_share
    ADD_DEFINITIONS(-DCPRIBM)
  ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL NAG)
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -kind=byte -wmismatch=mpi_send,mpi_recv,mpi_bcast,mpi_allreduce,mpi_reduce,mpi_isend,mpi_irecv,mpi_irsend,mpi_rsend,mpi_gatherv,mpi_gather,mpi_scatterv,mpi_allgather,mpi_alltoallv,mpi_file_read_all,mpi_file_write_all,mpi_file_read_at")
#    SET(OPT_FFLAGS "${OPT_FFLAGS} -ieee=full -O2")
    SET(DEBUG_FFLAGS "${DEBUG_FFLAGS} -g -time -f2003 -ieee=stop")
    ADD_DEFINITIONS(-DHAVE_F2003_PTR_BND_REMAP)
    # Needed by both PIO and csm_share
    ADD_DEFINITIONS(-DCPRNAG)
  ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL Cray)
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DHAVE_F2003_PTR_BND_REMAP")
    # Needed by csm_share
    ADD_DEFINITIONS(-DCPRCRAY)
 ENDIF ()
ENDIF ()

IF (${HOMME_USE_CXX})
  IF (DEFINED BASE_CPPFLAGS)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${BASE_CPPFLAGS}")
  ENDIF ()

  # C++ Flags

  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

  INCLUDE(CheckCXXCompilerFlag)
  # CHECK_CXX_COMPILER_FLAG("-cxxlib" CXXLIB_SUPPORTED)
  # IF (CXXLIB_SUPPORTED)
    # SET(CXXLIB_SUPPORTED_CACHE TRUE CACHE BOOL "")
  # ELSE()
    # SET(CXXLIB_SUPPORTED_CACHE FALSE CACHE BOOL "")
  # ENDIF ()

  # Handle Cuda.
  FIND_PACKAGE(CUDA QUIET)
  IF (${CUDA_FOUND})
    EXECUTE_PROCESS(COMMAND ${CMAKE_CXX_COMPILER} --version
                    OUTPUT_VARIABLE CXX_COMPILER_VERSION_OUT)
    STRING (FIND "${CXX_COMPILER_VERSION_OUT}" "nvcc" pos)
    IF (${pos} GREATER -1)
      SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --expt-extended-lambda")
      SET (CUDA_BUILD TRUE)
      MESSAGE (STATUS "Found CUDA, with nvcc as C++ backend compiler. Great.")
    ELSE ()
      MESSAGE ("\n ********************************** WARNING ********************************** ")
      MESSAGE ("  Cuda was found, but the backend C++ compiler is not nvcc_wrapper.")
      MESSAGE ("  If you are using mpicxx from OpenMPI, you must set the OMPI_CXX")
      MESSAGE ("  env variable to point to the nvcc_wrapper in your kokkos installation.")
      MESSAGE ("  If you are using mpicxx from MPICH, you must set the MPICH_CXX")
      MESSAGE ("  env variable to point to the nvcc_wrapper in your kokkos installation.")
      MESSAGE ("  If you are building kokkos inside E3SM, that is:\n")
      MESSAGE ("${HOMME_SOURCE_DIR}/../../externals/kokkos/bin/nvcc_wrapper\n\n")
      MESSAGE ("  If you set HOMMEXX_EXEC_SPACE=Cuda, or if you are using Kokkos")
      MESSAGE ("  default device and that happens to be Cuda, this build may not compile,")
      MESSAGE ("  or give wrong results.")
      MESSAGE ("  For the record, this was the command used:\n")
      MESSAGE ("${CMAKE_CXX_COMPILER} --version\n")
      MESSAGE ("  And this was the output:\n")
      MESSAGE ("${CXX_COMPILER_VERSION_OUT}")
      MESSAGE (" ***************************************************************************** \n")
    ENDIF ()
  ENDIF ()
ENDIF()

##############################################################################
# Optimization flags
# 1) OPT_FLAGS if specified sets the Fortran,C, and CXX optimization flags
# 2) OPT_FFLAGS if specified sets the Fortran optimization flags
# 3) OPT_CFLAGS if specified sets the C optimization flags
# 4) OPT_CXXFLAGS if specified sets the CXX optimization flags
##############################################################################
IF (OPT_FLAGS)
  # Flags for Fortran C and CXX
  SET (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OPT_FLAGS}")
  SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OPT_FLAGS}")
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPT_FLAGS}")
ELSE ()

  IF (OPT_FFLAGS)
    # User specified optimization flags
    SET (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OPT_FFLAGS}")
  ELSE ()
    # Defaults
    IF (CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
      SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3")
    ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL PGI)
      SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2")
    ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL PathScale)
    ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
      SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3")
      #SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mavx -DTEMP_INTEL_COMPILER_WORKAROUND_001")
    ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL XL)
      SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2 -qmaxmem=-1")
    ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL Cray)
      SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2")
    ENDIF ()
  ENDIF ()

  IF (OPT_CFLAGS)
    SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OPT_CFLAGS}")
  ELSE ()
    IF (CMAKE_C_COMPILER_ID STREQUAL GNU)
      SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2")
    ELSEIF (CMAKE_C_COMPILER_ID STREQUAL PGI)
    ELSEIF (CMAKE_C_COMPILER_ID STREQUAL PathScale)
    ELSEIF (CMAKE_C_COMPILER_ID STREQUAL Intel)
      SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
      #SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mavx -DTEMP_INTEL_COMPILER_WORKAROUND_001")
    ELSEIF (CMAKE_C_COMPILER_ID STREQUAL XL)
      SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2 -qmaxmem=-1")
    ELSEIF (CMAKE_C_COMPILER_ID STREQUAL Cray)
      SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2")
    ENDIF ()
  ENDIF ()

  IF (OPT_CXXFLAGS)
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPT_CXXFLAGS}")
  ELSE ()
    IF (CMAKE_CXX_COMPILER_ID STREQUAL GNU)
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG")
    ELSEIF (CMAKE_CXX_COMPILER_ID STREQUAL PGI)
    ELSEIF (CMAKE_CXX_COMPILER_ID STREQUAL PathScale)
    ELSEIF (CMAKE_CXX_COMPILER_ID STREQUAL Intel)
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG")
      #SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx -DTEMP_INTEL_COMPILER_WORKAROUND_001")
    ELSEIF (CMAKE_CXX_COMPILER_ID STREQUAL XL)
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -DNDEBUG -qmaxmem=-1")
    ELSEIF (CMAKE_CXX_COMPILER_ID STREQUAL Cray)
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -DNDEBUG")
    ENDIF ()
  ENDIF ()

ENDIF ()

##############################################################################
# DEBUG flags
# 1) DEBUG_FLAGS if specified sets the Fortran,C, and CXX debug flags
# 2) DEBUG_FFLAGS if specified sets the Fortran debugflags
# 3) DEBUG_CFLAGS if specified sets the C debug flags
# 4) DEBUG_CXXFLAGS if specified sets the CXX debug flags
##############################################################################
IF (DEBUG_FLAGS)
  SET (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${DEBUG_FLAGS}")
  SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${DEBUG_FLAGS}")
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DEBUG_FLAGS}")
ELSE ()
  IF (DEBUG_FFLAGS)
    SET (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${DEBUG_FFLAGS}")
  ELSE ()
    IF(${CMAKE_Fortran_COMPILER_ID} STREQUAL PGI)
      SET (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -gopt")
    ELSEIF(NOT ${CMAKE_Fortran_COMPILER_ID} STREQUAL Cray)
      SET (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g")
    ENDIF ()
  ENDIF ()

  IF (DEBUG_CFLAGS)
    SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${DEBUG_CFLAGS}")
  ELSE ()
    IF(${CMAKE_Fortran_COMPILER_ID} STREQUAL PGI)
      SET (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -gopt")
    ELSE()
      SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")
    ENDIF()
  ENDIF ()

  IF (DEBUG_CXXFLAGS)
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DEBUG_CXXFLAGS}")
  ELSE ()
    IF(${CMAKE_Fortran_COMPILER_ID} STREQUAL PGI)
      SET (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -gopt")
    ELSE()
      SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
    ENDIF ()
  ENDIF ()

ENDIF ()

OPTION(DEBUG_TRACE "Enables TRACE level debugging checks. Very slow" FALSE)
IF (${DEBUG_TRACE})
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDEBUG_TRACE")
ENDIF ()

##############################################################################
# OpenMP
# Two flavors:
#   1) HORIZ_OPENMP OpenMP over elements (standard OPENMP)
#   2) COLUMN_OPENMP OpenMP within an element (previously called ELEMENT_OPENMP)
# COLUMN_OPENMP will be disabled by the openACC exectuables.
#
# Kokkos does not distinguish between the two because it does not use
# nested OpenMP. Nested OpenMP is the reason the two are distinguished in the
# Fortran code.
##############################################################################

OPTION(ENABLE_OPENMP "OpenMP across elements" TRUE)
OPTION(ENABLE_HORIZ_OPENMP "OpenMP across elements" TRUE)
OPTION(ENABLE_COLUMN_OPENMP "OpenMP within an element" FALSE)

# If OpenMP is turned off also turn off ENABLE_HORIZ_OPENMP
IF (NOT ${ENABLE_OPENMP})
  SET(ENABLE_HORIZ_OPENMP FALSE)
  SET(ENABLE_COLUMN_OPENMP FALSE)
ENDIF ()

##############################################################################
IF (ENABLE_HORIZ_OPENMP OR ENABLE_COLUMN_OPENMP)
  IF(NOT ${CMAKE_Fortran_COMPILER_ID} STREQUAL Cray)
    FIND_PACKAGE(OpenMP)
    IF(OPENMP_FOUND)
      MESSAGE(STATUS "Found OpenMP Flags")
      IF (CMAKE_Fortran_COMPILER_ID STREQUAL XL)
        SET(OpenMP_C_FLAGS "-qsmp=omp")
        IF (ENABLE_COLUMN_OPENMP)
          SET(OpenMP_C_FLAGS "-qsmp=omp:nested_par -qsuppress=1520-045:1506-793")
        ENDIF ()
      ENDIF ()
      # This file is needed for the timing library - this is currently
      # inaccessible from the timing CMake script
      SET(OpenMP_Fortran_FLAGS "${OpenMP_C_FLAGS}")
      MESSAGE(STATUS "OpenMP_Fortran_FLAGS: ${OpenMP_Fortran_FLAGS}")
      MESSAGE(STATUS "OpenMP_C_FLAGS: ${OpenMP_C_FLAGS}")
      MESSAGE(STATUS "OpenMP_CXX_FLAGS: ${OpenMP_CXX_FLAGS}")
      MESSAGE(STATUS "OpenMP_EXE_LINKER_FLAGS: ${OpenMP_EXE_LINKER_FLAGS}")
      # The fortran openmp flag should be the same as the C Flag
      SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")
      SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    ELSE ()
      MESSAGE(FATAL_ERROR "Unable to find OpenMP")
    ENDIF()
  ENDIF()
 IF (${ENABLE_HORIZ_OPENMP})
   # Set this as global so it can be picked up by all executables
#   SET(HORIZ_OPENMP TRUE CACHE BOOL "Threading in the horizontal direction")
   SET(HORIZ_OPENMP TRUE BOOL "Threading in the horizontal direction")
   MESSAGE(STATUS "  Using HORIZ_OPENMP")
 ENDIF ()

 IF (${ENABLE_COLUMN_OPENMP})
   # Set this as global so it can be picked up by all executables
#   SET(COLUMN_OPENMP TRUE CACHE BOOL "Threading in the horizontal direction")
   SET(COLUMN_OPENMP TRUE BOOL "Threading in the horizontal direction")
   MESSAGE(STATUS "  Using COLUMN_OPENMP")
 ENDIF ()
ENDIF ()
##############################################################################

##############################################################################
# Intel Phi (MIC) specific flags - only supporting the Intel compiler
##############################################################################
OPTION(ENABLE_INTEL_PHI "Whether to build with Intel Xeon Phi (MIC) support" FALSE)

IF (ENABLE_INTEL_PHI)
  IF (NOT ${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
    MESSAGE(FATAL_ERROR "Intel Phi acceleration only supported through the Intel compiler")
  ELSE ()
    SET(INTEL_PHI_FLAGS "-mmic")
    # TODO: think about changing the line above with the following commented one
    # SET(INTEL_PHI_FLAGS "-xMIC-AVX512")
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${INTEL_PHI_FLAGS}")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  ${INTEL_PHI_FLAGS}")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${INTEL_PHI_FLAGS}")
    # CMake magic for cross-compilation
    SET(CMAKE_SYSTEM_NAME Linux)
    SET(CMAKE_SYSTEM_PROCESSOR k1om)
    SET(CMAKE_SYSTEM_VERSION 1)
    SET(_CMAKE_TOOLCHAIN_PREFIX  x86_64-k1om-linux-)
    # Specify the location of the target environment
    IF (TARGET_ROOT_PATH)
      SET(CMAKE_FIND_ROOT_PATH ${TARGET_ROOT_PATH})
    ELSE ()
      SET(CMAKE_FIND_ROOT_PATH /usr/linux-k1om-4.7)
    ENDIF ()
  ENDIF ()
ENDIF ()

##############################################################################
# Allow the option to add compiler flags to those provided
##############################################################################
SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${ADD_Fortran_FLAGS}")
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${ADD_C_FLAGS}")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ADD_CXX_FLAGS}")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${ADD_LINKER_FLAGS}")

##############################################################################
# Allow the option to override compiler flags
##############################################################################
IF (FORCE_Fortran_FLAGS)
  SET(CMAKE_Fortran_FLAGS ${FORCE_Fortran_FLAGS})
ENDIF ()
IF (FORCE_C_FLAGS)
  SET(CMAKE_C_FLAGS ${FORCE_C_FLAGS})
ENDIF ()
IF (FORCE_CXX_FLAGS)
  SET(CMAKE_CXX_FLAGS ${FORCE_CXX_FLAGS})
ENDIF ()
IF (FORCE_LINKER_FLAGS)
  SET(CMAKE_EXE_LINKER_FLAGS ${FORCE_LINKER_FLAGS})
ENDIF ()
