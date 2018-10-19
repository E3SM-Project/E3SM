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
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
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
  CHECK_CXX_COMPILER_FLAG("-std=c++11" CXX11_SUPPORTED)
  IF (${CXX11_SUPPORTED})
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  ELSEIF (${HOMME_USE_KOKKOS})
    MESSAGE (FATAL_ERROR "Kokkos needs C++11, but the C++ compiler does not support it.")
  ENDIF ()
  CHECK_CXX_COMPILER_FLAG("-cxxlib" CXXLIB_SUPPORTED)

  # Handle Cuda.
  FIND_PACKAGE(CUDA QUIET)
  IF (${CUDA_FOUND})
    STRING (FIND ${CMAKE_CXX_COMPILER} "nvcc" pos)
    IF (${pos} GREATER -1)
      SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --expt-extended-lambda -DCUDA_BUILD")
    ELSE ()
      MESSAGE ("Cuda was found, but the C++ compiler is not nvcc_wrapper, so building without Cuda support.")
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
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
    ELSEIF (CMAKE_CXX_COMPILER_ID STREQUAL PGI)
    ELSEIF (CMAKE_CXX_COMPILER_ID STREQUAL PathScale)
    ELSEIF (CMAKE_CXX_COMPILER_ID STREQUAL Intel)
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
      #SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx -DTEMP_INTEL_COMPILER_WORKAROUND_001")
    ELSEIF (CMAKE_CXX_COMPILER_ID STREQUAL XL)
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -qmaxmem=-1")
    ELSEIF (CMAKE_CXX_COMPILER_ID STREQUAL Cray)
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
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
    SET(AVX_VERSION "512")
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
# Compiler FLAGS for AVX1 and AVX2 (CXX compiler only)
##############################################################################
IF (${HOMME_USE_CXX})
  IF (NOT DEFINED AVX_VERSION)
    INCLUDE(FindAVX)
    FindAVX()
    IF (AVX512_FOUND)
      SET(AVX_VERSION "512")
    ELSEIF (AVX2_FOUND)
      SET(AVX_VERSION "2")
    ELSEIF (AVX_FOUND)
      SET(AVX_VERSION "1")
    ELSE ()
      SET(AVX_VERSION "0")
    ENDIF ()
  ENDIF ()

  IF (AVX_VERSION STREQUAL "512")
    IF (NOT ENABLE_INTEL_PHI)
      IF (CMAKE_CXX_COMPILER_ID STREQUAL Intel)
        SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xCORE-AVX512")
      ENDIF()
    ENDIF()
  ELSEIF (AVX_VERSION STREQUAL "2")
    IF (CMAKE_CXX_COMPILER_ID STREQUAL GNU)
      SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2")
    ELSEIF (CMAKE_CXX_COMPILER_ID STREQUAL Intel)
      SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xCORE-AVX2")
    ENDIF()
  ELSEIF (AVX_VERSION STREQUAL "1")
    IF (CMAKE_CXX_COMPILER_ID STREQUAL GNU)
      SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx")
    ELSEIF (CMAKE_CXX_COMPILER_ID STREQUAL Intel)
      SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xAVX")
    ENDIF()
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
