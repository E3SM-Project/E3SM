##############################################################################
# Compiler specific options
##############################################################################

# Small function to set the compiler flag '-fp model' given the name of the model
# Note: this is an Intel-only flag
function (SCREAM_set_fpmodel_flags fpmodel_string flags)
  string(TOLOWER "${fpmodel_string}" fpmodel_string_lower)
  if (("${fpmodel_string_lower}" STREQUAL "precise") OR
      ("${fpmodel_string_lower}" STREQUAL "strict") OR
      ("${fpmodel_string_lower}" STREQUAL "fast") OR
      ("${fpmodel_string_lower}" STREQUAL "fast=1") OR
      ("${fpmodel_string_lower}" STREQUAL "fast=2"))
    if (CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
      set (${flags} "-fp-model ${fpmodel_string_lower}" PARENT_SCOPE)
    elseif (CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
      if ("${fpmodel_string_lower}" STREQUAL "strict")
        set (${flags} "-ffp-contract=off" PARENT_SCOPE)
      endif ()
    endif ()
  elseif ("${fpmodel_string_lower}" STREQUAL "")
    set (${flags} "" PARENT_SCOPE)
  else()
    message(FATAL_ERROR "FP_MODEL_FLAG string '${fpmodel_string}' is not recognized.")
  endif()
endfunction()

set (FP_MODEL_FLAG "")
set (UT_FP_MODEL_FLAG "")
IF (DEFINED SCREAM_FPMODEL)
  SCREAM_set_fpmodel_flags("${SCREAM_FPMODEL}" FP_MODEL_FLAG)
  string(TOLOWER "${SCREAM_FPMODEL}" fpmodel_string_lower)
  if (fpmodel_string_lower STREQUAL "strict")
    set (SCREAM_STRICT_FP TRUE PARENT_SCOPE BOOL)
  endif()
ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  SET(${FP_MODEL_FLAG} "-fp-model precise")
  SET(${UT_FP_MODEL_FLAG} "-fp-model precise")
ENDIF ()

# enable all warning but disable vectorization remarks like "remark: simd loop has only one iteration"
# since we would get hit with 1000's of those anytime we set packsize to 1.
SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
IF (CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -diag-disable=remark")
  SET (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -warn all -diag-disable=remark")
ELSE()
  SET (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall")
ENDIF()

IF (DEFINED BASE_FFLAGS)
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${BASE_FFLAGS}")
ELSE ()
  IF (CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none ${FP_MODEL_FLAG}")
    add_definitions(-DCPRGNU)
  ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
    add_definitions(-DCPRINTEL)
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -assume byterecl")
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FP_MODEL_FLAG}")
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ftz")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${FP_MODEL_FLAG}")
  ENDIF ()
ENDIF ()

IF (DEFINED BASE_CPPFLAGS)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${BASE_CPPFLAGS}")
ELSE ()
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FP_MODEL_FLAG}")
ENDIF ()

# C++ Flags

IF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -restrict")
ENDIF()

STRING(TOUPPER "${PERFORMANCE_PROFILE}" PERF_PROF_UPPER)
IF ("${PERF_PROF_UPPER}" STREQUAL "VTUNE")
  ADD_DEFINITIONS(-DVTUNE_PROFILE)
ELSEIF ("${PERF_PROF_UPPER}" STREQUAL "CUDA")
  ADD_DEFINITIONS(-DCUDA_PROFILE)
ELSEIF ("${PERF_PROF_UPPER}" STREQUAL "GPROF")
  ADD_DEFINITIONS(-DGPROF_PROFILE -pg)
  SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
ENDIF ()

# Handle Cuda.
find_package(CUDA QUIET)
if (${CUDA_FOUND})
  # We found cuda, but we may be only interested in running on host.
  # Check if the compiler is not nvcc; if not, do not add cuda support
  execute_process(COMMAND ${CMAKE_CXX_COMPILER} "--nvcc-wrapper-show"
    RESULT_VARIABLE WRAPS_NVCC
    OUTPUT_VARIABLE WRAPS_NVCC_OUT1
    ERROR_QUIET)

  # Need to check OMPI_CXX if user is using mpicxx
  if (DEFINED ENV{OMPI_CXX})
    execute_process(COMMAND $ENV{OMPI_CXX} "--nvcc-wrapper-show"
      RESULT_VARIABLE WRAPS_NVCC
      OUTPUT_VARIABLE WRAPS_NVCC_OUT2
      ERROR_QUIET)
  endif()

  string (FIND "${WRAPS_NVCC_OUT1} ${WRAPS_NVCC_OUT2}" "nvcc" pos)
  if (${pos} GREATER -1)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --expt-extended-lambda")
    # Turn off fused multiply add for debug so we can stay BFB with host
    set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} --fmad=false")
    message (STATUS "Cuda enabled!")
  else ()
    message (STATUS "Cuda was found, but the C++ compiler is not nvcc_wrapper, so building without Cuda support.")
  endif ()
endif ()

##############################################################################
# Optimization flags
# 1) OPT_FLAGS if specified sets the Fortran,C, and CXX optimization flags
# 2) OPT_FFLAGS if specified sets the Fortran optimization flags
# 3) OPT_CFLAGS if specified sets the C optimization flags
# 4) OPT_CXXFLAGS if specified sets the CXX optimization flags
##############################################################################
IF (OPT_FLAGS)
  # Flags for Fortran C and CXX
  SET (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${OPT_FLAGS}")
  SET (CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${OPT_FLAGS}")
  SET (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OPT_FLAGS}")

ELSE ()

  IF (OPT_FFLAGS)
    # User specified optimization flags
    SET (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${OPT_FFLAGS}")
  ELSE ()
    # Defaults
    IF (CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
      SET(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O3")
    ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
      SET(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O3")
    ENDIF ()
  ENDIF ()

  IF (OPT_CFLAGS)
    SET (CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${OPT_CFLAGS}")
  ELSE ()
    IF (CMAKE_C_COMPILER_ID STREQUAL GNU)
      SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")
    ELSEIF (CMAKE_C_COMPILER_ID STREQUAL Intel)
      SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")
    ENDIF ()
  ENDIF ()

  IF (OPT_CXXFLAGS)
    SET (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OPT_CXXFLAGS}")
  ELSE ()
    IF (CMAKE_CXX_COMPILER_ID STREQUAL GNU)
      SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
    ELSEIF (CMAKE_CXX_COMPILER_ID STREQUAL Intel)
      SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
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
  SET (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${DEBUG_FLAGS}")
  SET (CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${DEBUG_FLAGS}")
  SET (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${DEBUG_FLAGS}")
ELSE ()
  IF (DEBUG_FFLAGS)
    SET (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${DEBUG_FFLAGS}")
  ELSE ()
    IF ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
      SET(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -no-vec")
    ENDIF()
    SET (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -g")
  ENDIF ()

  IF (DEBUG_CFLAGS)
    SET (CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${DEBUG_CFLAGS}")
  ELSE ()
    SET (CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g")
  ENDIF ()

  IF (DEBUG_CXXFLAGS)
    SET (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${DEBUG_CXXFLAGS}")
  ELSE ()
    IF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
      set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -no-vec")
    ENDIF()
    SET (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
  ENDIF ()

ENDIF ()

OPTION(DEBUG_TRACE "Enables TRACE level debugging checks. Very slow" FALSE)
IF (${DEBUG_TRACE})
  SET (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG_TRACE")
ENDIF ()

##############################################################################
# OpenMP
##############################################################################
string(FIND "${KOKKOS_GMAKE_DEVICES}" "OpenMP" openmp_str_pos)
if (${openmp_str_pos} GREATER -1)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    message(STATUS "Found OpenMP Flags")

    message(STATUS "OpenMP_Fortran_FLAGS: ${OpenMP_Fortran_FLAGS}")
    message(STATUS "OpenMP_C_FLAGS: ${OpenMP_C_FLAGS}")
    message(STATUS "OpenMP_CXX_FLAGS: ${OpenMP_CXX_FLAGS}")
    message(STATUS "OpenMP_EXE_LINKER_FLAGS: ${OpenMP_EXE_LINKER_FLAGS}")
    # The fortran openmp flag should be the same as the C Flag
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  else()
    message(FATAL_ERROR "Unable to find OpenMP")
  endif()
endif()

##############################################################################
# Intel Phi (MIC) specific flags - only supporting the Intel compiler
##############################################################################

# If kokkos thinks the archicture is KNL, we should probably have enable-phi on by default.
if ("${KOKKOS_GMAKE_ARCH}" STREQUAL "KNL")
  set(ENABLE_INTEL_PHI_DEFAULT TRUE)
else()
  set(ENABLE_INTEL_PHI_DEFAULT FALSE)
endif()

OPTION(ENABLE_INTEL_PHI "Whether to build with Intel Xeon Phi (MIC) support" ${ENABLE_INTEL_PHI_DEFAULT})

IF (ENABLE_INTEL_PHI)
  IF (NOT ${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
    MESSAGE(FATAL_ERROR "Intel Phi acceleration only supported through the Intel compiler")
  ELSE ()
    SET(INTEL_PHI_FLAGS "-xMIC-AVX512")
    SET(AVX_VERSION "512")
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${INTEL_PHI_FLAGS}")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  ${INTEL_PHI_FLAGS}")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${INTEL_PHI_FLAGS}")
    SET(IS_ACCELERATOR TRUE)
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

# NOTE: This won't work on batch machines where the architecture of the
# interactive node is different than the compute nodes.
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
