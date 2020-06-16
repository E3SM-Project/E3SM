##############################################################################
# Compiler specific options
##############################################################################

# Small function to set the compiler flag '-fp model' given the name of the model
# Note: this is an Intel-only flag
function (EKAT_set_fpmodel_flags fpmodel_string flags)
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

macro (SetCompilerFlags)
  set (FP_MODEL_FLAG "")
  set (UT_FP_MODEL_FLAG "")
  if (DEFINED SCREAM_FPMODEL)
    EKAT_set_fpmodel_flags("${SCREAM_FPMODEL}" FP_MODEL_FLAG)
    string(TOLOWER "${SCREAM_FPMODEL}" fpmodel_string_lower)
    if (fpmodel_string_lower STREQUAL "strict")
      set (SCREAM_STRICT_FP TRUE PARENT_SCOPE BOOL)
    endif()
  elseif (CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
    set(${FP_MODEL_FLAG} "-fp-model precise")
    set(${UT_FP_MODEL_FLAG} "-fp-model precise")
  endif ()

  # enable all warning but disable vectorization remarks like "remark: simd loop has only one iteration"
  # since we would get hit with 1000's of those anytime we set packsize to 1.
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
  if (CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -diag-disable=remark")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -warn all -diag-disable=remark -fpscomp logicals")
  else()
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall")
  endif()

  if (DEFINED BASE_FFLAGS)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${BASE_FFLAGS}")
  else ()
    if (CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none ${FP_MODEL_FLAG}")
      add_definitions(-DCPRGNU)
    elseif (CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
      add_definitions(-DCPRINTEL)
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -assume byterecl")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FP_MODEL_FLAG}")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ftz")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${FP_MODEL_FLAG}")
    endif ()
  endif ()

  if (DEFINED BASE_CPPFLAGS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${BASE_CPPFLAGS}")
  else ()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FP_MODEL_FLAG}")
  endif ()

  # C++ Flags

  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -restrict")
  endif()

  STRING(TOUPPER "${PERFORMANCE_PROFILE}" PERF_PROF_UPPER)
  if ("${PERF_PROF_UPPER}" STREQUAL "VTUNE")
    add_definitions(-DVTUNE_PROFILE)
  elseif ("${PERF_PROF_UPPER}" STREQUAL "CUDA")
    add_definitions(-DCUDA_PROFILE)
  elseif ("${PERF_PROF_UPPER}" STREQUAL "GPROF")
    add_definitions(-DGPROF_PROFILE -pg)
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
  endif ()

  # Handle Cuda.
  find_package(CUDA QUIET)
  if (${CUDA_FOUND})
    # We found cuda, but we may be only interested in running on host.
    # Check if the compiler is not nvcc; if not, do not add cuda support
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} "--nvcc-wrapper-show"
      RESULT_VARIABLE WRAPS_NVCC
      OUTPUT_VARIABLE WRAPS_NVCC_OUT1
      ERROR_QUIET)

    # Need to check OMPI_CXX/MPICH_CXX (if user is using mpicxx)
    set (mpi_distro_name)
    GetMpiDistributionName(mpi_distro_name)
    if (NOT "${mpi_distro_name}" STREQUAL "unknown")
      set (mpi_cxx_backend_var_name)
      SetMpiCxxBackendCompilerVarName(mpi_cxx_backend_var_name)
      if (DEFINED ENV{${mpi_cxx_backend_var_name}})
        execute_process(COMMAND $ENV{${mpi_cxx_backend_var_name}} "--nvcc-wrapper-show"
          RESULT_VARIABLE WRAPS_NVCC
          OUTPUT_VARIABLE WRAPS_NVCC_OUT2
          ERROR_QUIET)
      endif()
      unset (mpi_cxx_backend_var_name)
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
  if (OPT_FLAGS)
    # Flags for Fortran C and CXX
    set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${OPT_FLAGS}")
    set (CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${OPT_FLAGS}")
    set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OPT_FLAGS}")

  else ()

    if (OPT_FFLAGS)
      # User specified optimization flags
      set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${OPT_FFLAGS}")
    else ()
      # Defaults
      if (CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
        set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O3")
      elseif (CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
        set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O3")
      endif ()
    endif ()

    if (OPT_CFLAGS)
      set (CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${OPT_CFLAGS}")
    else ()
      if (CMAKE_C_COMPILER_ID STREQUAL GNU)
        set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")
      elseif (CMAKE_C_COMPILER_ID STREQUAL Intel)
        set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")
      endif ()
    endif ()

    if (OPT_CXXFLAGS)
      set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OPT_CXXFLAGS}")
    else ()
      if (CMAKE_CXX_COMPILER_ID STREQUAL GNU)
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
      elseif (CMAKE_CXX_COMPILER_ID STREQUAL Intel)
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
      endif ()
    endif ()
  endif ()

  ##############################################################################
  # DEBUG flags
  # 1) DEBUG_FLAGS if specified sets the Fortran,C, and CXX debug flags
  # 2) DEBUG_FFLAGS if specified sets the Fortran debugflags
  # 3) DEBUG_CFLAGS if specified sets the C debug flags
  # 4) DEBUG_CXXFLAGS if specified sets the CXX debug flags
  ##############################################################################
  if (DEBUG_FLAGS)
    set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${DEBUG_FLAGS}")
    set (CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${DEBUG_FLAGS}")
    set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${DEBUG_FLAGS}")
  else ()
    if (DEBUG_FFLAGS)
      set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${DEBUG_FFLAGS}")
    else ()
      if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
        set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -no-vec")
      endif()
      set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -g -O0")
    endif ()

    if (DEBUG_CFLAGS)
      set (CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${DEBUG_CFLAGS}")
    else ()
      set (CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g -O0")
    endif ()

    if (DEBUG_CXXFLAGS)
      set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${DEBUG_CXXFLAGS}")
    else ()
      if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -no-vec")
      endif()
      set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")
    endif ()

  endif ()

  option(DEBUG_TRACE "Enables TRACE level debugging checks. Very slow" FALSE)
  if (${DEBUG_TRACE})
    set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG_TRACE")
  endif ()

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

  option(ENABLE_INTEL_PHI "Whether to build with Intel Xeon Phi (MIC) support" ${ENABLE_INTEL_PHI_DEFAULT})

  if (ENABLE_INTEL_PHI)
    if (NOT ${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
      message(FATAL_ERROR "Intel Phi acceleration only supported through the Intel compiler")
    else ()
      set(INTEL_PHI_FLAGS "-xMIC-AVX512")
      set(AVX_VERSION "512")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${INTEL_PHI_FLAGS}")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  ${INTEL_PHI_FLAGS}")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${INTEL_PHI_FLAGS}")
      set(IS_ACCELERATOR TRUE)
      # CMake magic for cross-compilation
      set(CMAKE_SYSTEM_NAME Linux)
      set(CMAKE_SYSTEM_PROCESSOR k1om)
      set(CMAKE_SYSTEM_VERSION 1)
      set(_CMAKE_TOOLCHAIN_PREFIX  x86_64-k1om-linux-)
      # Specify the location of the target environment
      if (TARGET_ROOT_PATH)
        set(CMAKE_FIND_ROOT_PATH ${TARGET_ROOT_PATH})
      else ()
        set(CMAKE_FIND_ROOT_PATH /usr/linux-k1om-4.7)
      endif ()
    endif ()
  endif ()
  ##############################################################################
  # Compiler FLAGS for AVX1 and AVX2 (CXX compiler only)
  ##############################################################################

  # NOTE: This won't work on batch machines where the architecture of the
  # interactive node is different than the compute nodes.
  if (NOT DEFINED AVX_VERSION)
    include(FindAVX)
    FindAVX()
    if (AVX512_FOUND)
      set(AVX_VERSION "512")
    elseif (AVX2_FOUND)
      set(AVX_VERSION "2")
    elseif (AVX_FOUND)
      set(AVX_VERSION "1")
    else ()
      set(AVX_VERSION "0")
    endif ()
  endif ()

  if (AVX_VERSION STREQUAL "512")
    if (NOT ENABLE_INTEL_PHI)
      if (CMAKE_CXX_COMPILER_ID STREQUAL Intel)
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xCORE-AVX512")
      endif()
    endif()
  elseif (AVX_VERSION STREQUAL "2")
    if (CMAKE_CXX_COMPILER_ID STREQUAL GNU)
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2")
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL Intel)
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xCORE-AVX2")
    endif()
  elseif (AVX_VERSION STREQUAL "1")
    if (CMAKE_CXX_COMPILER_ID STREQUAL GNU)
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx")
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL Intel)
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xAVX")
    endif()
  endif ()

  ##############################################################################
  # Allow the option to add compiler flags to those provided
  ##############################################################################
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${ADD_Fortran_FLAGS}")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${ADD_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ADD_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${ADD_LINKER_FLAGS}")

  ##############################################################################
  # Allow the option to override compiler flags
  ##############################################################################
  if (FORCE_Fortran_FLAGS)
    set(CMAKE_Fortran_FLAGS ${FORCE_Fortran_FLAGS})
  endif ()
  if (FORCE_C_FLAGS)
    set(CMAKE_C_FLAGS ${FORCE_C_FLAGS})
  endif ()
  if (FORCE_CXX_FLAGS)
    set(CMAKE_CXX_FLAGS ${FORCE_CXX_FLAGS})
  endif ()
  if (FORCE_LINKER_FLAGS)
    set(CMAKE_EXE_LINKER_FLAGS ${FORCE_LINKER_FLAGS})
  endif ()
endmacro()
