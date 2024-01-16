###########################
# Build Control Variables #
###########################

set(OMEGA_PROJECT_NAME            "OmegaOceanModel")
set(OMEGA_EXE_NAME                "omega.exe")
set(OMEGA_LIB_NAME                "OmegaLib")

set(OMEGA_BUILD_MODES             "E3SM" "STANDALONE" "NOT_DEFINED")
set(OMEGA_BUILD_MODE              NOT_DEFINED CACHE STRING "Omega build mode")
set_property(CACHE OMEGA_BUILD_MODE PROPERTY STRINGS ${OMEGA_BUILD_MODES})

set(OMEGA_SOURCE_DIR              ${CMAKE_CURRENT_LIST_DIR})
set(OMEGA_DEFAULT_BUILD_TYPE      Release) # Debug or Release

set(E3SM_CIME_ROOT ${OMEGA_SOURCE_DIR}/../../cime)
set(E3SM_CIMECONFIG_ROOT ${OMEGA_SOURCE_DIR}/../../cime_config)

###########################
# Macros                  #
###########################

# set build control variables used for both e3sm build and standalone build
macro(setup_common_variables)

  option(OMEGA_DEBUG "Turn on error message throwing (default OFF)." OFF)

  if(NOT DEFINED OMEGA_CXX_FLAGS)
    set(OMEGA_CXX_FLAGS "")
  endif()

  if(NOT DEFINED OMEGA_LINK_OPTIONS)
    set(OMEGA_LINK_OPTIONS "")
  endif()

endmacro()

# Collect machine and compiler info from CIME
macro(preset)

  find_package (Python COMPONENTS Interpreter)

  if(NOT Python_FOUND)
    message(FATAL_ERROR "Python is not available, CMake will exit." )
  endif()

  set(_TMP_CMAKE_FILE ${CMAKE_CURRENT_BINARY_DIR}/_Omega.cmake)
  set(_PY_OPTS "-p;${E3SM_CIME_ROOT};-o;${_TMP_CMAKE_FILE}")

  if(DEFINED OMEGA_CIME_COMPILER)
    list(APPEND _PY_OPTS "-c" "${OMEGA_CIME_COMPILER}")
  endif()

  if(DEFINED OMEGA_CIME_MACHINE)
    list(APPEND _PY_OPTS "-m" "${OMEGA_CIME_MACHINE}")
  endif()

  if(OMEGA_BUILD_TYPE STREQUAL "Debug")

    list(APPEND _PY_OPTS "-d")

    execute_process(COMMAND ${Python_EXECUTABLE} create_scripts.py ${_PY_OPTS}
      WORKING_DIRECTORY ${OMEGA_SOURCE_DIR}
      OUTPUT_VARIABLE _MACHINE_INFO
      ERROR_VARIABLE _ERROR_INFO)

      message(STATUS "create_scripts.py output: ${_MACHINE_INFO}")
      message(STATUS "create_scripts.py error: ${_ERROR_INFO}")

  else()

    execute_process(COMMAND ${Python_EXECUTABLE} create_scripts.py ${_PY_OPTS}
      OUTPUT_QUIET ERROR_QUIET
      WORKING_DIRECTORY ${OMEGA_SOURCE_DIR}
      OUTPUT_VARIABLE _MACHINE_INFO)

  endif()

  include(${_TMP_CMAKE_FILE})
  if(OMEGA_BUILD_TYPE STREQUAL "Release")
    file(REMOVE ${_TMP_CMAKE_FILE})
  endif()

  if(OMEGA_C_COMPILER)
    find_program(_OMEGA_C_COMPILER ${OMEGA_C_COMPILER})
  else()
    if (MPILIB STREQUAL "mpi-serial")
      find_program(_OMEGA_C_COMPILER ${SCC})
    else()
      find_program(_OMEGA_C_COMPILER ${MPICC})
    endif()
  endif()

  if(OMEGA_CXX_COMPILER)
    find_program(_OMEGA_CXX_COMPILER ${OMEGA_CXX_COMPILER})
  else()
    if (MPILIB STREQUAL "mpi-serial")
      find_program(_OMEGA_CXX_COMPILER ${SCXX})
    else()
      find_program(_OMEGA_CXX_COMPILER ${MPICXX})
    endif()
  endif()

  if(OMEGA_Fortran_COMPILER)
    find_program(_OMEGA_Fortran_COMPILER ${OMEGA_Fortran_COMPILER})
  else()
    if (MPILIB STREQUAL "mpi-serial")
      find_program(_OMEGA_Fortran_COMPILER ${SFC})
    else()
      find_program(_OMEGA_Fortran_COMPILER ${MPIFC})
    endif()
  endif()

  set(OMEGA_C_COMPILER ${_OMEGA_C_COMPILER})
  set(CMAKE_C_COMPILER ${OMEGA_C_COMPILER})

  set(OMEGA_CXX_COMPILER ${_OMEGA_CXX_COMPILER})
  set(CMAKE_CXX_COMPILER ${OMEGA_CXX_COMPILER})

  set(OMEGA_Fortran_COMPILER ${_OMEGA_Fortran_COMPILER})
  set(CMAKE_Fortran_COMPILER ${OMEGA_Fortran_COMPILER})

  message(STATUS "OMEGA_C_COMPILER = ${OMEGA_C_COMPILER}")
  message(STATUS "OMEGA_CXX_COMPILER = ${OMEGA_CXX_COMPILER}")
  message(STATUS "OMEGA_Fortran_COMPILER = ${OMEGA_Fortran_COMPILER}")
  message(STATUS "MPI_EXEC = ${MPI_EXEC}")

endmacro()

# set build-control-variables for standalone build
macro(setup_standalone_build)

  setup_common_variables()

  if(NOT DEFINED OMEGA_BUILD_TYPE)
    set(OMEGA_BUILD_TYPE ${OMEGA_DEFAULT_BUILD_TYPE})
  endif()

  if( EXISTS ${OMEGA_SOURCE_DIR}/../../components AND
      EXISTS ${OMEGA_SOURCE_DIR}/../../cime AND
      EXISTS ${OMEGA_SOURCE_DIR}/../../cime_config AND
      EXISTS ${OMEGA_SOURCE_DIR}/../../externals)

    set(E3SM_SOURCE_DIR ${OMEGA_SOURCE_DIR}/../../components)
    set(E3SM_EXTERNALS_ROOT ${OMEGA_SOURCE_DIR}/../../externals)

  else()
    # so far, we assume that Omega exists inside of E3SM.
    # However, we leave this else part for later usage.

  endif()

  set(OMEGA_BUILD_MODE "STANDALONE")
  set(OMEGA_BUILD_EXECUTABLE ON)

endmacro()

# set build-control-variables for e3sm build
macro(setup_e3sm_build)

  setup_common_variables()

  set(OMEGA_BUILD_TYPE ${E3SM_DEFAULT_BUILD_TYPE})
  set(E3SM_EXTERNALS_ROOT ${E3SM_SOURCE_DIR}/../externals)

  set(OMEGA_CXX_COMPILER ${CMAKE_CXX_COMPILER})

  #TODO: set OMEGA_ARCH according to E3SM variables
  set(OMEGA_ARCH "NOT_DEFINED")
  set(OMEGA_BUILD_MODE "E3SM")

  message(STATUS "OMEGA_CXX_COMPILER = ${OMEGA_CXX_COMPILER}")

endmacro()


################################
# Set cmake and YAKL variables #
################################
macro(update_variables)

  # Set the build type
  set(CMAKE_BUILD_TYPE ${OMEGA_BUILD_TYPE})

  # Set compiler and linker flags
  if (CXXFLAGS)
    separate_arguments(_CXXFLAGS NATIVE_COMMAND ${CXXFLAGS})
    list(APPEND OMEGA_CXX_FLAGS ${_CXXFLAGS})
  endif()

  if (LDFLAGS)
    separate_arguments(_LDFLAGS NATIVE_COMMAND ${LDFLAGS})
    list(APPEND OMEGA_LINK_OPTIONS ${_LDFLAGS})
  endif()

  if (SLIBS)
    separate_arguments(_SLIBS NATIVE_COMMAND ${SLIBS})
    list(APPEND OMEGA_LINK_OPTIONS ${_SLIBS})
  endif()

  # check if MPI is supported
  string(CONCAT _TestMPISource
    "#include \"mpi.h\"\n"
    "int main(int argc, char* argv[])\n"
    "{MPI_Init(&argc, &argv)\; return 0\;}\n")
  set(_TestMPISrcFile ${CMAKE_CURRENT_BINARY_DIR}/_testMPI.cpp)
  set(_TestMPIObjFile ${CMAKE_CURRENT_BINARY_DIR}/_testMPI.o)
  file(WRITE ${_TestMPISrcFile}  ${_TestMPISource})

  execute_process(
    COMMAND ${OMEGA_CXX_COMPILER} -c ${_TestMPISrcFile} -o ${_TestMPIObjFile}
    OUTPUT_QUIET ERROR_QUIET
    RESULT_VARIABLE _MPI_TEST_RESULT
    OUTPUT_VARIABLE _MPI_TEST_OUTPUT
    ERROR_VARIABLE _MPI_TEST_ERROR)

  if(OMEGA_BUILD_TYPE EQUAL Release)
    file(REMOVE ${_TestMPISrcFile})
    file(REMOVE ${_TestMPIObjFile})
  endif()

  if (NOT _MPI_TEST_RESULT EQUAL 0)
    if (_MPI_TEST_RESULT MATCHES "^[-]?[0-9]+$")
      find_package(MPI)
      if(MPI_CXX_FOUND)
        list(APPEND OMEGA_CXX_FLAGS "-I${MPI_CXX_INCLUDE_DIRS}")
        list(APPEND OMEGA_LINK_OPTIONS
          "-L${MPI_CXX_INCLUDE_DIRS}/../lib" "-lmpi"
        )
      else()
        message(FATAL_ERROR "MPI is not found" )
      endif()
    else()
      message(FATAL_ERROR "MPI test failure: ${_MPI_TEST_RESULT}" )
    endif()
  endif()

  message(STATUS "OMEGA_CXX_FLAGS           = ${OMEGA_CXX_FLAGS}")
  message(STATUS "OMEGA_LINK_OPTIONS        = ${OMEGA_LINK_OPTIONS}")

  if(OMEGA_INSTALL_PREFIX)
    set(CMAKE_INSTALL_PREFIX ${OMEGA_INSTALL_PREFIX})
  endif()

  # Check if CUDA or HIP is supported
  if((NOT DEFINED OMEGA_ARCH) OR (OMEGA_ARCH STREQUAL "NOT_DEFINED"))

    if(USE_CUDA)
      set(OMEGA_ARCH "CUDA")

    elseif(USE_HIP)
      set(OMEGA_ARCH "HIP")

    else()

      execute_process(
        COMMAND ${OMEGA_CXX_COMPILER} --version
        RESULT_VARIABLE _CXX_VER_RESULT
        OUTPUT_VARIABLE _CXX_VER_OUTPUT)

      if (_CXX_VER_RESULT EQUAL 0)

        string(REGEX MATCH "HIP|hip"       _HIP_CHECK "${_CXX_VER_OUTPUT}")
        string(REGEX MATCH "AMD|amd"       _AMD_CHECK "${_CXX_VER_OUTPUT}")
        string(REGEX MATCH "NVCC|nvcc"     _NVCC_CHECK "${_CXX_VER_OUTPUT}")
        string(REGEX MATCH "NVIDIA|nvidia" _NVIDIA_CHECK "${_CXX_VER_OUTPUT}")

        if(_HIP_CHECK AND _AMD_CHECK)
          set(OMEGA_ARCH "HIP")

        elseif(_NVCC_CHECK AND _NVIDIA_CHECK)
          set(OMEGA_ARCH "CUDA")

        else()
          set(OMEGA_ARCH "")

        endif()
      else()
        set(OMEGA_ARCH "")

      endif()
    endif()
  endif()

  if(OMEGA_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "OMEGA_ARCH = ${OMEGA_ARCH}")
  endif()

  if(OMEGA_ARCH STREQUAL "CUDA")
    set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})
    find_program(CMAKE_CUDA_COMPILER "nvcc")

    if(OMEGA_BUILD_TYPE STREQUAL "Debug")
      message(STATUS "CMAKE_CUDA_COMPILER = ${CMAKE_CUDA_COMPILER}")
      message(STATUS "CMAKE_CUDA_HOST_COMPILER = ${CMAKE_CUDA_HOST_COMPILER}")
    endif()

  elseif(OMEGA_ARCH STREQUAL "HIP")
    set(CMAKE_HIP_HOST_COMPILER ${CMAKE_CXX_COMPILER})
    find_program(CMAKE_HIP_COMPILER "hipcc")

    if(OMEGA_BUILD_TYPE STREQUAL "Debug")
      message(STATUS "CMAKE_HIP_COMPILER = ${CMAKE_HIP_COMPILER}")
      message(STATUS "CMAKE_HIP_HOST_COMPILER = ${CMAKE_HIP_HOST_COMPILER}")
    endif()

  endif()

  set(YAKL_ARCH "${OMEGA_ARCH}")

  if(YAKL_ARCH)

      if(OMEGA_${YAKL_ARCH}_FLAGS)
        set(YAKL_${YAKL_ARCH}_FLAGS ${OMEGA_${YAKL_ARCH}_FLAGS})
      endif()

  endif()

  # Include the findParmetis script
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")
  find_package(Parmetis REQUIRED)

#  # prints generates all cmake variables
#  get_cmake_property(_variableNames VARIABLES)
#  list (SORT _variableNames)
#  foreach (_variableName ${_variableNames})
#      message(STATUS "${_variableName}=${${_variableName}}")
#  endforeach()

endmacro()



################################
# Verify variable integrity    #
################################
macro(check_setup)

  #message("OMEGA_BUILD_MODE = ${OMEGA_BUILD_MODE}")

  if(OMEGA_BUILD_MODE STREQUAL "E3SM")
    message(STATUS "*** Omega E3SM-component Build ***")

  elseif(${OMEGA_BUILD_MODE} STREQUAL "STANDALONE")
    message(STATUS "*** Omega Standalone Build ***")

  else()

    message(FATAL_ERROR "OMEGA_BUILD_MODE is neither E3SM nor STANDALONE.")

  endif()

  if (NOT DEFINED YAKL_ARCH)
    message(FATAL_ERROR "YAKL_ARCH is not defined.")
  endif()

endmacro()


################################
# Prepare output               #
################################
macro(wrap_outputs)

  if(OMEGA_INSTALL_PREFIX)

    install(TARGETS ${OMEGA_LIB_NAME}
      LIBRARY DESTINATION "${OMEGA_INSTALL_PREFIX}/lib"
    )

    if(OMEGA_BUILD_EXECUTABLE)
      install(TARGETS ${OMEGA_EXE_NAME}
        RUNTIME DESTINATION "${OMEGA_INSTALL_PREFIX}/bin"
      )
    endif()

  endif()

endmacro()
