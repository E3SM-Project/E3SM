###########################
# Internal variables      #
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
# Public variables        #
###########################
macro(setup_common_variables)

  option(OMEGA_DEBUG "Turn on error message throwing (default OFF)." OFF)

  if(NOT DEFINED OMEGA_CXX_FLAGS)
    set(OMEGA_CXX_FLAGS "")
  endif()

  if(NOT DEFINED OMEGA_INCLUDE_DIRECTORIES)
    set(OMEGA_INCLUDE_DIRECTORIES "")
  endif()

  if(NOT DEFINED OMEGA_LINK_DIRECTORIES)
    set(OMEGA_LINK_DIRECTORIES "")
  endif()

  if(NOT DEFINED OMEGA_LINK_OPTIONS)
    set(OMEGA_LINK_OPTIONS "")
  endif()

endmacro()

###########################
# Preset Standalone build #
###########################
macro(preset)

  find_package (Python COMPONENTS Interpreter)

  if(NOT Python_FOUND)
    message(FATAL_ERROR "Python is not available, CMake will exit." )
  endif()

  set(_TMP_CMAKE_FILE ${CMAKE_CURRENT_BINARY_DIR}/_Omega.cmake)

  if(DEFINED OMEGA_CIME_COMPILER)
    set(_TMP_COMPILER "-c ${OMEGA_CIME_COMPILER}")
  else()
    set(_TMP_COMPILER "")
  endif()

  execute_process(COMMAND ${Python_EXECUTABLE} gen_machine_info.py
    -p ${E3SM_CIME_ROOT} -o ${_TMP_CMAKE_FILE} ${_TMP_COMPILER}
    OUTPUT_QUIET ERROR_QUIET
    WORKING_DIRECTORY ${OMEGA_SOURCE_DIR}
    OUTPUT_VARIABLE _MACHINE_INFO)

  include(${_TMP_CMAKE_FILE})
  file(REMOVE ${_TMP_CMAKE_FILE})

  if(NOT OMEGA_CXX_COMPILER)
    if (MPILIB STREQUAL "mpi-serial")
      set(OMEGA_CXX_COMPILER ${SCXX})
    else()
      set(OMEGA_CXX_COMPILER ${MPICXX})
    endif()
  endif()

  set(CMAKE_CXX_COMPILER ${OMEGA_CXX_COMPILER})
  message(STATUS "OMEGA_CXX_COMPILER = ${OMEGA_CXX_COMPILER}")

endmacro()

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

macro(setup_e3sm_build)

  setup_common_variables()

  set(OMEGA_BUILD_TYPE ${E3SM_DEFAULT_BUILD_TYPE})
  set(E3SM_EXTERNALS_ROOT ${E3SM_SOURCE_DIR}/../externals)


  #TODO: set OMEGA_ARCH according to E3SM variables
  set(OMEGA_ARCH "NOT_DEFINED")
  set(OMEGA_BUILD_MODE "E3SM")

endmacro()


################################
# Set cmake and YAKL variables #
################################
macro(update_variables)

  # Set the build type
  set(CMAKE_BUILD_TYPE ${OMEGA_BUILD_TYPE})

  if (NETCDF_PATH)
    list(APPEND OMEGA_INCLUDE_DIRECTORIES "${NETCDF_PATH}/include")
    list(APPEND OMEGA_LINK_DIRECTORIES "${NETCDF_PATH}/lib")
    list(APPEND OMEGA_LINK_OPTIONS "-lnetcdf")
  endif()

  if (PNETCDF_PATH)
    list(APPEND OMEGA_INCLUDE_DIRECTORIES "${PNETCDF_PATH}/include")
    list(APPEND OMEGA_LINK_DIRECTORIES "${PNETCDF_PATH}/lib")
    list(APPEND OMEGA_LINK_OPTIONS "-lpnetcdf")
  endif()

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
    RESULT_VARIABLE _MPI_TEST_RESULT)

  if(OMEGA_BUILD_TYPE EQUAL Release)
	file(REMOVE ${_TestMPISrcFile})
	file(REMOVE ${_TestMPIObjFile})
  endif()

  if (NOT _MPI_TEST_RESULT EQUAL 0)
    if (_MPI_TEST_RESULT MATCHES "^[-]?[0-9]+$")
      find_package(MPI)
      if(MPI_CXX_FOUND)
        list(APPEND OMEGA_INCLUDE_DIRECTORIES "${MPI_CXX_INCLUDE_DIRS}")
        list(APPEND OMEGA_LINK_DIRECTORIES "${MPI_CXX_INCLUDE_DIRS}/../lib")
        list(APPEND OMEGA_LINK_OPTIONS "-lmpi")
      else()
        message(FATAL_ERROR "MPI is not found" )
      endif()
    else()
      message(FATAL_ERROR "MPI test failure: ${_MPI_TEST_RESULT}" )
    endif()
  endif()

  if(USE_CUDA)
    list(APPEND OMEGA_CXX_FLAGS ${CUDA_FLAGS})

  elseif(USE_HIP)
    list(APPEND OMEGA_CXX_FLAGS ${HIP_FLAGS})

  elseif(USE_SYCL)
    list(APPEND OMEGA_CXX_FLAGS ${SYCL_FLAGS})

  else()
    list(APPEND OMEGA_CXX_FLAGS ${CXXFLAGS})

    if(CXX_LIBS)
      list(APPEND OMEGA_LINK_OPTIONS  ${CXX_LIBS})
    endif()

  endif()

  message(STATUS "OMEGA_CXX_FLAGS           = ${CMAKE_CXX_FLAGS}")
  message(STATUS "OMEGA_INCLUDE_DIRECTORIES = ${OMEGA_INCLUDE_DIRECTORIES}")
  message(STATUS "OMEGA_LINK_DIRECTORIES    = ${OMEGA_LINK_DIRECTORIES}")
  message(STATUS "OMEGA_LINK_OPTIONS        = ${OMEGA_LINK_OPTIONS}")

  if(OMEGA_INSTALL_PREFIX)
    set(CMAKE_INSTALL_PREFIX ${OMEGA_INSTALL_PREFIX})
  endif()

  if(DEFINED OMEGA_ARCH)

    if(OMEGA_ARCH STREQUAL "NOT_DEFINED")
      set(YAKL_ARCH "")

    else()
      set(YAKL_ARCH ${OMEGA_ARCH})

      if(OMEGA_${YAKL_ARCH}_FLAGS)
        set(YAKL_${YAKL_ARCH}_FLAGS ${OMEGA_${YAKL_ARCH}_FLAGS})
      endif()

    endif()

  else()

    set(YAKL_ARCH "")

  endif()

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

    install(TARGETS ${OMEGA_LIB_NAME} LIBRARY DESTINATION "${OMEGA_INSTALL_PREFIX}/lib")

    if(OMEGA_BUILD_EXECUTABLE)
      install(TARGETS ${OMEGA_EXE_NAME} RUNTIME DESTINATION "${OMEGA_INSTALL_PREFIX}/bin")
    endif()

  endif()

endmacro()
