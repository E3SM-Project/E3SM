# Module used for CIME testing.
#
# This module does some initial setup that must be done BEFORE the 'project'
# line in the main CMakeLists.txt file.

include(${CMAKE_BINARY_DIR}/Macros.cmake RESULT_VARIABLE FOUND)
if(NOT FOUND)
  message(FATAL_ERROR "You must generate a Macros.cmake file using CIME's configure")
endif()
if(MPILIB STREQUAL "mpi-serial")
  set(CMAKE_C_COMPILER ${SCC})
  set(CMAKE_Fortran_COMPILER ${SFC})
else()
  set(CMAKE_C_COMPILER ${MPICC})
  set(CMAKE_Fortran_COMPILER ${MPIFC})
endif()
