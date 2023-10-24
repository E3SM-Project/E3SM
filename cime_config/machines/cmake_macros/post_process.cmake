# Error checking?

if (MPILIB STREQUAL "mpi-serial")
  set(CMAKE_Fortran_COMPILER ${SFC})
  set(CMAKE_C_COMPILER ${SCC})
  set(CMAKE_CXX_COMPILER ${SCXX})
else()
  set(CMAKE_Fortran_COMPILER ${MPIFC})
  set(CMAKE_C_COMPILER ${MPICC})
  set(CMAKE_CXX_COMPILER ${MPICXX})
endif()
