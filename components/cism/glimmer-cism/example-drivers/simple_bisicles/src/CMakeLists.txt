# simple_bisicles build
# Get libraries for link line from Trilinos build information
# 

IF (NOT ${NO_TRILINOS})
  set(CISM_TRILINOS_LIBS glimmercismcpp ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}
     ${Trilinos_EXTRA_LD_FLAGS} ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
ELSE()
  set(CISM_TRILINOS_LIBS ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES} )
ENDIF()

#message("")
#message("  CMake detected the following libraries for linking Fortran with C++ compiler:")
#message("     ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES} ")

# Need include directories from Trilinos but also mod files from glimmer
include_directories (${GLIMMER_BINARY_DIR}/include  ${PYTHON_INC_DIR}
                     ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})

link_directories (${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS}
                  ${GLIMMER_DYCORE_DIR} ${GLIMMER_BISICLES_DIR} ${BISICLES_INTERFACE_DIR}
                  ${CISM_HDF5_LIB_DIR} ${PYTHON_LIB_DIR} )

# These are local source files needed to make  simple_bisicles  executable
add_executable(simple_bisicles  simple_bisicles.F90 simple_forcing.F90 testsfg.F90)

MESSAGE("GLIMMER_BISICLES_DIR: " ${GLIMMER_BISICLES_DIR})
MESSAGE("CISM_HDF5_LIBS: " ${CISM_HDF5_LIBS})

# Executable depends on several glimmer libraries and Trilinos,
# and potentially an f90main.o file (on jaguar with PGI at least)
target_link_libraries(simple_bisicles
      ${GLIMMER_FMAIN}
      glimmercismfortran
      DyCoreToGlimmer
      libBisicles.a
      libChomboLibs.a
      ${PYTHON_LIBS}
      ${GLIMMER_NETCDF_LIBS}
      ${CISM_HDF5_LIBS}
      ${GLIMMER_MPI_LIBS}
      ${CISM_TRILINOS_LIBS}
      ${GLIMMER_EXTRA_LIBS}
      )



# Helpful(?) message near end of configuration step
MESSAGE("\n   Execute 'make -j 8'")
MESSAGE("   Executable  simple_bisicles  should appear in dir: build/example-drivers/simple_bisicles/src")
MESSAGE("")
