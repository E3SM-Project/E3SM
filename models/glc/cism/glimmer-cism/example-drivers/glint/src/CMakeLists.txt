# glint_example build
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
include_directories (${GLIMMER_BINARY_DIR}/include)

# These are local source files needed to make  simple_bisicles  executable
add_executable(glint_example glint_example.F90)

# Executable depends on several glimmer libraries and Trilinos,
# and potentially an f90main.o file (on jaguar with PGI at least)
target_link_libraries(glint_example
      glimmercismfortran
      ${GLIMMER_NETCDF_LIBS}
      ${GLIMMER_GPTL_LIB}
      ${CISM_TRILINOS_LIBS}
      ${GLIMMER_EXTRA_LIBS})

# Helpful(?) message near end of configuration step
MESSAGE("\n   Execute 'make -j 8'")
MESSAGE("   Executable eis_glide  should appear in dir: build/example-drivers/eis/src")
MESSAGE("")
