# Get libraries for link line from Trilinos build information
set(TRILINOS_LIBS ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}
   ${Trilinos_EXTRA_LD_FLAGS} ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES} )

message("")
message("  CMake detected the following libraries for linking Fortran with C++ compiler:")
message("     ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES} ")
message("  Compare this to the implicit libraries for C++:")
message("     ${CMAKE_CXX_IMPLICIT_LINK_LIBRARIES} ")
message("  CMake detected the following link directories for linking Fortran with C++ compiler:")
message("     ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES} ")
message("  Compare this to the implicit link directories for C++:")
message("     ${CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES} ")
message("  Bill S. and Jim E.: try this list in CESM builds when linking with C++")

# Need include directories from Trilinos but also mod files from glimmer
include_directories (${GLIMMER_BINARY_DIR}
                     ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})

link_directories (${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS})

# These are local source files needed to make  simple_glide  executable
add_executable(simple_glide  simple_glide.F90 simple_forcing.F90 testsfg.F90)

# Executable depends on several glimmer libraries and Trilinos,
# and potentially an f90main.o file (on jaguar with PGI at least)

IF (NOT ${NO_TRILINOS})
   SET(CISM_USE_LIBS ${GLIMMER_FMAIN}
                     glimmercismfortran
                     glimmercismcpp 
                     ${GLIMMER_NETCDF_LIBS}
                     ${GLIMMER_GPTL_LIB}
		     ${TRILINOS_LIBS}
		     ${GLIMMER_EXTRA_LIBS}) 
ELSE()
   SET(CISM_USE_LIBS glimmercismfortran
                     ${GLIMMER_NETCDF_LIBS}
                     ${GLIMMER_GPTL_LIB}
		     ${GLIMMER_EXTRA_LIBS})
ENDIF()
    
target_link_libraries(simple_glide ${CISM_USE_LIBS})
 
# Helpful(?) message near end of configuration step
MESSAGE("\n   Execute 'make -j 8'")
MESSAGE("   Executable  simple_glide  should appear in dir: build/example-drivers/simple_glide/src")
MESSAGE("")
