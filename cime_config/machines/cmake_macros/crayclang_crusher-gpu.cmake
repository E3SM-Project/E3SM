if (compile_threaded)
  if (NOT COMP_NAME STREQUAL pio2)
    string(APPEND CFLAGS " -fopenmp")
    string(APPEND FFLAGS " -fopenmp")
    string(APPEND CXXFLAGS " -fopenmp")
  endif ()
  string(APPEND LDFLAGS " -fopenmp")
endif()

string(APPEND SLIBS " -L$ENV{PNETCDF_PATH}/lib -lpnetcdf")
set(NETCDF_PATH "$ENV{NETCDF_DIR}")
set(PNETCDF_PATH "$ENV{PNETCDF_DIR}")
set(PIO_FILESYSTEM_HINTS "gpfs")
string(APPEND CXX_LIBS " -lstdc++")

SET(CMAKE_C_COMPILER "cc" CACHE STRING "")
SET(CMAKE_Fortran_COMPILER "ftn" CACHE STRING "")
SET(CMAKE_CXX_COMPILER "hipcc" CACHE STRING "")

#recommended I and L options acc. to crusher man pages
#need to be sorted, append for fflags and cflags prob not needed
string(APPEND CXXFLAGS " -I${MPICH_DIR}/include -L${MPICH_DIR}/lib -lmpi -L/opt/cray/pe/mpich/8.1.16/gtl/lib -lmpi_gtl_hsa")
string(APPEND CFLAGS " -I${MPICH_DIR}/include -L${MPICH_DIR}/lib -lmpi -L/opt/cray/pe/mpich/8.1.16/gtl/lib -lmpi_gtl_hsa")
string(APPEND FFLAGS " -hnoacc -I${MPICH_DIR}/include -L${MPICH_DIR}/lib -lmpi -L/opt/cray/pe/mpich/8.1.16/gtl/lib -lmpi_gtl_hsa")

#this resolves a crash in mct in docn init
if (NOT DEBUG)
  string(APPEND CFLAGS " -O2 -hnoacc -hfp0 -hipa0")
  string(APPEND FFLAGS " -O2 -hnoacc -hfp0 -hipa0")
endif()

string(APPEND CPPDEFS " -DCPRCRAY")

