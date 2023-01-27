string(APPEND FFLAGS " -Wno-implicit-interface")

if (NOT DEBUG)
  string(APPEND CFLAGS   " -O2")
  string(APPEND CXXFLAGS " -O2")
  string(APPEND FFLAGS   " -O2")
endif()
if (COMP_NAME STREQUAL gptl)
  string(APPEND CPPDEFS " -DHAVE_SLASHPROC")
endif()
set(NETCDF_PATH "$ENV{NETCDF_DIR}")
set(PNETCDF_PATH "$ENV{PNETCDF_DIR}")
set(PIO_FILESYSTEM_HINTS "gpfs")
string(APPEND CMAKE_OPTS " -DPIO_ENABLE_TOOLS:BOOL=OFF")

set(MPICC  "cc")
set(MPICXX "hipcc")
set(MPIFC  "ftn")
set(SCC  "cc")
set(SCXX "hipcc")
set(SFC  "ftn")

string(APPEND CXXFLAGS " -I$ENV{MPICH_DIR}/include --amdgpu-target=gfx90a")
string(APPEND SLIBS    " -Wl,--copy-dt-needed-entries -L/opt/cray/pe/gcc-libs -lgfortran -L$ENV{MPICH_DIR}/lib -lmpi -L$ENV{CRAY_MPICH_ROOTDIR}/gtl/lib -lmpi_gtl_hsa ")

string(APPEND KOKKOS_OPTIONS " -DKokkos_ENABLE_HIP=On -DKokkos_ARCH_ZEN3=On -DKokkos_ARCH_VEGA90A=On")

