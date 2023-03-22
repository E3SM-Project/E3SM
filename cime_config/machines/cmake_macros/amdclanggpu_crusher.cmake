if (NOT DEBUG)
  string(APPEND CFLAGS   " -O2")
  string(APPEND CXXFLAGS " -O2")
  string(APPEND FFLAGS   " -O2")
endif()

set(PIO_FILESYSTEM_HINTS "gpfs")
set(NETCDF_PATH "$ENV{NETCDF_DIR}")
set(PNETCDF_PATH "$ENV{PNETCDF_DIR}")
string(APPEND CMAKE_OPTS " -DPIO_ENABLE_TOOLS:BOOL=OFF")
string(APPEND SLIBS " -L$ENV{CRAY_LIBSCI_PREFIX_DIR}/lib -lsci_amd")

set(MPICXX "hipcc")
set(SCXX "hipcc")
string(APPEND CXXFLAGS " -I$ENV{MPICH_DIR}/include --amdgpu-target=gfx90a")
string(APPEND SLIBS    " -L$ENV{MPICH_DIR}/lib -lmpi -L$ENV{CRAY_MPICH_ROOTDIR}/gtl/lib -lmpi_gtl_hsa")

string(APPEND KOKKOS_OPTIONS " -DKokkos_ENABLE_HIP=On -DKokkos_ARCH_ZEN3=On -DKokkos_ARCH_VEGA90A=On")

