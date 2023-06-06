if (compile_threaded)
  string(APPEND FFLAGS   " -fopenmp")
  string(APPEND CFLAGS   " -fopenmp")
  string(APPEND CXXFLAGS " -fopenmp")
  string(APPEND LDFLAGS  " -fopenmp")
endif()
if (DEBUG)
  string(APPEND CFLAGS   " -O0 -g")
  string(APPEND FFLAGS   " -O0 -g")
  string(APPEND CXXFLAGS " -O0 -g")
  string(APPEND CPPDEFS " -DYAKL_DEBUG")
endif()
string(APPEND CPPDEFS " -DFORTRANUNDERSCORE -DNO_R16 -DCPRCRAY")
string(APPEND FC_AUTO_R8 " -s real64")
string(APPEND FFLAGS " -f free  -em")
if (NOT compile_threaded)
  # -M1077 flag used to suppress message about OpenMP directives
  # that are ignored for non-threaded builds. (-h omp inactive)
  # Details: `explain ftn-1077`
  string(APPEND FFLAGS " -M1077")
endif()
string(APPEND FFLAGS_NOOPT " -O0")
set(HAS_F2008_CONTIGUOUS "TRUE")

# -Wl,--allow-shlib-undefined was added to address rocm 5.4.3 Fortran linker issue:
# /opt/rocm-5.4.3/lib/libhsa-runtime64.so.1: undefined reference to `std::condition_variable::wait(std::unique_lock<std::mutex>&)@GLIBCXX_3.4.30'
# AMD started building with GCC 12.2.0, which brings in a GLIBCXX symbol that isn't in CCE's default GCC toolchain.
string(APPEND LDFLAGS " -Wl,--allow-multiple-definition -Wl,--allow-shlib-undefined")

set(SUPPORTS_CXX "TRUE")
set(CXX_LINKER "FORTRAN")
set(MPICC "cc")
set(MPICXX "hipcc")
set(MPIFC "ftn")
set(SCC "cc")
set(SCXX "hipcc")
set(SFC "ftn")

# Switch to O3 for better performance
# Using O2 to ensure passing tests
if (NOT DEBUG)
  string(APPEND CFLAGS   " -O2")
  string(APPEND CXXFLAGS " -O2")
  string(APPEND FFLAGS   " -O2")
endif()

if (COMP_NAME STREQUAL elm)
  # See Land NaNs in conditionals: https://github.com/E3SM-Project/E3SM/issues/4996
  string(APPEND FFLAGS " -hfp0")
endif()
# -em -ef generates modulename.mod (lowercase files) to support
# Scorpio installs
# Disable ipa and zero initialization are for other NaN isues:
# https://github.com/E3SM-Project/E3SM/pull/5208
string(APPEND FFLAGS " -hipa0 -hzero -em -ef -hnoacc")

set(NETCDF_PATH "$ENV{NETCDF_DIR}")
set(PNETCDF_PATH "$ENV{PNETCDF_DIR}")
string(APPEND CMAKE_OPTS " -DPIO_ENABLE_TOOLS:BOOL=OFF")
string(APPEND CXX_LIBS " -lstdc++")

string(APPEND CXXFLAGS " -I$ENV{MPICH_DIR}/include --offload-arch=gfx90a")
string(APPEND SLIBS    " -L$ENV{MPICH_DIR}/lib -lmpi -L$ENV{CRAY_MPICH_ROOTDIR}/gtl/lib -lmpi_gtl_hsa")
if (NOT MPILIB STREQUAL mpi-serial)
  string(APPEND SLIBS " -L$ENV{ADIOS2_DIR}/lib64 -ladios2_c_mpi -ladios2_c -ladios2_core_mpi -ladios2_core -ladios2_evpath -ladios2_ffs -ladios2_dill -ladios2_atl -ladios2_enet")
endif()

string(APPEND KOKKOS_OPTIONS " -DKokkos_ENABLE_HIP=On -DKokkos_ARCH_ZEN3=On -DKokkos_ARCH_VEGA90A=On")

set(USE_HIP "TRUE")
string(APPEND HIP_FLAGS "${CXXFLAGS} -munsafe-fp-atomics -x hip")
