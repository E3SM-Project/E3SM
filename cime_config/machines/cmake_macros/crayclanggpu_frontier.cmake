set(MPICC "cc")
set(MPICXX "mpicxx")
#set(MPICXX "CC")
set(MPIFC "ftn")
set(SCC "cc")
set(SCXX "hipcc")
set(SFC "ftn")

string(APPEND CPPDEFS " -DLINUX")
if (COMP_NAME STREQUAL gptl)
    string(APPEND CPPDEFS " -DHAVE_NANOTIME -DBIT64 -DHAVE_SLASHPROC -DHAVE_COMM_F2C -DHAVE_TIMES -DHAVE_GETTIMEOFDAY")
endif()

if (compile_threaded)
  string(APPEND CMAKE_Fortran_FLAGS   " -fopenmp")
  string(APPEND CMAKE_C_FLAGS   " -fopenmp")
  string(APPEND CMAKE_CXX_FLAGS " -fopenmp")
  string(APPEND CMAKE_EXE_LINKER_FLAGS  " -fopenmp")
endif()
string(APPEND CMAKE_C_FLAGS_DEBUG   " -O0 -g")
string(APPEND CMAKE_Fortran_FLAGS_DEBUG   " -O0 -g")
string(APPEND CMAKE_CXX_FLAGS_DEBUG " -O0 -g")
string(APPEND CPPDEFS_DEBUG " -DYAKL_DEBUG")
string(APPEND CPPDEFS " -DFORTRANUNDERSCORE -DNO_R16 -DCPRCRAY")
string(APPEND CMAKE_Fortran_FLAGS " -f free  -em")
if (NOT compile_threaded)
  # -M1077 flag used to suppress message about OpenMP directives
  # that are ignored for non-threaded builds. (-h omp inactive)
  # Details: `explain ftn-1077`
  string(APPEND CMAKE_Fortran_FLAGS " -M1077")
endif()
set(HAS_F2008_CONTIGUOUS "TRUE")

# -Wl,--allow-shlib-undefined was added to address rocm 5.4.3 Fortran linker issue:
# /opt/rocm-5.4.3/lib/libhsa-runtime64.so.1: undefined reference to `std::condition_variable::wait(std::unique_lock<std::mutex>&)@GLIBCXX_3.4.30'
# AMD started building with GCC 12.2.0, which brings in a GLIBCXX symbol that isn't in CCE's default GCC toolchain.
#string(APPEND CMAKE_EXE_LINKER_FLAGS " -Wl,--allow-multiple-definition -Wl,--allow-shlib-undefined")

# Switching to O3 for performance benchmarking
# Will revisit any failing tests
string(APPEND CMAKE_C_FLAGS_RELEASE   " -O3")
string(APPEND CMAKE_CXX_FLAGS_RELEASE " -O3")
string(APPEND CMAKE_Fortran_FLAGS_RELEASE   " -O3")

if (COMP_NAME STREQUAL elm)
  # See Land NaNs in conditionals: https://github.com/E3SM-Project/E3SM/issues/4996
  string(APPEND CMAKE_Fortran_FLAGS " -hfp0")
endif()
# -em -ef generates modulename.mod (lowercase files) to support
# Scorpio installs
# Disable ipa and zero initialization are for other NaN isues:
# https://github.com/E3SM-Project/E3SM/pull/5208
string(APPEND CMAKE_Fortran_FLAGS " -hipa0 -hzero -em -ef -hnoacc")

string(APPEND CMAKE_CXX_FLAGS " --offload-arch=gfx90a")
string(APPEND CMAKE_EXE_LINKER_FLAGS " -L$ENV{CRAY_MPICH_ROOTDIR}/gtl/lib -lmpi_gtl_hsa")
string(APPEND CMAKE_EXE_LINKER_FLAGS " -L$ENV{ROCM_PATH}/lib -lamdhip64")
string(APPEND CMAKE_EXE_LINKER_FLAGS " -L/opt/cray/pe/gcc/12.2.0/snos/lib64 -lgfortran -lstdc++")

string(APPEND KOKKOS_OPTIONS " -DKokkos_ENABLE_HIP=On -DKokkos_ARCH_ZEN3=On -DKokkos_ARCH_VEGA90A=On")

set(USE_HIP "TRUE")
string(APPEND CMAKE_HIP_FLAGS "${CXXFLAGS} -munsafe-fp-atomics -x hip")
set(E3SM_LINK_WITH_FORTRAN "TRUE")
