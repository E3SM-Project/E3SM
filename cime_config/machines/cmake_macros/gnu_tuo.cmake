string(APPEND CONFIG_ARGS " --host=cray")
# string(APPEND CPPDEFS " -DGPU") # ndk try add --dont think its needed? or not sure where
string(APPEND CPPDEFS " -DTHRUST_IGNORE_CUB_VERSION_CHECK")

set(MPICC "cc")
#set(MPICXX "hipcc") # Needs MPICH_CXX to use hipcc
set(MPICXX "mpicxx") # Needs MPICH_CXX=hipcc
#set(MPICXX "CC") # Needs MPICH_CXX=hipcc  does not like -fno-gpu-rdc flag 
set(MPIFC "ftn") # Linker needs to be the Cray wrapper ftn, not mpif90
set(SCC "cc")
set(SCXX "mpicxx") # kokkos build fail at one point, maybe ok now
#set(SCXX "CC") # kokkos build fail
#set(SCXX "hipcc")
set(SFC "ftn")

#string(APPEND CPPDEFS " -DLINUX -DSCREAM_SYSTEM_WORKAROUND=1") # frontier
if (COMP_NAME STREQUAL gptl)
  #string(APPEND CPPDEFS " -DHAVE_NANOTIME -DBIT64 -DHAVE_SLASHPROC -DHAVE_GETTIMEOFDAY") # worked
  string(APPEND CPPDEFS " -DHAVE_NANOTIME -DBIT64 -DHAVE_VPRINTF -DHAVE_BACKTRACE -DHAVE_SLASHPROC -DHAVE_COMM_F2C -DHAVE_TIMES -DHAVE_GETTIMEOFDAY")
endif()

string(APPEND CMAKE_Fortran_FLAGS " -Wno-implicit-interface")


#string(APPEND SLIBS " -L$ENV{CRAY_HDF5_PARALLEL_PREFIX}/lib -lhdf5_hl -lhdf5 -L$ENV{CRAY_NETCDF_HDF5PARALLEL_PREFIX} -L$ENV{CRAY_PARALLEL_NETCDF_PREFIX}/lib -lpnetcdf -lnetcdf -lnetcdff")
string(APPEND SLIBS " -L$ENV{CRAY_NETCDF_HDF5PARALLEL_PREFIX} -L$ENV{CRAY_PARALLEL_NETCDF_PREFIX}/lib -lpnetcdf -lnetcdf -lnetcdff")
#string(APPEND SLIBS " -L/opt/rocm-5.4.3/lib -lblas -llapack")
string(APPEND SLIBS " -lblas -llapack")
# https://hpc.llnl.gov/documentation/user-guides/using-el-capitan-systems/known-issues#XPMEM_and_GTL
string(APPEND SLIBS " -lxpmem  -L$ENV{CRAY_MPICH_ROOTDIR}/gtl/lib -lmpi_gtl_hsa -Wl,-rpath,$ENV{CRAY_MPICH_ROOTDIR}/gtl/lib") 

# https://hpc.llnl.gov/documentation/user-guides/using-el-capitan-systems/mpi-overview
#LIBS="$PE_MPICH_GTL_DIR_amd_gfx942 $PE_MPICH_GTL_LIBS_amd_gfx942"
#LDFLAGS="-Wl,-rpath,${PE_MPICH_GTL_DIR_amd_gfx942:2}"

string(APPEND SPIO_CMAKE_OPTS " -DPIO_ENABLE_TOOLS:BOOL=OFF")

#string(APPEND KOKKOS_OPTIONS " -DKokkos_ENABLE_HIP=ON -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_OPENMP=OFF")
# gets us around a build error with newer kokkos -- not sure other consequences
string(APPEND KOKKOS_OPTIONS " -DKokkos_ENABLE_HIP=ON -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_OPENMP=OFF -DKokkos_ARCH_AMD_GFX942=ON -DKokkos_ARCH_AMD_GFX942_APU=OFF")
#string(APPEND KOKKOS_OPTIONS " -DKokkos_ENABLE_HIP=ON -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_OPENMP=OFF -DCMAKE_CXX_FLAGS='-std=gnu++14'")
set(USE_HIP "TRUE")

#ndk is this being used?  i don't see "unsafe" flag in build logs
# it may only be used for YAKL code
#string(APPEND CMAKE_HIP_FLAGS "${CXXFLAGS} -munsafe-fp-atomics -x hip")
string(APPEND CMAKE_HIP_FLAGS "${CXXFLAGS} -O2 -munsafe-fp-atomics -x hip")
#string(APPEND CMAKE_HIP_FLAGS " -munsafe-fp-atomics -x hip")

string(APPEND CMAKE_C_FLAGS_RELEASE " -O2")
#string(APPEND CMAKE_CXX_FLAGS_RELEASE " -O2 -g -famd-opt")
#string(APPEND CMAKE_CXX_FLAGS_RELEASE " -O2 -g")
#string(APPEND CMAKE_Fortran_FLAGS_RELEASE   " -O2 -g")
string(APPEND CMAKE_CXX_FLAGS_RELEASE " -O2 -munsafe-fp-atomics -x hip")
string(APPEND CMAKE_Fortran_FLAGS_RELEASE   " -O2")
#string(APPEND CMAKE_C_FLAGS_RELEASE " -O3")
#string(APPEND CMAKE_CXX_FLAGS_RELEASE " -O2")
#string(APPEND CMAKE_Fortran_FLAGS_RELEASE   " -O3")
#string(APPEND CMAKE_C_FLAGS_RELEASE " -g")
#string(APPEND CMAKE_CXX_FLAGS_RELEASE " -O0 -g")
#string(APPEND CMAKE_Fortran_FLAGS_RELEASE   " -O0 -g")

string(APPEND CXX_LIBS " -lstdc++")

set(E3SM_LINK_WITH_FORTRAN "TRUE")
#string(APPEND CMAKE_CXX_FLAGS " -I$ENV{MPICH_DIR}/include --offload-arch=gfx90a") # frontier
#string(APPEND CMAKE_EXE_LINKER_FLAGS    " -lgfortran -L$ENV{CRAY_ROCM_DIR}/lib -lhsa-runtime64 -L$ENV{MPICH_DIR}/lib -lmpi -L$ENV{CRAY_MPICH_ROOTDIR}/gtl/lib -lmpi_gtl_hsa ") # frontier
#string(APPEND CMAKE_EXE_LINKER_FLAGS " -L$ENV{ROCM_PATH}/lib -lamdhip64 -L/opt/gcc/12.2.0/snos/lib64") # frontier
#string(APPEND KOKKOS_OPTIONS " -DKokkos_ENABLE_HIP=On -DKokkos_ARCH_VEGA90A=On -DCMAKE_CXX_FLAGS='-std=gnu++14' -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_OPENMP=OFF") # frontier

set(PIO_FILESYSTEM_HINTS "lustre")
set(SCREAM_ENABLE_MAM OFF) # ndk
