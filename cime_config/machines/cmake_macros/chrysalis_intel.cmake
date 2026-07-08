string(APPEND CMAKE_EXE_LINKER_FLAGS " -L/gpfs/fs1/soft/chrysalis/spack-latest/opt/spack/linux-rhel8-x86_64/gcc-8.5.0/gcc-11.3.0-jkpmtgq/lib64 -lstdc++")

# Workaround for oneapi ifx v2025.2.0 (and earlier than 2025.3.0):
# use -mllvm -disable-hir-temp-cleanup to avoid ICE
# See: https://github.com/argonne-lcf/AuroraBugTracking/issues/64
if (CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM"
    AND CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "2025.3.0")
  string(APPEND CMAKE_Fortran_FLAGS " -mllvm -disable-hir-temp-cleanup")
endif()

string(APPEND CPPDEFS " -DHAVE_SLASHPROC -DHIDE_MPI")

# Hack similar to pm-cpu_intel.cmake. Reset the following flags and load
# with oneapi supported options.
set(CMAKE_CXX_FLAGS " ")
set(CMAKE_Fortran_FLAGS " ")
set(CMAKE_Fortran_FLAGS_DEBUG " ")

if (compile_threaded)
  string(APPEND CMAKE_CXX_FLAGS " -qopenmp")
  string(APPEND CMAKE_Fortran_FLAGS   " -qopenmp")
endif()
string(APPEND CMAKE_CXX_FLAGS " -fp-model precise")
string(APPEND CMAKE_Fortran_FLAGS " -convert big_endian -assume byterecl -traceback -assume realloc_lhs -fpscomp logicals -fp-model precise")
string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -O0 -g ")

set(MPIFC "mpifort")
set(SCC "icx")
set(SCXX "icpx")
set(SFC "ifx")
