
string(APPEND CMAKE_EXE_LINKER_FLAGS " -L/gpfs/fs1/soft/chrysalis/spack-latest/opt/spack/linux-rhel8-x86_64/gcc-8.5.0/gcc-11.3.0-jkpmtgq/lib64 -lstdc++")

# -mllvm -disable-hir-temp-cleanup for oneapi v2025.2.0 https://github.com/argonne-lcf/AuroraBugTracking/issues/64
string(APPEND CMAKE_Fortran_FLAGS " -mllvm -disable-hir-temp-cleanup")

