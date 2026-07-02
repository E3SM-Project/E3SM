string(APPEND CMAKE_EXE_LINKER_FLAGS " -L/gpfs/fs1/soft/chrysalis/spack-latest/opt/spack/linux-rhel8-x86_64/gcc-8.5.0/gcc-11.3.0-jkpmtgq/lib64 -lstdc++")

if (CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
  if (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "2025.3.0")
    # Workaround for oneapi ifx v2025.2.0 (and earlier than 2025.3.0):
    # use -mllvm -disable-hir-temp-cleanup to avoid ICE
    # See: https://github.com/argonne-lcf/AuroraBugTracking/issues/64
    string(APPEND CMAKE_Fortran_FLAGS " -mllvm -disable-hir-temp-cleanup")
  endif()
  if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "2025.0")
      # Sanitization workaround for Intel 2025.0 and later
      string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -check nouninit")
  endif()
endif()
