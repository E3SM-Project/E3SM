# Common settings for our ghci images
include(${CMAKE_CURRENT_LIST_DIR}/ghci.cmake)

# Set Fortran flags
set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch" CACHE STRING "Fortran compiler flags" FORCE)

# Set the path to BLAS/LAPACK libraries
set(BLAS_LIBRARIES "$ENV{BLAS_ROOT}/lib/libopenblas.so" CACHE STRING "Path to BLAS library" FORCE)
set(LAPACK_LIBRARIES "$ENV{BLAS_ROOT}/lib/libopenblas.so" CACHE STRING "Path to LAPACK library" FORCE)

# Set SCREAM_MACHINE
set(SCREAM_MACHINE ghci-hip CACHE STRING "")

# EAMxx enables the HIP language, whose compiler (clang from ROCm) must find the GCC install used by the host
# compiler (libstdc++ headers). The runners' environment (HIPFLAGS) may provide '--gcc-install-dir' with no value,
# which clang rejects. So set CMAKE_HIP_FLAGS explicitly, with the install dir of the gcc behind the MPI wrappers
# (the dir of its libgcc, as reported by the compiler itself).
execute_process (COMMAND gcc -print-libgcc-file-name OUTPUT_VARIABLE GCC_LIBGCC_FILE
                 OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET RESULT_VARIABLE GCC_RC)
if (GCC_RC EQUAL 0 AND GCC_LIBGCC_FILE)
  get_filename_component (GCC_INSTALL_DIR "${GCC_LIBGCC_FILE}" DIRECTORY)
  set (CMAKE_HIP_FLAGS "--gcc-install-dir=${GCC_INSTALL_DIR}" CACHE STRING "" FORCE)
else()
  set (CMAKE_HIP_FLAGS "" CACHE STRING "" FORCE)
endif()

# Enable HIP in kokkos (the arch is passed in by the CI workflow, via Kokkos_ARCH_AMD_GFX* cmake config)
set (EKAT_MACH_FILES_PATH ${CMAKE_CURRENT_LIST_DIR}/../../../../externals/ekat/cmake/machine-files)
include (${EKAT_MACH_FILES_PATH}/kokkos/hip.cmake)

# Our runners do not (yet) feature GPU-aware MPI
option(SCREAM_MPI_ON_DEVICE "Whether to use device pointers for MPI calls" OFF)

# Python tests rely on the venv in ghci.cmake, which is available on the hip image too
option (EAMXX_ENABLE_PYTHON "Whether to enable python interface from eamxx" ON)
