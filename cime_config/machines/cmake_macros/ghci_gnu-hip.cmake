# Pick up the generic gnu settings (CIME would not load them automatically for this compiler name)
include (${CMAKE_CURRENT_LIST_DIR}/gnu.cmake)

string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -fcheck=bounds")

# Detect the AMD GPU architecture of the runner (via rocminfo), and enable HIP in kokkos
# NOTE: the MPI wrapper (mpicxx) must use hipcc as backend compiler (set OMPI_CXX/MPICH_CXX before configuring/building)
find_program (ROCMINFO_EXE rocminfo REQUIRED)
execute_process (COMMAND ${ROCMINFO_EXE} OUTPUT_VARIABLE ROCMINFO_OUT OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
if (ROCMINFO_OUT MATCHES "gfx942")
  set (HIP_ARCH_FLAG " -DKokkos_ARCH_AMD_GFX942=ON")
elseif (ROCMINFO_OUT MATCHES "gfx90a")
  set (HIP_ARCH_FLAG " -DKokkos_ARCH_AMD_GFX90A=ON")
elseif (ROCMINFO_OUT MATCHES "gfx908")
  set (HIP_ARCH_FLAG " -DKokkos_ARCH_AMD_GFX908=ON")
else()
  message (FATAL_ERROR "Unsupported AMD GPU architecture (or rocminfo failed). rocminfo output:\n${ROCMINFO_OUT}")
endif()
string (APPEND KOKKOS_OPTIONS " -DKokkos_ENABLE_HIP=ON${HIP_ARCH_FLAG}")

# Change to ON if/when our runners feature GPU-aware MPI
set(SCREAM_MPI_ON_DEVICE OFF CACHE STRING "")
