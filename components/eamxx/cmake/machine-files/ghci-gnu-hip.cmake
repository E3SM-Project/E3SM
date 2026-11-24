# Machine file used by CIME builds of EAMxx on ghci with the gnu-hip compiler
# (build_eamxx.cmake looks for <mach>-<compiler>.cmake, and unsets KOKKOS_OPTIONS,
# so the HIP backend and the AMD arch must be set here, as direct Kokkos settings).
include(${CMAKE_CURRENT_LIST_DIR}/ghci-hip.cmake)

# Detect the AMD GPU arch of the runner (the standalone CI passes it in via -c instead)
find_program (ROCMINFO_EXE rocminfo REQUIRED)
execute_process (COMMAND ${ROCMINFO_EXE} OUTPUT_VARIABLE ROCMINFO_OUT OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
if (ROCMINFO_OUT MATCHES "gfx942")
  set (Kokkos_ARCH_AMD_GFX942 ON CACHE BOOL "")
elseif (ROCMINFO_OUT MATCHES "gfx90a")
  set (Kokkos_ARCH_AMD_GFX90A ON CACHE BOOL "")
elseif (ROCMINFO_OUT MATCHES "gfx908")
  set (Kokkos_ARCH_AMD_GFX908 ON CACHE BOOL "")
else()
  message (FATAL_ERROR "Unsupported AMD GPU architecture (or rocminfo failed). rocminfo output:\n${ROCMINFO_OUT}")
endif()
