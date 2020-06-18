#ifndef INCLUDE_SCREAM_KOKKOS
#define INCLUDE_SCREAM_KOKKOS

// Funnel all Kokkos includes through this file.

// We do not want to see warnings coming from Kokkos. This is currently
// handled through CMAKE, using the SYSTEM qualifier when setting up
// target_include_directories. Future Kokkos changes may complicate
// this approach.
//#pragma GCC system_header

#include <Kokkos_Core.hpp>

#endif
