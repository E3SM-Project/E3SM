#ifndef SCREAM_TYPES_HPP
#define SCREAM_TYPES_HPP

#include <Kokkos_Core.hpp>

#include "share/scream_config.hpp"

namespace scream {

#ifdef SCREAM_DOUBLE_PRECISION
using Real = double;
#else
using Real = float;
#endif

typedef int Int;

// Detect ExecSpace using some ifdef logic. I put Serial here as a placeholder
using ExecSpace = Kokkos::Serial;
using HostSpace = Kokkos::HostSpace;

// Memory spaces
using ExecMemSpace = ExecSpace::memory_space;
using HostMemSpace = HostSpace::memory_space;

// Memory traits
using MemoryManaged   = Kokkos::MemoryTraits<Kokkos::Restrict>;
using MemoryUnmanaged = Kokkos::MemoryTraits<Kokkos::Restrict | Kokkos::Unmanaged>;

// A general view in Scream has layout right ('row-major'-like)
template<typename DataType, typename... Properties>
using ViewType = Kokkos::View<DataType, Kokkos::LayoutRight, Properties...>;

// ======== Some particular view's aliases ========== //

// Managed/Unmanaged view
template <typename DataType, typename... Properties>
using ViewManaged = ViewType<DataType, Properties..., MemoryManaged>;
template <typename DataType, typename... Properties>
using ViewUnmanaged = ViewType<DataType, Properties..., MemoryUnmanaged>;

// Host/Deviceviews
template <typename DataType, typename... Properties>
using HostView = ViewType<DataType, HostMemSpace, Properties...>;
template <typename DataType, typename... Properties>
using ExecView = ViewType<DataType, ExecMemSpace, Properties...>;

// Further specializations for execution space and managed/unmanaged memory
template <typename DataType>
using ExecViewManaged = ExecView<DataType, MemoryManaged>;
template <typename DataType>
using ExecViewUnmanaged = ExecView<DataType, MemoryUnmanaged>;

// Further specializations for host space.
template <typename DataType>
using HostViewManaged = HostView<DataType, MemoryManaged>;
template <typename DataType>
using HostViewUnmanaged = HostView<DataType, MemoryUnmanaged>;

namespace util {
// Helper structure template on a type T. It establishes
//  1) if T is a pack
//  2) what's the core type
//  3) what's the size of the pack
// The default impl says T is not a pack, the core type is T itself, and the size is 1.

template<typename T>
struct is_pack {
  using scalar_type = T;
  static constexpr bool value = false;
  static constexpr int  size  = 1;
};
} // namespace util

} // namespace scream

#endif // SCREAM_TYPES_HPP
