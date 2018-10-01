#ifndef INCLUDE_SCREAM_KOKKOS_META
#define INCLUDE_SCREAM_KOKKOS_META

#include <Kokkos_Core.hpp>

namespace scream {
namespace ko {

// Kokkos-related types not specific to the app. Thus, do not include
// app-specific array types here; rather, include only those things that one
// could imagine could go into Kokkos itself.

// Make the input View Unmanaged, whether or not it already is. One might
// imagine that View::unmanaged_type would provide this.
//   Use: Unmanged<ViewType>
template <typename View>
using Unmanaged = 
  // Provide a full View type specification, augmented with Unmanaged.
  Kokkos::View<typename View::traits::scalar_array_type,
               typename View::traits::array_layout,
               typename View::traits::device_type,
               Kokkos::MemoryTraits<
                 // All the current values...
                 (View::traits::memory_traits::RandomAccess ? Kokkos::RandomAccess : 0) |
                 (View::traits::memory_traits::Atomic ? Kokkos::Atomic : 0) |
                 (View::traits::memory_traits::Restrict ? Kokkos::Restrict : 0) |
                 (View::traits::memory_traits::Aligned ? Kokkos::Aligned : 0) |
                 // ... |ed with the one we want, whether or not it's
                 // already there.
                 Kokkos::Unmanaged> >;

} // namespace ko
} // namespace scream

#endif // INCLUDE_SCREAM_KOKKOS_META
