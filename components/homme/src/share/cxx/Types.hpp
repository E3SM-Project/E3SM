/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_TYPES_HPP
#define HOMMEXX_TYPES_HPP

#include <Kokkos_Core.hpp>

#include "PackTraits.hpp"
#include "Config.hpp"
#include "ExecSpaceDefs.hpp"
#include "Dimensions.hpp"

#include <vector/KokkosKernels_Vector.hpp>

#define __MACRO_STRING(MacroVal) #MacroVal
#define MACRO_STRING(MacroVal) __MACRO_STRING(MacroVal)

namespace Homme {

// Usual typedef for real scalar type
using Real = double;
using RCPtr = Real *const;
using CRCPtr = const Real *const;
using F90Ptr = Real *const; // Using this in a function signature emphasizes
                            // that the ordering is Fortran
using CF90Ptr = const Real *const; // Using this in a function signature
                                   // emphasizes that the ordering is Fortran

using VectorTagType = KokkosKernels::Batched::Experimental::SIMD<Real, ExecSpace>;

using VectorType = KokkosKernels::Batched::Experimental::VectorTag<VectorTagType, VECTOR_SIZE>;

using Scalar = KokkosKernels::Batched::Experimental::Vector<VectorType>;

// Specialize PackTraits for Scalar
template<>
struct PackTraits<Scalar> {
  static constexpr int pack_length = Scalar::vector_length;
  using value_type = Real;
};

static_assert(sizeof(Scalar) > 0, "Vector type has 0 size");
static_assert(sizeof(Scalar) == sizeof(Real[VECTOR_SIZE]), "Vector type is not correctly defined");
static_assert(Scalar::vector_length>0, "Vector type is not correctly defined (vector_length=0)");

using MemoryManaged   = Kokkos::MemoryTraits<Kokkos::Restrict>;
using MemoryUnmanaged = Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::Restrict>;

// The memory spaces
using ExecMemSpace    = ExecSpace::memory_space;
using ScratchMemSpace = ExecSpace::scratch_memory_space;
using HostMemSpace    = Kokkos::HostSpace;
// Selecting the memory space for the MPI buffers, that is, where the MPI
// buffers will be allocated. In a CPU build, this is always going to be
// the same as ExecMemSpace, i.e., HostSpace. In a GPU build, one can choose
// whether the MPI is done on host or on device. If on device, then all MPI
// calls will be issued from device pointers.
// NOTE: this has nothing to do with pack/unpack, which ALWAYS happen on
//       device views (to be done in parallel). The difference is ONLY in
//       the location of the MPI buffer for send/receive.

#if HOMMEXX_MPI_ON_DEVICE
  using MPIMemSpace = ExecMemSpace;
#else
  using MPIMemSpace = HostMemSpace;
#endif

// A team member type
using TeamMember     = Kokkos::TeamPolicy<ExecSpace>::member_type;

// Short name for views
template <typename DataType, typename... Properties>
using ViewType = Kokkos::View<DataType, Kokkos::LayoutRight, Properties...>;

// Managed/Unmanaged view
template <typename DataType, typename... Properties>
using ViewManaged = ViewType<DataType, Properties..., MemoryManaged>;
template <typename DataType, typename... Properties>
using ViewUnmanaged = ViewType<DataType, Properties..., MemoryUnmanaged>;

// Host/Device/MPI views
template <typename DataType, typename... Properties>
using HostView = ViewType<DataType, HostMemSpace, Properties...>;
template <typename DataType, typename... Properties>
using ExecView = ViewType<DataType, ExecMemSpace, Properties...>;
// Work around Cuda 9.1.85 parse bugs.
#if CUDA_PARSE_BUG_FIXED
template <typename DataType, typename... Properties>
using MPIView = typename std::conditional<std::is_same<MPIMemSpace,ExecMemSpace>::value,
                                          ExecView<DataType,Properties...>,
                                          typename ExecView<DataType,Properties...>::HostMirror>::type;
#else
# if HOMMEXX_MPI_ON_DEVICE
template <typename DataType, typename... Properties>
using MPIView = ExecView<DataType,Properties...>;
# else
/// A Cuda 9.1.85 (and probably other 9.1.* versions) parse bug requires us to
/// spell out this type rather than rely on our layers of abstraction.
template <typename DataType, typename... Properties>
using MPIView = typename Kokkos::View<DataType, Kokkos::LayoutRight, ExecMemSpace, Properties...>::HostMirror;
# endif
#endif

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

// Further specializations for MPI memory space.
template <typename DataType>
using MPIViewManaged = MPIView<DataType, MemoryManaged>;
template <typename DataType>
using MPIViewUnmanaged = MPIView<DataType, MemoryUnmanaged>;

// The scratch view type: always unmanaged, and always with c pointers
template <typename DataType>
using ScratchView = ViewType<DataType, ScratchMemSpace, MemoryUnmanaged>;

namespace Impl {
// Turn a View's MemoryTraits (traits::memory_traits) into the equivalent
// unsigned int mask. This is an implementation detail for Unmanaged; see next.
template <typename View>
  struct MemoryTraitsMask {
    enum : unsigned int {
      value = ((View::traits::memory_traits::is_random_access ? Kokkos::RandomAccess : 0) |
               (View::traits::memory_traits::is_atomic ? Kokkos::Atomic : 0) |
               (View::traits::memory_traits::is_restrict ? Kokkos::Restrict : 0) |
               (View::traits::memory_traits::is_aligned ? Kokkos::Aligned : 0) |
               (View::traits::memory_traits::is_unmanaged ? Kokkos::Unmanaged : 0))
        };
  };
}

template <typename View>
using Unmanaged =
  // Provide a full View type specification, augmented with Unmanaged.
  Kokkos::View<typename View::traits::scalar_array_type,
               typename View::traits::array_layout,
               typename View::traits::device_type,
               Kokkos::MemoryTraits<
                 // All the current values...
                 Impl::MemoryTraitsMask<View>::value |
                 // ... |ed with the one we want, whether or not it's
                 // already there.
                 Kokkos::Unmanaged> >;


// To view the fully expanded name of a complicated template type T,
// just try to access some non-existent field of MyDebug<T>. E.g.:
// MyDebug<T>::type i;
template <typename T> struct MyDebug {};

namespace Remap {
// All VertRemapAlg types must provide the following methods:
// compute_grids_phase, and compute_remap_phase
//
// compute_grids_phase is expected to have less parallelism available and to
// compute quantities which are independent of the tracers,
// based on the computed partitions
//
// compute_remap_phase remaps each of the tracers based on the quantities
// previously computed in compute_grids_phase.
// It is also expected to have a large amount of parallelism, specifically
// qsize * num_elems
struct VertRemapAlg {};
} // namespace Remap

} // Homme

// A kokkos-compatible implementation of a static array of 2 Real's
namespace Kokkos {

struct Real2 {
  Homme::Real v[2];
  KOKKOS_FORCEINLINE_FUNCTION
  Real2 () : v{0,0} {}

  KOKKOS_FORCEINLINE_FUNCTION
  Real2& operator+= (const Real2& o) {
    v[0] += o.v[0];
    v[1] += o.v[1];
    return *this;
  }
};

template<>
struct reduction_identity<Homme::Scalar> {
  KOKKOS_FORCEINLINE_FUNCTION
  static Homme::Scalar sum()  {return Homme::Scalar(reduction_identity<Homme::Real>::sum());}
};

// Specialization of a Kokkos structure, needed in the initialization of reduction operations.
template<> struct reduction_identity<Real2> {
  KOKKOS_FORCEINLINE_FUNCTION
  static Real2 sum() { return Real2(); }
};

} // namespace Kokkos

#endif // HOMMEXX_TYPES_HPP
