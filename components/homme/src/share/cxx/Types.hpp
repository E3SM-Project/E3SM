/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_TYPES_HPP
#define HOMMEXX_TYPES_HPP

#include <Kokkos_Core.hpp>

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

#if (HOMMEXX_AVX_VERSION > 0)
using VectorTagType =
    KokkosKernels::Batched::Experimental::AVX<Real, ExecSpace>;
#else
using VectorTagType =
    KokkosKernels::Batched::Experimental::SIMD<Real, ExecSpace>;
#endif // HOMMEXX_AVX_VERSION

using VectorType =
    KokkosKernels::Batched::Experimental::VectorTag<VectorTagType, VECTOR_SIZE>;

using Scalar = KokkosKernels::Batched::Experimental::Vector<VectorType>;

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
template <typename DataType, typename... Properties>
using MPIView = typename ExecView<DataType,Properties...>::HostMirror>::type;
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

// To view the fully expanded name of a complicated template type T,
// just try to access some non-existent field of MyDebug<T>. E.g.:
// MyDebug<T>::type i;
template <typename T> struct MyDebug {};

} // Homme

#endif // HOMMEXX_TYPES_HPP
