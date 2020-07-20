#ifndef EKAT_TYPES_HPP
#define EKAT_TYPES_HPP

#include "ekat_config.h"
#include "ekat/scream_kokkos_meta.hpp"

/*
 * Header contains globally useful types for Scream.
 * The global Int and Real types are defined here along
 * with a type dictionary for accessing commonly-used
 * Kokkos types.
 */

namespace scream {

#ifdef EKAT_DOUBLE_PRECISION
using Real = double;
#else
using Real = float;
#endif

typedef int Int;

#if defined KOKKOS_COMPILER_GNU
// See https://github.com/kokkos/kokkos-kernels/issues/129
# define ConstExceptGnu
#else
# define ConstExceptGnu const
#endif

// The default device we use.
using DefaultDevice = Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;

// A device type to force host execution
using HostDevice = Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::DefaultHostExecutionSpace::memory_space>;

// Struct for getting useful Kokkos types based on the device
template <typename DeviceType>
struct KokkosTypes
{
  using Device = DeviceType;
  using Layout = Kokkos::LayoutRight;
  using MemSpace = typename Device::memory_space;
  using ExeSpace = typename Device::execution_space;
  using TeamPolicy = Kokkos::TeamPolicy<ExeSpace>;
  using MemberType = typename TeamPolicy::member_type;
  using RangePolicy = Kokkos::RangePolicy<ExeSpace>;
  template<typename TagType>
  using RangeTagPolicy = Kokkos::RangePolicy<ExeSpace,TagType>;
  template<typename TagType>
  using TeamTagPolicy = Kokkos::TeamPolicy<ExeSpace,TagType>;

  template <typename DataType>
  using view = Kokkos::View<DataType, Layout, Device>;

  // left-layout views, may be useful for interacting with fortran
  template <typename DataType>
  using lview = Kokkos::View<DataType, Kokkos::LayoutLeft, Device>;

  template <typename Scalar>
  using view_1d = view<Scalar*>;

  template <typename Scalar>
  using view_2d = view<Scalar**>;

  template <typename Scalar>
  using view_3d = view<Scalar***>;

  template <typename Scalar, int X>
  using view_1d_table = view<const Scalar[X]>;

  template <typename Scalar, int X, int Y>
  using view_2d_table = view<const Scalar[X][Y]>;

  // Our workspace implementation makes this a useful type
  template <typename Scalar, int N>
  using view_1d_ptr_array = Kokkos::Array<ko::Unmanaged<view_1d<Scalar> >*, N>;

  template <typename Scalar, int N>
  using view_1d_ptr_carray = Kokkos::Array<const ko::Unmanaged<view_1d<Scalar> >*, N>;
};

// Memory traits
using MemoryManaged   = Kokkos::MemoryTraits<Kokkos::Restrict>;
using MemoryUnmanaged = Kokkos::MemoryTraits<Kokkos::Restrict | Kokkos::Unmanaged>;

namespace util {
// Helper structure templated on a type T. It establishes
//  1) if T is a pack
//  2) what's the underlying scalar type
// The default impl says T is not a pack, and that the scalar type is T itself.

template<typename T>
struct ScalarProperties {
  using scalar_type = T;
  static constexpr bool is_pack = false;
};

} // namespace util

// An enum to be used with object that have 'repository'-like behavior
enum class RepoState {
  Clean,
  Open,
  Closed
};

} // namespace scream

#endif // EKAT_TYPES_HPP
