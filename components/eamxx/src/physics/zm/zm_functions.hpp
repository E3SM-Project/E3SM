#ifndef ZM_FUNCTIONS_HPP
#define ZM_FUNCTIONS_HPP

#include "physics/share/physics_constants.hpp"

#include "share/eamxx_types.hpp"

#include <ekat_pack_kokkos.hpp>
#include <ekat_workspace.hpp>

namespace scream {
namespace zm {

/*
 * Functions is a stateless struct used to encapsulate a number of functions for ZM.
 */

template <typename ScalarT, typename DeviceT>
struct Functions
{
  // ---------------------------------------------------------------------------
  // Types

  using Scalar = ScalarT;
  using Device = DeviceT;

  template <typename S> using BigPack = ekat::Pack<S,SCREAM_PACK_SIZE>;
  template <typename S> using SmallPack = ekat::Pack<S,SCREAM_SMALL_PACK_SIZE>;
  
  using SPackInt = SmallPack<Int>;
  using BPack    = BigPack<Scalar>;
  using Spack    = SmallPack<Scalar>;

  using KT = ekat::KokkosTypes<Device>;

  template <typename S> using view_1d           = typename KT::template view_1d<S>;
  template <typename S> using view_2d           = typename KT::template view_2d<S>;
  template <typename S> using view_3d           = typename KT::template view_3d<S>;
  template <typename S> using view_2d_strided   = typename KT::template sview<S**>;
  template <typename S> using view_3d_strided   = typename KT::template sview<S***>;
  template <typename S> using uview_1d          = typename ekat::template Unmanaged<view_1d<S> >;
  template <typename S> using uview_2d          = typename ekat::template Unmanaged<view_2d<S> >;
  template <typename S> using uview_2d_strided  = typename ekat::template Unmanaged<view_2d_strided<S> >;

  // ---------------------------------------------------------------------------
  // Structs
  struct ZMRuntime {
    ZMRuntime() = default;
  };

  // This struct stores input views for ZM_main.
  struct ZMInput {
    ZMInput() = default;
  };

  // This struct stores input/outputs views for ZM_main.
  struct ZMInputOutput {
    ZMInputOutput() = default;
  };

  // This struct stores output only views for ZM_main.
  struct ZMOutput {
    ZMOutput() = default;
  };

  // This struct stores output views for ZM diagnostics for ZM_main.
  struct ZMHistoryOutput {
    ZMHistoryOutput() = default;
  };

  // ---------------------------------------------------------------------------
  // Functions


}; // struct Functions

} // namespace zm
} // namespace scream

#endif // ZM_FUNCTIONS_HPP
