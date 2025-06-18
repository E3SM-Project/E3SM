#ifndef ZM_FUNCTIONS_HPP
#define ZM_FUNCTIONS_HPP

#include "physics/share/physics_constants.hpp"

#include "share/eamxx_types.hpp"

#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/ekat_workspace.hpp"

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
    // Runtime options
    // Scalar ????;
  };

  // This struct stores input views for ZM_main.
  struct ZMInput {
    ZMInput() = default;
    // // heights, for thermo grid [m]
    // view_2d<const Spack>  zmid;
    // // heights, for interface grid [m]
    // view_2d<const Spack>  zint;
    // // pressure levels on thermo grid [Pa]
    // view_2d<const Spack>  pmid;
    // // pressure levels on interface grid [Pa]
    // view_2d<const Spack>  pint;
    // // Differences in pressure levels [Pa]
    // view_2d<const Spack>  pdel;
    // //  temperature [K]
    // view_2d<const Spack>  T_mid;
    // // large scale vertical velocity [m/s]
    // view_2d<const Spack>  w_field;
    // // Host model surface geopotential height
    // view_1d<const Scalar> phis;
  };

  // This struct stores input/outputs views for ZM_main.
  struct ZMInputOutput {
    ZMInputOutput() = default;
    // // prognostic temp variable of host model
    // // dry static energy [J/kg]
    // // dse = Cp*T + g*z + phis
    // view_2d<Spack>  host_dse;
  };

  // This struct stores output only views for ZM_main.
  struct ZMOutput {
    ZMOutput() = default;
    // // planetary boundary layer depth [m]
    // view_1d<Scalar> pblh;
    // // cloud liquid mixing ratio variance [kg^2/kg^2]
    // view_2d<Spack>  ZM_ql2;
  };

  // This struct stores output views for ZM diagnostics for ZM_main.
  struct ZMHistoryOutput {
    ZMHistoryOutput() = default;
    // // ??? [???]
    // view_2d<Spack>  ZM_???;
  };

}; // struct Functions

} // namespace zm
} // namespace scream

#endif // ZM_FUNCTIONS_HPP
