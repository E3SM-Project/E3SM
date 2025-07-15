#ifndef ZM_FUNCTIONS_HPP
#define ZM_FUNCTIONS_HPP

#include "physics/share/physics_constants.hpp"

#include "share/core/eamxx_types.hpp"

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

  struct zm_runtime_opt {
    zm_runtime_opt() = default;
  };

  struct zm_input_state {
    zm_input_state() = default;
    // -------------------------------------------------------------------------
    // number of columns for current task/chunk
    Int ncol;
    // flag for first call
    bool is_first_step;
    // Temperature [K]
    view_2d<Spack> T_mid;
    // Water vapor mixing ratio [kg kg-1]
    view_2d<Spack> qv;
    // -------------------------------------------------------------------------
    // transpose method for fortran bridging
    template <ekat::TransposeDirection::Enum D>
    void transpose()
    {
      std::vector<view_2d<Spack>> transposed_views(2);
      ekat::host_to_device<D>({ekat::scalarize(T_mid).data(),
                               ekat::scalarize(qv).data()},
                               T_mid.extent(0), T_mid.extent(1) * Spack::n,
                               transposed_views, true);
      if (D == ekat::TransposeDirection::c2f) {
        T_mid = transposed_views[0];
        qv    = transposed_views[1];
      // else {
        // ???
      }
    };
    // -------------------------------------------------------------------------
  };

    // // surface geopotential height [m2/s]
    // view_1d<Spack> phis;
    // // mid-point level altitude [m]
    // view_2d<Spack> zmid;
    // // interface level altitude [m]
    // view_2d<Spack> zint;
    // // mid-point level pressure [Pa]
    // view_2d<Spack> pmid;
    // // interface level pressure [Pa]
    // view_2d<Spack> pint;
    // // pressure thickness [Pa]
    // view_2d<Spack> pdel;
    // // PBL height [m]
    // view_1d<Spack> pblh;
    // // Temperature [K]
    // view_2d<Spack> T_mid;
    // // Water vapor mixing ratio [kg kg-1]
    // view_2d<Spack> qv;
    // // Cloud mass mixing ratio [kg kg-1]
    // view_2d<Spack> qc;
    // // Ice total mass mixing ratio [kg kg-1]
    // view_2d<Spack> qi;
    // // vertical pressure velocity [Pa/s]
    // view_2d<Spack> omega;

  struct zm_output_tend {
    zm_output_tend() = default;
  };

  struct zm_output_diag {
    zm_output_diag() = default;
  };

  // struct zm_runtime_opts {
  //   zm_runtime_opts() = default;
  // };

  // // This struct stores input views for ZM_main.
  // struct zm_input_state {
    // zm_input_state() = default;
  // };

  // // This struct stores input/outputs views for ZM_main.
  // struct zm_inout {
  //   zm_inout() = default;
  // };

  // // This struct stores output only views for ZM_main.
  // struct zm_output {
  //   zm_output() = default;
  // };

  // // This struct stores output views for ZM diagnostics for ZM_main.
  // struct ZMHistoryOutput {
  //   ZMHistoryOutput() = default;
  // };

  // ---------------------------------------------------------------------------
  // Functions

  // static Int zm_main()

}; // struct Functions

} // namespace zm
} // namespace scream

#endif // ZM_FUNCTIONS_HPP
