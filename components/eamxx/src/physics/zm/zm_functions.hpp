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
  template <typename S> using view_2dl          = typename KT::template lview<S**>;
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
    Int             ncol;           // number of columns for current task/chunk
    bool            is_first_step;  // flag for first call

    // view_2d<Spack>  zmid;           // mid-point level altitude [m]
    // view_2d<Spack>  zint;           // interface level altitude [m]
    // view_2d<Spack>  pmid;           // mid-point level pressure [Pa]
    // view_2d<Spack>  pint;           // interface level pressure [Pa]
    // view_2d<Spack>  pdel;           // pressure thickness [Pa]
    view_1d<Scalar>  phis;           // surface geopotential height [m2/s]
    // view_1d<Scalar>  pblh;           // PBL height [m]

    view_2d<Spack>  T_mid;          // Temperature [K]
    view_2d<Spack>  qv;             // Water vapor mixing ratio [kg kg-1]

    // view_2d<Spack>  qc;             // Cloud mass mixing ratio [kg kg-1]
    // view_2d<Spack>  qi;             // Ice total mass mixing ratio [kg kg-1]
    // view_2d<Spack>  omega;          // vertical pressure velocity [Pa/s]

    // -------------------------------------------------------------------------
    // vectors for alternate transpose method
    std::vector<Real> phis_v;
    std::vector<Real> T_mid_v;
    std::vector<Real> qv_v;
    // -------------------------------------------------------------------------
    // // transpose method for fortran bridging
    // template <ekat::TransposeDirection::Enum D>
    // void transpose()
    // {
    //   std::vector<view_2d<Spack>> transposed_views(2);
    //   ekat::host_to_device<D>({ekat::scalarize(T_mid).data(),
    //                            ekat::scalarize(qv).data()},
    //                            T_mid.extent(0), T_mid.extent(1) * Spack::n, transposed_views, true);
    //   if (D == ekat::TransposeDirection::c2f) {
    //     T_mid = transposed_views[0];
    //     qv    = transposed_views[1];
    //   // else {
    //     // ???
    //   }
    // };
    // -------------------------------------------------------------------------
    // alternate transpose method
    template <ekat::TransposeDirection::Enum D>
    void transpose(int pver)
    {
      std::vector<view_2d<Spack>> tviews = { T_mid, qv };
      if (D == ekat::TransposeDirection::c2f) {
        auto rsz = T_mid.size() * Spack::n;
        // auto rsz = T_mid.extent(1) * Spack::n;
        // ---------------------------------------------------------------------
        T_mid_v .resize(rsz);
        qv_v    .resize(rsz);
        ekat::device_to_host( { T_mid_v .data(),
                                qv_v    .data() },
                              T_mid.extent(0), T_mid.extent(1)*Spack::n, tviews, true);
        T_mid = tviews[0];
        qv    = tviews[1];
        // ---------------------------------------------------------------------
      }
      // else {
      //   ekat::host_to_device({T_mid_v .data(),
      //                         qv_v    .data()},
      //                        T_mid.extent(0), T_mid.extent(1) * Spack::n, transposed_views, true);
      // }
    }
    // -------------------------------------------------------------------------
  };

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
