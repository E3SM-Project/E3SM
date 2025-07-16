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
    Int             pcol;           // max number of columns across tasks/chunks
    bool            is_first_step;  // flag for first call

    view_1d<const Scalar>  phis;           // surface geopotential height [m2/s]
    // view_2d<Spack>   z_mid;           // mid-point level altitude [m]
    // view_2d<Spack>   z_int;           // interface level altitude [m]
    view_2d<const Spack>   p_mid;           // mid-point level pressure [Pa]
    view_2d<const Spack>   p_int;           // interface level pressure [Pa]
    view_2d<const Spack>   p_del;           // pressure thickness [Pa]

    // view_1d<Scalar>  pblh;           // PBL height [m]

    view_2d<Spack>  T_mid;          // Temperature [K]
    view_2d<Spack>  qv;             // Water vapor mixing ratio [kg kg-1]
    view_2d<Spack>  qc;             // Cloud mass mixing ratio [kg kg-1]

    // view_2d<Spack>  qi;             // Ice total mass mixing ratio [kg kg-1]
    // view_2d<Spack>  omega;          // vertical pressure velocity [Pa/s]

    // LayoutLeft views for fortran bridging
    // view_2dl<Real>  f_z_mid;
    // view_2dl<Real>  f_z_int;
    view_2dl<Real>  f_p_mid;
    view_2dl<Real>  f_p_int;
    view_2dl<Real>  f_p_del;

    view_2dl<Real>  f_T_mid;
    view_2dl<Real>  f_qv;
    view_2dl<Real>  f_qc;

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
      if (D == ekat::TransposeDirection::c2f) {
        // f_z_mid = view_2dl<Real>("f_z_mid", ncol, pver);
        // f_z_int = view_2dl<Real>("f_z_int", ncol, pver);
        f_p_mid = view_2dl<Real>("f_p_mid", pcol, pver);
        f_p_int = view_2dl<Real>("f_p_int", pcol, pver);
        f_p_del = view_2dl<Real>("f_p_del", pcol, pver);
        f_T_mid = view_2dl<Real>("f_T_mid", pcol, pver);
        f_qv    = view_2dl<Real>("f_qv",    pcol, pver);
        f_qc    = view_2dl<Real>("f_qc",    pcol, pver);

        for (int i=0; i<ncol; ++i) {
          for (int j=0; j<pver; ++j) {
            // f_z_mid(i,j) = z_mid(i, j / Spack::n)[j % Spack::n];
            // f_z_int(i,j) = z_int(i, j / Spack::n)[j % Spack::n];
            f_p_mid(i,j) = p_mid(i, j / Spack::n)[j % Spack::n];
            f_p_int(i,j) = p_int(i, j / Spack::n)[j % Spack::n];
            f_p_del(i,j) = p_del(i, j / Spack::n)[j % Spack::n];
            f_T_mid(i,j) = T_mid(i, j / Spack::n)[j % Spack::n];
            f_qv   (i,j) = qv   (i, j / Spack::n)[j % Spack::n];
            f_qc   (i,j) = qc   (i, j / Spack::n)[j % Spack::n];
          }
        }
      }
      // else {
      //   std::vector<view_2d<Spack>> dev_views = { T_mid,
      //                                             qv,
      //                                             qc,
      //                                           };
      //   ekat::host_to_device( { f_T_mid.data(),
      //                           f_qv   .data(),
      //                           f_qc   .data(),
      //                         }, ncol, pver, dev_views);
      // }
    }
    // -------------------------------------------------------------------------
  };

  struct zm_output_tend {
    zm_output_tend() = default;

    view_2d<const Spack>   tend_s;           // output tendency of water vapor
    view_2d<const Spack>   tend_q;           // output tendency of dry statis energy
    view_1d<const Scalar>  precip;           // surface precipitation [m/s]

    // LayoutLeft views for fortran bridging
    view_2dl<Real>  f_tend_s;
    view_2dl<Real>  f_tend_q;
    // -------------------------------------------------------------------------
    // transpose method for fortran bridging
    template <ekat::TransposeDirection::Enum D>
    void transpose(int pver)
    {
      // if (D == ekat::TransposeDirection::c2f) {
      //   f_p_mid = view_2dl<Real>("f_p_mid", pcol, pver);
      //   for (int i=0; i<ncol; ++i) {
      //     for (int j=0; j<pver; ++j) {
      //       f_p_mid(i,j) = p_mid(i, j / Spack::n)[j % Spack::n];
      //     }
      //   }
      // }
      if (D == ekat::TransposeDirection::f2c) {
        std::vector<view_2d<Spack>> dev_views = { tend_s,
                                                  tend_q,
                                                };
        ekat::host_to_device( { f_tend_s.data(),
                                f_tend_q.data(),
                              }, ncol, pver, dev_views);
      }
    }
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
