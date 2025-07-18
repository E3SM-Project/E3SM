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
struct Functions {
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
  template <typename S> using uview_2dl         = typename ekat::template Unmanaged<view_2dl<S> >;
  template <typename S> using uview_2d_strided  = typename ekat::template Unmanaged<view_2d_strided<S> >;

  // ---------------------------------------------------------------------------
  // Structs

  struct zm_runtime_opt {
    zm_runtime_opt() = default;
  };

  struct zm_input_state {
    zm_input_state() = default;
    // -------------------------------------------------------------------------
    Int ncol;                       // number of columns for current task/chunk
    Int pcol;                       // max number of columns across tasks/chunks
    Real dtime;                     // model phsyics time step [s]
    bool is_first_step;             // flag for first call

    view_1d<const Scalar>  phis;    // surface geopotential height [m2/s]
    // view_2d<Spack>   z_mid;         // mid-point level altitude [m]
    // view_2d<Spack>   z_int;         // interface level altitude [m]
    view_2d<const Spack>  p_mid;    // mid-point level pressure [Pa]
    view_2d<const Spack>  p_int;    // interface level pressure [Pa]
    view_2d<const Spack>  p_del;    // pressure thickness [Pa]

    view_2d<      Spack>  T_mid;    // Temperature [K]
    view_2d<      Spack>  qv;       // Water vapor mixing ratio [kg kg-1]
    view_2d<      Spack>  qc;       // Cloud mass mixing ratio [kg kg-1]
    view_2d<const Spack>  omega;    // vertical pressure velocity [Pa/s]
    view_1d<const Scalar> pblh;     // PBL height [m]

    // LayoutLeft views for fortran bridging
    // view_2dl<Real>  f_z_mid;
    // view_2dl<Real>  f_z_int;
    view_2dl<Real>  f_p_mid;
    view_2dl<Real>  f_p_int;
    view_2dl<Real>  f_p_del;

    view_2dl<Real>  f_T_mid;
    view_2dl<Real>  f_qv;
    view_2dl<Real>  f_qc;
    view_2dl<Real>  f_omega;

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
    void transpose(int pver){
      auto pverp = pver+1;
      if (D == ekat::TransposeDirection::c2f) {
        // f_z_mid = view_2dl<Real>("f_z_mid", ncol, pver);
        // f_z_int = view_2dl<Real>("f_z_int", ncol, pverp);
        f_p_mid = view_2dl<Real>("f_p_mid", pcol, pver);
        f_p_int = view_2dl<Real>("f_p_int", pcol, pverp);
        f_p_del = view_2dl<Real>("f_p_del", pcol, pver);
        f_T_mid = view_2dl<Real>("f_T_mid", pcol, pver);
        f_qv    = view_2dl<Real>("f_qv",    pcol, pver);
        f_qc    = view_2dl<Real>("f_qc",    pcol, pver);
        f_omega = view_2dl<Real>("f_omega", pcol, pver);
        for (int i=0; i<ncol; ++i) {
          for (int j=0; j<pver; ++j) {
            // f_z_mid(i,j) = z_mid(i, j / Spack::n)[j % Spack::n];
            f_p_mid(i,j) = p_mid(i, j / Spack::n)[j % Spack::n];
            f_p_del(i,j) = p_del(i, j / Spack::n)[j % Spack::n];
            f_T_mid(i,j) = T_mid(i, j / Spack::n)[j % Spack::n];
            f_qv   (i,j) = qv   (i, j / Spack::n)[j % Spack::n];
            f_qc   (i,j) = qc   (i, j / Spack::n)[j % Spack::n];
            f_omega(i,j) = omega(i, j / Spack::n)[j % Spack::n];
          }
          for (int j=0; j<pverp; ++j) {
            // f_z_int(i,j) = z_int(i, j / Spack::n)[j % Spack::n];
            f_p_int(i,j) = p_int(i, j / Spack::n)[j % Spack::n];
          }
        }
      }
    }
    // -------------------------------------------------------------------------
  };

  struct zm_output_tend {
    zm_output_tend() = default;

    Int ncol;                       // number of columns for current task/chunk
    Int pcol;                       // max number of columns across tasks/chunks

    view_1d<Scalar> precip;         // surface precipitation [m/s]
    view_2d<Spack>  tend_s;         // output tendency of water vapor
    view_2d<Spack>  tend_q;         // output tendency of dry statis energy
    view_2d<Spack>  prec_flux;      // output convective precipitation flux
    view_2d<Spack>  mass_flux;      // output convective mass flux

    // LayoutLeft views for fortran bridging
    view_2dl<Real>  f_tend_s;
    view_2dl<Real>  f_tend_q;
    view_2dl<Real>  f_prec_flux;
    view_2dl<Real>  f_mass_flux;

    // -------------------------------------------------------------------------
    // transpose method for fortran bridging
    template <ekat::TransposeDirection::Enum D>
    void transpose(int pver) {
      auto pverp = pver+1;
      if (D == ekat::TransposeDirection::c2f) {
        f_tend_s    = view_2dl<Real>("f_tend_s",    pcol, pver);
        f_tend_q    = view_2dl<Real>("f_tend_q",    pcol, pver);
        f_prec_flux = view_2dl<Real>("f_prec_flux", pcol, pverp);
        f_mass_flux = view_2dl<Real>("f_mass_flux", pcol, pverp);
        for (int i=0; i<pcol; ++i) {
          // mid-point level variables
          for (int j=0; j<pver; ++j) {
            f_tend_s   (i,j) = tend_s   (i, j / Spack::n)[j % Spack::n];
            f_tend_q   (i,j) = tend_q   (i, j / Spack::n)[j % Spack::n];
          }
          // interface level variables
          for (int j=0; j<pverp; ++j) {
            f_prec_flux(i,j) = prec_flux(i, j / Spack::n)[j % Spack::n];
            f_mass_flux(i,j) = mass_flux(i, j / Spack::n)[j % Spack::n];
          }
        }
        // sync_to_host?
      }
      if (D == ekat::TransposeDirection::f2c) {
        // sync_to_device?
        for (int i=0; i<pcol; ++i) {
          // mid-point level variables
          for (int j=0; j<pver; ++j) {
            tend_s   (i, j / Spack::n)[j % Spack::n] = f_tend_s   (i,j);
            tend_q   (i, j / Spack::n)[j % Spack::n] = f_tend_q   (i,j);
          }
          // interface level variables
          for (int j=0; j<pverp; ++j) {
            prec_flux(i, j / Spack::n)[j % Spack::n] = f_prec_flux(i,j);
            mass_flux(i, j / Spack::n)[j % Spack::n] = f_mass_flux(i,j);
          }
        }
      }
    }
  };

  struct zm_output_diag {
    zm_output_diag() = default;
  };

  // Structure for storing local variables initialized using the ATMBufferManager
  struct zm_buffer_data {

    static constexpr int num_1d_scalr_views = 1; // number of 1D variables
    static constexpr int num_2d_midlv_views = 2; // number of 2D variables on mid-point levels
    static constexpr int num_2d_intfc_views = 2; // number of 2D variables on interface levels

    uview_1d<Scalar>  precip;       // surface precipitation [m/s]
    uview_2d<Spack>   tend_s;       // output tendency of water vapor
    uview_2d<Spack>   tend_q;       // output tendency of dry statis energy
    uview_2d<Spack>   prec_flux;    // output convective precipitation flux
    uview_2d<Spack>   mass_flux;    // output convective mass flux

    uview_2dl<Real>   f_tend_s;       // output tendency of water vapor
    uview_2dl<Real>   f_tend_q;       // output tendency of dry statis energy

    void init(int ncol,int pver) {
      auto pverp = pver+1;
      // 1D scalar variables
      for (int i=0; i<ncol; ++i) {
        precip(i) = 0;
      }
      // mid-point level variables
      for (int i=0; i<ncol; ++i) {
        for (int j=0; j<pver; ++j) {
          tend_s(i,j/Spack::n)[j%Spack::n] = 0;
          tend_q(i,j/Spack::n)[j%Spack::n] = 0;
          f_tend_s(i,j) = 0;
          f_tend_q(i,j) = 0;
        }
      }
      // interface level variables
      for (int i=0; i<ncol; ++i) {
        for (int j=0; j<pverp; ++j) {
          prec_flux(i,j) = 0;
          mass_flux(i,j) = 0;
        }
      }
    };
  };

  // ---------------------------------------------------------------------------
  // Functions

  // static Int zm_main()

}; // struct Functions

} // namespace zm
} // namespace scream

#endif // ZM_FUNCTIONS_HPP
