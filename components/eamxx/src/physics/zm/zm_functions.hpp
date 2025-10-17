#ifndef ZM_FUNCTIONS_HPP
#define ZM_FUNCTIONS_HPP

#include "share/physics/physics_constants.hpp"
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
  template <typename S> using uview_2dh         = typename ekat::template Unmanaged<view_2dl<S>>::HostMirror;
  template <typename S> using uview_1dh         = typename ekat::template Unmanaged<view_1d<S>>::HostMirror;

  // ---------------------------------------------------------------------------
  // Structs

  struct zm_runtime_opt {
    zm_runtime_opt() = default;
    bool apply_tendencies = false;

    void load_runtime_options(ekat::ParameterList& params) {
      apply_tendencies = params.get<bool>("apply_tendencies", apply_tendencies);
    }
  };

  struct zm_input_state {
    zm_input_state() = default;
    // -------------------------------------------------------------------------
    Real dtime;                     // model phsyics time step [s]
    bool is_first_step;             // flag for first call

    static constexpr int num_1d_scalr_views   = 1;  // number of 1D scalar variables
    static constexpr int num_1d_intgr_views   = 0;  // number of 1D integer variables
    static constexpr int num_2d_midlv_c_views = 2;  // number of 2D variables on mid-point levels
    static constexpr int num_2d_intfc_c_views = 1;  // number of 2D variables on interface levels
    static constexpr int num_2d_midlv_f_views = 9;  // number of 2D variables on mid-point levels
    static constexpr int num_2d_intfc_f_views = 2;  // number of 2D variables on interface levels

    uview_1d<     Scalar> tpert;    // temperature perturbation at top of PBL

    uview_2d<     Spack>  z_mid;    // mid-point level altitude     [m]
    uview_2d<     Spack>  z_del;    // altitude thickness           [m]

    uview_2d<     Spack>  z_int;    // interface level altitude     [m]

    // variables we get from the field manager
    view_1d<const Scalar> phis;     // surface geopotential height  [m2/s]
    view_2d<const Spack>  p_mid;    // mid-point level pressure     [Pa]
    view_2d<const Spack>  p_int;    // interface level pressure     [Pa]
    view_2d<const Spack>  p_del;    // pressure thickness           [Pa]
    view_2d<      Spack>  T_mid;    // temperature                  [K]
    view_2d<      Spack>  qv;       // water vapor mixing ratio     [kg kg-1]
    view_2d<      Spack>  uwind;    // zonal wind                   [m/s]
    view_2d<      Spack>  vwind;    // meridional wind              [m/s]
    view_2d<const Spack>  omega;    // vertical pressure velocity   [Pa/s]
    view_2d<const Spack>  cldfrac;  // total cloud fraction
    view_1d<const Scalar> pblh;     // PBL height                   [m]
    view_1d<const Scalar> landfrac; // land area fraction

    // unmanaged LayoutLeft views for fortran bridging
    uview_2dl<Real> f_z_mid;
    uview_2dl<Real> f_p_mid;
    uview_2dl<Real> f_p_del;
    uview_2dl<Real> f_T_mid;
    uview_2dl<Real> f_qv;
    uview_2dl<Real> f_uwind;
    uview_2dl<Real> f_vwind;
    uview_2dl<Real> f_omega;
    uview_2dl<Real> f_cldfrac;
    uview_2dl<Real> f_z_int;
    uview_2dl<Real> f_p_int;

    // host mirror versions of ZM interface variables
    uview_2dh<Real> h_z_mid;
    uview_2dh<Real> h_p_mid;
    uview_2dh<Real> h_p_del;
    uview_2dh<Real> h_T_mid;
    uview_2dh<Real> h_qv;
    uview_2dh<Real> h_uwind;
    uview_2dh<Real> h_vwind;
    uview_2dh<Real> h_omega;
    uview_2dh<Real> h_cldfrac;
    uview_2dh<Real> h_z_int;
    uview_2dh<Real> h_p_int;

    uview_1dh<Scalar> h_phis;
    uview_1dh<Scalar> h_pblh;
    uview_1dh<Scalar> h_tpert;
    uview_1dh<Scalar> h_landfrac;

    // -------------------------------------------------------------------------
    // transpose method for fortran bridging
    template <ekat::TransposeDirection::Enum D>
    void transpose(int ncol_in, int pver_in) {
      auto pverp = pver_in+1;
      if (D == ekat::TransposeDirection::c2f) {
        for (int i=0; i<ncol_in; ++i) {
          for (int j=0; j<pver_in; ++j) {
            f_z_mid   (i,j) = z_mid   (i,j/Spack::n)[j%Spack::n];
            f_p_mid   (i,j) = p_mid   (i,j/Spack::n)[j%Spack::n];
            f_p_del   (i,j) = p_del   (i,j/Spack::n)[j%Spack::n];
            f_T_mid   (i,j) = T_mid   (i,j/Spack::n)[j%Spack::n];
            f_qv      (i,j) = qv      (i,j/Spack::n)[j%Spack::n];
            f_uwind   (i,j) = uwind   (i,j/Spack::n)[j%Spack::n];
            f_vwind   (i,j) = vwind   (i,j/Spack::n)[j%Spack::n];
            f_omega   (i,j) = omega   (i,j/Spack::n)[j%Spack::n];
            f_cldfrac (i,j) = cldfrac (i,j/Spack::n)[j%Spack::n];
          }
          for (int j=0; j<pverp; ++j) {
            f_z_int   (i,j) = z_int   (i,j/Spack::n)[j%Spack::n];
            f_p_int   (i,j) = p_int   (i,j/Spack::n)[j%Spack::n];
          }
        }
        //----------------------------------------------------------------------
        // copy to host mirrors
        Kokkos::deep_copy(h_z_mid,    f_z_mid);
        Kokkos::deep_copy(h_p_mid,    f_p_mid);
        Kokkos::deep_copy(h_p_del,    f_p_del);
        Kokkos::deep_copy(h_T_mid,    f_T_mid);
        Kokkos::deep_copy(h_qv,       f_qv);
        Kokkos::deep_copy(h_uwind,    f_uwind);
        Kokkos::deep_copy(h_vwind,    f_vwind);
        Kokkos::deep_copy(h_omega,    f_omega);
        Kokkos::deep_copy(h_cldfrac,  f_cldfrac);
        Kokkos::deep_copy(h_z_int,    f_z_int);
        Kokkos::deep_copy(h_p_int,    f_p_int);
        Kokkos::deep_copy(h_phis,     phis);
        Kokkos::deep_copy(h_pblh,     pblh);
        Kokkos::deep_copy(h_tpert,    tpert);
        Kokkos::deep_copy(h_landfrac, landfrac);
      }
    }
    // -------------------------------------------------------------------------
  };

  struct zm_output_tend {
    zm_output_tend() = default;

    static constexpr int num_1d_intgr_views   = 1;  // number of 1D integer variables
    static constexpr int num_1d_scalr_views   = 3;  // number of 1D scalar variables
    static constexpr int num_2d_midlv_c_views = 6;  // number of 2D variables on mid-point levels
    static constexpr int num_2d_intfc_c_views = 3;  // number of 2D variables on interface levels
    static constexpr int num_2d_midlv_f_views = 6;  // number of 2D variables on mid-point levels
    static constexpr int num_2d_intfc_f_views = 3;  // number of 2D variables on interface levels

    uview_1d<Int>    activity;       // integer deep convection activity flag

    uview_1d<Scalar> prec;           // surface precipitation                   [m/s]
    uview_1d<Scalar> snow;           // surface snow                            [m/s]
    uview_1d<Scalar> cape;           // convective available potential energy   [J]

    uview_2d<Spack>  tend_s;         // output tendency of dry static energy    []
    uview_2d<Spack>  tend_qv;        // output tendency of water vapor          []
    uview_2d<Spack>  tend_u;         // output tendency of zonal wind           []
    uview_2d<Spack>  tend_v;         // output tendency of meridional wind      []
    uview_2d<Spack>  rain_prod;      // rain production rate
    uview_2d<Spack>  snow_prod;      // snow production rate

    uview_2d<Spack>  prec_flux;      // output convective precipitation flux    []
    uview_2d<Spack>  snow_flux;      // output convective precipitation flux    []
    uview_2d<Spack>  mass_flux;      // output convective mass flux             []

    // LayoutLeft views for fortran bridging
    uview_2dl<Real>  f_tend_s;
    uview_2dl<Real>  f_tend_qv;
    uview_2dl<Real>  f_tend_u;
    uview_2dl<Real>  f_tend_v;
    uview_2dl<Real>  f_rain_prod;
    uview_2dl<Real>  f_snow_prod;
    uview_2dl<Real>  f_prec_flux;
    uview_2dl<Real>  f_snow_flux;
    uview_2dl<Real>  f_mass_flux;

    // host mirror versions of ZM interface variables
    uview_2dh<Real> h_tend_s;
    uview_2dh<Real> h_tend_qv;
    uview_2dh<Real> h_tend_u;
    uview_2dh<Real> h_tend_v;
    uview_2dh<Real> h_rain_prod;
    uview_2dh<Real> h_snow_prod;
    uview_2dh<Real> h_prec_flux;
    uview_2dh<Real> h_snow_flux;
    uview_2dh<Real> h_mass_flux;

    uview_1dh<Scalar> h_prec;
    uview_1dh<Scalar> h_snow;
    uview_1dh<Scalar> h_cape;
    uview_1dh<Int>    h_activity;

    // -------------------------------------------------------------------------
    // transpose method for fortran bridging
    template <ekat::TransposeDirection::Enum D>
    void transpose(int ncol_in, int pver_in) {
      auto pverp = pver_in+1;
      if (D == ekat::TransposeDirection::c2f) {
        for (int i=0; i<ncol_in; ++i) {
          // mid-point level variables
          for (int j=0; j<pver_in; ++j) {
            f_tend_s   (i,j) = tend_s   (i,j/Spack::n)[j%Spack::n];
            f_tend_qv  (i,j) = tend_qv  (i,j/Spack::n)[j%Spack::n];
            f_tend_u   (i,j) = tend_u   (i,j/Spack::n)[j%Spack::n];
            f_tend_v   (i,j) = tend_v   (i,j/Spack::n)[j%Spack::n];
            f_rain_prod(i,j) = rain_prod(i,j/Spack::n)[j%Spack::n];
            f_snow_prod(i,j) = snow_prod(i,j/Spack::n)[j%Spack::n];
          }
          // interface level variables
          for (int j=0; j<pverp; ++j) {
            f_prec_flux(i,j) = prec_flux(i,j/Spack::n)[j%Spack::n];
            f_snow_flux(i,j) = snow_flux(i,j/Spack::n)[j%Spack::n];
            f_mass_flux(i,j) = mass_flux(i,j/Spack::n)[j%Spack::n];
          }
        }
        //----------------------------------------------------------------------
        // copy to host mirror
        Kokkos::deep_copy(h_tend_s,   f_tend_s);
        Kokkos::deep_copy(h_tend_qv,  f_tend_qv);
        Kokkos::deep_copy(h_tend_u,   f_tend_u);
        Kokkos::deep_copy(h_tend_v,   f_tend_v);
        Kokkos::deep_copy(h_rain_prod,f_rain_prod);
        Kokkos::deep_copy(h_snow_prod,f_snow_prod);
        Kokkos::deep_copy(h_prec_flux,f_prec_flux);
        Kokkos::deep_copy(h_snow_flux,f_snow_flux);
        Kokkos::deep_copy(h_mass_flux,f_mass_flux);
        Kokkos::deep_copy(h_prec,     prec);
        Kokkos::deep_copy(h_snow,     snow);
        Kokkos::deep_copy(h_cape,     cape);
        Kokkos::deep_copy(h_activity, activity);
      }
      if (D == ekat::TransposeDirection::f2c) {
        //----------------------------------------------------------------------
        // copy to host mirror
        Kokkos::deep_copy(f_tend_s,   h_tend_s);
        Kokkos::deep_copy(f_tend_qv,  h_tend_qv);
        Kokkos::deep_copy(f_tend_u,   h_tend_u);
        Kokkos::deep_copy(f_tend_v,   h_tend_v);
        Kokkos::deep_copy(f_rain_prod,h_rain_prod);
        Kokkos::deep_copy(f_snow_prod,h_snow_prod);
        Kokkos::deep_copy(f_prec_flux,h_prec_flux);
        Kokkos::deep_copy(f_snow_flux,h_snow_flux);
        Kokkos::deep_copy(f_mass_flux,h_mass_flux);
        Kokkos::deep_copy(prec,       h_prec);
        Kokkos::deep_copy(snow,       h_snow);
        Kokkos::deep_copy(cape,       h_cape);
        Kokkos::deep_copy(activity,   h_activity);
        //----------------------------------------------------------------------
        for (int i=0; i<ncol_in; ++i) {
          // mid-point level variables
          for (int j=0; j<pver_in; ++j) {
            tend_s   (i,j/Spack::n)[j%Spack::n] = f_tend_s   (i,j);
            tend_qv  (i,j/Spack::n)[j%Spack::n] = f_tend_qv  (i,j);
            tend_u   (i,j/Spack::n)[j%Spack::n] = f_tend_u   (i,j);
            tend_v   (i,j/Spack::n)[j%Spack::n] = f_tend_v   (i,j);
            rain_prod(i,j/Spack::n)[j%Spack::n] = f_rain_prod(i,j);
            snow_prod(i,j/Spack::n)[j%Spack::n] = f_snow_prod(i,j);
          }
          // interface level variables
          for (int j=0; j<pverp; ++j) {
            prec_flux(i,j/Spack::n)[j%Spack::n] = f_prec_flux(i,j);
            snow_flux(i,j/Spack::n)[j%Spack::n] = f_snow_flux(i,j);
            mass_flux(i,j/Spack::n)[j%Spack::n] = f_mass_flux(i,j);
          }
        }
      }
    };
    // -------------------------------------------------------------------------
    void init(int ncol_in, int pver_in) {
      Real init_fill_value = -999;
      // 1D scalar variables
      for (int i=0; i<ncol_in; ++i) {
        prec(i) = init_fill_value;
        snow(i) = init_fill_value;
        cape(i) = init_fill_value;
        activity(i) = -1;
      }
      // mid-point level variables
      for (int i=0; i<ncol_in; ++i) {
        for (int j=0; j<pver_in; ++j) {
          tend_s   (i,j/Spack::n)[j%Spack::n] = init_fill_value;
          tend_qv  (i,j/Spack::n)[j%Spack::n] = init_fill_value;
          tend_u   (i,j/Spack::n)[j%Spack::n] = init_fill_value;
          tend_v   (i,j/Spack::n)[j%Spack::n] = init_fill_value;
          rain_prod(i,j/Spack::n)[j%Spack::n] = init_fill_value;
          snow_prod(i,j/Spack::n)[j%Spack::n] = init_fill_value;
          f_tend_s   (i,j) = init_fill_value;
          f_tend_qv  (i,j) = init_fill_value;
          f_tend_u   (i,j) = init_fill_value;
          f_tend_v   (i,j) = init_fill_value;
          f_rain_prod(i,j) = init_fill_value;
          f_snow_prod(i,j) = init_fill_value;
        }
      }
      auto pverp = pver_in+1;
      // interface level variables
      for (int i=0; i<ncol_in; ++i) {
        for (int j=0; j<pverp; ++j) {
          prec_flux(i,j/Spack::n)[j%Spack::n] = init_fill_value;
          snow_flux(i,j/Spack::n)[j%Spack::n] = init_fill_value;
          mass_flux(i,j/Spack::n)[j%Spack::n] = init_fill_value;
          f_prec_flux(i,j) = init_fill_value;
          f_snow_flux(i,j) = init_fill_value;
          f_mass_flux(i,j) = init_fill_value;
        }
      }
    };
    // -------------------------------------------------------------------------
  };

  struct zm_output_diag {
    zm_output_diag() = default;
  };

  // ---------------------------------------------------------------------------
  // Functions

  // static Int zm_main()

}; // struct Functions

} // namespace zm
} // namespace scream

#endif // ZM_FUNCTIONS_HPP
