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
  // template <typename S> using uview_2d_strided  = typename ekat::template Unmanaged<view_2d_strided<S> >;
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

    static constexpr int num_dvc_1d_intgr = 0;  // number of device 1D integer views
    static constexpr int num_dvc_1d_scalr = 1;  // number of device 1D scalar views
    static constexpr int num_dvc_2d_midlv = 2;  // number of device 2D views on mid-point levels
    static constexpr int num_dvc_2d_intfc = 1;  // number of device 2D views on interface levels

    // static constexpr int num_lol_1d_intgr = 0;  // number of layout-left 1D integer views
    // static constexpr int num_lol_1d_scalr = 0;  // number of layout-left 1D scalar views
    // static constexpr int num_lol_2d_midlv = 9;  // number of layout-left 2D views on mid-point levels
    // static constexpr int num_lol_2d_intfc = 2;  // number of layout-left 2D views on interface levels

    static constexpr int num_hst_1d_intgr = 0;  // number of host 1D integer views
    static constexpr int num_hst_1d_scalr = 4;  // number of host 1D scalar views
    static constexpr int num_hst_2d_midlv = 9;  // number of host 2D views on mid-point levels
    static constexpr int num_hst_2d_intfc = 2;  // number of host 2D views on interface levels

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

    // LayoutLeft views for fortran bridging
    // uview_2dl<Real> f_z_mid;
    // uview_2dl<Real> f_p_mid;
    // uview_2dl<Real> f_p_del;
    // uview_2dl<Real> f_T_mid;
    // uview_2dl<Real> f_qv;
    // uview_2dl<Real> f_uwind;
    // uview_2dl<Real> f_vwind;
    // uview_2dl<Real> f_omega;
    // uview_2dl<Real> f_cldfrac;
    // uview_2dl<Real> f_z_int;
    // uview_2dl<Real> f_p_int;

    // host mirror versions of ZM interface variables
    uview_1dh<Scalar> h_phis;
    uview_1dh<Scalar> h_pblh;
    uview_1dh<Scalar> h_tpert;
    uview_1dh<Scalar> h_landfrac;

    uview_2dh<Real>   h_z_mid;
    uview_2dh<Real>   h_p_mid;
    uview_2dh<Real>   h_p_del;
    uview_2dh<Real>   h_T_mid;
    uview_2dh<Real>   h_qv;
    uview_2dh<Real>   h_uwind;
    uview_2dh<Real>   h_vwind;
    uview_2dh<Real>   h_omega;
    uview_2dh<Real>   h_cldfrac;

    uview_2dh<Real>   h_z_int;
    uview_2dh<Real>   h_p_int;

    // -------------------------------------------------------------------------
    // transpose method for fortran bridging
    template <ekat::TransposeDirection::Enum D>
    void transpose(int ncol, int nlev_mid) {

      // auto pverp = pver_in+1;
      // if (D == ekat::TransposeDirection::c2f) {
      //   for (int i=0; i<ncol_in; ++i) {
      //     for (int j=0; j<pver_in; ++j) {
      //       f_z_mid   (i,j) = z_mid   (i,j/Spack::n)[j%Spack::n];
      //       f_p_mid   (i,j) = p_mid   (i,j/Spack::n)[j%Spack::n];
      //       f_p_del   (i,j) = p_del   (i,j/Spack::n)[j%Spack::n];
      //       f_T_mid   (i,j) = T_mid   (i,j/Spack::n)[j%Spack::n];
      //       f_qv      (i,j) = qv      (i,j/Spack::n)[j%Spack::n];
      //       f_uwind   (i,j) = uwind   (i,j/Spack::n)[j%Spack::n];
      //       f_vwind   (i,j) = vwind   (i,j/Spack::n)[j%Spack::n];
      //       f_omega   (i,j) = omega   (i,j/Spack::n)[j%Spack::n];
      //       f_cldfrac (i,j) = cldfrac (i,j/Spack::n)[j%Spack::n];
      //     }
      //     for (int j=0; j<pverp; ++j) {
      //       f_z_int   (i,j) = z_int   (i,j/Spack::n)[j%Spack::n];
      //       f_p_int   (i,j) = p_int   (i,j/Spack::n)[j%Spack::n];
      //     }
      //   }
      //   //----------------------------------------------------------------------
      //   // copy to host mirrors
      //   Kokkos::deep_copy(h_z_mid,    f_z_mid);
      //   Kokkos::deep_copy(h_p_mid,    f_p_mid);
      //   Kokkos::deep_copy(h_p_del,    f_p_del);
      //   Kokkos::deep_copy(h_T_mid,    f_T_mid);
      //   Kokkos::deep_copy(h_qv,       f_qv);
      //   Kokkos::deep_copy(h_uwind,    f_uwind);
      //   Kokkos::deep_copy(h_vwind,    f_vwind);
      //   Kokkos::deep_copy(h_omega,    f_omega);
      //   Kokkos::deep_copy(h_cldfrac,  f_cldfrac);
      //   Kokkos::deep_copy(h_z_int,    f_z_int);
      //   Kokkos::deep_copy(h_p_int,    f_p_int);
      //   Kokkos::deep_copy(h_phis,     phis);
      //   Kokkos::deep_copy(h_pblh,     pblh);
      //   Kokkos::deep_copy(h_tpert,    tpert);
      //   Kokkos::deep_copy(h_landfrac, landfrac);
      // }

      auto nlev_int = nlev_mid+1;
      if (D == ekat::TransposeDirection::c2f) {
        ekat::device_to_host({h_phis.data()},     ncol,           std::vector< view_1d<const Scalar>>{phis});
        ekat::device_to_host({h_pblh.data()},     ncol,           std::vector< view_1d<const Scalar>>{pblh});
        ekat::device_to_host({h_tpert.data()},    ncol,           std::vector<uview_1d<      Scalar>>{tpert});
        ekat::device_to_host({h_landfrac.data()}, ncol,           std::vector< view_1d<const Scalar>>{landfrac});
        ekat::device_to_host({h_z_mid.data()},    ncol, nlev_mid, std::vector<uview_2d<const Spack >>{z_mid});
        ekat::device_to_host({h_p_mid.data()},    ncol, nlev_mid, std::vector< view_2d<const Spack >>{p_mid});
        ekat::device_to_host({h_p_del.data()},    ncol, nlev_mid, std::vector< view_2d<const Spack >>{p_del});
        ekat::device_to_host({h_T_mid.data()},    ncol, nlev_mid, std::vector< view_2d<      Spack >>{T_mid});
        ekat::device_to_host({h_qv.data()},       ncol, nlev_mid, std::vector< view_2d<      Spack >>{qv});
        ekat::device_to_host({h_uwind.data()},    ncol, nlev_mid, std::vector< view_2d<      Spack >>{uwind});
        ekat::device_to_host({h_vwind.data()},    ncol, nlev_mid, std::vector< view_2d<      Spack >>{vwind});
        ekat::device_to_host({h_omega.data()},    ncol, nlev_mid, std::vector< view_2d<const Spack >>{omega});
        ekat::device_to_host({h_cldfrac.data()},  ncol, nlev_mid, std::vector< view_2d<const Spack >>{cldfrac});
        ekat::device_to_host({h_z_int.data()},    ncol, nlev_int, std::vector<uview_2d<      Spack >>{z_int});
        ekat::device_to_host({h_p_int.data()},    ncol, nlev_int, std::vector< view_2d<const Spack >>{p_int});
      }

    }
    // -------------------------------------------------------------------------
  };

  struct zm_output_tend {
    zm_output_tend() = default;

    static constexpr int num_dvc_1d_intgr = 1;  // number of device 1D integer views
    static constexpr int num_dvc_1d_scalr = 3;  // number of device 1D scalar views
    static constexpr int num_dvc_2d_midlv = 6;  // number of device 2D views on mid-point levels
    static constexpr int num_dvc_2d_intfc = 3;  // number of device 2D views on interface levels

    // static constexpr int num_lol_1d_intgr = 0;  // number of layout-left 1D integer views
    // static constexpr int num_lol_1d_scalr = 0;  // number of layout-left 1D scalar views
    // static constexpr int num_lol_2d_midlv = 6;  // number of layout-left 2D views on mid-point levels
    // static constexpr int num_lol_2d_intfc = 3;  // number of layout-left 2D views on interface levels

    static constexpr int num_hst_1d_intgr = 1;  // number of host 1D integer views
    static constexpr int num_hst_1d_scalr = 3;  // number of host 1D scalar views
    static constexpr int num_hst_2d_midlv = 6;  // number of host 2D views on mid-point levels
    static constexpr int num_hst_2d_intfc = 3;  // number of host 2D views on interface levels

    uview_1d<Int>    activity;       // integer deep convection activity flag

    uview_1d<Scalar> prec;           // surface precipitation                   [m/s]
    uview_1d<Scalar> snow;           // surface snow                            [m/s]
    uview_1d<Scalar> cape;           // convective available potential energy   [J]

    uview_2d<Spack>  tend_t;         // output tendency of dry static energy    []
    uview_2d<Spack>  tend_qv;        // output tendency of water vapor          []
    uview_2d<Spack>  tend_u;         // output tendency of zonal wind           []
    uview_2d<Spack>  tend_v;         // output tendency of meridional wind      []
    uview_2d<Spack>  rain_prod;      // rain production rate
    uview_2d<Spack>  snow_prod;      // snow production rate

    uview_2d<Spack>  prec_flux;      // output convective precipitation flux    []
    uview_2d<Spack>  snow_flux;      // output convective precipitation flux    []
    uview_2d<Spack>  mass_flux;      // output convective mass flux             []

    // LayoutLeft views for fortran bridging
    // uview_2dl<Real>  f_tend_t;
    // uview_2dl<Real>  f_tend_qv;
    // uview_2dl<Real>  f_tend_u;
    // uview_2dl<Real>  f_tend_v;
    // uview_2dl<Real>  f_rain_prod;
    // uview_2dl<Real>  f_snow_prod;

    // uview_2dl<Real>  f_prec_flux;
    // uview_2dl<Real>  f_snow_flux;
    // uview_2dl<Real>  f_mass_flux;

    // host versions of ZM interface variables
    uview_1dh<Int>    h_activity;

    uview_1dh<Scalar> h_prec;
    uview_1dh<Scalar> h_snow;
    uview_1dh<Scalar> h_cape;

    uview_2dh<Real>   h_tend_t;
    uview_2dh<Real>   h_tend_qv;
    uview_2dh<Real>   h_tend_u;
    uview_2dh<Real>   h_tend_v;
    uview_2dh<Real>   h_rain_prod;
    uview_2dh<Real>   h_snow_prod;

    uview_2dh<Real>   h_prec_flux;
    uview_2dh<Real>   h_snow_flux;
    uview_2dh<Real>   h_mass_flux;

    // -------------------------------------------------------------------------
    // transpose method for fortran bridging
    template <ekat::TransposeDirection::Enum D>
    void transpose(int ncol, int nlev_mid) {

      // if (D == ekat::TransposeDirection::c2f) {
      //   for (int i=0; i<ncol; ++i) {
      //     // mid-point level variables
      //     for (int j=0; j<nlev_mid; ++j) {
      //       f_tend_t   (i,j) = tend_t   (i,j/Spack::n)[j%Spack::n];
      //       f_tend_qv  (i,j) = tend_qv  (i,j/Spack::n)[j%Spack::n];
      //       f_tend_u   (i,j) = tend_u   (i,j/Spack::n)[j%Spack::n];
      //       f_tend_v   (i,j) = tend_v   (i,j/Spack::n)[j%Spack::n];
      //       f_rain_prod(i,j) = rain_prod(i,j/Spack::n)[j%Spack::n];
      //       f_snow_prod(i,j) = snow_prod(i,j/Spack::n)[j%Spack::n];
      //     }
      //     // interface level variables
      //     for (int j=0; j<pverp; ++j) {
      //       f_prec_flux(i,j) = prec_flux(i,j/Spack::n)[j%Spack::n];
      //       f_snow_flux(i,j) = snow_flux(i,j/Spack::n)[j%Spack::n];
      //       f_mass_flux(i,j) = mass_flux(i,j/Spack::n)[j%Spack::n];
      //     }
      //   }
      //   //----------------------------------------------------------------------
      //   // copy to host mirror
      //   Kokkos::deep_copy(h_tend_t,   f_tend_t);
      //   Kokkos::deep_copy(h_tend_qv,  f_tend_qv);
      //   Kokkos::deep_copy(h_tend_u,   f_tend_u);
      //   Kokkos::deep_copy(h_tend_v,   f_tend_v);
      //   Kokkos::deep_copy(h_rain_prod,f_rain_prod);
      //   Kokkos::deep_copy(h_snow_prod,f_snow_prod);
      //   Kokkos::deep_copy(h_prec_flux,f_prec_flux);
      //   Kokkos::deep_copy(h_snow_flux,f_snow_flux);
      //   Kokkos::deep_copy(h_mass_flux,f_mass_flux);
      //   Kokkos::deep_copy(h_prec,     prec);
      //   Kokkos::deep_copy(h_snow,     snow);
      //   Kokkos::deep_copy(h_cape,     cape);
      //   Kokkos::deep_copy(h_activity, activity);
      //   // ----------------------------------------------------------------------
      // }

      // if (D == ekat::TransposeDirection::f2c) {
        // // copy to host mirror
        // Kokkos::deep_copy(f_tend_t,   h_tend_t);
        // Kokkos::deep_copy(f_tend_qv,  h_tend_qv);
        // Kokkos::deep_copy(f_tend_u,   h_tend_u);
        // Kokkos::deep_copy(f_tend_v,   h_tend_v);
        // Kokkos::deep_copy(f_rain_prod,h_rain_prod);
        // Kokkos::deep_copy(f_snow_prod,h_snow_prod);
        // Kokkos::deep_copy(f_prec_flux,h_prec_flux);
        // Kokkos::deep_copy(f_snow_flux,h_snow_flux);
        // Kokkos::deep_copy(f_mass_flux,h_mass_flux);
        // Kokkos::deep_copy(prec,       h_prec);
        // Kokkos::deep_copy(snow,       h_snow);
        // Kokkos::deep_copy(cape,       h_cape);
        // Kokkos::deep_copy(activity,   h_activity);
        // //----------------------------------------------------------------------
        // for (int i=0; i<ncol_in; ++i) {
        //   // mid-point level variables
        //   for (int j=0; j<pver_in; ++j) {
        //     tend_t   (i,j/Spack::n)[j%Spack::n] = f_tend_t   (i,j);
        //     tend_qv  (i,j/Spack::n)[j%Spack::n] = f_tend_qv  (i,j);
        //     tend_u   (i,j/Spack::n)[j%Spack::n] = f_tend_u   (i,j);
        //     tend_v   (i,j/Spack::n)[j%Spack::n] = f_tend_v   (i,j);
        //     rain_prod(i,j/Spack::n)[j%Spack::n] = f_rain_prod(i,j);
        //     snow_prod(i,j/Spack::n)[j%Spack::n] = f_snow_prod(i,j);
        //   }
        //   // interface level variables
        //   for (int j=0; j<pverp; ++j) {
        //     prec_flux(i,j/Spack::n)[j%Spack::n] = f_prec_flux(i,j);
        //     snow_flux(i,j/Spack::n)[j%Spack::n] = f_snow_flux(i,j);
        //     mass_flux(i,j/Spack::n)[j%Spack::n] = f_mass_flux(i,j);
        //   }
        // }
      // }


      auto nlev_int = nlev_mid+1;
      if (D == ekat::TransposeDirection::f2c) {

        // std::vector<uview_1dh<Scalar>> tmp_vector_d[1];

        // std::vector<uview_1dh<Scalar>> tmp_vector_d = {h_prec};
        // std::vector< view_1d <Scalar>> tmp_vector_h = {prec.data()};

        // std::vector<uview_1d<Scalar>> tmp(1);
        // ekat::host_to_device({h_prec.data()}, ncol, tmp);
        // Kokkos::deep_copy(prec,tmp[0]);

        ekat::host_to_device({h_prec.data()}, ncol, std::vector<uview_1d<Scalar>> {prec});


        // prec = tmp_vector[0];
        // ekat::host_to_device({prec},     ncol,           std::vector<uview_1dh<Scalar>>{h_prec});
        // ekat::host_to_device({h_prec.data()},     ncol,           std::vector< view_1d <Scalar>>{prec});

        // ekat::host_to_device({h_activity}, ncol,           std::vector<uview_1d <Int>>   {activity});
        // ekat::host_to_device({h_prec},     ncol,           std::vector<uview_1d <Scalar>>{prec});
        // ekat::host_to_device({h_snow},     ncol,           std::vector<uview_1d <Scalar>>{snow});
        // ekat::host_to_device({h_cape},     ncol,           std::vector<uview_1d <Scalar>>{cape});
        // ekat::host_to_device({h_tend_t},   ncol, nlev_mid, std::vector<uview_2d <Spack>> {tend_t});
        // ekat::host_to_device({h_tend_qv},  ncol, nlev_mid, std::vector<uview_2d <Spack>> {tend_qv});
        // ekat::host_to_device({h_tend_u},   ncol, nlev_mid, std::vector<uview_2d <Spack>> {tend_u});
        // ekat::host_to_device({h_tend_v},   ncol, nlev_mid, std::vector<uview_2d <Spack>> {tend_v});
        // ekat::host_to_device({h_rain_prod},ncol, nlev_mid, std::vector<uview_2d <Spack>> {rain_prod});
        // ekat::host_to_device({h_snow_prod},ncol, nlev_mid, std::vector<uview_2d <Spack>> {snow_prod});
        // ekat::host_to_device({h_prec_flux},ncol, nlev_int, std::vector<uview_2d <Spack>> {prec_flux});
        // ekat::host_to_device({h_snow_flux},ncol, nlev_int, std::vector<uview_2d <Spack>> {snow_flux});
        // ekat::host_to_device({h_mass_flux},ncol, nlev_int, std::vector<uview_2d <Spack>> {mass_flux});
      }


      // auto nlev_int = nlev_mid+1;
      // if (D == ekat::TransposeDirection::f2c) {
      //   ekat::host_to_device({h_activity}, ncol,           std::vector<uview_1d <Int>>   {activity});
      //   ekat::host_to_device({h_prec},     ncol,           std::vector< view_1d <Scalar>>{prec});
      //   ekat::host_to_device({h_snow},     ncol,           std::vector< view_1d <Scalar>>{snow});
      //   ekat::host_to_device({h_cape},     ncol,           std::vector< view_1d <Scalar>>{cape});
      //   ekat::host_to_device({h_tend_t},   ncol, nlev_mid, std::vector<uview_2d <Spack>> {tend_t});
      //   ekat::host_to_device({h_tend_qv},  ncol, nlev_mid, std::vector<uview_2d <Spack>> {tend_qv});
      //   ekat::host_to_device({h_tend_u},   ncol, nlev_mid, std::vector<uview_2d <Spack>> {tend_u});
      //   ekat::host_to_device({h_tend_v},   ncol, nlev_mid, std::vector<uview_2d <Spack>> {tend_v});
      //   ekat::host_to_device({h_rain_prod},ncol, nlev_mid, std::vector<uview_2d <Spack>> {rain_prod});
      //   ekat::host_to_device({h_snow_prod},ncol, nlev_mid, std::vector<uview_2d <Spack>> {snow_prod});
      //   ekat::host_to_device({h_prec_flux},ncol, nlev_int, std::vector<uview_2d <Spack>> {prec_flux});
      //   ekat::host_to_device({h_snow_flux},ncol, nlev_int, std::vector<uview_2d <Spack>> {snow_flux});
      //   ekat::host_to_device({h_mass_flux},ncol, nlev_int, std::vector<uview_2d <Spack>> {mass_flux});
      // }

        // std::vector<uview_1dh<Int>>    h_int_vector_1d = {h_activity};
        // std::vector<uview_1d <Int>>    d_int_vector_1d = {activity};
        // std::vector<uview_1dh<Scalar>> h_scl_vector_1d = {h_prec,
        //                                                   h_snow,
        //                                                   h_cape};
        // std::vector< view_1d <Scalar>> d_scl_vector_1d = {prec,
        //                                                   snow,
        //                                                   cape};
        // std::vector<uview_2dh<Real>>   h_vector_2d_mid = {h_tend_t,
        //                                                   h_tend_qv,
        //                                                   h_tend_u,
        //                                                   h_tend_v,
        //                                                   h_rain_prod,
        //                                                   h_snow_prod};
        // std::vector<uview_2d <Spack>>   d_vector_2d_mid = {tend_t,
        //                                                   tend_qv,
        //                                                   tend_u,
        //                                                   tend_v,
        //                                                   rain_prod,
        //                                                   snow_prod};
        // std::vector<uview_2dh<Real>>   h_vector_2d_int = {h_prec_flux,
        //                                                   h_snow_flux,
        //                                                   h_mass_flux};
        // std::vector< view_2d <Real>>   d_vector_2d_int = {prec_flux,
        //                                                   snow_flux,
        //                                                   mass_flux};
        
        // ekat::device_to_host(h_int_vector_1d, ncol,           d_int_vector_1d);
        // ekat::device_to_host(h_scl_vector_1d, ncol,           d_scl_vector_1d);
        // ekat::device_to_host(h_vector_2d_mid, ncol, nlev_mid, d_vector_2d_mid);
        // ekat::device_to_host(h_vector_2d_int, ncol, nlev_int, d_vector_2d_int);
      // }

    };
    // // -------------------------------------------------------------------------
    // void init(int ncol_in, int nlev_mid_in) {
    //   auto nlev_int_in = nlev_mid_in+1;
    //   auto nlev_int_packs = ekat::npack<Spack>(pverp);
    //   auto nlev_int_packs = ekat::npack<Spack>(pverp);
    //   Kokkos::parallel_for("zm_output_init", KT::RangePolicy(0, m_ncol), KOKKOS_LAMBDA (const int i) {
    //     Real init_fill_value = -999;
    //     // 1D scalar variables
    //     prec(i) = init_fill_value;
    //     snow(i) = init_fill_value;
    //     cape(i) = init_fill_value;
    //     activity(i) = -1;
    //   }
    //   // mid-point level variables
    //   Kokkos::parallel_for("zm_update_precip",KT::RangePolicy(0, m_ncol*nlevm_packs), KOKKOS_LAMBDA (const int i) {
    //     const int icol = i/nlevm_packs;
    //     const int ilev = i%nlevm_packs;
    //     tend_t   (i,k) = init_fill_value;
    //     tend_qv  (i,k) = init_fill_value;
    //     tend_u   (i,k) = init_fill_value;
    //     tend_v   (i,k) = init_fill_value;
    //     rain_prod(i,k) = init_fill_value;
    //     snow_prod(i,k) = init_fill_value;
    //   }
    //   Kokkos::parallel_for("zm_update_precip",KT::RangePolicy(0, m_ncol*pver_in), KOKKOS_LAMBDA (const int i) {
    //     const int icol = i/pver_in;
    //     const int ilev = i%pver_in;
    //     f_tend_t   (i,k) = init_fill_value;
    //     f_tend_qv  (i,k) = init_fill_value;
    //     f_tend_u   (i,k) = init_fill_value;
    //     f_tend_v   (i,k) = init_fill_value;
    //     f_rain_prod(i,k) = init_fill_value;
    //     f_snow_prod(i,k) = init_fill_value;
    //   }
    //   // interface level variables
    //   Kokkos::parallel_for("zm_update_precip",KT::RangePolicy(0, m_ncol*nlevi_packs), KOKKOS_LAMBDA (const int i) {
    //     const int icol = i/nlevi_packs;
    //     const int ilev = i%nlevi_packs;
    //     prec_flux(i,k) = init_fill_value;
    //     snow_flux(i,k) = init_fill_value;
    //     mass_flux(i,k) = init_fill_value;
    //   }
    //   Kokkos::parallel_for("zm_update_precip",KT::RangePolicy(0, m_ncol*pverp), KOKKOS_LAMBDA (const int i) {
    //     const int icol = i/pverp;
    //     const int ilev = i%pverp;
    //     f_prec_flux(i,k) = init_fill_value;
    //     f_snow_flux(i,k) = init_fill_value;
    //     f_mass_flux(i,k) = init_fill_value;
    //   }

      // // mid-point level variables
      // for (int i=0; i<ncol_in; ++i) {
      //   for (int j=0; j<pver_in; ++j) {
      //     tend_t   (i,j/Spack::n)[j%Spack::n] = init_fill_value;
      //     tend_qv  (i,j/Spack::n)[j%Spack::n] = init_fill_value;
      //     tend_u   (i,j/Spack::n)[j%Spack::n] = init_fill_value;
      //     tend_v   (i,j/Spack::n)[j%Spack::n] = init_fill_value;
      //     rain_prod(i,j/Spack::n)[j%Spack::n] = init_fill_value;
      //     snow_prod(i,j/Spack::n)[j%Spack::n] = init_fill_value;
      //     f_tend_t   (i,j) = init_fill_value;
      //     f_tend_qv  (i,j) = init_fill_value;
      //     f_tend_u   (i,j) = init_fill_value;
      //     f_tend_v   (i,j) = init_fill_value;
      //     f_rain_prod(i,j) = init_fill_value;
      //     f_snow_prod(i,j) = init_fill_value;
      //   }
      // }
      // auto pverp = pver_in+1;
      // // interface level variables
      // for (int i=0; i<ncol_in; ++i) {
      //   for (int j=0; j<pverp; ++j) {
      //     prec_flux(i,j/Spack::n)[j%Spack::n] = init_fill_value;
      //     snow_flux(i,j/Spack::n)[j%Spack::n] = init_fill_value;
      //     mass_flux(i,j/Spack::n)[j%Spack::n] = init_fill_value;
      //     f_prec_flux(i,j) = init_fill_value;
      //     f_snow_flux(i,j) = init_fill_value;
      //     f_mass_flux(i,j) = init_fill_value;
      //   }
      // }
    // };
    // // -------------------------------------------------------------------------
    // void init(int ncol_in, int pver_in) {
    //   Real init_fill_value = -999;
    //   // 1D scalar variables
    //   Kokkos::deep_copy(prec, init_fill_value);
    //   Kokkos::deep_copy(snow, init_fill_value);
    //   Kokkos::deep_copy(cape, init_fill_value);
    //   Kokkos::deep_copy(activity, -1);
    //   Kokkos::deep_copy(tend_t   , init_fill_value);
    //   Kokkos::deep_copy(tend_qv  , init_fill_value);
    //   Kokkos::deep_copy(tend_u   , init_fill_value);
    //   Kokkos::deep_copy(tend_v   , init_fill_value);
    //   Kokkos::deep_copy(rain_prod, init_fill_value);
    //   Kokkos::deep_copy(snow_prod, init_fill_value);
    //   // Kokkos::deep_copy(f_tend_t,    init_fill_value);
    //   // Kokkos::deep_copy(f_tend_qv,   init_fill_value);
    //   // Kokkos::deep_copy(f_tend_u,    init_fill_value);
    //   // Kokkos::deep_copy(f_tend_v,    init_fill_value);
    //   // Kokkos::deep_copy(f_rain_prod, init_fill_value);
    //   // Kokkos::deep_copy(f_snow_prod, init_fill_value);
    //   Kokkos::deep_copy(prec_flux, init_fill_value);
    //   Kokkos::deep_copy(snow_flux, init_fill_value);
    //   Kokkos::deep_copy(mass_flux, init_fill_value);
    //   // Kokkos::deep_copy(f_prec_flux, init_fill_value);
    //   // Kokkos::deep_copy(f_snow_flux, init_fill_value);
    //   // Kokkos::deep_copy(f_mass_flux, init_fill_value);
    //   // mid-point level variables
    //   for (int i=0; i<ncol_in; ++i) {
    //     for (int j=0; j<pver_in; ++j) {
    //       f_tend_t   (i,j) = init_fill_value;
    //       f_tend_qv  (i,j) = init_fill_value;
    //       f_tend_u   (i,j) = init_fill_value;
    //       f_tend_v   (i,j) = init_fill_value;
    //       f_rain_prod(i,j) = init_fill_value;
    //       f_snow_prod(i,j) = init_fill_value;
    //     }
    //   }
    //   auto pverp = pver_in+1;
    //   // interface level variables
    //   for (int i=0; i<ncol_in; ++i) {
    //     for (int j=0; j<pverp; ++j) {
    //       f_prec_flux(i,j) = init_fill_value;
    //       f_snow_flux(i,j) = init_fill_value;
    //       f_mass_flux(i,j) = init_fill_value;
    //     }
    //   }
    // };
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
