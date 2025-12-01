#ifndef ZM_FUNCTIONS_HPP
#define ZM_FUNCTIONS_HPP

#include "share/physics/physics_constants.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"
#include "share/core/eamxx_types.hpp"

#include <ekat_pack_kokkos.hpp>
#include <ekat_workspace.hpp>
#include <ekat_reduction_utils.hpp>

namespace scream {
namespace zm {

/*
 * Functions is a stateless struct used to encapsulate a number of functions for ZM.
 */

template <typename ScalarT, typename DeviceT>
struct Functions {
  // -----------------------------------------------------------------------------------------------
  // Types

  using Scalar = ScalarT;
  using Device = DeviceT;

  template <typename S> using BigPack = ekat::Pack<S,SCREAM_PACK_SIZE>;
  template <typename S> using SmallPack = ekat::Pack<S,SCREAM_SMALL_PACK_SIZE>;

  using SPackInt = SmallPack<Int>;
  using BPack    = BigPack<Scalar>;
  using Spack    = SmallPack<Scalar>;

  using PF  = scream::PhysicsFunctions<DefaultDevice>;
  using PC  = scream::physics::Constants<Real>;

  using KT  = ekat::KokkosTypes<DefaultDevice>;

  template <typename S> using view_1d           = typename KT::template view_1d<S>;
  template <typename S> using view_2d           = typename KT::template view_2d<S>;
  template <typename S> using view_2dl          = typename KT::template lview<S**>;
  template <typename S> using view_3d           = typename KT::template view_3d<S>;
  template <typename S> using view_2d_strided   = typename KT::template sview<S**>;
  template <typename S> using view_3d_strided   = typename KT::template sview<S***>;
  template <typename S> using uview_1d          = typename ekat::template Unmanaged<view_1d<S> >;
  template <typename S> using uview_2d          = typename ekat::template Unmanaged<view_2d<S> >;
  template <typename S> using uview_2dl         = typename ekat::template Unmanaged<view_2dl<S> >;
  template <typename S> using view_2dh          = typename view_2dl<S>::HostMirror;
  template <typename S> using view_1dh          = typename view_1d<S>::HostMirror;

  // -----------------------------------------------------------------------------------------------
  // Structs

  struct zm_runtime_opt {
    zm_runtime_opt() = default;
    bool apply_tendencies = false;

    void load_runtime_options(ekat::ParameterList& params) {
      apply_tendencies = params.get<bool>("apply_tendencies", apply_tendencies);
    }
  };

  // -----------------------------------------------------------------------------------------------

  struct zm_input_state {
    zm_input_state() = default;
    // -------------------------------------------------------------------------
    Real dtime;                     // model phsyics time step [s]
    bool is_first_step;             // flag for first call

    // variable counters for device-side only
    static constexpr int num_1d_intgr = 0;  // number of 1D integer views
    static constexpr int num_1d_scalr = 1;  // number of 1D scalar views
    static constexpr int num_2d_midlv = 2;  // number of 2D mid-point views
    static constexpr int num_2d_intfc = 1;  // number of 2D interface views

    uview_1d<     Scalar> tpert;    // PBL top temperature perturb. [K]
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
    view_2d<const Spack>  qc;       // cloud liquid water           [kg kg-1]
    view_2d<      Spack>  uwind;    // zonal wind                   [m/s]
    view_2d<      Spack>  vwind;    // meridional wind              [m/s]
    view_2d<const Spack>  omega;    // vertical pressure velocity   [Pa/s]
    view_2d<const Spack>  cldfrac;  // total cloud fraction         [frac]
    view_1d<const Scalar> pblh;     // PBL height                   [m]
    view_1d<const Scalar> landfrac; // land area fraction           [frac]
    view_2d<const Spack>  thl_sec;  // thetal variance from SHOC    [K^2]

    // *************************************************************************
    // TEMPORARY
    // *************************************************************************
    // LayoutLeft views for fortran bridging
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
    // *************************************************************************
    // TEMPORARY
    // *************************************************************************

    // host mirror versions of ZM interface variables
    view_1dh<Scalar> h_phis;
    view_1dh<Scalar> h_pblh;
    view_1dh<Scalar> h_tpert;
    view_1dh<Scalar> h_landfrac;
    view_2dh<Real>   h_z_mid;
    view_2dh<Real>   h_p_mid;
    view_2dh<Real>   h_p_del;
    view_2dh<Real>   h_T_mid;
    view_2dh<Real>   h_qv;
    view_2dh<Real>   h_uwind;
    view_2dh<Real>   h_vwind;
    view_2dh<Real>   h_omega;
    view_2dh<Real>   h_cldfrac;
    view_2dh<Real>   h_z_int;
    view_2dh<Real>   h_p_int;

    // -------------------------------------------------------------------------
    // transpose method for fortran bridging
    template <ekat::TransposeDirection::Enum D>
    void transpose(int ncol, int nlev_mid) {
      auto nlev_int = nlev_mid+1;

      // ***********************************************************************
      // TEMPORARY
      // ***********************************************************************
      auto nlev_mid_packs = ekat::npack<Spack>(nlev_mid);
      auto nlev_int_packs = ekat::npack<Spack>(nlev_int);
      if (D == ekat::TransposeDirection::c2f) {
        //----------------------------------------------------------------------
        // mid-point level variables
        Kokkos::parallel_for("zm_output_tx_mid",KT::RangePolicy(0, ncol*nlev_mid_packs), KOKKOS_LAMBDA (const int i) {
          const int icol = i/nlev_mid_packs;
          const int klev = i%nlev_mid_packs;
          f_z_mid   (icol,klev) = z_mid   (icol,klev/Spack::n)[klev%Spack::n];
          f_p_mid   (icol,klev) = p_mid   (icol,klev/Spack::n)[klev%Spack::n];
          f_p_del   (icol,klev) = p_del   (icol,klev/Spack::n)[klev%Spack::n];
          f_T_mid   (icol,klev) = T_mid   (icol,klev/Spack::n)[klev%Spack::n];
          f_qv      (icol,klev) = qv      (icol,klev/Spack::n)[klev%Spack::n];
          f_uwind   (icol,klev) = uwind   (icol,klev/Spack::n)[klev%Spack::n];
          f_vwind   (icol,klev) = vwind   (icol,klev/Spack::n)[klev%Spack::n];
          f_omega   (icol,klev) = omega   (icol,klev/Spack::n)[klev%Spack::n];
          f_cldfrac (icol,klev) = cldfrac (icol,klev/Spack::n)[klev%Spack::n];
        });
        // interface level variables
        Kokkos::parallel_for("zm_output_tx_mid",KT::RangePolicy(0, ncol*nlev_int_packs), KOKKOS_LAMBDA (const int i) {
          const int icol = i/nlev_int_packs;
          const int klev = i%nlev_int_packs;
          f_z_int   (icol,klev) = z_int   (icol,klev/Spack::n)[klev%Spack::n];
          f_p_int   (icol,klev) = p_int   (icol,klev/Spack::n)[klev%Spack::n];
        });
        //----------------------------------------------------------------------
        // copy to host mirrors
        Kokkos::deep_copy(h_phis,     phis );
        Kokkos::deep_copy(h_pblh,     pblh );
        Kokkos::deep_copy(h_tpert,    tpert );
        Kokkos::deep_copy(h_landfrac, landfrac );
        Kokkos::deep_copy(h_z_mid,    f_z_mid );
        Kokkos::deep_copy(h_p_mid,    f_p_mid );
        Kokkos::deep_copy(h_p_del,    f_p_del );
        Kokkos::deep_copy(h_T_mid,    f_T_mid );
        Kokkos::deep_copy(h_qv,       f_qv );
        Kokkos::deep_copy(h_uwind,    f_uwind );
        Kokkos::deep_copy(h_vwind,    f_vwind );
        Kokkos::deep_copy(h_omega,    f_omega );
        Kokkos::deep_copy(h_cldfrac,  f_cldfrac );
        Kokkos::deep_copy(h_z_int,    f_z_int );
        Kokkos::deep_copy(h_p_int,    f_p_int );
      }
      // ***********************************************************************
      // TEMPORARY
      // ***********************************************************************

      // if (D == ekat::TransposeDirection::c2f) {
      //   ekat::device_to_host({h_phis.data()},     ncol,           std::vector< view_1d<const Scalar>>{phis});
      //   ekat::device_to_host({h_pblh.data()},     ncol,           std::vector< view_1d<const Scalar>>{pblh});
      //   ekat::device_to_host({h_tpert.data()},    ncol,           std::vector<uview_1d<      Scalar>>{tpert});
      //   ekat::device_to_host({h_landfrac.data()}, ncol,           std::vector< view_1d<const Scalar>>{landfrac});
      //   ekat::device_to_host({h_z_mid.data()},    ncol, nlev_mid, std::vector<uview_2d<const Spack >>{z_mid});
      //   ekat::device_to_host({h_p_mid.data()},    ncol, nlev_mid, std::vector< view_2d<const Spack >>{p_mid});
      //   ekat::device_to_host({h_p_del.data()},    ncol, nlev_mid, std::vector< view_2d<const Spack >>{p_del});
      //   ekat::device_to_host({h_T_mid.data()},    ncol, nlev_mid, std::vector< view_2d<      Spack >>{T_mid});
      //   ekat::device_to_host({h_qv.data()},       ncol, nlev_mid, std::vector< view_2d<      Spack >>{qv});
      //   ekat::device_to_host({h_uwind.data()},    ncol, nlev_mid, std::vector< view_2d<      Spack >>{uwind});
      //   ekat::device_to_host({h_vwind.data()},    ncol, nlev_mid, std::vector< view_2d<      Spack >>{vwind});
      //   ekat::device_to_host({h_omega.data()},    ncol, nlev_mid, std::vector< view_2d<const Spack >>{omega});
      //   ekat::device_to_host({h_cldfrac.data()},  ncol, nlev_mid, std::vector< view_2d<const Spack >>{cldfrac});
      //   ekat::device_to_host({h_z_int.data()},    ncol, nlev_int, std::vector<uview_2d<      Spack >>{z_int});
      //   ekat::device_to_host({h_p_int.data()},    ncol, nlev_int, std::vector< view_2d<const Spack >>{p_int});
      // }
    }
    // -------------------------------------------------------------------------
    void calculate_tpert(int ncol,int nlev,bool is_first_step) {
      const Real cpair  = PC::Cpair;
      const Real latvap = PC::LatVap;

      // create temporaries to avoid "Implicit capture" warning
      auto loc_tpert    = tpert;
      auto loc_pblh     = pblh;
      auto loc_z_int    = z_int;
      auto loc_p_mid    = p_mid;
      auto loc_qc       = qc;
      auto loc_thl_sec  = thl_sec;

      Kokkos::parallel_for("zm_calculate_tpert",ncol, KOKKOS_LAMBDA (const int i) {
        if (is_first_step) {
          loc_tpert(i) = 0.0;
        } else {
          // identify interface index for top of PBL
          int pblh_k_ind = -1;
          for (int k=0; k<nlev; ++k) {
            auto z_int_tmp_k   = loc_z_int(i,k/Spack::n)[k%Spack::n];
            auto z_int_tmp_kp1 = loc_z_int(i,k/Spack::n)[k%Spack::n];
            if ( z_int_tmp_k>loc_pblh(i) && z_int_tmp_kp1<=loc_pblh(i) ) {
              pblh_k_ind = k;
            }
          }
          if (pblh_k_ind==-1) {
            // PBL top index not found, so just set the perturbation to zero
            loc_tpert(i) = 0.0;
          } else {
            // calculate tpert as std deviation of temperature from SHOC's theta-l variance
            auto exner_pbl    = PF::exner_function( loc_p_mid(i,pblh_k_ind/Spack::n)[pblh_k_ind%Spack::n] );
            auto qc_pbl       = loc_qc(i,pblh_k_ind/Spack::n)[pblh_k_ind%Spack::n];
            auto thl_sec_pbl  = loc_thl_sec(i,pblh_k_ind/Spack::n)[pblh_k_ind%Spack::n];
            auto thl_std_pbl  = sqrt( thl_sec_pbl ); // std deviation of thetal;
            loc_tpert(i) = ( thl_std_pbl + (latvap/cpair)*qc_pbl ) / exner_pbl;
            loc_tpert(i) = ekat::impl::min(2.0,loc_tpert(i)); // apply limiter
          }
        }
      });
    }
  };

  // -----------------------------------------------------------------------------------------------

  struct zm_output_tend {
    zm_output_tend() = default;

    // variable counters for device-side only
    static constexpr int num_1d_intgr = 1;  // number of 1D integer views
    static constexpr int num_1d_scalr = 3;  // number of 1D scalar views
    static constexpr int num_2d_midlv = 6;  // number of 2D mid-point views
    static constexpr int num_2d_intfc = 3;  // number of 2D interface views

    uview_1d<Int>    activity;       // integer deep convection activity flag
    uview_1d<Scalar> prec;           // surface precipitation                   [m/s]
    uview_1d<Scalar> snow;           // surface snow                            [m/s]
    uview_1d<Scalar> cape;           // convective available potential energy   [J]
    uview_2d<Spack>  tend_t;         // output tendency of temperature          [K/s]
    uview_2d<Spack>  tend_qv;        // output tendency of water vapor          [kg/kg/s]
    uview_2d<Spack>  tend_u;         // output tendency of zonal wind           [m/s/s]
    uview_2d<Spack>  tend_v;         // output tendency of meridional wind      [m/s/s]
    uview_2d<Spack>  rain_prod;      // rain production rate                    [?]
    uview_2d<Spack>  snow_prod;      // snow production rate                    [?]
    uview_2d<Spack>  prec_flux;      // output convective precipitation flux    [?]
    uview_2d<Spack>  snow_flux;      // output convective precipitation flux    [?]
    uview_2d<Spack>  mass_flux;      // output convective mass flux             [?]

    // *************************************************************************
    // TEMPORARY
    // *************************************************************************
    // LayoutLeft views for fortran bridging
    uview_2dl<Real>  f_tend_t;
    uview_2dl<Real>  f_tend_qv;
    uview_2dl<Real>  f_tend_u;
    uview_2dl<Real>  f_tend_v;
    uview_2dl<Real>  f_rain_prod;
    uview_2dl<Real>  f_snow_prod;
    uview_2dl<Real>  f_prec_flux;
    uview_2dl<Real>  f_snow_flux;
    uview_2dl<Real>  f_mass_flux;
    // *************************************************************************
    // TEMPORARY
    // *************************************************************************

    // host versions of ZM interface variables
    view_1dh<Int>    h_activity;
    view_1dh<Scalar> h_prec;
    view_1dh<Scalar> h_snow;
    view_1dh<Scalar> h_cape;
    view_2dh<Real>   h_tend_t;
    view_2dh<Real>   h_tend_qv;
    view_2dh<Real>   h_tend_u;
    view_2dh<Real>   h_tend_v;
    view_2dh<Real>   h_rain_prod;
    view_2dh<Real>   h_snow_prod;
    view_2dh<Real>   h_prec_flux;
    view_2dh<Real>   h_snow_flux;
    view_2dh<Real>   h_mass_flux;

    // -------------------------------------------------------------------------
    // transpose method for fortran bridging
    template <ekat::TransposeDirection::Enum D>
    void transpose(int ncol, int nlev_mid) {
      auto nlev_int = nlev_mid+1;

      // ***********************************************************************
      // TEMPORARY
      // ***********************************************************************
      auto nlev_mid_packs = ekat::npack<Spack>(nlev_mid);
      auto nlev_int_packs = ekat::npack<Spack>(nlev_int);
      if (D == ekat::TransposeDirection::f2c) {
        // copy back to device
        Kokkos::deep_copy(f_tend_t,   h_tend_t);
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
        // mid-point level variables
        Kokkos::parallel_for("zm_output_tx_mid",KT::RangePolicy(0, ncol*nlev_mid_packs), KOKKOS_LAMBDA (const int i) {
          const int icol = i/nlev_mid_packs;
          const int klev = i%nlev_mid_packs;
          tend_t   (icol,klev/Spack::n)[klev%Spack::n] = f_tend_t   (icol,klev);
          tend_qv  (icol,klev/Spack::n)[klev%Spack::n] = f_tend_qv  (icol,klev);
          tend_u   (icol,klev/Spack::n)[klev%Spack::n] = f_tend_u   (icol,klev);
          tend_v   (icol,klev/Spack::n)[klev%Spack::n] = f_tend_v   (icol,klev);
          rain_prod(icol,klev/Spack::n)[klev%Spack::n] = f_rain_prod(icol,klev);
          snow_prod(icol,klev/Spack::n)[klev%Spack::n] = f_snow_prod(icol,klev);
        });
        // interface level variables
        Kokkos::parallel_for("zm_output_tx_mid",KT::RangePolicy(0, ncol*nlev_int_packs), KOKKOS_LAMBDA (const int i) {
          const int icol = i/nlev_int_packs;
          const int klev = i%nlev_int_packs;
          prec_flux(icol,klev/Spack::n)[klev%Spack::n] = f_prec_flux(icol,klev);
          snow_flux(icol,klev/Spack::n)[klev%Spack::n] = f_snow_flux(icol,klev);
          mass_flux(icol,klev/Spack::n)[klev%Spack::n] = f_mass_flux(icol,klev);
        });
      }
      // ***********************************************************************
      // TEMPORARY
      // ***********************************************************************

      // if (D == ekat::TransposeDirection::f2c) {
      //   ekat::host_to_device({h_prec.data()},     ncol,           std::vector<uview_1d<Scalar>>{prec});
      //   ekat::host_to_device({h_activity.data()}, ncol,           std::vector<uview_1d<Int>>   {activity});
      //   ekat::host_to_device({h_prec.data()},     ncol,           std::vector<uview_1d<Scalar>>{prec});
      //   ekat::host_to_device({h_snow.data()},     ncol,           std::vector<uview_1d<Scalar>>{snow});
      //   ekat::host_to_device({h_cape.data()},     ncol,           std::vector<uview_1d<Scalar>>{cape});
      //   ekat::host_to_device({h_tend_t.data()},   ncol, nlev_mid, std::vector<uview_2d<Spack>> {tend_t},    true);
      //   ekat::host_to_device({h_tend_qv.data()},  ncol, nlev_mid, std::vector<uview_2d<Spack>> {tend_qv},   true);
      //   ekat::host_to_device({h_tend_u.data()},   ncol, nlev_mid, std::vector<uview_2d<Spack>> {tend_u},    true);
      //   ekat::host_to_device({h_tend_v.data()},   ncol, nlev_mid, std::vector<uview_2d<Spack>> {tend_v},    true);
      //   ekat::host_to_device({h_rain_prod.data()},ncol, nlev_mid, std::vector<uview_2d<Spack>> {rain_prod}, true);
      //   ekat::host_to_device({h_snow_prod.data()},ncol, nlev_mid, std::vector<uview_2d<Spack>> {snow_prod}, true);
      //   ekat::host_to_device({h_prec_flux.data()},ncol, nlev_int, std::vector<uview_2d<Spack>> {prec_flux}, true);
      //   ekat::host_to_device({h_snow_flux.data()},ncol, nlev_int, std::vector<uview_2d<Spack>> {snow_flux}, true);
      //   ekat::host_to_device({h_mass_flux.data()},ncol, nlev_int, std::vector<uview_2d<Spack>> {mass_flux}, true);
      // }
    };

    // -------------------------------------------------------------------------
    void init(int ncol, int nlev_mid) {
      auto nlev_int = nlev_mid+1;
      auto nlev_mid_packs = ekat::npack<Spack>(nlev_mid);
      auto nlev_int_packs = ekat::npack<Spack>(nlev_int);
      Real init_fill_value = 0;
      // create temporaries to avoid "Implicit capture" warning
      auto loc_prec       = prec;
      auto loc_snow       = snow;
      auto loc_cape       = cape;
      auto loc_activity   = activity;
      auto loc_tend_t     = tend_t;
      auto loc_tend_qv    = tend_qv;
      auto loc_tend_u     = tend_u;
      auto loc_tend_v     = tend_v;
      auto loc_rain_prod  = rain_prod;
      auto loc_snow_prod  = snow_prod;
      auto loc_prec_flux  = prec_flux;
      auto loc_snow_flux  = snow_flux;
      auto loc_mass_flux  = mass_flux;
      // 1D scalar variables
      Kokkos::parallel_for("zm_output_init_s", KT::RangePolicy(0, ncol), KOKKOS_LAMBDA (const int i) {
        loc_prec(i)     = init_fill_value;
        loc_snow(i)     = init_fill_value;
        loc_cape(i)     = init_fill_value;
        loc_activity(i) = -1;
      });
      // mid-point level variables
      Kokkos::parallel_for("zm_output_init_m",KT::RangePolicy(0, ncol*nlev_mid_packs), KOKKOS_LAMBDA (const int i) {
        const int icol = i/nlev_mid_packs;
        const int klev = i%nlev_mid_packs;
        loc_tend_t   (icol,klev) = init_fill_value;
        loc_tend_qv  (icol,klev) = init_fill_value;
        loc_tend_u   (icol,klev) = init_fill_value;
        loc_tend_v   (icol,klev) = init_fill_value;
        loc_rain_prod(icol,klev) = init_fill_value;
        loc_snow_prod(icol,klev) = init_fill_value;
      });
      // interface level variables
      Kokkos::parallel_for("zm_output_init_i",KT::RangePolicy(0, ncol*nlev_int_packs), KOKKOS_LAMBDA (const int i) {
        const int icol = i/nlev_int_packs;
        const int klev = i%nlev_int_packs;
        loc_prec_flux(icol,klev) = init_fill_value;
        loc_snow_flux(icol,klev) = init_fill_value;
        loc_mass_flux(icol,klev) = init_fill_value;
      });
    };
    // -------------------------------------------------------------------------
  };

  // -----------------------------------------------------------------------------------------------

  struct zm_output_diag {
    zm_output_diag() = default;
  };

  // -----------------------------------------------------------------------------------------------
  // Functions

  // static Int zm_main()

}; // struct Functions

} // namespace zm
} // namespace scream

#endif // ZM_FUNCTIONS_HPP
