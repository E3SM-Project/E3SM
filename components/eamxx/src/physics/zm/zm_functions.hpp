#ifndef ZM_FUNCTIONS_HPP
#define ZM_FUNCTIONS_HPP

#include "share/physics/physics_constants.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"
#include "share/core/eamxx_types.hpp"

#include <ekat_pack_kokkos.hpp>
#include <ekat_workspace.hpp>
#include <ekat_reduction_utils.hpp>
#include <ekat_math_utils.hpp>

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

  template <typename S> using view_0d           = typename KT::template view<S>;
  template <typename S> using view_1d           = typename KT::template view_1d<S>;
  template <typename S> using view_2d           = typename KT::template view_2d<S>;
  template <typename S> using view_2dl          = typename KT::template lview<S**>;
  template <typename S> using view_3d           = typename KT::template view_3d<S>;
  template <typename S> using view_2d_strided   = typename KT::template sview<S**>;
  template <typename S> using view_3d_strided   = typename KT::template sview<S***>;
  template <typename S> using uview_1d          = typename ekat::template Unmanaged<view_1d<S> >;
  template <typename S> using uview_2d          = typename ekat::template Unmanaged<view_2d<S> >;
  template <typename S> using uview_2dl         = typename ekat::template Unmanaged<view_2dl<S> >;
  template <typename S> using uview_3d          = typename ekat::template Unmanaged<view_3d<S> >;
  template <typename S> using view_2dh          = typename view_2dl<S>::HostMirror;
  template <typename S> using view_1dh          = typename view_1d<S>::HostMirror;

  using MemberType = typename KT::MemberType;

  using WorkspaceManager = typename ekat::WorkspaceManager<Scalar, Device>;
  using Workspace        = typename WorkspaceManager::Workspace;

  // -----------------------------------------------------------------------------------------------
  // Structs

  //
  // ---------- ZM constants ---------
  //
  struct ZMC {
    // This value is slightly high, but it seems to be the value for the
    // steam point of water originally (and most frequently) used in the
    // Goff & Gratch scheme.
    static inline constexpr Real tboil = 373.16;

    static inline constexpr Real omeps = 1 - PC::ep_2.value;

    static inline constexpr Real pref = 1000;

    static inline constexpr Real cpwv = 1.810e3; // specific heat of water vapor (J/K/kg)

    static inline constexpr Int LOOPMAX = 100; // Max number of iteration loops for invert_entropy

    static inline constexpr Real tol_coeff = 0.001; // tolerance coefficient

    static inline constexpr Real tol_eps   = 3.e-8; // small value for tolerance calculation

    static inline constexpr Real half = 0.5; // Useful for bfb with fortran

    static inline constexpr Real small = 1.e-36; // a small number to avoid division by zero

    static inline constexpr Real cdifr_min = 1.e-6; // minimum layer difference for geometric averaging

    static inline constexpr Real maxc_factor = 1.e-12; // Small numerical regularization constant used in the maximum cloud fraction

    static inline constexpr Real flux_factor = 1.e-12; // Small numerical regularization constant used in convective flux related calculations

    static inline constexpr Real mbsth = 1.e-15; // threshold below which we treat the mass fluxes as zero (in mb/s)

    static inline constexpr Real lcl_pressure_threshold = 600.0; // if LCL pressure is lower => no convection and cape is zero
  };

  //----------------------------------------------------------------------------
  // Purpose: derived type to hold ZM tunable parameters
  //----------------------------------------------------------------------------
  struct ZmRuntimeOpt {
    ZmRuntimeOpt() = default;

    void load_runtime_options(ekat::ParameterList& params) {
      apply_tendencies = params.get<bool>("apply_tendencies", apply_tendencies);
    }

    Real tau;           // convective adjustment time scale
    Real alfa;          // max downdraft mass flux fraction
    Real ke;            // evaporation efficiency
    Real dmpdz;         // fractional mass entrainment rate [1/m]
    bool tpert_fix;     // flag to disable using applying tpert to PBL-rooted convection
    Real tpert_fac;     // tunable temperature perturbation factor
    Real tiedke_add;    // tunable temperature perturbation
    Real c0_lnd;        // autoconversion coefficient over land
    Real c0_ocn;        // autoconversion coefficient over ocean
    int num_cin;        // num of neg buoyancy regions allowed before the conv top and CAPE calc are completed
    int limcnv;         // upper pressure interface level to limit deep convection
    int mx_bot_lyr_adj; // bot layer index adjustment for launch level search
    bool trig_dcape;    // true if to using DCAPE trigger - based on CAPE generation by the dycor
    bool trig_ull;      // true if to using the "unrestricted launch level" (ULL) mode
    bool clos_dyn_adj;  // flag for mass flux adjustment to CAPE closure
    bool no_deep_pbl;   // flag to eliminate deep convection within PBL
    bool apply_tendencies;
    // ZM micro parameters
    bool zm_microp;     // switch for convective microphysics
    bool old_snow;      // switch to calculate snow prod in zm_conv_evap() (old treatment before zm_microp was implemented)
    Real auto_fac;      // ZM microphysics enhancement factor for droplet-rain autoconversion
    Real accr_fac;      // ZM microphysics enhancement factor for droplet-rain accretion
    Real micro_dcs;     // ZM microphysics size threshold for cloud ice to snow autoconversion [m]
    // MCSP parameters
    bool mcsp_enabled;  // flag for mesoscale coherent structure parameterization (MSCP)
    Real mcsp_t_coeff;  // MCSP coefficient for temperature tendencies
    Real mcsp_q_coeff;  // MCSP coefficient for specific humidity tendencies
    Real mcsp_u_coeff;  // MCSP coefficient for zonal momentum tendencies
    Real mcsp_v_coeff;  // MCSP coefficient for meridional momentum tendencies
  };

  // -----------------------------------------------------------------------------------------------

  struct ZmInputState {
    ZmInputState() = default;
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
    void transpose(int ncol, int nlev_mid);

    // -------------------------------------------------------------------------
    void calculate_tpert(int ncol,int nlev,bool is_first_step);
  };

  // -----------------------------------------------------------------------------------------------

  struct ZmOutputTend {
    ZmOutputTend() = default;

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
    void transpose(int ncol, int nlev_mid);

    // -------------------------------------------------------------------------
    void init(int ncol, int nlev_mid);
  };

  // -----------------------------------------------------------------------------------------------

  struct ZmOutputDiag {
    ZmOutputDiag() = default;
  };


  //
  // --------- Init/Finalize Functions ---------
  //
  static void zm_common_init();

  static void zm_finalize() {}

  // static Int zm_main()

  //
  // --------- Functions ---------
  //

  KOKKOS_FUNCTION
  static void invert_entropy(
    // Inputs
    const MemberType& team,
    const Real& s,    // entropy                           [J/kg]
    const Real& p,    // pressure                          [mb]
    const Real& qt,   // total water mixing ratio          [kg/kg]
    const Real& tfg,  // input temperature for first guess [K]
    // Outputs
    Real& t,          // temperature                       [k]
    Real& qst);       // saturation vapor mixing ratio     [kg/kg]

  KOKKOS_FUNCTION
  static Real entropy(
    // Inputs
    const Real& tk,    // temperature              [K]
    const Real& p,     // pressure                 [mb]
    const Real& qtot); // total water mixing ratio [kg/kg]

  KOKKOS_INLINE_FUNCTION
  static void qsat_hPa(
    // Inputs
    const Real& t,    // Temperature                  [K]
    const Real& p,    // Pressure                     [hPa]
    Real& es,         // Saturation vapor pressure    [hPa]
    Real& qm)         // Saturation mass mixing ratio [kg/kg] (vapor mass over dry mass)
  {
    // GoffGratch_svp_water
    static constexpr Real magic1 = -7.90298;
    static constexpr Real magic2 = 5.02808;
    static constexpr Real magic3 = 1.3816e-7;
    static constexpr Real magic4 = 11.344;
    static constexpr Real magic5 = 8.1328e-3;
    static constexpr Real magic6 = -3.49149;
    static constexpr Real magic7 = 1013.246;
    es = std::pow(10, magic1*(ZMC::tboil/t-1) +
                  magic2*std::log10(ZMC::tboil/t) -
                  magic3*(std::pow(10, magic4*(1 - t/ZMC::tboil)) - 1) +
                  magic5*(std::pow(10, magic6*(ZMC::tboil/t-1)) - 1) +
                  std::log10(magic7))*100;

    // If pressure is less than SVP, set qs to maximum of 1.
    if ( (p*100 - es) <= 0 ) {
      qm = 1;
    }
    else {
      qm = PC::ep_2.value*es / (p*100 - ZMC::omeps*es);
    }

    // Ensures returned es is consistent with limiters on qs.
    es = ekat::impl::min(es, p*100);

    es = es*0.01;
  }

  KOKKOS_FUNCTION
  static void zm_transport_tracer(
    // Inputs
    const MemberType& team,
    const Workspace& workspace,
    const ZmRuntimeOpt& runtime_opt,
    const Int& pver,                        // number of mid-point levels
    const uview_1d<const bool>& doconvtran, // flag for doing convective transport
    const uview_2d<const Real>& q,          // tracer array (including water vapor)
    const Int& ncnst,                       // number of tracers to transport
    const uview_1d<const Real>& mu,         // mass flux up
    const uview_1d<const Real>& md,         // mass flux down
    const uview_1d<const Real>& du,         // mass detraining from updraft
    const uview_1d<const Real>& eu,         // mass entraining from updraft
    const uview_1d<const Real>& ed,         // mass entraining from downdraft
    const uview_1d<const Real>& dp,         // delta pressure between interfaces
    const Int& jt,                          // index of cloud top for this column
    const Int& mx,                          // index of cloud bottom for this column
    const Int& ktm,                         // Highest top level for any column
    const Int& kbm,                         // Highest bottom level for any column
    const uview_2d<const Real>& fracis,     // fraction of tracer that is insoluble
    const uview_1d<const Real>& dpdry,      // delta pressure between interfaces
    const Real& dt,                         // model time increment)
    // Outputs
    const uview_2d<Real>& dqdt);            // output tendency array

    KOKKOS_FUNCTION
  static void zm_transport_momentum(
    // Inputs
    const MemberType& team,
    const Workspace& workspace,
    const Int& pver, // number of mid-point levels
    const Int& pverp, // number of interface levels
    const uview_2d<const Real>& wind_in, // input Momentum array
    const Int& nwind, // number of tracers to transport
    const uview_1d<const Real>& mu, // mass flux up
    const uview_1d<const Real>& md, // mass flux down
    const uview_1d<const Real>& du, // mass detraining from updraft
    const uview_1d<const Real>& eu, // mass entraining from updraft
    const uview_1d<const Real>& ed, // mass entraining from downdraft
    const uview_1d<const Real>& dp, // gathered pressure delta between interfaces
    const Int& jt, // index of cloud top for each column
    const Int& mx, // index of cloud top for each column
    const Int& ideep, // gathering array
    const Int& il1g, // gathered min ncol index
    const Int& il2g, // gathered max ncol index
    const Real& dt, // time step in seconds : 2*delta_t
    const Int& ktm, // Highest top level for any column
    const Int& kbm, // Highest bottom level for any column
    // Outputs
    const uview_2d<Real>& wind_tend, // output momentum tendency
    const uview_2d<Real>& pguall, // apparent force from  updraft PG
    const uview_2d<Real>& pgdall, // apparent force from  downdraft PG
    const uview_2d<Real>& icwu, // in-cloud winds in updraft
    const uview_2d<Real>& icwd, // in-cloud winds in downdraft
    const uview_1d<Real>& seten); // dry static energy tendency);

  KOKKOS_FUNCTION
  static void compute_dilute_cape(
    // Inputs
    const MemberType& team,
    const Workspace& workspace,
    const Int& pver, // number of mid-point vertical levels
    const Int& pverp, // number of interface vertical levels
    const Int& num_cin, // num of negative buoyancy regions that are allowed before the conv. top and CAPE calc are completed
    const Int& num_msg, // index of highest level convection is allowed
    const uview_1d<const Real>& sp_humidity_in, // specific humidity
    const uview_1d<const Real>& temperature_in, // temperature
    const uview_1d<const Real>& zmid, // altitude/height at mid-levels
    const uview_1d<const Real>& pmid, // pressure at mid-levels
    const uview_1d<const Real>& pint, // pressure at interfaces
    const Int& pblt, // index of pbl top used as upper limit index of max MSE search
    const Real& tpert, // perturbation temperature by pbl processes
    const bool& calc_msemax_klev, // true for normal procedure, otherwise use prev_msemax_klev from 1st call
    const Int& prev_msemax_klev, // values of msemax_klev from previous call for dcape closure
    const bool& use_input_tq_mx, // if .true., use input values of prev_msemax_klev, q_mx, t_mx in the CAPE calculation
    // Inputs/Outputs
    const uview_1d<Real>& parcel_qsat, // parcel saturation mixing ratio
    Int& msemax_klev, // index of max MSE at parcel launch level
    Int& lcl_klev, // index of lifting condensation level (i.e. cloud bottom)
    Int& eql_klev, // index of equilibrium level (i.e. cloud top)
    Real& cape, // convective available potential energy
    Real& q_mx, // specified sp humidity to apply at level of max MSE if use_input_tq_mx=.true.
    Real& t_mx, // specified temperature to apply at level of max MSE if use_input_tq_mx=.true.)
    // Outputs
    const uview_1d<Real>& parcel_temp, // parcel temperature
    Real& lcl_temperature); // lifting condensation level (LCL) temperature

  KOKKOS_FUNCTION
  static void find_mse_max(
    // Inputs
    const MemberType& team,
    const Int& pver, // number of mid-point vertical levels
    const Int& num_msg, // number of missing moisture levels at the top of model
    const Int& msemax_top_k, // upper limit index of max MSE search
    const bool& pergro_active, // flag for perturbation growth test (pergro)
    const uview_1d<const Real>& temperature, // environement temperature
    const uview_1d<const Real>& zmid, // height/altitude at mid-levels
    const uview_1d<const Real>& sp_humidity, // specific humidity
    // Inputs/Outputs
    Int& msemax_klev, // index of max MSE at parcel launch level
    Real& mse_max_val); // value of max MSE at parcel launch level

  KOKKOS_FUNCTION
  static void compute_dilute_parcel(
    // Inputs
    const MemberType& team,
    const Int& pver, // number of mid-point vertical levels
    const Int& num_msg, // number of missing moisture levels at the top of model
    const Int& klaunch, // index of parcel launch level based on max MSE
    const uview_1d<const Real>& pmid, // ambient env pressure at cell center
    const uview_1d<const Real>& temperature, // ambient env temperature at cell center
    const uview_1d<const Real>& sp_humidity, // ambient env specific humidity at cell center
    const Real& tpert, // PBL temperature perturbation
    const Int& pblt, // index of pbl depth
    // Inputs/Outputs
    const uview_1d<Real>& parcel_temp, // Parcel temperature
    const uview_1d<Real>& parcel_vtemp, // Parcel virtual temperature
    const uview_1d<Real>& parcel_qsat, // Parcel water vapour (sat value above lcl)
    Real& lcl_pmid, // lifting condensation level (LCL) pressure
    Real& lcl_temperature, // lifting condensation level (LCL) temperature
    Int& lcl_klev); // lifting condensation level (LCL) vertical index

  KOKKOS_FUNCTION
  static void compute_cape_from_parcel(
    // Inputs
    const MemberType& team,
    const Workspace& workspace,
    const ZmRuntimeOpt& runtime_opt,
    const Int& pver, // number of mid-point vertical levels
    const Int& pverp, // number of interface vertical levels
    const Int& num_cin, // num of negative buoyancy regions that are allowed before the conv. top and CAPE calc are completed
    const Int& num_msg, // number of missing moisture levels at the top of model
    const uview_1d<const Real>& temperature, // temperature
    const uview_1d<const Real>& tv, // virtual temperature
    const uview_1d<const Real>& zmid, // height/altitude at mid-levels
    const uview_1d<const Real>& sp_humidity, // specific humidity
    const uview_1d<const Real>& pint, // pressure at interfaces
    const Int& msemax_klev, // index of max MSE at parcel launch level
    const Real& lcl_pmid, // lifting condensation level (LCL) pressure
    const Int& lcl_klev, // lifting condensation level (LCL) index
    // Inputs/Outputs
    const uview_1d<Real>& parcel_qsat, // parcel saturation mixing ratio
    const uview_1d<Real>& parcel_temp, // parcel temperature
    const uview_1d<Real>& parcel_vtemp, // parcel virtual temperature
    Int& eql_klev, // index of equilibrium level (i.e. cloud top)
    Real& cape); // convective available potential energy

  //
  // --------- Members ---------
  //
  inline static ZmRuntimeOpt s_common_init;
}; // struct Functions

} // namespace zm
} // namespace scream

#if defined(EAMXX_ENABLE_GPU) && !defined(KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE) \
                                && !defined(KOKKOS_ENABLE_HIP_RELOCATABLE_DEVICE_CODE)
# include "impl/zm_input_state_impl.hpp"
# include "impl/zm_output_tend_impl.hpp"
# include "impl/zm_common_init_impl.hpp"
# include "impl/zm_invert_entropy_impl.hpp"
# include "impl/zm_entropy_impl.hpp"
# include "impl/zm_zm_transport_tracer_impl.hpp"
# include "impl/zm_zm_transport_momentum_impl.hpp"
# include "impl/zm_compute_dilute_cape_impl.hpp"
# include "impl/zm_find_mse_max_impl.hpp"
# include "impl/zm_compute_dilute_parcel_impl.hpp"
# include "impl/zm_compute_cape_from_parcel_impl.hpp"
#endif // GPU && !KOKKOS_ENABLE_*_RELOCATABLE_DEVICE_CODE
#endif // ZM_FUNCTIONS_HPP
