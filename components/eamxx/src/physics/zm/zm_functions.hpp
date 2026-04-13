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

  using Pack     = ekat::Pack<Scalar,SCREAM_PACK_SIZE>;
  using IntPack = ekat::Pack<Int,SCREAM_PACK_SIZE>;

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

    static inline constexpr Int nit_lheat = 2; // Number of iterations for condensation/freezing loop

    static inline constexpr Real lwmax = 1.e-3; // maximum condensate that can be held in cloud before rainout

    static inline constexpr Real ull_upper_launch_pressure = 600.0; // upper search limit for unrestricted launch level (ULL)

    static inline constexpr Real MCSP_storm_speed_pref = 600e2; // pressure level for winds in MCSP calculation [Pa]

    static inline constexpr Real MCSP_conv_depth_min = 700e2; // pressure thickness of convective heating [Pa]

    static inline constexpr Real MCSP_shear_min = 3.0;   // min shear value for MCSP to be active

    static inline constexpr Real MCSP_shear_max = 200.0; // max shear value for MCSP to be active

    static inline constexpr Real cape_threshold_old = 70.;     // threshold value of cape for deep convection (old value before DCAPE)

    static inline constexpr Real cape_threshold_new = 0.;      // threshold value of cape for deep convection

    static inline constexpr Real dcape_threshold    = 0.;      // threshold value of dcape for deep convection

    static inline constexpr Real interp_diff_min    = 1.E-6;   // minimum threshold for interpolation method - see eq (4.109), (4.118), (4.119)

    static inline constexpr Real omsm               = 0.99999; // to prevent problems due to round off error

    static inline constexpr Real small_conv         = 1.e-20;  // small number to limit blowup when normalizing by mass flux

    static inline constexpr Real beta               = 0;       // proportion of liquid water from layer below used in closure
    static inline constexpr Real mu_min             = 0.02;    // minimum updraft mass flux threshold [mb/s]
    static inline constexpr Real hu_diff_min        = -2000;   // updraft MSE undershoot threshold for cloud top determination [J/kg]
    static inline constexpr Real lambda_limit_min   = 0;       // minimum fractional entrainment limiter [1/m]
    static inline constexpr Real lambda_limit_max   = 0.0002;  // maximum fractional entrainment limiter [1/m]
    static inline constexpr Real lambda_threshold   = 1.e-6;   // threshold for moving detrainment level downward
    static inline constexpr Real pergro_rhd_threshold  = -1.e-4;  // MSE relative difference threshold for perturbation growth test
    static inline constexpr Real pergro_perturbation   = 8.64e-11; // perturbation magnitude added to avoid div-by-zero in pergro test
    static inline constexpr Real momcu             = 0.4;      // pressure gradient term constant for updrafts
    static inline constexpr Real momcd             = 0.4;      // pressure gradient term constant for downdrafts

    // Table of saturation vapor pressure values (estbl) from tmin to
    // tmax+1 Kelvin, in one degree increments.  ttrice defines the
    // transition region, estbl contains a combination of ice & water
    // values.
    static inline constexpr Real tmin = 127.16;
    static inline constexpr Real tmax = 375.16;

    static inline constexpr Real h2otrip = 273.16;

    static inline constexpr Real ttrice = 20.0;  // transition range from es over H2O to es over ice
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
    view_1d<Real> estbl; // table values of saturation vapor pressure
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
    uview_2d<     Pack>  z_mid;    // mid-point level altitude     [m]
    uview_2d<     Pack>  z_del;    // altitude thickness           [m]
    uview_2d<     Pack>  z_int;    // interface level altitude     [m]

    // variables we get from the field manager
    view_1d<const Scalar> phis;     // surface geopotential height  [m2/s]
    view_2d<const Pack>  p_mid;    // mid-point level pressure     [Pa]
    view_2d<const Pack>  p_int;    // interface level pressure     [Pa]
    view_2d<const Pack>  p_del;    // pressure thickness           [Pa]
    view_2d<      Pack>  T_mid;    // temperature                  [K]
    view_2d<      Pack>  qv;       // water vapor mixing ratio     [kg kg-1]
    view_2d<const Pack>  qc;       // cloud liquid water           [kg kg-1]
    view_2d<      Pack>  uwind;    // zonal wind                   [m/s]
    view_2d<      Pack>  vwind;    // meridional wind              [m/s]
    view_2d<const Pack>  omega;    // vertical pressure velocity   [Pa/s]
    view_2d<const Pack>  cldfrac;  // total cloud fraction         [frac]
    view_1d<const Scalar> pblh;     // PBL height                   [m]
    view_1d<const Scalar> landfrac; // land area fraction           [frac]
    view_2d<const Pack>  thl_sec;  // thetal variance from SHOC    [K^2]

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
    uview_2d<Pack>  tend_t;         // output tendency of temperature          [K/s]
    uview_2d<Pack>  tend_qv;        // output tendency of water vapor          [kg/kg/s]
    uview_2d<Pack>  tend_u;         // output tendency of zonal wind           [m/s/s]
    uview_2d<Pack>  tend_v;         // output tendency of meridional wind      [m/s/s]
    uview_2d<Pack>  rain_prod;      // rain production rate                    [?]
    uview_2d<Pack>  snow_prod;      // snow production rate                    [?]
    uview_2d<Pack>  prec_flux;      // output convective precipitation flux    [?]
    uview_2d<Pack>  snow_flux;      // output convective precipitation flux    [?]
    uview_2d<Pack>  mass_flux;      // output convective mass flux             [?]

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

  static void zm_finalize() {
    s_common_init.estbl = view_1d<Real>();
  }

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
  static Real goffgratch_svp_water(
    // Inputs
    const Real& t)    // Temperature                  [K]
  {
    // GoffGratch_svp_water
    static constexpr Real magic1 = -7.90298;
    static constexpr Real magic2 = 5.02808;
    static constexpr Real magic3 = 1.3816e-7;
    static constexpr Real magic4 = 11.344;
    static constexpr Real magic5 = 8.1328e-3;
    static constexpr Real magic6 = -3.49149;
    static constexpr Real magic7 = 1013.246;
    return std::pow(10, magic1*(ZMC::tboil/t-1) +
                    magic2*std::log10(ZMC::tboil/t) -
                    magic3*(std::pow(10, magic4*(1 - t/ZMC::tboil)) - 1) +
                    magic5*(std::pow(10, magic6*(ZMC::tboil/t-1)) - 1) +
                    std::log10(magic7))*100;
  }

  KOKKOS_INLINE_FUNCTION
  static Real goffgratch_svp_ice(
    // Inputs
    const Real& t)    // Temperature                  [K]
  {
    // good down to -100 C
    static constexpr Real magic1 = -9.09718;
    static constexpr Real magic2 = 3.56654;
    static constexpr Real magic3 = 0.876793;
    static constexpr Real magic4 = 6.1071;

    return std::pow(10, magic1*(ZMC::h2otrip/t - 1) -
                    magic2*std::log10(ZMC::h2otrip/t) +
                    magic3*(1 - t/ZMC::h2otrip) +
                    std::log10(magic4))*100;
  }

  KOKKOS_INLINE_FUNCTION
  static Real svp_trans(
    // Inputs
    const Real& t)    // Temperature                  [K]
  {
    Real es;

    //
    // Water
    //
    if (t >= (PC::Tmelt.value - ZMC::ttrice)) {
      es = goffgratch_svp_water(t);
    }
    else {
      es = 0;
    }

    //
    // Ice
    //
    if (t < PC::Tmelt.value) {

      const Real esice = goffgratch_svp_ice(t);

      Real weight;
      if ( (PC::Tmelt.value - t) > ZMC::ttrice ) {
        weight = 1;
      }
      else {
        weight = (PC::Tmelt.value - t)/ZMC::ttrice;
      }

      es = weight*esice + (1 - weight)*es;
    }
    return es;
  }

  KOKKOS_INLINE_FUNCTION
  static void qsat_hPa(
    // Inputs
    const Real& t,    // Temperature                  [K]
    const Real& p,    // Pressure                     [hPa]
    // Outputs
    Real& es,         // Saturation vapor pressure    [hPa]
    Real& qm)         // Saturation mass mixing ratio [kg/kg] (vapor mass over dry mass)
  {
    es = goffgratch_svp_water(t);

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

  KOKKOS_INLINE_FUNCTION
  static void qsat(
    // Inputs
    const Real& t,    // Temperature                  [K]
    const Real& p,    // Pressure                     [hPa]
    const ZmRuntimeOpt& runtime_opt,
    // Outputs
    Real& es,         // Saturation vapor pressure    [hPa]
    Real& qs)         // Saturation mass mixing ratio [kg/kg] (vapor mass over dry mass)
  {
    constexpr Real mmin = 0.0;
    const Real t_tmp = Kokkos::max(Kokkos::min(t,ZMC::tmax)-ZMC::tmin, mmin);     // Number of table entries above tmin
    const Int i = int(t_tmp);                        // Corresponding index.
    const Real weight = t_tmp - Kokkos::trunc(t_tmp);// Fractional part of t_tmp (for interpolation).
    es = (1 - weight)*runtime_opt.estbl(i) + weight*runtime_opt.estbl(i+1);

    // If pressure is less than SVP, set qs to maximum of 1.
    if ( (p - es) < 0 ) {
      qs = 1;
    }
    else {
      qs = PC::ep_2.value*es / (p - ZMC::omeps*es);
    }

    // Ensures returned es is consistent with limiters on qs.
    es = Kokkos::min(es, p);
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
    const ZmRuntimeOpt& runtime_opt,
    const Int& pver, // number of mid-point vertical levels
    const Int& pverp, // number of interface vertical levels
    const Int& num_cin, // num of negative buoyancy regions that are allowed before the conv. top and CAPE calc are completed
    const Int& num_msg, // index of highest level convection is allowed
    const uview_1d<const Real>& sp_humidity_in, // specific humidity [kg/kg]
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
    const ZmRuntimeOpt& runtime_opt,
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
    const Workspace& workspace,
    const ZmRuntimeOpt& runtime_opt,
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

  KOKKOS_FUNCTION
  static void zm_conv_mcsp_calculate_shear(
    // Inputs
    const MemberType& team,
    const Int& pver, // number of mid-point vertical levels
    const uview_1d<const Real>& state_pmid, // physics state mid-point pressure
    const uview_1d<const Real>& state_u, // physics state u momentum
    // Outputs
    Real& mcsp_shear);

  KOKKOS_FUNCTION
  static void zm_conv_mcsp_tend(
    // Inputs
    const MemberType& team,
    const Workspace& workspace,
    const ZmRuntimeOpt& runtime_opt,
    const Int& pver, // number of mid-point vertical levels
    const Int& pverp, // number of interface vertical levels
    const Real& ztodt, // 2x physics time step
    const Int& jctop, // cloud top level indices
    const uview_1d<const Real>& state_pmid, // physics state mid-point pressure
    const uview_1d<const Real>& state_pint, // physics state interface pressure
    const uview_1d<const Real>& state_pdel, // physics state pressure thickness
    const uview_1d<const Real>& state_s, // physics state dry static energy
    const uview_1d<const Real>& state_q, // physics state specific humidity
    const uview_1d<const Real>& state_u, // physics state u momentum
    const uview_1d<const Real>& state_v, // physics state v momentum
    const uview_1d<const Real>& ptend_zm_s, // input ZM tendency for dry static energy (DSE)
    const uview_1d<const Real>& ptend_zm_q, // input ZM tendency for specific humidity (qv)
    // Inputs/Outputs
    const uview_1d<Real>& ptend_s, // output tendency of DSE
    const uview_1d<Real>& ptend_q, // output tendency of qv
    const uview_1d<Real>& ptend_u, // output tendency of u-wind
    const uview_1d<Real>& ptend_v, // output tendency of v-wind
    // Outputs
    const uview_1d<Real>& mcsp_dt_out, // final MCSP tendency for DSE
    const uview_1d<Real>& mcsp_dq_out, // final MCSP tendency for qv
    const uview_1d<Real>& mcsp_du_out, // final MCSP tendency for u wind
    const uview_1d<Real>& mcsp_dv_out, // final MCSP tendency for v wind
    Real& mcsp_freq, // MSCP frequency for output
    Real& mcsp_shear, // shear used to check against threshold
    Real& zm_depth); // pressure depth of ZM heating


  KOKKOS_INLINE_FUNCTION
  static bool is_conv_active(
    // Inputs
    const ZmRuntimeOpt& runtime_opt,
    const bool& is_first_step,       // flag for first step of run
    const Real& cape,                // conv. avail. potential energy     [J]
    const Real& dcape,               // CAPE generated by dycor (dCAPE)   [J]
    // Outputs
    Real& cape_threshold_loc)        // cape threshold
  {
    // set local threshold to be used for zm_closure()
    if ( runtime_opt.trig_dcape && !is_first_step ) {
      cape_threshold_loc = ZMC::cape_threshold_new;
    }
    else {
      cape_threshold_loc = ZMC::cape_threshold_old;
    }

    //--------------------------------------------------------------------------
    // determine if column is active
    if ( runtime_opt.trig_dcape && ! is_first_step ) {
      if ( cape > cape_threshold_loc && dcape > ZMC::dcape_threshold ) {
        return true;
      }
    }
    else {
      if (cape > cape_threshold_loc) {
        return true;
      }
    }

    return false;
  }

  KOKKOS_FUNCTION
  static bool zm_conv_main(
    // Inputs
    const MemberType& team,
    const Workspace& workspace,
    const ZmRuntimeOpt& runtime_opt,
    const Int& pver, // number of mid-point levels
    const Int& pverp, // number of interface levels
    const bool& is_first_step, // flag for first step of run
    const Real& time_step, // model time-step                         [s]
    const uview_1d<const Real>& t_mid, // temperature                             [K]
    const uview_1d<const Real>& q_mid_in, // specific humidity                       [kg/kg]
    const uview_1d<const Real>& omega, // vertical pressure velocity              [Pa/s]
    const uview_1d<const Real>& p_mid_in, // mid-point pressure                      [Pa]
    const uview_1d<const Real>& p_int_in, // interface pressure                      [Pa]
    const uview_1d<const Real>& p_del_in, // pressure thickness                      [Pa]
    const Real& geos, // surface geopotential                    [m2/s2]
    const uview_1d<const Real>& z_mid_in, // mid-point geopotential                  [m2/s2]
    const uview_1d<const Real>& z_int_in, // interface geopotential                  [m2/s2]
    const Real& pbl_hgt, // boundary layer height                   [m]
    const Real& tpert, // parcel temperature perturbation         [K]
    const Real& landfrac, // land fraction                           []
    const uview_1d<const Real>& t_star, // for DCAPE - prev temperature            [K]
    const uview_1d<const Real>& q_star, // for DCAPE - prev sp. humidity           [kg/kg]
    // Outputs
    Int& msemax_klev, // level indices of max MSE
    Int& jctop, // top-of-deep-convection indices
    Int& jcbot, // base of cloud indices
    Int& jt, // gathered top level index of convection
    Real& prec, // output precipitation                    [m/s]
    const uview_1d<Real>& heat, // dry static energy tendency              [W/kg]
    const uview_1d<Real>& qtnd, // specific humidity tendency              [kg/kg/s]
    Real& cape, // conv. avail. potential energy           [J]
    Real& dcape, // CAPE generated by dycor (dCAPE)         [J]
    const uview_1d<Real>& mcon, // convective mass flux                    [mb/s]
    const uview_1d<Real>& pflx, // precip flux at each level               [kg/m2/s]
    const uview_1d<Real>& zdu, // detraining mass flux                    [1/s]
    const uview_1d<Real>& mflx_up, // updraft mass flux                       [mb/s]
    const uview_1d<Real>& entr_up, // updraft entrainment                     [1/s]
    const uview_1d<Real>& detr_up, // updraft detrainment                     [1/s]
    const uview_1d<Real>& mflx_dn, // downdraft mass flux                     [mb/s]
    const uview_1d<Real>& entr_dn, // downdraft entrainment                   [1/s]
    const uview_1d<Real>& p_del, // layer thickness                         [mb]
    Real& dsubcld, // thickness between lcl and msemax_klev   [mb]
    const uview_1d<Real>& ql, // cloud liquid water for chem/wetdep      [?]
    Real& rliq, // reserved liquid (not yet in cldliq) for energy integrals
    const uview_1d<Real>& rprd, // rain production rate                    [kg/kg/s]
    const uview_1d<Real>& dlf); // detrainment rate of cloud liquid water  [kg/kg/s]);

  KOKKOS_FUNCTION
  static void zm_conv_evap(
    // Inputs
    const MemberType& team,
    const ZmRuntimeOpt& runtime_opt,
    const Int& pver, // number of mid-point vertical levels
    const Int& pverp, // number of interface vertical levels
    const Real& time_step, // model time step                         [s]
    const uview_1d<const Real>& p_mid, // midpoint pressure                       [Pa]
    const uview_1d<const Real>& p_del, // layer thickness                         [Pa]
    const uview_1d<const Real>& t_mid, // temperature                             [K]
    const uview_1d<const Real>& q_mid, // water vapor                             [kg/kg]
    const uview_1d<const Real>& prdprec, // precipitation production                [kg/kg/s]
    const uview_1d<const Real>& cldfrc, // cloud fraction
    // Inputs/Outputs
    const uview_1d<Real>& tend_s, // heating rate                            [J/kg/s]
    const uview_1d<Real>& tend_q, // water vapor tendency                    [kg/kg/s]
    // Outputs
    const uview_1d<Real>& tend_s_snwprd, // Heating rate of snow production         [J/kg/s]
    const uview_1d<Real>& tend_s_snwevmlt, // Heating rate of snow evap/melt          [J/kg/s]
    // Inputs/Outputs
    Real& prec, // Convective-scale prec rate              [m/s]
    // Outputs
    Real& snow, // Convective-scale snow rate              [m/s]
    const uview_1d<Real>& ntprprd, // net precip production in layer          [kg/kg/s]
    const uview_1d<Real>& ntsnprd, // net snow production in layer            [kg/kg/s]
    const uview_1d<Real>& flxprec, // Convective flux of prec at interfaces   [kg/m2/s]
    const uview_1d<Real>& flxsnow); // Convective flux of snow at interfaces   [kg/m2/s]);

  KOKKOS_FUNCTION
  static void zm_calc_fractional_entrainment(
    // Inputs
    const MemberType& team,
    const Workspace& workspace,
    const Int& pver, // number of mid-point vertical levels
    const Int& pverp, // number of interface vertical levels
    const Int& msg, // number of levels to ignore at model top
    const Int& jb, // updraft base level
    const Int& jt, // updraft top level
    // Inputs/Outputs
    Int& j0, // level where updraft begins detraining
    // Inputs
    const uview_1d<const Real>& z_mid, // env altitude at mid-point
    const uview_1d<const Real>& z_int, // env altitude at interface
    const uview_1d<const Real>& dz, // layer thickness
    const uview_1d<const Real>& h_env, // env moist stat energy
    const uview_1d<const Real>& h_env_sat, // env saturated moist stat energy
    // Inputs/Outputs
    Real& h_env_min, // mid-tropospheric MSE minimum
    // Outputs
    const uview_1d<Real>& lambda, // fractional entrainment
    Real& lambda_max); // fractional entrainment maximum

  KOKKOS_FUNCTION
  static void zm_downdraft_properties(
    // Inputs
    const MemberType& team,
    const ZmRuntimeOpt& runtime_opt,
    const Int& pver, // number of mid-point vertical levels
    const Int& pverp, // number of interface vertical levels
    const Int& msg, // number of levels to ignore at model top
    const Int& jb, // updraft base level
    // Inputs/Outputs
    Int& jt, // updraft top level
    // Inputs
    const Int& j0, // level where updraft begins detraining
    // Inputs/Outputs
    Int& jd, // level of downdraft
    // Inputs
    const uview_1d<const Real>& z_int, // env altitude at interface
    const uview_1d<const Real>& dz, // layer thickness
    const uview_1d<const Real>& s_mid, // env dry static energy of env [K] (normalized)
    const uview_1d<const Real>& q_mid, // env specific humidity
    const uview_1d<const Real>& h_env, // ambient env moist stat energy
    const uview_1d<const Real>& lambda, // fractional entrainment
    const Real& lambda_max, // fractional entrainment max
    const uview_1d<const Real>& qsthat, // interface interpolated qst
    const uview_1d<const Real>& hsthat, // interface interpolated hst
    const uview_1d<const Real>& gamhat, // interface interpolated gamma
    const uview_1d<const Real>& rprd, // rate of production of precip at that layer
    const uview_1d<const Real>& mflx_up, // updraft mass flux
    // Inputs/Outputs
    const uview_1d<Real>& mflx_dn, // downdraft mass flux
    const uview_1d<Real>& entr_dn, // downdraft entrainment rate
    const uview_1d<Real>& s_dnd, // dndraft dry static energy [K] (normalized)
    const uview_1d<Real>& q_dnd, // dndraft specific humidity [kg/kg]
    const uview_1d<Real>& h_dnd, // dndraft moist static energy
    const uview_1d<Real>& q_dnd_sat, // dndraft saturation specific humdity
    const uview_1d<Real>& evp, // evaporation rate
    Real& totevp); // total evap   for dndraft proportionality factor - see eq (4.106)

  KOKKOS_FUNCTION
  static void zm_cloud_properties(
    // Inputs
    const MemberType& team,
    const Workspace& workspace,
    const ZmRuntimeOpt& runtime_opt,
    const Int& pver, // number of mid-point vertical levels
    const Int& pverp, // number of interface vertical levels
    const Int& msg, // number of levels to ignore at model top
    const Int& limcnv, // convection limiting level
    const uview_1d<const Real>& p_mid, // env pressure at mid-point
    const uview_1d<const Real>& z_mid, // env altitude at mid-point
    const uview_1d<const Real>& z_int, // env altitude at interface
    const uview_1d<const Real>& t_mid, // env temperature
    const uview_1d<const Real>& s_mid, // env dry static energy of env [K] (normalized)
    const uview_1d<const Real>& s_int, // interface values of dry stat energy
    const uview_1d<const Real>& q_mid, // env specific humidity
    const Real& landfrac, // Land fraction
    const Real& tpert_g, // PBL temperature perturbation
    const Int& jb, // updraft base level
    const Int& lel, // updraft parcel launch level
    // Outputs
    Int& jt, // updraft plume top
    Int& jlcl, // updraft lifting cond level
    Int& j0, // level where detrainment begins (starting at h_env_min)
    Int& jd, // level of downdraft
    const uview_1d<Real>& mflx_up, // updraft mass flux
    const uview_1d<Real>& entr_up, // entrainment rate of updraft
    const uview_1d<Real>& detr_up, // detrainement rate of updraft
    const uview_1d<Real>& mflx_dn, // downdraft mass flux
    const uview_1d<Real>& entr_dn, // downdraft entrainment rate
    const uview_1d<Real>& mflx_net, // net mass flux
    const uview_1d<Real>& s_upd, // updraft dry static energy [K] (normalized)
    const uview_1d<Real>& q_upd, // updraft specific humidity [kg/kg]
    const uview_1d<Real>& ql, // updraft liq water
    const uview_1d<Real>& s_dnd, // dndraft dry static energy [K] (normalized)
    const uview_1d<Real>& q_dnd, // dndraft specific humidity [kg/kg]
    const uview_1d<Real>& qst, // env saturation mixing ratio
    const uview_1d<Real>& cu, // condensation rate
    const uview_1d<Real>& evp, // evaporation rate
    const uview_1d<Real>& pflx, // precipitation flux thru layer
    const uview_1d<Real>& rprd); // rate of production of precip at that layer

  KOKKOS_FUNCTION
  static void zm_closure(
    // Inputs
    const MemberType& team,
    const Workspace& workspace,
    const ZmRuntimeOpt& runtime_opt,
    const Int& pver, // number of mid-point vertical levels
    const Int& pverp, // number of interface vertical levels
    const Int& msg, // number of levels to ignore at model top
    const Real& cape_threshold_in, // CAPE threshold for "cloud work function" (i.e. A)
    const Int& lcl, // index of lcl
    const Int& lel, // index of launch leve
    const Int& jt, // top of updraft
    const Int& mx, // base of updraft
    const Real& dsubcld, // thickness of subcloud layer
    const uview_1d<const Real>& z_mid, // altitude (m)
    const uview_1d<const Real>& z_int, // height of interface levels
    const uview_1d<const Real>& p_mid, // ambient pressure (mb)
    const uview_1d<const Real>& p_del, // pressure thickness of layers
    const uview_1d<const Real>& t_mid, // ambient temperature
    const uview_1d<const Real>& s_mid, // ambient dry static energy (normalized)
    const uview_1d<const Real>& q_mid, // ambient specific humidity
    const uview_1d<const Real>& qs, // ambient saturation specific humidity
    const uview_1d<const Real>& ql, // ambient liquid water mixing ratio
    const uview_1d<const Real>& s_int, // env. normalized dry static energy at intrfcs
    const uview_1d<const Real>& q_int, // environment specific humidity at interfaces
    const Real& t_pcl_lcl, // parcel temperature at LCL
    const uview_1d<const Real>& t_pcl, // parcel temperature
    const uview_1d<const Real>& q_pcl_sat, // parcel specific humidity
    const uview_1d<const Real>& s_upd, // updraft dry static energy (normalized)
    const uview_1d<const Real>& q_upd, // updraft specific humidity
    const uview_1d<const Real>& mflx_net, // net convective mass flux
    const uview_1d<const Real>& detr_up, // detrainment from updraft
    const uview_1d<const Real>& mflx_up, // updraft mass flux
    const uview_1d<const Real>& mflx_dn, // dndraft mass flux
    const uview_1d<const Real>& q_dnd, // dndraft specific humidity
    const uview_1d<const Real>& s_dnd, // dndraft dry static energy
    const Real& cape, // convective available potential energy
    // Outputs
    Real& cld_base_mass_flux); // cloud base mass flux

  KOKKOS_FUNCTION
  static void zm_calc_output_tend(
    // Inputs
    const MemberType& team,
    const Int& pver, // number of mid-point vertical levels
    const Int& pverp, // number of interface vertical levels
    const Int& msg, // number of levels to ignore at model top
    const Int& jt, // level index of updraft top
    const Int& mx, // level index of updraft base
    const Int& ktm, // min jt over all columns
    const Int& kbm, // min mx over all columns
    const Real& dsubcld, // sub-cloud layer thickness
    const uview_1d<const Real>& p_del, // pressure thickness
    const uview_1d<const Real>& s_int, // ambient interface dry static energy
    const uview_1d<const Real>& q_int, // ambient interface specific humidity
    const uview_1d<const Real>& s_upd, // updraft dry static energy
    const uview_1d<const Real>& q_upd, // updraft specific humidity
    const uview_1d<const Real>& mflx_up, // updraft mass flux
    const uview_1d<const Real>& detr_up, // updraft detrainment
    const uview_1d<const Real>& mflx_dn, // downdraft mass flux
    const uview_1d<const Real>& s_dnd, // downdraft dry static energy
    const uview_1d<const Real>& q_dnd, // downdraft specific humidity
    const uview_1d<const Real>& ql, // cloud liquid water
    const uview_1d<const Real>& evp, // evaporation
    const uview_1d<const Real>& cu, // updraft condensation
    // Outputs
    const uview_1d<Real>& dsdt, // output tendency for dry static energy
    const uview_1d<Real>& dqdt, // output tendency for specific humidity
    const uview_1d<Real>& dl); // output tendency for cloud liquid water

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
# include "impl/zm_zm_conv_mcsp_calculate_shear_impl.hpp"
# include "impl/zm_zm_conv_mcsp_tend_impl.hpp"
# include "impl/zm_zm_conv_main_impl.hpp"
# include "impl/zm_zm_conv_evap_impl.hpp"
# include "impl/zm_zm_calc_fractional_entrainment_impl.hpp"
# include "impl/zm_zm_downdraft_properties_impl.hpp"
# include "impl/zm_zm_cloud_properties_impl.hpp"
# include "impl/zm_zm_closure_impl.hpp"
# include "impl/zm_zm_calc_output_tend_impl.hpp"
#endif // GPU && !KOKKOS_ENABLE_*_RELOCATABLE_DEVICE_CODE
#endif // ZM_FUNCTIONS_HPP
