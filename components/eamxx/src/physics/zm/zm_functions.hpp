#ifndef ZM_FUNCTIONS_HPP
#define ZM_FUNCTIONS_HPP

#include "share/physics/physics_constants.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"
#include "share/core/eamxx_types.hpp"
#include "share/grid/abstract_grid.hpp"

#include <ekat_pack_kokkos.hpp>
#include <ekat_workspace.hpp>
#include <ekat_reduction_utils.hpp>
#include <ekat_math_utils.hpp>

#include <iostream>
#include <iomanip>

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

  using Pack    = ekat::Pack<Scalar,SCREAM_PACK_SIZE>;
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
  template <typename S> using view_2dh          = typename view_2dl<S>::host_mirror_type;
  template <typename S> using view_1dh          = typename view_1d<S>::host_mirror_type;

  using MemberType = typename KT::MemberType;

  using WorkspaceManager = typename ekat::WorkspaceManager<Scalar, Device>;
  using Workspace        = typename WorkspaceManager::Workspace;

  // -----------------------------------------------------------------------------------------------
  // Structs

  //
  // ---------- ZM constants ---------
  //
  struct ZMC {
    static constexpr int nwind = 2;  // number of wind components in tmp_winds
    // This value is slightly high, but it seems to be the value for the
    // steam point of water originally (and most frequently) used in the
    // Goff & Gratch scheme.
    static inline constexpr Real tboil            = 373.16;
    static inline constexpr Real pref             = 1000;     // reference pressure [hPa]
    static inline constexpr Real cpwv             = 1.810e3;  // specific heat of water vapor (J/K/kg)
    // iteration limits
    static inline constexpr Int LOOPMAX           = 100;      // Max number of iteration loops for invert_entropy
    static inline constexpr Int nit_lheat         = 2;        // Number of iterations for condensation/freezing loop
    // miscellaneous tolerance values
    static inline constexpr Real omeps            = 1 - PC::ep_2.value;
    static inline constexpr Real omsm             = 0.99999;  // to prevent problems due to round off error
    static inline constexpr Real tol_coeff        = 0.001;    // tolerance coefficient
    static inline constexpr Real tol_eps          = 3.e-8;    // small value for tolerance calculation
    static inline constexpr Real half             = 0.5;      // Useful for bfb with fortran
    static inline constexpr Real small            = 1.e-36;   // a small number to avoid division by zero
    static inline constexpr Real cdifr_min        = 1.e-6;    // minimum layer difference for geometric averaging
    static inline constexpr Real maxc_factor      = 1.e-12;   // Small numerical regularization constant used in the maximum cloud fraction
    static inline constexpr Real flux_factor      = 1.e-12;   // Small numerical regularization constant used in convective flux related calculations
    static inline constexpr Real mbsth            = 1.e-15;   // threshold below which we treat the mass fluxes as zero (in mb/s)
    static inline constexpr Real small_conv       = 1.e-20;   // small number to limit blowup when normalizing by mass flux
    static inline constexpr Real interp_diff_min  = 1.E-6;    // minimum threshold for interpolation method - see eq (4.109), (4.118), (4.119)
    // conversion factors
    static inline constexpr Real zvir             = 0.608;    // virtual temperature factor (Rv/Rd - 1)
    static inline constexpr Real lcl_coeff_a      = 2840;     // Bolton (1980) LCL temperature formula coefficient A
    static inline constexpr Real lcl_coeff_b      = 3.5;      // Bolton (1980) LCL temperature formula coefficient B
    static inline constexpr Real lcl_coeff_c      = 4.805;    // Bolton (1980) LCL temperature formula coefficient C
    static inline constexpr Real pa_to_mb         = 0.01;     // Pa to mb conversion factor
    static inline constexpr Real mb_to_pa         = 100;      // mb to Pa conversion factor
    // other ZM parameters and thresholds
    static inline constexpr Real lcl_pressure_threshold    = 600.0;    // if LCL pressure is lower => no convection and cape is zero
    static inline constexpr Real lwmax                     = 1.e-3;    // maximum condensate that can be held in cloud before rainout
    static inline constexpr Real ull_upper_launch_pressure = 600.0;    // upper search limit for unrestricted launch level (ULL)
    static inline constexpr Real cape_threshold_old        = 70.;      // threshold value of cape for deep convection (old value before DCAPE)
    static inline constexpr Real cape_threshold_new        = 0.;       // threshold value of cape for deep convection
    static inline constexpr Real dcape_threshold           = 0.;       // threshold value of dcape for deep convection
    static inline constexpr Real beta                      = 0;        // proportion of liquid water from layer below used in closure
    static inline constexpr Real mu_min                    = 0.02;     // minimum updraft mass flux threshold [mb/s]
    static inline constexpr Real hu_diff_min               = -2000;    // updraft MSE undershoot threshold for cloud top determination [J/kg]
    static inline constexpr Real lambda_limit_min          = 0;        // minimum fractional entrainment limiter [1/m]
    static inline constexpr Real lambda_limit_max          = 0.0002;   // maximum fractional entrainment limiter [1/m]
    static inline constexpr Real lambda_threshold          = 1.e-6;    // threshold for moving detrainment level downward
    static inline constexpr Real pergro_rhd_threshold      = -1.e-4;   // MSE relative difference threshold for perturbation growth test
    static inline constexpr Real pergro_perturbation       = 8.64e-11; // perturbation magnitude added to avoid div-by-zero in pergro test
    static inline constexpr Real momcu                     = 0.4;      // pressure gradient term constant for updrafts
    static inline constexpr Real momcd                     = 0.4;      // pressure gradient term constant for downdrafts
    static inline constexpr Real mse_min_diff              = 100;      // min MSE buoyancy difference for Taylor series in entrainment [J/kg]
    static inline constexpr Real tpert_limiter             = 2;        // upper limit on temperature perturbation in input state [K]
    // MCSP parameters
    static inline constexpr Real MCSP_storm_speed_pref     = 600e2;    // pressure level for winds in MCSP calculation Pa]
    static inline constexpr Real MCSP_conv_depth_min       = 700e2;    // pressure thickness of convective heating [Pa]
    static inline constexpr Real MCSP_shear_min            = 3.0;      // min shear value for MCSP to be active
    static inline constexpr Real MCSP_shear_max            = 200.0;    // max shear value for MCSP to be active
    static inline constexpr Real MCSP_t_coeff              = 0.3;      // default MCSP temperature coefficient
    static inline constexpr Real MCSP_q_coeff              = 0.0;      // default MCSP sp. humidity coefficient
    static inline constexpr Real MCSP_u_coeff              = 0.0;      // default MCSP U-wind coefficient
    static inline constexpr Real MCSP_v_coeff              = 0.0;      // default MCSP V-wind coefficient
    // Default values for ZmRuntimeOpt fields
    static inline constexpr Real alfa                      = 0.14;     // default downdraft proportionality factor
    static inline constexpr Real ke                        = 2.5E-6;   // default evaporation efficiency
    static inline constexpr Real dmpdz                     = -0.7e-3;  // default convective entrainment parameter [1/m]
    static inline constexpr Real tiedke_add                = 0.8;      // default Tiedke temperature perturbation addition [K]
    static inline constexpr Real c0                        = 0.0020;   // default autoconversion coefficient
    static inline constexpr Real auto_fac                  = 7;        // enhancement factor for droplet-rain autoconversion
    static inline constexpr Real accr_fac                  = 1.5;      // default accretion factor
    static inline constexpr Real micro_dcs                 = 150.E-6;  // default size threshold for cloud ice to snow autoconversion

    // Table of saturation vapor pressure values (estbl) from tmin to
    // tmax+1 Kelvin, in one degree increments. ttrice defines the transition
    // region, estbl contains a combination of ice & water values.
    static inline constexpr Real tmin    = 127.16;
    static inline constexpr Real tmax    = 375.16;
    static inline constexpr Real h2otrip = 273.16;
    static inline constexpr Real ttrice  = 20.0;  // transition range from es over H2O to es over ice

    static inline constexpr Real dlf_tk1 = 268.15; // T at/above which detrained condensate (dlf) is assumed to be liquid
    static inline constexpr Real dlf_tk2 = 238.15; // T at/below which detrained condensate (dlf) is assumed to be ice
  };

  //----------------------------------------------------------------------------
  // Purpose: derived type to hold ZM tunable parameters
  //----------------------------------------------------------------------------
  struct ZmRuntimeOpt {
    ZmRuntimeOpt() = default;

    void load_runtime_options(ekat::ParameterList& params) {
      apply_detr_tend     = params.get<bool>("apply_detr_tend",     true);
      use_fortran_bridge  = params.get<bool>("use_fortran_bridge",  true);
      upper_limit_pref    = params.get<Real>("upper_limit_pref",    40e2);
      tau                 = params.get<Real>("tau",                 3600);
      alfa                = params.get<Real>("alfa",                ZMC::alfa);
      ke                  = params.get<Real>("ke",                  ZMC::ke);
      dmpdz               = params.get<Real>("dmpdz",               ZMC::dmpdz);
      tpert_fix           = params.get<bool>("tpert_fix",           true);
      tpert_fac           = params.get<Real>("tpert_fac",           2);
      tiedke_add          = params.get<Real>("tiedke_add",          ZMC::tiedke_add);
      c0_lnd              = params.get<Real>("c0_lnd",              ZMC::c0);
      c0_ocn              = params.get<Real>("c0_ocn",              ZMC::c0);
      num_cin             = params.get<int>("num_cin",              1);
      mx_bot_lyr_adj      = params.get<int>("mx_bot_lyr_adj",       1);
      trig_dcape          = params.get<bool>("trig_dcape",          true);
      trig_ull            = params.get<bool>("trig_ull",            true);
      clos_dyn_adj        = params.get<bool>("clos_dyn_adj",        true);
      no_deep_pbl         = params.get<bool>("no_deep_pbl",         false);
      // ZM micro parameters
      zm_microp           = params.get<bool>("zm_microp",           false);
      old_snow            = params.get<bool>("old_snow",            false);
      auto_fac            = params.get<Real>("auto_fac",            ZMC::auto_fac);
      accr_fac            = params.get<Real>("accr_fac",            ZMC::accr_fac);
      micro_dcs           = params.get<Real>("micro_dcs",           ZMC::micro_dcs);
      // MCSP parameters
      mcsp_enabled        = params.get<bool>("mcsp_enabled",        true);
      mcsp_t_coeff        = params.get<Real>("mcsp_t_coeff",        ZMC::MCSP_t_coeff);
      mcsp_q_coeff        = params.get<Real>("mcsp_q_coeff",        ZMC::MCSP_q_coeff);
      mcsp_u_coeff        = params.get<Real>("mcsp_u_coeff",        ZMC::MCSP_u_coeff);
      mcsp_v_coeff        = params.get<Real>("mcsp_v_coeff",        ZMC::MCSP_v_coeff);

      // determine SVP table size (add two to make the table slightly too big, just in case)
      plenest = static_cast<Int>(ZMC::tmax-ZMC::tmin) + 3;

      // Build the table in a local view first. Referencing the data member
      // 'estbl' directly inside the lambda would capture 'this' (a host
      // pointer), causing an illegal device-memory access on GPU.
      view_1d<Real> estbl_tmp("estbl",plenest);
      Kokkos::parallel_for(Kokkos::RangePolicy<typename KT::ExeSpace>(0, plenest), KOKKOS_LAMBDA(const int i) {
        estbl_tmp(i) = svp_trans(ZMC::tmin + i);
      });
      estbl = estbl_tmp;
    }

    void set_limcnv(std::shared_ptr<const AbstractGrid> grid) {
      EKAT_REQUIRE_MSG(upper_limit_pref > 0,
      "Error! ZmRuntimeOpt::set_limcnv: upper_limit_pref must be positive [Pa], but got "
      << upper_limit_pref << "\n");
      EKAT_REQUIRE_MSG(upper_limit_pref < PC::P0.value,
      "Error! ZmRuntimeOpt::set_limcnv: upper_limit_pref must be less than the reference surface pressure, but got "
      << upper_limit_pref << "\n");
      // Determine upper limit level index of deep convection based on the reference pressure profile
      const auto nlev   = grid->get_num_vertical_levels();
      const auto hyai_h = grid->get_geometry_data("hyai").get_view<const Real*, Host>();
      const auto hybi_h = grid->get_geometry_data("hybi").get_view<const Real*, Host>();
      const auto ps0 = PC::P0.value;
      limcnv = -1;
      if (ps0*hyai_h(0) + ps0*hybi_h(0) >= upper_limit_pref) {
        limcnv = 0;
      } else {
        for (int k = 0; k < nlev; ++k) {
          Real pk0 = ps0*hyai_h(k)   + ps0*hybi_h(k);
          Real pk1 = ps0*hyai_h(k+1) + ps0*hybi_h(k+1);
          if (pk0 < upper_limit_pref && pk1 >= upper_limit_pref) {
            limcnv = k;
            break;
          }
        }
        if (limcnv == -1) { limcnv = nlev+1; }
      }
    }

    // -------------------------------------------------------------------------
    // print parameter values for the log file (C++ analog of zm_param_print)
    void print(std::ostream& os = std::cout) const {
      const std::string indent = "  ";
      // preserve and restore the stream's formatting flags
      const std::ios::fmtflags saved_flags = os.flags();
      os << std::boolalpha;
      os << "\n";
      os << "ZM deep convection parameter values:\n";
      os << indent << "tau             : " << tau            << "\n";
      os << indent << "alfa            : " << alfa           << "\n";
      os << indent << "ke              : " << ke             << "\n";
      os << indent << "dmpdz           : " << dmpdz          << "\n";
      os << indent << "tpert_fix       : " << tpert_fix      << "\n";
      os << indent << "tpert_fac       : " << tpert_fac      << "\n";
      os << indent << "tiedke_add      : " << tiedke_add     << "\n";
      os << indent << "c0_lnd          : " << c0_lnd         << "\n";
      os << indent << "c0_ocn          : " << c0_ocn         << "\n";
      os << indent << "num_cin         : " << num_cin        << "\n";
      os << indent << "limcnv          : " << limcnv         << "\n";
      os << indent << "mx_bot_lyr_adj  : " << mx_bot_lyr_adj << "\n";
      os << indent << "trig_dcape      : " << trig_dcape     << "\n";
      os << indent << "trig_ull        : " << trig_ull       << "\n";
      os << indent << "clos_dyn_adj    : " << clos_dyn_adj   << "\n";
      os << indent << "no_deep_pbl     : " << no_deep_pbl    << "\n";
      // ZM micro parameters
      os << indent << "zm_microp       : " << zm_microp      << "\n";
      os << indent << "old_snow        : " << old_snow       << "\n";
      os << indent << "auto_fac        : " << auto_fac       << "\n";
      os << indent << "accr_fac        : " << accr_fac       << "\n";
      os << indent << "micro_dcs       : " << micro_dcs      << "\n";
      // MCSP parameters
      os << indent << "mcsp_enabled    : " << mcsp_enabled   << "\n";
      os << indent << "mcsp_t_coeff    : " << mcsp_t_coeff   << "\n";
      os << indent << "mcsp_q_coeff    : " << mcsp_q_coeff   << "\n";
      os << indent << "mcsp_u_coeff    : " << mcsp_u_coeff   << "\n";
      os << indent << "mcsp_v_coeff    : " << mcsp_v_coeff   << "\n";
      os << std::endl;
      os.flags(saved_flags);
    }

    Real tau;               // convective adjustment time scale
    Real alfa;              // max downdraft mass flux fraction
    Real ke;                // evaporation efficiency
    Real dmpdz;             // fractional mass entrainment rate [1/m]
    bool tpert_fix;         // flag to disable using applying tpert to PBL-rooted convection
    Real tpert_fac;         // tunable temperature perturbation factor
    Real tiedke_add;        // tunable temperature perturbation
    Real c0_lnd;            // autoconversion coefficient over land
    Real c0_ocn;            // autoconversion coefficient over ocean
    int num_cin;            // num of neg buoyancy regions allowed before the conv top and CAPE calc are completed
    Real upper_limit_pref;  // pressure limit above which deep convection is not allowed [Pa] (used to set limcnv)
    int limcnv;             // upper pressure interface level to limit deep convection
    int mx_bot_lyr_adj;     // bot layer index adjustment for launch level search
    bool trig_dcape;        // true if to using DCAPE trigger - based on CAPE generation by the dycor
    bool trig_ull;          // true if to using the "unrestricted launch level" (ULL) mode
    bool clos_dyn_adj;      // flag for mass flux adjustment to CAPE closure
    bool no_deep_pbl;       // flag to eliminate deep convection within PBL
    bool apply_detr_tend;
    bool use_fortran_bridge;
    // ZM micro parameters
    bool zm_microp;         // switch for convective microphysics
    bool old_snow;          // switch to calculate snow prod in zm_conv_evap() (old treatment before zm_microp was implemented)
    Real auto_fac;          // ZM microphysics enhancement factor for droplet-rain autoconversion
    Real accr_fac;          // ZM microphysics enhancement factor for droplet-rain accretion
    Real micro_dcs;         // ZM microphysics size threshold for cloud ice to snow autoconversion [m]
    // MCSP parameters
    bool mcsp_enabled;      // flag for mesoscale coherent structure parameterization (MSCP)
    Real mcsp_t_coeff;      // MCSP coefficient for temperature tendencies
    Real mcsp_q_coeff;      // MCSP coefficient for specific humidity tendencies
    Real mcsp_u_coeff;      // MCSP coefficient for zonal momentum tendencies
    Real mcsp_v_coeff;      // MCSP coefficient for meridional momentum tendencies
    // saturation vapor pressure (svp) table
    Int plenest;            // saturation vapor pressure table size
    view_1d<Real> estbl;    // saturation vapor pressure table values

  };

  // -----------------------------------------------------------------------------------------------

  struct ZmInputState {
    ZmInputState() = default;
    // -------------------------------------------------------------------------
    // variable counters for device-side only
    static constexpr int num_1d_intgr = 0;  // number of 1D integer views
    static constexpr int num_1d_scalr = 1;  // number of 1D scalar views
    static constexpr int num_2d_midlv = 5;  // number of 2D mid-point views
    static constexpr int num_2d_intfc = 1;  // number of 2D interface views
    static constexpr int num_3d_midlv = 1;  // number of 3D mid-point views (ncol,nwind,nlev)

    uview_1d<     Scalar> tpert;    // PBL top temperature perturb. [K]
    uview_2d<     Real>   z_mid;    // mid-point level altitude     [m]
    uview_2d<     Real>   z_del;    // altitude thickness           [m]
    uview_2d<     Real>   z_int;    // interface level altitude     [m]
    // variables we get from the field manager
    view_1d<const Scalar> phis;     // surface geopotential height  [m2/s]
    view_1d<const Scalar> pblh;     // PBL height                   [m]
    view_1d<const Scalar> landfrac; // land area fraction           [frac]
    view_2d<const Real>   p_mid;    // mid-point level pressure     [Pa]
    view_2d<const Real>   p_int;    // interface level pressure     [Pa]
    view_2d<const Real>   p_del;    // pressure thickness           [Pa]
    view_2d<      Real>   T_mid;    // temperature                  [K]
    view_2d<      Real>   qv;       // water vapor mixing ratio     [kg/kg]
    view_2d<const Real>   qc;       // cloud liquid water           [kg/kg]
    view_2d<      Real>   uwind;    // zonal wind                   [m/s]
    view_2d<      Real>   vwind;    // meridional wind              [m/s]
    view_2d<const Real>   omega;    // vertical pressure velocity   [Pa/s]
    view_2d<const Real>   cldfrac;  // total cloud fraction         [frac]
    view_2d<const Real>   thl_sec;  // thetal variance from SHOC    [K^2]
    // intermediate state variables updated with intermediate tendencies
    uview_2d<     Real>   tmp_s_mid;// dry static energy            [J]
    uview_2d<     Real>   tmp_T_mid;// temperature                  [K]
    uview_2d<     Real>   tmp_qv;   // water vapor mixing ratio     [kg/kg]
    uview_3d<     Real>   tmp_winds;// horizontal winds (ncol,nwind,nlev) [m/s]
    // variables only needed for calling the C++ version of ZM
    view_2d<      Real>   t_prev;   // DCAPE T from previous time step [K]
    view_2d<      Real>   q_prev;   // DCAPE q from previous time step [kg/kg]

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
    uview_2dl<Real> f_t_prev;
    uview_2dl<Real> f_q_prev;
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
    view_2dh<Real>   h_t_prev;
    view_2dh<Real>   h_q_prev;
    view_2dh<Real>   h_z_int;
    view_2dh<Real>   h_p_int;

    // -------------------------------------------------------------------------
    // allocate host mirror input variables (only needed for fortran bridge)
    void init_host_mirrors(int ncol, int nlev) {
      h_phis     = view_1dh<Scalar>("zm_input.h_phis",     ncol);
      h_pblh     = view_1dh<Scalar>("zm_input.h_pblh",     ncol);
      h_tpert    = view_1dh<Scalar>("zm_input.h_tpert",    ncol);
      h_landfrac = view_1dh<Scalar>("zm_input.h_landfrac", ncol);
      h_z_mid    = view_2dh<Real>  ("zm_input.h_z_mid",    ncol, nlev);
      h_p_mid    = view_2dh<Real>  ("zm_input.h_p_mid",    ncol, nlev);
      h_p_del    = view_2dh<Real>  ("zm_input.h_p_del",    ncol, nlev);
      h_T_mid    = view_2dh<Real>  ("zm_input.h_T_mid",    ncol, nlev);
      h_qv       = view_2dh<Real>  ("zm_input.h_qv",       ncol, nlev);
      h_uwind    = view_2dh<Real>  ("zm_input.h_uwind",    ncol, nlev);
      h_vwind    = view_2dh<Real>  ("zm_input.h_vwind",    ncol, nlev);
      h_omega    = view_2dh<Real>  ("zm_input.h_omega",    ncol, nlev);
      h_cldfrac  = view_2dh<Real>  ("zm_input.h_cldfrac",  ncol, nlev);
      h_t_prev   = view_2dh<Real>  ("zm_input.h_t_prev",   ncol, nlev);
      h_q_prev   = view_2dh<Real>  ("zm_input.h_q_prev",   ncol, nlev);
      h_z_int    = view_2dh<Real>  ("zm_input.h_z_int",    ncol, nlev+1);
      h_p_int    = view_2dh<Real>  ("zm_input.h_p_int",    ncol, nlev+1);
    }

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
    static constexpr int num_1d_intgr =  5; // number of 1D integer views
    static constexpr int num_1d_scalr =  6; // number of 1D scalar views
    static constexpr int num_2d_midlv = 24; // number of 2D mid-point views
    static constexpr int num_2d_intfc =  3; // number of 2D interface views
    static constexpr int num_3d_midlv =  1; // number of 3D mid-point views (ncol,nwind,nlev)

    uview_1d<Int>    activity;       // integer deep convection activity flag
    uview_1d<Scalar> prec;           // surface precipitation                   [m/s]
    uview_1d<Scalar> snow;           // surface snow                            [m/s]
    uview_1d<Scalar> cape;           // convective available potential energy   [J]
    uview_2d<Real>   tend_out_t;     // output tendency of temperature          [K/s]
    uview_2d<Real>   tend_out_s;     // output tendency of dry static energy    [J/s]
    uview_2d<Real>   tend_out_qv;    // output tendency of water vapor          [kg/kg/s]
    uview_2d<Real>   tend_out_u;     // output tendency of zonal wind           [m/s/s]
    uview_2d<Real>   tend_out_v;     // output tendency of meridional wind      [m/s/s]
    uview_2d<Real>   tend_tmp_s;     // temporary tendency of dry static energy [J/s]
    uview_2d<Real>   tend_tmp_qv;    // temporary tendency of water vapor       [kg/kg/s]
    uview_3d<Real>   tend_tmp_winds; // temporary wind tendency (ncol,nwind,nlev) - used by MCSP and momentum transport [m/s/s]
    uview_2d<Real>   tend_s_snwprd;  // Heating rate of snow production         [J/kg/s]
    uview_2d<Real>   tend_s_snwevmlt;// Heating rate of snow evap/melt          [J/kg/s]
    uview_2d<Real>   rain_prod;      // rain production rate                    [kg/kg/s]
    uview_2d<Real>   snow_prod;      // snow production rate                    [kg/kg/s]
    uview_2d<Real>   ntprprd;        // net precip production in layer          [kg/kg/s]
    uview_2d<Real>   ntsnprd;        // net snow production in layer            [kg/kg/s]
    uview_2d<Real>   flxprec;        // Convective flux of prec at interfaces   [kg/m2/s]
    uview_2d<Real>   flxsnow;        // Convective flux of snow at interfaces   [kg/m2/s]
    uview_2d<Real>   prec_flux;      // output convective prec flux             [kg/m2/s]
    uview_2d<Real>   snow_flux;      // output convective snow flux             [kg/m2/s]
    uview_2d<Real>   mass_flux;      // output convective mass flux             [mb/s]
    // variables only needed for calling the C++ version of ZM
    uview_1d<Int>    msemax_klev;    // level indices of max MSE
    uview_1d<Int>    jctop;          // top-of-deep-convection indices
    uview_1d<Int>    jcbot;          // base of cloud indices
    uview_1d<Int>    jt;             // top level index of convection
    Int              ktm = 0;        // highest cloud-top level index over active columns
    Int              kbm = 0;        // highest cloud-base level index over active columns
    uview_1d<Scalar> dcape;          // CAPE generated by dycor (dCAPE)         [J]
    uview_2d<Real>   zdu;            // detraining mass flux                    [1/s]
    uview_2d<Real>   mflx_up;        // updraft mass flux                       [mb/s]
    uview_2d<Real>   entr_up;        // updraft entrainment                     [1/s]
    uview_2d<Real>   detr_up;        // updraft detrainment                     [1/s]
    uview_2d<Real>   mflx_dn;        // downdraft mass flux                     [mb/s]
    uview_2d<Real>   entr_dn;        // downdraft entrainment                   [1/s]
    uview_2d<Real>   p_del_mb;       // layer thickness                         [mb]
    uview_1d<Scalar> dsubcld;        // thickness between lcl and msemax_klev   [mb]
    uview_2d<Real>   ql;             // cloud liquid water for chem/wetdep      [?]
    uview_1d<Scalar> rliq;           // reserved liquid (not yet in cldliq) for energy integrals
    uview_2d<Real>   dlf;            // detrainment rate of cloud liquid water  [kg/kg/s]
    // MCSP diagnostic outputs
    view_2d<Real>    mcsp_dt_out;    // MCSP tendency for DSE
    view_2d<Real>    mcsp_dq_out;    // MCSP tendency for qv
    view_2d<Real>    mcsp_du_out;    // MCSP tendency for u wind
    view_2d<Real>    mcsp_dv_out;    // MCSP tendency for v wind
    view_1d<Real>    mcsp_freq;      // MSCP frequency for output
    view_1d<Real>    mcsp_shear;     // MCSP shear used to check against threshold
    view_1d<Real>    zm_depth;       // MCSP pressure depth of ZM heating

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
    uview_2dl<Real>  f_dlf;
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
    view_2dh<Real>   h_dlf;
    view_2dh<Real>   h_prec_flux;
    view_2dh<Real>   h_snow_flux;
    view_2dh<Real>   h_mass_flux;

    // -------------------------------------------------------------------------
    // allocate host mirror output variables (only needed for fortran bridge)
    void init_host_mirrors(int ncol, int nlev) {
      h_activity  = view_1dh<Int>   ("zm_output.h_activity",  ncol);
      h_prec      = view_1dh<Scalar>("zm_output.h_prec",      ncol);
      h_snow      = view_1dh<Scalar>("zm_output.h_snow",      ncol);
      h_cape      = view_1dh<Scalar>("zm_output.h_cape",      ncol);
      h_tend_t    = view_2dh<Real>  ("zm_output.h_tend_t",    ncol, nlev);
      h_tend_qv   = view_2dh<Real>  ("zm_output.h_tend_qv",   ncol, nlev);
      h_tend_u    = view_2dh<Real>  ("zm_output.h_tend_u",    ncol, nlev);
      h_tend_v    = view_2dh<Real>  ("zm_output.h_tend_v",    ncol, nlev);
      h_rain_prod = view_2dh<Real>  ("zm_output.h_rain_prod", ncol, nlev);
      h_snow_prod = view_2dh<Real>  ("zm_output.h_snow_prod", ncol, nlev);
      h_dlf       = view_2dh<Real>  ("zm_output.h_dlf",       ncol, nlev);
      h_prec_flux = view_2dh<Real>  ("zm_output.h_prec_flux", ncol, nlev+1);
      h_snow_flux = view_2dh<Real>  ("zm_output.h_snow_flux", ncol, nlev+1);
      h_mass_flux = view_2dh<Real>  ("zm_output.h_mass_flux", ncol, nlev+1);
    }

    // -------------------------------------------------------------------------
    // transpose method for fortran bridging
    template <ekat::TransposeDirection::Enum D>
    void transpose(int ncol, int nlev_mid);

    // -------------------------------------------------------------------------
    void init_all(int ncol, int nlev_mid); // initialize all variables in struct
    void init_tmp(int ncol, int nlev_mid); // initialize temporary tendencies variables only
  };

  // -----------------------------------------------------------------------------------------------

  //
  // --------- Init/Finalize Functions ---------
  //
  static void zm_opts_init();

  static void zm_finalize() {
    // release Kokkos Views held by the static structs
    s_zm_opts.estbl = view_1d<Real>();
  }

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
    const Real& dt, // time step in seconds : 2*delta_t
    const uview_2d<const Real>& wind_mid, // input Momentum array (nwind, pver)
    const Int& nwind, // number of tracers to transport
    const uview_1d<const Real>& mu, // mass flux up
    const uview_1d<const Real>& md, // mass flux down
    const uview_1d<const Real>& du, // mass detraining from updraft
    const uview_1d<const Real>& eu, // mass entraining from updraft
    const uview_1d<const Real>& ed, // mass entraining from downdraft
    const uview_1d<const Real>& dp, // gathered pressure delta between interfaces
    const Int& jt, // index of cloud top for each column
    const Int& mx, // index of cloud top for each column
    const Int& ktm, // Highest top level for any column
    const Int& kbm, // Highest bottom level for any column
    // Outputs (all wind arrays are (nwind, pver))
    const uview_2d<Real>& wind_tend, // output momentum tendency
    // unused diagnostics - commented out but kept for easy restoration
    // const uview_2d<Real>& pguall, // apparent force from  updraft PG
    // const uview_2d<Real>& pgdall, // apparent force from  downdraft PG
    // const uview_2d<Real>& icwu, // in-cloud winds in updraft
    // const uview_2d<Real>& icwd, // in-cloud winds in downdraft
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

  static void zm_conv_main(
    // Inputs
    const ZmRuntimeOpt& runtime_opt,
    const Int& ncol, // number of columns
    const Int& pver, // number of mid-point levels
    const Int& pverp, // number of interface levels
    const bool& is_first_step, // flag for first step of run
    const Real& time_step, // model time-step                         [s]
    const uview_2d<const Real>& t_mid, // temperature                             [K]       [ncol,pver]
    const uview_2d<const Real>& q_mid_in, // specific humidity                       [kg/kg]  [ncol,pver]
    const uview_2d<const Real>& omega, // vertical pressure velocity              [Pa/s]   [ncol,pver]
    const uview_2d<const Real>& p_mid_in, // mid-point pressure                      [Pa]     [ncol,pver]
    const uview_2d<const Real>& p_int_in, // interface pressure                      [Pa]     [ncol,pverp]
    const uview_2d<const Real>& p_del_in, // pressure thickness                      [Pa]     [ncol,pver]
    const uview_1d<const Real>& geos, // surface geopotential                    [m2/s2]  [ncol]
    const uview_2d<const Real>& z_mid_in, // mid-point geopotential                  [m]      [ncol,pver]
    const uview_2d<const Real>& z_int_in, // interface geopotential                  [m]      [ncol,pverp]
    const uview_1d<const Real>& pbl_hgt, // boundary layer height                   [m]      [ncol]
    const uview_1d<const Real>& tpert, // parcel temperature perturbation         [K]      [ncol]
    const uview_1d<const Real>& landfrac, // land fraction                           []       [ncol]
    const uview_2d<const Real>& t_star, // for DCAPE - prev temperature            [K]      [ncol,pver]
    const uview_2d<const Real>& q_star, // for DCAPE - prev sp. humidity           [kg/kg]  [ncol,pver]
    // Outputs
    const uview_1d<Int>& msemax_klev, // level index of max MSE                  [ncol]
    const uview_1d<Int>& jctop, // top-of-deep-convection index            [ncol]
    const uview_1d<Int>& jcbot, // base of cloud index                     [ncol]
    const uview_1d<Int>& jt, // top level index of convection           [ncol]
    const uview_1d<Int>& active, // deep convection activity flag (1/0)     [ncol]
    const uview_1d<Real>& prec, // output precipitation                    [m/s]    [ncol]
    const uview_2d<Real>& heat, // dry static energy tendency              [W/kg]   [ncol,pver]
    const uview_2d<Real>& qtnd, // specific humidity tendency              [kg/kg/s][ncol,pver]
    const uview_1d<Real>& cape, // conv. avail. potential energy           [J]      [ncol]
    const uview_1d<Real>& dcape, // CAPE generated by dycor (dCAPE)         [J]      [ncol]
    const uview_2d<Real>& mcon, // convective mass flux                    [mb/s]   [ncol,pverp]
    const uview_2d<Real>& pflx, // precip flux at each level               [kg/m2/s][ncol,pverp]
    const uview_2d<Real>& zdu, // detraining mass flux                    [1/s]    [ncol,pver]
    const uview_2d<Real>& mflx_up, // updraft mass flux                       [mb/s]   [ncol,pver]
    const uview_2d<Real>& entr_up, // updraft entrainment                     [1/s]    [ncol,pver]
    const uview_2d<Real>& detr_up, // updraft detrainment                     [1/s]    [ncol,pver]
    const uview_2d<Real>& mflx_dn, // downdraft mass flux                     [mb/s]   [ncol,pver]
    const uview_2d<Real>& entr_dn, // downdraft entrainment                   [1/s]    [ncol,pver]
    const uview_2d<Real>& p_del, // layer thickness                         [mb]     [ncol,pver]
    const uview_1d<Real>& dsubcld, // thickness between lcl and msemax_klev   [mb]     [ncol]
    const uview_2d<Real>& ql, // cloud liquid water for chem/wetdep      [?]      [ncol,pver]
    const uview_1d<Real>& rliq, // reserved liquid (not yet in cldliq)     []       [ncol]
    const uview_2d<Real>& rprd, // rain production rate                    [kg/kg/s][ncol,pver]
    const uview_2d<Real>& dlf, // detrainment rate of cloud liquid water  [kg/kg/s][ncol,pver]
    Int& ktm, // highest cloud-top level index over active columns
    Int& kbm); // highest cloud-base level index over active columns

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
  inline static ZmRuntimeOpt s_zm_opts;

  //
  // inline functions were moved to an impl, so include them here
  //
  #include "impl/zm_inline_functions.hpp"

}; // struct Functions

} // namespace zm
} // namespace scream

#if defined(EAMXX_ENABLE_GPU) && !defined(KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE) \
                              && !defined(KOKKOS_ENABLE_HIP_RELOCATABLE_DEVICE_CODE)
# include "impl/zm_input_state_impl.hpp"
# include "impl/zm_output_tend_impl.hpp"
# include "impl/zm_opts_impl.hpp"
# include "impl/zm_invert_entropy_impl.hpp"
# include "impl/zm_entropy_impl.hpp"
# include "impl/zm_transport_tracer_impl.hpp"
# include "impl/zm_transport_momentum_impl.hpp"
# include "impl/zm_compute_dilute_cape_impl.hpp"
# include "impl/zm_find_mse_max_impl.hpp"
# include "impl/zm_compute_dilute_parcel_impl.hpp"
# include "impl/zm_compute_cape_from_parcel_impl.hpp"
# include "impl/zm_conv_mcsp_calculate_shear_impl.hpp"
# include "impl/zm_conv_mcsp_tend_impl.hpp"
# include "impl/zm_conv_main_impl.hpp"
# include "impl/zm_conv_evap_impl.hpp"
# include "impl/zm_calc_fractional_entrainment_impl.hpp"
# include "impl/zm_downdraft_properties_impl.hpp"
# include "impl/zm_cloud_properties_impl.hpp"
# include "impl/zm_closure_impl.hpp"
# include "impl/zm_calc_output_tend_impl.hpp"
#endif // GPU && !KOKKOS_ENABLE_*_RELOCATABLE_DEVICE_CODE
#endif // ZM_FUNCTIONS_HPP
