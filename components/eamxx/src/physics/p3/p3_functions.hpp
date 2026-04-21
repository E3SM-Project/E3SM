#ifndef P3_FUNCTIONS_HPP
#define P3_FUNCTIONS_HPP

#include "share/physics/physics_constants.hpp"

#include "share/core/eamxx_types.hpp"

#include <ekat_pack_kokkos.hpp>
#include <ekat_parameter_list.hpp>
#include <ekat_workspace.hpp>

namespace scream
{
namespace p3
{

/*
 * Functions is a stateless struct used to encapsulate a
 * number of functions for p3. We use the ETI pattern for
 * these functions.
 *
 * P3 assumptions:
 *  - Kokkos team policies have a vector length of 1
 */

template <typename ScalarT, typename DeviceT> struct Functions {
  //
  // ---------- P3 constants ---------
  //
  struct P3C {
    // Constants for ice lookup tables
    enum {
      densize            = 5,
      rimsize            = 4,
      isize              = 50,
      ice_table_size     = 12, // number of quantities used from lookup table
      rcollsize          = 30,
      collect_table_size = 2, // number of ice-rain collection  quantities used from lookup table

      // switch for warm-rain parameterization
      // 1 => Seifert and Beheng 2001
      // 2 => Beheng 1994
      // 3 => Khairoutdinov and Kogan 2000
      iparam  = 3,
      dnusize = 16,
    };

    static constexpr ScalarT lookup_table_1a_dum1_c =
        4.135985029041767e+00; // 1.0/(0.1*log10(261.7))
    static constexpr const char *p3_lookup_base = SCREAM_DATA_DIR "/tables/p3_lookup_table_1.dat-v";

    static constexpr const char *p3_version =
        "4.1.1"; // TODO: Change this so that the table version and table path is a runtime option.
  };

  //
  // ------- Types --------
  //

  using Scalar = ScalarT;
  using Device = DeviceT;

  using Pack         = ekat::Pack<Scalar, SCREAM_PACK_SIZE>;
  using IntPack = ekat::Pack<Int, SCREAM_PACK_SIZE>;

  using Mask = ekat::Mask<Pack::n>;

  using KT = KokkosTypes<Device>;

  using C = scream::physics::Constants<Scalar>;

  template <typename S> using view_1d = typename KT::template view_1d<S>;
  template <typename S> using view_2d = typename KT::template view_2d<S>;

  // lookup table values for rain shape parameter mu_r
  using view_1d_table = typename KT::template view_1d_table<Scalar, C::MU_R_TABLE_DIM>;

  // lookup table values for rain number- and mass-weighted fallspeeds and ventilation parameters
  using view_2d_table = typename KT::template view_2d_table<Scalar, C::VTABLE_DIM0, C::VTABLE_DIM1>;

  // ice lookup table values
  using view_ice_table = typename KT::template view<
      const Scalar[P3C::densize][P3C::rimsize][P3C::isize][P3C::ice_table_size]>;

  // ice lookup table values for ice-rain collision/collection
  using view_collect_table =
      typename KT::template view<const Scalar[P3C::densize][P3C::rimsize][P3C::isize]
                                             [P3C::rcollsize][P3C::collect_table_size]>;

  // droplet spectral shape parameter for mass spectra, used for Seifert and Beheng (2001)
  // warm rain autoconversion/accretion option only (iparam = 1)
  using view_dnu_table = typename KT::template view_1d_table<Scalar, P3C::dnusize>;

  template <typename S, int N>
  using view_1d_ptr_array = typename KT::template view_1d_ptr_carray<S, N>;

  template <typename S> using uview_1d = typename ekat::template Unmanaged<view_1d<S>>;
  template <typename S> using uview_2d = typename ekat::template Unmanaged<view_2d<S>>;

  using MemberType = typename KT::MemberType;

  using WorkspaceManager = typename ekat::WorkspaceManager<Pack, Device>;
  using Workspace        = typename WorkspaceManager::Workspace;

  // Structure to store p3 runtime options
  struct P3Runtime {

    Scalar max_total_ni                         = 740.0e3;
    Scalar autoconversion_prefactor             = 1350.0;
    Scalar autoconversion_qc_exponent           = 2.47;
    Scalar autoconversion_nc_exponent           = 1.79;
    Scalar autoconversion_radius                = 25.0e-6;
    Scalar accretion_prefactor                  = 67.0;
    Scalar accretion_qc_exponent                = 1.15;
    Scalar accretion_qr_exponent                = 1.15;
    Scalar rain_selfcollection_prefactor        = 5.78;
    Scalar rain_selfcollection_breakup_diameter = 0.00028;
    Scalar constant_mu_rain                     = 1.0;
    Scalar spa_ccn_to_nc_factor                 = 1.0;
    Scalar spa_ccn_to_nc_exponent               = 1.0;
    Scalar cldliq_to_ice_collection_factor      = 0.5;
    Scalar rain_to_ice_collection_factor        = 1.0;
    Scalar min_rime_rho                         = 50.0;
    Scalar max_rime_rho                         = 900.0;
    Scalar immersion_freezing_exponent          = 0.65;
    Scalar deposition_nucleation_exponent       = 0.304;
    Scalar ice_sedimentation_factor             = 1.0;
    bool do_ice_production                      = true;
    bool set_cld_frac_l_to_one                  = false;
    bool set_cld_frac_i_to_one                  = false;
    bool set_cld_frac_r_to_one                  = false;
    bool use_hetfrz_classnuc                    = false;
    bool use_separate_ice_liq_frac              = false;
    bool extra_p3_diags                         = false;

    void
    load_runtime_options_from_file(ekat::ParameterList &params)
    {
      max_total_ni = params.get<double>("max_total_ni", max_total_ni);
      autoconversion_prefactor =
          params.get<double>("autoconversion_prefactor", autoconversion_prefactor);
      autoconversion_qc_exponent =
          params.get<double>("autoconversion_qc_exponent", autoconversion_qc_exponent);
      autoconversion_nc_exponent =
          params.get<double>("autoconversion_nc_exponent", autoconversion_nc_exponent);
      autoconversion_radius = params.get<double>("autoconversion_radius", autoconversion_radius);
      accretion_prefactor   = params.get<double>("accretion_prefactor", accretion_prefactor);
      accretion_qc_exponent = params.get<double>("accretion_qc_exponent", accretion_qc_exponent);
      accretion_qr_exponent = params.get<double>("accretion_qr_exponent", accretion_qr_exponent);
      rain_selfcollection_prefactor =
          params.get<double>("rain_selfcollection_prefactor", rain_selfcollection_prefactor);
      rain_selfcollection_breakup_diameter = params.get<double>(
          "rain_selfcollection_breakup_diameter", rain_selfcollection_breakup_diameter);
      constant_mu_rain       = params.get<double>("constant_mu_rain", constant_mu_rain);
      spa_ccn_to_nc_factor   = params.get<double>("spa_ccn_to_nc_factor", spa_ccn_to_nc_factor);
      spa_ccn_to_nc_exponent = params.get<double>("spa_ccn_to_nc_exponent", spa_ccn_to_nc_exponent);
      cldliq_to_ice_collection_factor =
          params.get<double>("cldliq_to_ice_collection_factor", cldliq_to_ice_collection_factor);
      rain_to_ice_collection_factor =
          params.get<double>("rain_to_ice_collection_factor", rain_to_ice_collection_factor);
      min_rime_rho = params.get<double>("min_rime_rho", min_rime_rho);
      max_rime_rho = params.get<double>("max_rime_rho", max_rime_rho);
      immersion_freezing_exponent =
          params.get<double>("immersion_freezing_exponent", immersion_freezing_exponent);
      deposition_nucleation_exponent =
          params.get<double>("deposition_nucleation_exponent", deposition_nucleation_exponent);
      ice_sedimentation_factor =
          params.get<double>("ice_sedimentation_factor", ice_sedimentation_factor);
      do_ice_production     = params.get<bool>("do_ice_production", do_ice_production);
      set_cld_frac_l_to_one = params.get<bool>("set_cld_frac_l_to_one", set_cld_frac_l_to_one);
      set_cld_frac_i_to_one = params.get<bool>("set_cld_frac_i_to_one", set_cld_frac_i_to_one);
      set_cld_frac_r_to_one = params.get<bool>("set_cld_frac_r_to_one", set_cld_frac_r_to_one);
      use_hetfrz_classnuc   = params.get<bool>("use_hetfrz_classnuc", use_hetfrz_classnuc);
      use_separate_ice_liq_frac =
          params.get<bool>("use_separate_ice_liq_frac", use_separate_ice_liq_frac);
      extra_p3_diags = params.get<bool>("extra_p3_diags", extra_p3_diags);
    }
  };

  // This struct stores prognostic variables evolved by P3.
  struct P3PrognosticState {
    // Cloud mass mixing ratio [kg kg-1]
    view_2d<Pack> qc;
    // Cloud number mixing ratio [# kg-1]
    view_2d<Pack> nc;
    // Rain mass mixing ratio [kg kg-1]
    view_2d<Pack> qr;
    // Rain number mixing ratio [# kg-1]
    view_2d<Pack> nr;
    // Ice total mass mixing ratio [kg kg-1]
    view_2d<Pack> qi;
    // Ice rime mass mixing ratio [kg kg-1]
    view_2d<Pack> qm;
    // Ice total number mixing ratio [# kg-1]
    view_2d<Pack> ni;
    // Ice rime volume mixing ratio [m3 kg-1]
    view_2d<Pack> bm;
    // Water vapor mixing ratio [kg kg-1]
    view_2d<Pack> qv;
    // Potential temperature [K]
    view_2d<Pack> th;
  };

  // This struct stores diagnostic variables used by P3.
  struct P3DiagnosticInputs {
    // CCN activated number tendency [kg-1 s-1]
    view_2d<const Pack> nc_nuceat_tend;
    // CCN prescribed number density [kg-1 s-1]
    view_2d<const Pack> nccn;
    // Activated ice nuclei concentration [1/kg]
    view_2d<const Pack> ni_activated;
    // Assumed SGS 1/(var(qc)/mean(qc)) [kg2/kg2]
    view_2d<const Pack> inv_qc_relvar;
    // Ice cloud fraction
    view_2d<const Pack> cld_frac_i;
    // Liquid cloud fraction
    view_2d<const Pack> cld_frac_l;
    // Rain cloud fraction
    view_2d<const Pack> cld_frac_r;
    // Pressure [Pa]
    view_2d<const Pack> pres;
    // Vertical grid spacing [m]
    view_2d<const Pack> dz;
    // Pressure thickness [Pa]
    view_2d<const Pack> dpres;
    // Exner expression
    view_2d<const Pack> inv_exner;
    // qv from previous step [kg/kg]
    view_2d<const Pack> qv_prev;
    // T from previous step [K]
    view_2d<const Pack> t_prev;
    // Heterogeneous freezing by immersion nucleation [cm^-3 s^-1]
    view_2d<const Pack> hetfrz_immersion_nucleation_tend;
    // Heterogeneous freezing by contact nucleation [cm^-3 s^-1]
    view_2d<const Pack> hetfrz_contact_nucleation_tend;
    // Heterogeneous freezing by deposition nucleation [cm^-3 s^-1]
    view_2d<const Pack> hetfrz_deposition_nucleation_tend;
  };

  // This struct stores diagnostic outputs computed by P3.
  struct P3DiagnosticOutputs {
    // qitend due to deposition/sublimation
    view_2d<Pack> qv2qi_depos_tend;
    // Precipitation rate, liquid [m s-1]
    view_1d<Scalar> precip_liq_surf;
    // Precipitation rate, solid [m s-1]
    view_1d<Scalar> precip_ice_surf;
    // Effective cloud radius [m]
    view_2d<Pack> diag_eff_radius_qc;
    // Effective ice radius [m]
    view_2d<Pack> diag_eff_radius_qi;
    // Effective rain radius [m]
    view_2d<Pack> diag_eff_radius_qr;
    // Bulk density of ice [kg m-3]
    view_2d<Pack> rho_qi;
    // Grid-box average rain flux [kg m^-2 s^-1] pverp
    view_2d<Pack> precip_liq_flux;
    // Grid-box average ice/snow flux [kg m^-2 s^-1] pverp
    view_2d<Pack> precip_ice_flux;
    // Total precipitation (rain + snow) [kg/kg/s]
    view_2d<Pack> precip_total_tend;
    // Evaporation of total precipitation (rain + snow) [kg/kg/s]
    view_2d<Pack> nevapr;
    // Equivalent radar reflectivity [dBz]
    view_2d<Pack> diag_equiv_reflectivity;
  };

  // This struct stores time stepping and grid-index-related information.
  struct P3Infrastructure {
    // Model time step [s]
    Real dt;
    // Time step counter (1-based)
    Int it;
    // Lower bound for horizontal column indices.
    Int its;
    // Upper bound for horizontal column indices.
    Int ite;
    // Lower bound for vertical level indices.
    Int kts;
    // Upper bound for vertical level indices.
    Int kte;
    // Set to true to have P3 predict Nc, false to have Nc specified.
    bool predictNc;
    // Set to true to use prescribed CCN
    bool prescribedCCN;
    // Coordinates of columns, nj x 3
    view_2d<const Scalar> col_location;
  };

  // This struct stores tendencies computed by P3 and used by other
  // parameterizations.
  struct P3HistoryOnly {
    // Sum of liq-ice phase change tendencies
    view_2d<Pack> liq_ice_exchange;
    // Sum of vap-liq phase change tendencies
    view_2d<Pack> vap_liq_exchange;
    // Sum of vap-ice phase change tendencies
    view_2d<Pack> vap_ice_exchange;
    // extra_p3_diags
    view_2d<Pack> qr2qv_evap;
    view_2d<Pack> qi2qv_sublim;
    view_2d<Pack> qc2qr_accret;
    view_2d<Pack> qc2qr_autoconv;
    view_2d<Pack> qv2qi_vapdep;
    view_2d<Pack> qc2qi_berg;
    view_2d<Pack> qc2qr_ice_shed;
    view_2d<Pack> qc2qi_collect;
    view_2d<Pack> qr2qi_collect;
    view_2d<Pack> qc2qi_hetero_freeze;
    view_2d<Pack> qr2qi_immers_freeze;
    view_2d<Pack> qi2qr_melt;
    view_2d<Pack> qr_sed;
    view_2d<Pack> qc_sed;
    view_2d<Pack> qi_sed;
  };

  // This struct stores kokkos views for the lookup tables needed in p3_main()
  struct P3LookupTables {
    // lookup table values for rain shape parameter mu_r
    view_1d_table mu_r_table_vals;
    // lookup table values for rain number- and mass-weighted fallspeeds and ventilation parameters
    view_2d_table vn_table_vals, vm_table_vals, revap_table_vals;
    // ice lookup table values
    view_ice_table ice_table_vals;
    // ice lookup table values for ice-rain collision/collection
    view_collect_table collect_table_vals;
    // droplet spectral shape parameter for mass spectra
    view_dnu_table dnu_table_vals;
  };

#ifdef SCREAM_P3_SMALL_KERNELS
  struct P3Temporaries {
    // shape parameter of rain
    view_2d<Pack> mu_r;
    // temperature at the beginning of the microphysics step [K]
    view_2d<Pack> T_atm;
    // 2D size distribution and fallspeed parameters
    view_2d<Pack> lamr, logn0r, nu;
    view_2d<Pack> cdist, cdist1, cdistr;
    // Variables needed for in-cloud calculations
    // Inverse cloud fractions (1/cld)
    view_2d<Pack> inv_cld_frac_i, inv_cld_frac_l, inv_cld_frac_r;
    // In cloud mass-mixing ratios
    view_2d<Pack> qc_incld, qr_incld, qi_incld, qm_incld;
    // In cloud number concentrations
    view_2d<Pack> nc_incld, nr_incld, ni_incld, bm_incld;
    // Other
    view_2d<Pack> inv_dz, inv_rho, ze_ice, ze_rain;
    view_2d<Pack> prec, rho, rhofacr, rhofaci;
    view_2d<Pack> acn, qv_sat_l, qv_sat_i, sup;
    view_2d<Pack> qv_supersat_i, tmparr2, exner;
    view_2d<Pack> diag_equiv_reflectivity, diag_vm_qi, diag_diam_qi;
    view_2d<Pack> pratot, prctot;
    // p3_tend_out, may not need these
    view_2d<Pack> qtend_ignore, ntend_ignore;
    // Variables still used in F90 but removed from C++ interface
    view_2d<Pack> mu_c, lamc;
    view_2d<Pack> qr_evap_tend;
    // cloud sedimentation
    view_2d<Pack> v_qc, v_nc, flux_qx, flux_nx;
    // ice sedimentation
    view_2d<Pack> v_qit, v_nit, flux_nit, flux_bir;
    view_2d<Pack> flux_qir, flux_qit;
    // rain sedimentation
    view_2d<Pack> v_qr, v_nr;
  };
#endif

  // -- Table3 --

  struct Table3 {
    IntPack dumii, dumjj;
    Pack rdumii, rdumjj;
  };

  struct TableIce {
    IntPack dumi, dumjj, dumii, dumzz;
    Pack dum1, dum4, dum5, dum6;
  };

  struct TableRain {
    IntPack dumj;
    Pack dum3;
  };

  //
  // --------- Functions ---------
  //

  // Call to get global tables
  static void get_global_tables(view_2d_table &vn_table_vals, view_2d_table &vm_table_vals,
                                view_2d_table &revap_table_vals, view_1d_table &mu_r_table_vals,
                                view_dnu_table &dnu);

  static void get_global_ice_lookup_tables(view_ice_table &ice_table_vals,
                                           view_collect_table &collect_table_vals);

  static P3LookupTables p3_init(const bool write_tables = false, const bool masterproc = false);

  // Map (mu_r, lamr) to Table3 data.
  KOKKOS_FUNCTION
  static void lookup(const Pack &mu_r, const Pack &lamr, Table3 &tab,
                     const Mask &context = Mask(true));

  // Converts quantities to cell averages
  KOKKOS_FUNCTION
  static void back_to_cell_average(
      const Pack &cld_frac_l, const Pack &cld_frac_r, const Pack &cld_frac_i,
      Pack &qc2qr_accret_tend, Pack &qr2qv_evap_tend, Pack &qc2qr_autoconv_tend,
      Pack &nc_accret_tend, Pack &nc_selfcollect_tend, Pack &nc2nr_autoconv_tend,
//[shanyp 20260220
//      Pack &nr_selfcollect_tend, Pack &nr_evap_tend, Pack &ncautr, Pack &qi2qv_sublim_tend,
      Pack &nr_selfcollect_tend, Pack &nr_breakup_tend, Pack &nr_evap_tend, Pack &ncautr, Pack &qi2qv_sublim_tend,
//shanyp 20260220]
      Pack &nr_ice_shed_tend, Pack &qc2qi_hetero_freeze_tend, Pack &qr2qi_collect_tend,
      Pack &qc2qr_ice_shed_tend, Pack &qi2qr_melt_tend, Pack &qc2qi_collect_tend,
      Pack &qr2qi_immers_freeze_tend, Pack &ni2nr_melt_tend, Pack &nc_collect_tend,
      Pack &ncshdc, Pack &nc2ni_immers_freeze_tend, Pack &nr_collect_tend,
      Pack &ni_selfcollect_tend, Pack &qv2qi_vapdep_tend, Pack &nr2ni_immers_freeze_tend,
      Pack &ni_sublim_tend, Pack &qv2qi_nucleat_tend, Pack &ni_nucleat_tend,
      Pack &qc2qi_berg_tend, Pack &ncheti_cnt, Pack &qcheti_cnt, Pack &nicnt, Pack &qicnt,
      Pack &ninuc_cnt, Pack &qinuc_cnt, const Mask &context = Mask(true),
      const P3Runtime &runtime_options = {});

  //------------------------------------------------------------------------------------------!
  // Finds indices in 3D ice (only) lookup table
  // ------------------------------------------------------------------------------------------!
  KOKKOS_FUNCTION
  static void lookup_ice(const Pack &qi, const Pack &ni, const Pack &qm, const Pack &rhop,
                         TableIce &tab, const Mask &context = Mask(true));

  //------------------------------------------------------------------------------------------!
  // Finds indices in 3D rain lookup table
  //------------------------------------------------------------------------------------------!
  KOKKOS_FUNCTION
  static void lookup_rain(const Pack &qr, const Pack &nr, TableRain &tab,
                          const Mask &context = Mask(true));

  // Apply Table3 data to the table to return a value. This performs bilinear
  // interpolation within the quad given by {tab.dumii, tab.dumjj} x {t.dumii+1,
  // tab.dumjj+1}.
  KOKKOS_FUNCTION
  static Pack apply_table(const view_2d_table &table, const Table3 &t);

  // Apply TableIce data to the ice tables to return a value.
  KOKKOS_FUNCTION
  static Pack apply_table_ice(const int &index, const view_ice_table &ice_table_vals,
                               const TableIce &tab, const Mask &context = Mask(true));

  // Interpolates lookup table values for rain/ice collection processes
  KOKKOS_FUNCTION
  static Pack apply_table_coll(const int &index, const view_collect_table &collect_table_vals,
                                const TableIce &ti, const TableRain &tr,
                                const Mask &context = Mask(true));

  // -- Sedimentation time step

  // Calculate the first-order upwind step in the region [k_bot,
  // k_top]. Velocity V is input, and flux is workspace and need not be
  // initialized. On input, r contains mixing ratio data at the time step start;
  // on output, it contains mixing ratio data at the time step end.
  // kdir = 1 -> vertical columns are processed from bottom to top, opposite for kdir = -1
  //
  // A subtlety is that this procedure does not do exact upwind of a mixing
  // ratio. That is because the background density rho is assumed to be static;
  // rho does not get advected. Thus, there is an inconsistency between rho and
  // r*rho at the level of |r|.

  // Evolve nfield mixing ratios simultaneously. nfield is a compile-time
  // parameter so the loops over nfield are compiled efficiently. So far the use
  // cases have no need of a runtime version.
  template <int nfield>
  KOKKOS_FUNCTION static void
  calc_first_order_upwind_step(const uview_1d<const Pack> &rho,
                               const uview_1d<const Pack> &inv_rho, // 1/rho
                               const uview_1d<const Pack> &inv_dz, const MemberType &team,
                               const Int &nk, const Int &k_bot, const Int &k_top, const Int &kdir,
                               const Scalar &dt_sub,
                               const view_1d_ptr_array<Pack, nfield> &flux, // workspace
                               const view_1d_ptr_array<Pack, nfield> &V,    // (behaviorally const)
                               const view_1d_ptr_array<Pack, nfield> &r);   // in/out

  // Evolve 1 mixing ratio. This is a syntax-convenience version of the above.
  KOKKOS_FUNCTION
  static void calc_first_order_upwind_step(const uview_1d<const Pack> &rho,
                                           const uview_1d<const Pack> &inv_rho, // 1/rho
                                           const uview_1d<const Pack> &inv_dz,
                                           const MemberType &team, const Int &nk, const Int &k_bot,
                                           const Int &k_top, const Int &kdir, const Scalar &dt_sub,
                                           const uview_1d<Pack> &flux,
                                           const uview_1d<const Pack> &V,
                                           const uview_1d<Pack> &r);

  // This is the main routine. It can be called by the user if kdir is known at
  // compile time. So far it is not, so the above versions are called instead.
  template <Int kdir, int nfield>
  KOKKOS_FUNCTION static void calc_first_order_upwind_step(
      const uview_1d<const Pack> &rho, const uview_1d<const Pack> &inv_rho,
      const uview_1d<const Pack> &inv_dz, const MemberType &team, const Int &nk, const Int &k_bot,
      const Int &k_top, const Scalar &dt_sub, const view_1d_ptr_array<Pack, nfield> &flux,
      const view_1d_ptr_array<Pack, nfield> &V, // (behaviorally const)
      const view_1d_ptr_array<Pack, nfield> &r);

  template <int nfield>
  KOKKOS_FUNCTION static void
  generalized_sedimentation(const uview_1d<const Pack> &rho, const uview_1d<const Pack> &inv_rho,
                            const uview_1d<const Pack> &inv_dz, const MemberType &team,
                            const Int &nk, const Int &k_qxtop, Int &k_qxbot, const Int &kbot,
                            const Int &kdir, const Scalar &Co_max, Scalar &dt_left,
                            Scalar &prt_accum, const view_1d_ptr_array<Pack, nfield> &fluxes,
                            const view_1d_ptr_array<Pack, nfield> &Vs, // (behaviorally const)
                            const view_1d_ptr_array<Pack, nfield> &rs);

  // Cloud sedimentation
  KOKKOS_FUNCTION
  static void cloud_sedimentation(
      const uview_1d<Pack> &qc_incld, const uview_1d<const Pack> &rho,
      const uview_1d<const Pack> &inv_rho, const uview_1d<const Pack> &cld_frac_l,
      const uview_1d<const Pack> &acn, const uview_1d<const Pack> &inv_dz,
      const view_dnu_table &dnu, const MemberType &team, const Workspace &workspace, const Int &nk,
      const Int &ktop, const Int &kbot, const Int &kdir, const Scalar &dt, const Scalar &inv_dt,
      const bool &do_predict_nc, const uview_1d<Pack> &qc, const uview_1d<Pack> &nc,
      const uview_1d<Pack> &nc_incld, const uview_1d<Pack> &mu_c, const uview_1d<Pack> &lamc,
      const uview_1d<Pack> &qc_tend, const uview_1d<Pack> &nc_tend, Scalar &precip_liq_surf);

#ifdef SCREAM_P3_SMALL_KERNELS
  static void cloud_sedimentation_disp(
      const uview_2d<Pack> &qc_incld, const uview_2d<const Pack> &rho,
      const uview_2d<const Pack> &inv_rho, const uview_2d<const Pack> &cld_frac_l,
      const uview_2d<const Pack> &acn, const uview_2d<const Pack> &inv_dz,
      const view_dnu_table &dnu, const WorkspaceManager &workspace_mgr, const Int &nj,
      const Int &nk, const Int &ktop, const Int &kbot, const Int &kdir, const Scalar &dt,
      const Scalar &inv_dt, const bool &do_predict_nc, const uview_2d<Pack> &qc,
      const uview_2d<Pack> &nc, const uview_2d<Pack> &nc_incld, const uview_2d<Pack> &mu_c,
      const uview_2d<Pack> &lamc, const uview_2d<Pack> &qc_tend, const uview_2d<Pack> &nc_tend,
      const uview_1d<Scalar> &precip_liq_surf, const uview_1d<bool> &is_nucleat_possible,
      const uview_1d<bool> &is_hydromet_present);
#endif

  // TODO: comment
  KOKKOS_FUNCTION
  static void rain_sedimentation(
      const uview_1d<const Pack> &rho, const uview_1d<const Pack> &inv_rho,
      const uview_1d<const Pack> &rhofacr, const uview_1d<const Pack> &cld_frac_r,
      const uview_1d<const Pack> &inv_dz, const uview_1d<Pack> &qr_incld, const MemberType &team,
      const Workspace &workspace, const view_2d_table &vn_table_vals,
      const view_2d_table &vm_table_vals, const Int &nk, const Int &ktop, const Int &kbot,
      const Int &kdir, const Scalar &dt, const Scalar &inv_dt, const uview_1d<Pack> &qr,
      const uview_1d<Pack> &nr, const uview_1d<Pack> &nr_incld, const uview_1d<Pack> &mu_r,
      const uview_1d<Pack> &lamr, const uview_1d<Pack> &precip_liq_flux,
      const uview_1d<Pack> &qr_tend, const uview_1d<Pack> &nr_tend, Scalar &precip_liq_surf,
      const P3Runtime &runtime_options);

#ifdef SCREAM_P3_SMALL_KERNELS
  static void rain_sedimentation_disp(
      const uview_2d<const Pack> &rho, const uview_2d<const Pack> &inv_rho,
      const uview_2d<const Pack> &rhofacr, const uview_2d<const Pack> &cld_frac_r,
      const uview_2d<const Pack> &inv_dz, const uview_2d<Pack> &qr_incld,
      const WorkspaceManager &workspace_mgr, const view_2d_table &vn_table_vals,
      const view_2d_table &vm_table_vals, const Int &nj, const Int &nk, const Int &ktop,
      const Int &kbot, const Int &kdir, const Scalar &dt, const Scalar &inv_dt,
      const uview_2d<Pack> &qr, const uview_2d<Pack> &nr, const uview_2d<Pack> &nr_incld,
      const uview_2d<Pack> &mu_r, const uview_2d<Pack> &lamr,
      const uview_2d<Pack> &precip_liq_flux, const uview_2d<Pack> &qr_tend,
      const uview_2d<Pack> &nr_tend, const uview_1d<Scalar> &precip_liq_surf,
      const uview_1d<bool> &is_nucleat_possible, const uview_1d<bool> &is_hydromet_present,
      const P3Runtime &runtime_options);
#endif

  // TODO: comment
  KOKKOS_FUNCTION
  static void ice_sedimentation(
      const uview_1d<const Pack> &rho, const uview_1d<const Pack> &inv_rho,
      const uview_1d<const Pack> &rhofaci, const uview_1d<const Pack> &cld_frac_i,
      const uview_1d<const Pack> &inv_dz, const MemberType &team, const Workspace &workspace,
      const Int &nk, const Int &ktop, const Int &kbot, const Int &kdir, const Scalar &dt,
      const Scalar &inv_dt, const uview_1d<Pack> &qi, const uview_1d<Pack> &qi_incld,
      const uview_1d<Pack> &ni, const uview_1d<Pack> &ni_incld, const uview_1d<Pack> &qm,
      const uview_1d<Pack> &qm_incld, const uview_1d<Pack> &bm, const uview_1d<Pack> &bm_incld,
      const uview_1d<Pack> &qi_tend, const uview_1d<Pack> &ni_tend,
      const view_ice_table &ice_table_vals, Scalar &precip_ice_surf,
      const P3Runtime &runtime_options);

#ifdef SCREAM_P3_SMALL_KERNELS
  static void ice_sedimentation_disp(
      const uview_2d<const Pack> &rho, const uview_2d<const Pack> &inv_rho,
      const uview_2d<const Pack> &rhofaci, const uview_2d<const Pack> &cld_frac_i,
      const uview_2d<const Pack> &inv_dz, const WorkspaceManager &workspace_mgr, const Int &nj,
      const Int &nk, const Int &ktop, const Int &kbot, const Int &kdir, const Scalar &dt,
      const Scalar &inv_dt, const uview_2d<Pack> &qi, const uview_2d<Pack> &qi_incld,
      const uview_2d<Pack> &ni, const uview_2d<Pack> &ni_incld, const uview_2d<Pack> &qm,
      const uview_2d<Pack> &qm_incld, const uview_2d<Pack> &bm, const uview_2d<Pack> &bm_incld,
      const uview_2d<Pack> &qi_tend, const uview_2d<Pack> &ni_tend,
      const view_ice_table &ice_table_vals, const uview_1d<Scalar> &precip_ice_surf,
      const uview_1d<bool> &is_nucleat_possible, const uview_1d<bool> &is_hydromet_present,
      const P3Runtime &runtime_options);
#endif

  // homogeneous freezing of cloud and rain
  KOKKOS_FUNCTION
  static void homogeneous_freezing(const uview_1d<const Pack> &T_atm,
                                   const uview_1d<const Pack> &inv_exner, const MemberType &team,
                                   const Int &nk, const Int &ktop, const Int &kbot, const Int &kdir,
                                   const uview_1d<Pack> &qc, const uview_1d<Pack> &nc,
                                   const uview_1d<Pack> &qr, const uview_1d<Pack> &nr,
                                   const uview_1d<Pack> &qi, const uview_1d<Pack> &ni,
                                   const uview_1d<Pack> &qm, const uview_1d<Pack> &bm,
                                   const uview_1d<Pack> &th_atm);

#ifdef SCREAM_P3_SMALL_KERNELS
  static void homogeneous_freezing_disp(
      const uview_2d<const Pack> &T_atm, const uview_2d<const Pack> &inv_exner, const Int &nj,
      const Int &nk, const Int &ktop, const Int &kbot, const Int &kdir, const uview_2d<Pack> &qc,
      const uview_2d<Pack> &nc, const uview_2d<Pack> &qr, const uview_2d<Pack> &nr,
      const uview_2d<Pack> &qi, const uview_2d<Pack> &ni, const uview_2d<Pack> &qm,
      const uview_2d<Pack> &bm, const uview_2d<Pack> &th_atm,
      const uview_1d<bool> &is_nucleat_possible, const uview_1d<bool> &is_hydromet_present);
#endif

  // -- Find layers

  // Find the bottom and top of the mixing ratio, e.g., qr. It's worth casing
  // these out in two ways: 1 thread/column vs many, and by kdir.
  KOKKOS_FUNCTION
  static Int find_bottom(const MemberType &team, const uview_1d<const Scalar> &v,
                         const Scalar &small, const Int &kbot, const Int &ktop, const Int &kdir,
                         bool &log_present);

  KOKKOS_FUNCTION
  static Int find_top(const MemberType &team, const uview_1d<const Scalar> &v, const Scalar &small,
                      const Int &kbot, const Int &ktop, const Int &kdir, bool &log_present);

  KOKKOS_FUNCTION
  static void cloud_water_conservation(
      const Pack &qc, const Scalar dt, Pack &qc2qr_autoconv_tend, Pack &qc2qr_accret_tend,
      Pack &qc2qi_collect_tend, Pack &qc2qi_hetero_freeze_tend, Pack &qc2qr_ice_shed_tend,
      Pack &qc2qi_berg_tend, Pack &qi2qv_sublim_tend, Pack &qv2qi_vapdep_tend, Pack &qcheti_cnt,
      Pack &qicnt, const bool &use_hetfrz_classnuc, const Mask &context = Mask(true),
      const Pack &cld_frac_l = Pack(), const Pack &cld_frac_i = Pack(),
      const P3Runtime &runtime_options = {});

  KOKKOS_FUNCTION
  static void rain_water_conservation(const Pack &qr, const Pack &qc2qr_autoconv_tend,
                                      const Pack &qc2qr_accret_tend, const Pack &qi2qr_melt_tend,
                                      const Pack &qc2qr_ice_shed_tend, const Scalar dt,
                                      Pack &qr2qv_evap_tend, Pack &qr2qi_collect_tend,
                                      Pack &qr2qi_immers_freeze_tend,
                                      const Mask &context = Mask(true));

  KOKKOS_FUNCTION
  static void ice_water_conservation(
      const Pack &qi, const Pack &qv2qi_vapdep_tend, const Pack &qv2qi_nucleat_tend,
      const Pack &qc2qi_berg_tend, const Pack &qr2qi_collect_tend,
      const Pack &qc2qi_collect_tend, const Pack &qr2qi_immers_freeze_tend,
      const Pack &qc2qi_hetero_freeze_tend, const Scalar dt, Pack &qinuc_cnt, Pack &qcheti_cnt,
      Pack &qicnt, Pack &qi2qv_sublim_tend, Pack &qi2qr_melt_tend,
      const bool &use_hetfrz_classnuc, const Mask &context = Mask(true));

  // TODO: comment
  KOKKOS_FUNCTION
  static void get_cloud_dsd2(const Pack &qc, Pack &nc, Pack &mu_c, const Pack &rho, Pack &nu,
                             const view_dnu_table &dnu, Pack &lamc, Pack &cdist, Pack &cdist1,
                             const Mask &context = Mask(true));

  // Computes and returns rain size distribution parameters
  KOKKOS_FUNCTION
  static void get_rain_dsd2(const Pack &qr, Pack &nr, Pack &mu_r, Pack &lamr,
                            const P3Runtime &runtime_options, const Mask &context = Mask(true));

  // Computes and returns additional rain size distribution parameters
  KOKKOS_FUNCTION
  static void get_cdistr_logn0r(const Pack &qr, const Pack &nr, const Pack &mu_r,
                                const Pack &lamr, Pack &cdistr, Pack &logn0r,
                                const Mask &context = Mask(true));

  // Calculates rime density
  KOKKOS_FUNCTION
  static void calc_rime_density(const Pack &T_atm, const Pack &rhofaci,
                                const Pack &table_val_qi_fallspd, const Pack &acn,
                                const Pack &lamc, const Pack &mu_c, const Pack &qc_incld,
                                const Pack &qc2qi_collect_tend, Pack &vtrmi1, Pack &rho_qm_cloud,
                                const Mask &context = Mask(true));

  // Computes contact and immersion freezing droplets
  KOKKOS_FUNCTION
  static void cldliq_immersion_freezing(const Pack &T_atm, const Pack &lamc, const Pack &mu_c,
                                        const Pack &cdist1, const Pack &qc_incld,
                                        const Pack &inv_qc_relvar, Pack &qc2qi_hetero_freeze_tend,
                                        Pack &nc2ni_immers_freeze_tend,
                                        const P3Runtime &runtime_options,
                                        const Mask &context = Mask(true));

  // Computes the immersion freezing of rain
  KOKKOS_FUNCTION
  static void rain_immersion_freezing(const Pack &T_atm, const Pack &lamr, const Pack &mu_r,
                                      const Pack &cdistr, const Pack &qr_incld,
                                      Pack &qr2qi_immers_freeze_tend,
                                      Pack &nr2ni_immers_freeze_tend,
                                      const P3Runtime &runtime_options,
                                      const Mask &context = Mask(true));

  //
  KOKKOS_FUNCTION
  static void ice_classical_nucleation(const Pack &frzimm, const Pack &frzcnt,
                                       const Pack &frzdep, const Pack &rho, const Pack &qc_incld,
                                       const Pack &nc_incld, const int Iflag, Pack &ncheti_cnt,
                                       Pack &qcheti_cnt, Pack &nicnt, Pack &qucnt,
                                       Pack &ninuc_cnt, Pack &qinuc_cnt);

  // Computes droplet self collection
  KOKKOS_FUNCTION
  static void droplet_self_collection(const Pack &rho, const Pack &inv_rho, const Pack &qc_incld,
                                      const Pack &mu_c, const Pack &nu,
                                      const Pack &nc2nr_autoconv_tend, Pack &nc_selfcollect_tend,
                                      const Mask &context = Mask(true));

  // Computes the accretion of clouds by rain
  KOKKOS_FUNCTION
  static void cloud_rain_accretion(const Pack &rho, const Pack &inv_rho, const Pack &qc_incld,
                                   const Pack &nc_incld, const Pack &qr_incld,
                                   const Pack &inv_qc_relvar, Pack &qc2qr_accret_tend,
                                   Pack &nc_accret_tend, const P3Runtime &runtime_options,
                                   const Mask &context = Mask(true));

  // Computes cloud water autoconversion process rate
  KOKKOS_FUNCTION
  static void cloud_water_autoconversion(const Pack &rho, const Pack &qc_incld,
                                         const Pack &nc_incld, const Pack &inv_qc_relvar,
                                         Pack &qc2qr_autoconv_tend, Pack &nc2nr_autoconv_tend,
                                         Pack &ncautr, const P3Runtime &runtime_options,
                                         const Mask &context = Mask(true));

  // Computes rain self collection process rate
  KOKKOS_FUNCTION
  static void rain_self_collection(const Pack &rho, const Pack &qr_incld, const Pack &nr_incld,
//[shanyp 20260220
//		  Pack &nr_selfcollect_tend, const P3Runtime &runtime_options,
		Pack &nr_selfcollect_tend, Pack &nr_breakup_tend, const P3Runtime &runtime_options,
//shanyp 20260220]
                                   const Mask &context = Mask(true));

  // Impose maximum ice number
  KOKKOS_FUNCTION
  static void impose_max_total_ni(Pack &ni_local, const Scalar &max_total_ni,
                                  const Pack &inv_rho_local, const Mask &context = Mask(true));

  //--------------------------------------------------------------------------------
  //  Calculates and returns the bulk rime density from the prognostic ice variables
  //  and adjusts qm and bm appropriately.
  //--------------------------------------------------------------------------------
  KOKKOS_FUNCTION
  static Pack calc_bulk_rho_rime(const Pack &qi_tot, Pack &qi_rim, Pack &bi_rim,
                                  const P3Runtime &runtime_options,
                                  const Mask &context = Mask(true));

  // TODO - comment
  KOKKOS_FUNCTION
  static void compute_rain_fall_velocity(const view_2d_table &vn_table_vals,
                                         const view_2d_table &vm_table_vals, const Pack &qr_incld,
                                         const Pack &rhofacr, Pack &nr_incld, Pack &mu_r,
                                         Pack &lamr, Pack &V_qr, Pack &V_nr,
                                         const P3Runtime &runtime_options,
                                         const Mask &context = Mask(true));

  //---------------------------------------------------------------------------------
  // update prognostic microphysics and thermodynamics variables
  //---------------------------------------------------------------------------------
  //-- ice-phase dependent processes:
  KOKKOS_FUNCTION
  static void update_prognostic_ice(
      const Pack &qc2qi_hetero_freeze_tend, const Pack &qc2qi_collect_tend,
      const Pack &qc2qr_ice_shed_tend, const Pack &nc_collect_tend,
      const Pack &nc2ni_immers_freeze_tend, const Pack &ncshdc, const Pack &qr2qi_collect_tend,
      const Pack &nr_collect_tend, const Pack &qr2qi_immers_freeze_tend,
      const Pack &nr2ni_immers_freeze_tend, const Pack &nr_ice_shed_tend,
      const Pack &qi2qr_melt_tend, const Pack &ni2nr_melt_tend, const Pack &qi2qv_sublim_tend,
      const Pack &qv2qi_vapdep_tend, const Pack &qv2qi_nucleat_tend, const Pack &ni_nucleat_tend,
      const Pack &ni_selfcollect_tend, const Pack &ni_sublim_tend, const Pack &qc2qi_berg_tend,
      const Pack &inv_exner, const bool do_predict_nc, const Mask &log_wetgrowth, const Scalar dt,
      const Scalar &nmltratio, const Pack &rho_qm_cloud, Pack &ncheti_cnt, Pack &nicnt,
      Pack &ninuc_cnt, Pack &qcheti_cnt, Pack &qicnt, Pack &qinuc_cnt, Pack &th_atm, Pack &qv,
      Pack &qi, Pack &ni, Pack &qm, Pack &bm, Pack &qc, Pack &nc, Pack &qr, Pack &nr,
      const bool &use_hetfrz_classnuc, const Mask &context = Mask(true));

  // TODO (comments)
  KOKKOS_FUNCTION
  static void get_time_space_phys_variables(const Pack &T_atm, const Pack &pres, const Pack &rho,
                                            const Pack &qv_sat_l, const Pack &qv_sat_i, Pack &mu,
                                            Pack &dv, Pack &sc, Pack &dqsdt, Pack &dqsidt,
                                            Pack &ab, Pack &abi, Pack &kap, Pack &eii,
                                            const Mask &context = Mask(true));

  // TODO (comments)
  KOKKOS_FUNCTION
  static void ice_cldliq_collection(const Pack &rho, const Pack &temp, const Pack &rhofaci,
                                    const Pack &table_val_qc2qi_collect, const Pack &qi_incld,
                                    const Pack &qc_incld, const Pack &ni_incld,
                                    const Pack &nc_incld, Pack &qc2qi_collect_tend,
                                    Pack &nc_collect_tend, Pack &qc2qr_ice_shed_tend,
                                    Pack &ncshdc, const P3Runtime &runtime_options,
                                    const Mask &context = Mask(true));

  // TODO (comments)
  KOKKOS_FUNCTION
  static void ice_rain_collection(const Pack &rho, const Pack &temp, const Pack &rhofaci,
                                  const Pack &logn0r, const Pack &table_val_nr_collect,
                                  const Pack &table_val_qr2qi_collect, const Pack &qi_incld,
                                  const Pack &ni_incld, const Pack &qr_incld,
                                  Pack &qr2qi_collect_tend, Pack &nr_collect_tend,
                                  const P3Runtime &runtime_options,
                                  const Mask &context = Mask(true));

  // TODO (comments)
  KOKKOS_FUNCTION
  static void ice_self_collection(const Pack &rho, const Pack &rhofaci,
                                  const Pack &table_val_ni_self_collect, const Pack &eii,
                                  const Pack &qm_incld, const Pack &qi_incld,
                                  const Pack &ni_incld, Pack &ni_selfcollect_tend,
                                  const Mask &context = Mask(true));

  // helper fn for evaporate_rain
  KOKKOS_FUNCTION
  static void rain_evap_tscale_weight(const Pack &dt_over_tau, Pack &weight,
                                      const Mask &context = Mask(true));

  // helper fn for evaporate_rain
  KOKKOS_FUNCTION
  static void rain_evap_equilib_tend(const Pack &A_c, const Pack &ab, const Pack &tau_eff,
                                     const Pack &tau_r, Pack &tend,
                                     const Mask &context = Mask(true));

  // helper fn for evaporate_rain
  KOKKOS_FUNCTION
  static void rain_evap_instant_tend(const Pack &ssat_r, const Pack &ab, const Pack &tau_r,
                                     Pack &tend, const Mask &context = Mask(true));

  // TODO (comments)
  KOKKOS_FUNCTION
  static void evaporate_rain(const Pack &qr_incld, const Pack &qc_incld, const Pack &nr_incld,
                             const Pack &qi_incld, const Pack &cld_frac_l,
                             const Pack &cld_frac_r, const Pack &qv, const Pack &qv_prev,
                             const Pack &qv_sat_l, const Pack &qv_sat_i, const Pack &ab,
                             const Pack &abi, const Pack &epsr, const Pack &epsi_tot,
                             const Pack &t, const Pack &t_prev, const Pack &dqsdt,
                             const Scalar &dt, Pack &qr2qv_evap_tend, Pack &nr_evap_tend,
                             const Mask &context = Mask(true));

  // get number and mass tendencies due to melting ice
  KOKKOS_FUNCTION
  static void ice_melting(const Pack &rho, const Pack &T_atm, const Pack &pres,
                          const Pack &rhofaci, const Pack &table_val_qi2qr_melting,
                          const Pack &table_val_qi2qr_vent_melt, const Pack &dv, const Pack &sc,
                          const Pack &mu, const Pack &kap, const Pack &qv, const Pack &qi_incld,
                          const Pack &ni_incld, Pack &qi2qr_melt_tend, Pack &ni2nr_melt_tend,
                          const Mask &context = Mask(true));

  // liquid-phase dependent processes:
  KOKKOS_FUNCTION
  static void update_prognostic_liquid(
      const Pack &qc2qr_accret_tend, const Pack &nc_accret_tend, const Pack &qc2qr_autoconv_tend,
      const Pack &nc2nr_autoconv_tend, const Pack &ncautr, const Pack &nc_selfcollect_tend,
//[shanyp 20260220
//      const Pack &qr2qv_evap_tend, const Pack &nr_evap_tend, const Pack &nr_selfcollect_tend,
      const Pack &qr2qv_evap_tend, const Pack &nr_evap_tend, const Pack &nr_selfcollect_tend, const Pack &nr_breakup_tend,
//shanyp 20260220]
      const bool do_predict_nc, const bool do_prescribed_CCN, const Pack &inv_rho,
      const Pack &inv_exner, const Scalar dt, Pack &th_atm, Pack &qv, Pack &qc, Pack &nc,
      Pack &qr, Pack &nr, const Mask &context = Mask(true));

  // compute deposition onto ice or sublimation from ice
  KOKKOS_FUNCTION
  static void ice_deposition_sublimation(const Pack &qi_incld, const Pack &ni_incld,
                                         const Pack &t_atm, const Pack &qv_sat_l,
                                         const Pack &qv_sat_i, const Pack &epsi, const Pack &abi,
                                         const Pack &qv, const Scalar &inv_dt, Pack &qidep,
                                         Pack &qi2qv_sublim_tend, Pack &ni_sublim_tend,
                                         Pack &qiberg, const Mask &context = Mask(true));

  KOKKOS_FUNCTION
  static void ice_relaxation_timescale(const Pack &rho, const Pack &temp, const Pack &rhofaci,
                                       const Pack &table_val_qi2qr_melting,
                                       const Pack &table_val_qi2qr_vent_melt, const Pack &dv,
                                       const Pack &mu, const Pack &sc, const Pack &qi_incld,
                                       const Pack &ni_incld, Pack &epsi, Pack &epsi_tot,
                                       const Mask &context = Mask(true));

  KOKKOS_FUNCTION
  static void calc_liq_relaxation_timescale(const view_2d_table &revap_table_vals, const Pack &rho,
                                            const Scalar &f1r, const Scalar &f2r, const Pack &dv,
                                            const Pack &mu, const Pack &sc, const Pack &mu_r,
                                            const Pack &lamr, const Pack &cdistr,
                                            const Pack &cdist, const Pack &qr_incld,
                                            const Pack &qc_incld, Pack &epsr, Pack &epsc,
                                            const Mask &context = Mask(true));

  // ice nucleation
  KOKKOS_FUNCTION
  static void ice_nucleation(const Pack &temp, const Pack &inv_rho, const Pack &ni,
                             const Pack &ni_activated, const Pack &qv_supersat_i,
                             const Scalar &inv_dt, const bool &do_predict_nc,
                             const bool &do_prescribed_CCN, Pack &qv2qi_nucleat_tend,
                             Pack &ni_nucleat_tend, const P3Runtime &runtime_options,
                             const Mask &context = Mask(true));

  KOKKOS_FUNCTION
  static Pack subgrid_variance_scaling(const Pack &relvar, const Scalar &expon);

  KOKKOS_FUNCTION
  static void ice_cldliq_wet_growth(
      const Pack &rho, const Pack &temp, const Pack &pres, const Pack &rhofaci,
      const Pack &table_val_qi2qr_melting, const Pack &table_val_qi2qr_vent_melt, const Pack &dv,
      const Pack &kap, const Pack &mu, const Pack &sc, const Pack &qv, const Pack &qc_incld,
      const Pack &qi_incld, const Pack &ni_incld, const Pack &qr_incld, Mask &log_wetgrowth,
      Pack &qr2qi_collect_tend, Pack &qc2qi_collect_tend, Pack &qc_growth_rate,
      Pack &nr_ice_shed_tend, Pack &qc2qr_ice_shed_tend, const Mask &context = Mask(true));

  KOKKOS_FUNCTION
  static void check_values(const uview_1d<const Pack> &qv, const uview_1d<const Pack> &temp,
                           const Int &ktop, const Int &kbot, const Int &timestepcount,
                           const bool &force_abort, const Int &source_ind, const MemberType &team,
                           const uview_1d<const Scalar> &col_loc);

#ifdef SCREAM_P3_SMALL_KERNELS
  static void check_values_disp(const uview_2d<const Pack> &qv, const uview_2d<const Pack> &temp,
                                const Int &ktop, const Int &kbot, const Int &timestepcount,
                                const bool &force_abort, const Int &source_ind,
                                const uview_2d<const Scalar> &col_loc, const Int &nj,
                                const Int &nk);
#endif

  KOKKOS_FUNCTION
  static void calculate_incloud_mixingratios(
      const Pack &qc, const Pack &qr, const Pack &qi, const Pack &qm, const Pack &nc,
      const Pack &nr, const Pack &ni, const Pack &bm, const Pack &inv_cld_frac_l,
      const Pack &inv_cld_frac_i, const Pack &inv_cld_frac_r, Pack &qc_incld, Pack &qr_incld,
      Pack &qi_incld, Pack &qm_incld, Pack &nc_incld, Pack &nr_incld, Pack &ni_incld,
      Pack &bm_incld, const Mask &context = Mask(true));

  //
  // main P3 functions
  //

  KOKKOS_FUNCTION
  static void
  p3_main_init(const MemberType &team, const Int &nk_pack, const uview_1d<const Pack> &cld_frac_i,
               const uview_1d<const Pack> &cld_frac_l, const uview_1d<const Pack> &cld_frac_r,
               const uview_1d<const Pack> &inv_exner, const uview_1d<const Pack> &th_atm,
               const uview_1d<const Pack> &dz, const uview_1d<Pack> &diag_equiv_reflectivity,
               const uview_1d<Pack> &ze_ice, const uview_1d<Pack> &ze_rain,
               const uview_1d<Pack> &diag_eff_radius_qc, const uview_1d<Pack> &diag_eff_radius_qi,
               const uview_1d<Pack> &diag_eff_radius_qr, const uview_1d<Pack> &inv_cld_frac_i,
               const uview_1d<Pack> &inv_cld_frac_l, const uview_1d<Pack> &inv_cld_frac_r,
               const uview_1d<Pack> &exner, const uview_1d<Pack> &T_atm,
               const uview_1d<Pack> &qv, const uview_1d<Pack> &inv_dz, Scalar &precip_liq_surf,
               Scalar &precip_ice_surf, view_1d_ptr_array<Pack, 36> &zero_init);

#ifdef SCREAM_P3_SMALL_KERNELS
  static void p3_main_init_disp(
      const Int &nj, const Int &nk_pack, const uview_2d<const Pack> &cld_frac_i,
      const uview_2d<const Pack> &cld_frac_l, const uview_2d<const Pack> &cld_frac_r,
      const uview_2d<const Pack> &inv_exner, const uview_2d<const Pack> &th_atm,
      const uview_2d<const Pack> &dz, const uview_2d<Pack> &diag_equiv_reflectivity,
      const uview_2d<Pack> &ze_ice, const uview_2d<Pack> &ze_rain,
      const uview_2d<Pack> &diag_eff_radius_qc, const uview_2d<Pack> &diag_eff_radius_qi,
      const uview_2d<Pack> &diag_eff_radius_qr, const uview_2d<Pack> &inv_cld_frac_i,
      const uview_2d<Pack> &inv_cld_frac_l, const uview_2d<Pack> &inv_cld_frac_r,
      const uview_2d<Pack> &exner, const uview_2d<Pack> &T_atm, const uview_2d<Pack> &qv,
      const uview_2d<Pack> &inv_dz, const uview_1d<Scalar> &precip_liq_surf,
      const uview_1d<Scalar> &precip_ice_surf, const uview_2d<Pack> &mu_r,
      const uview_2d<Pack> &lamr, const uview_2d<Pack> &logn0r, const uview_2d<Pack> &nu,
      const uview_2d<Pack> &cdist, const uview_2d<Pack> &cdist1, const uview_2d<Pack> &cdistr,
      const uview_2d<Pack> &qc_incld, const uview_2d<Pack> &qr_incld,
      const uview_2d<Pack> &qi_incld, const uview_2d<Pack> &qm_incld,
      const uview_2d<Pack> &nc_incld, const uview_2d<Pack> &nr_incld,
      const uview_2d<Pack> &ni_incld, const uview_2d<Pack> &bm_incld,
      const uview_2d<Pack> &inv_rho, const uview_2d<Pack> &prec, const uview_2d<Pack> &rho,
      const uview_2d<Pack> &rhofacr, const uview_2d<Pack> &rhofaci, const uview_2d<Pack> &acn,
      const uview_2d<Pack> &qv_sat_l, const uview_2d<Pack> &qv_sat_i, const uview_2d<Pack> &sup,
      const uview_2d<Pack> &qv_supersat_i, const uview_2d<Pack> &qtend_ignore,
      const uview_2d<Pack> &ntend_ignore, const uview_2d<Pack> &mu_c, const uview_2d<Pack> &lamc,
      const uview_2d<Pack> &rho_qi, const uview_2d<Pack> &qv2qi_depos_tend,
      const uview_2d<Pack> &precip_total_tend, const uview_2d<Pack> &nevapr,
      const uview_2d<Pack> &precip_liq_flux, const uview_2d<Pack> &precip_ice_flux);
#endif

  KOKKOS_FUNCTION
  static void p3_main_part1(
      const MemberType &team, const Int &nk, const bool &do_predict_nc,
      const bool &do_prescribed_CCN, const Scalar &dt, const uview_1d<const Pack> &pres,
      const uview_1d<const Pack> &dpres, const uview_1d<const Pack> &dz,
      const uview_1d<const Pack> &nc_nuceat_tend, const uview_1d<const Pack> &nccn_prescribed,
      const uview_1d<const Pack> &inv_exner, const uview_1d<const Pack> &exner,
      const uview_1d<const Pack> &inv_cld_frac_l, const uview_1d<const Pack> &inv_cld_frac_i,
      const uview_1d<const Pack> &inv_cld_frac_r, const uview_1d<Pack> &T_atm,
      const uview_1d<Pack> &rho, const uview_1d<Pack> &inv_rho, const uview_1d<Pack> &qv_sat_l,
      const uview_1d<Pack> &qv_sat_i, const uview_1d<Pack> &qv_supersat_i,
      const uview_1d<Pack> &rhofacr, const uview_1d<Pack> &rhofaci, const uview_1d<Pack> &acn,
      const uview_1d<Pack> &qv, const uview_1d<Pack> &th_atm, const uview_1d<Pack> &qc,
      const uview_1d<Pack> &nc, const uview_1d<Pack> &qr, const uview_1d<Pack> &nr,
      const uview_1d<Pack> &qi, const uview_1d<Pack> &ni, const uview_1d<Pack> &qm,
      const uview_1d<Pack> &bm, const uview_1d<Pack> &qc_incld, const uview_1d<Pack> &qr_incld,
      const uview_1d<Pack> &qi_incld, const uview_1d<Pack> &qm_incld,
      const uview_1d<Pack> &nc_incld, const uview_1d<Pack> &nr_incld,
      const uview_1d<Pack> &ni_incld, const uview_1d<Pack> &bm_incld, bool &is_nucleat_possible,
      bool &is_hydromet_present, const P3Runtime &runtime_options);

#ifdef SCREAM_P3_SMALL_KERNELS
  static void p3_main_part1_disp(
      const Int &nj, const Int &nk, const bool &do_predict_nc, const bool &do_prescribed_CCN,
      const Scalar &dt, const uview_2d<const Pack> &pres, const uview_2d<const Pack> &dpres,
      const uview_2d<const Pack> &dz, const uview_2d<const Pack> &nc_nuceat_tend,
      const uview_2d<const Pack> &nccn_prescribed, const uview_2d<const Pack> &inv_exner,
      const uview_2d<const Pack> &exner, const uview_2d<const Pack> &inv_cld_frac_l,
      const uview_2d<const Pack> &inv_cld_frac_i, const uview_2d<const Pack> &inv_cld_frac_r,
      const uview_2d<Pack> &T_atm, const uview_2d<Pack> &rho, const uview_2d<Pack> &inv_rho,
      const uview_2d<Pack> &qv_sat_l, const uview_2d<Pack> &qv_sat_i,
      const uview_2d<Pack> &qv_supersat_i, const uview_2d<Pack> &rhofacr,
      const uview_2d<Pack> &rhofaci, const uview_2d<Pack> &acn, const uview_2d<Pack> &qv,
      const uview_2d<Pack> &th_atm, const uview_2d<Pack> &qc, const uview_2d<Pack> &nc,
      const uview_2d<Pack> &qr, const uview_2d<Pack> &nr, const uview_2d<Pack> &qi,
      const uview_2d<Pack> &ni, const uview_2d<Pack> &qm, const uview_2d<Pack> &bm,
      const uview_2d<Pack> &qc_incld, const uview_2d<Pack> &qr_incld,
      const uview_2d<Pack> &qi_incld, const uview_2d<Pack> &qm_incld,
      const uview_2d<Pack> &nc_incld, const uview_2d<Pack> &nr_incld,
      const uview_2d<Pack> &ni_incld, const uview_2d<Pack> &bm_incld,
      const uview_1d<bool> &is_nucleat_possible, const uview_1d<bool> &is_hydromet_present,
      const P3Runtime &runtime_options);
#endif

  KOKKOS_FUNCTION
  static void p3_main_part2(
      const MemberType &team, const Int &nk_pack, const Scalar &max_total_ni,
      const bool &do_predict_nc, const bool &do_prescribed_CCN, const Scalar &dt,
      const Scalar &inv_dt, const uview_1d<const Pack> &ohetfrz_immersion_nucleation_tend,
      const uview_1d<const Pack> &ohetfrz_contact_nucleation_tend,
      const uview_1d<const Pack> &ohetfrz_deposition_nucleation_tend, const view_dnu_table &dnu,
      const view_ice_table &ice_table_vals, const view_collect_table &collect_table_vals,
      const view_2d_table &revap_table_vals, const uview_1d<const Pack> &pres,
      const uview_1d<const Pack> &dpres, const uview_1d<const Pack> &dz,
      const uview_1d<const Pack> &nc_nuceat_tend, const uview_1d<const Pack> &inv_exner,
      const uview_1d<const Pack> &exner, const uview_1d<const Pack> &inv_cld_frac_l,
      const uview_1d<const Pack> &inv_cld_frac_i, const uview_1d<const Pack> &inv_cld_frac_r,
      const uview_1d<const Pack> &ni_activated, const uview_1d<const Pack> &inv_qc_relvar,
      const uview_1d<const Pack> &cld_frac_i, const uview_1d<const Pack> &cld_frac_l,
      const uview_1d<const Pack> &cld_frac_r, const uview_1d<const Pack> &qv_prev,
      const uview_1d<const Pack> &t_prev, const uview_1d<Pack> &T_atm, const uview_1d<Pack> &rho,
      const uview_1d<Pack> &inv_rho, const uview_1d<Pack> &qv_sat_l,
      const uview_1d<Pack> &qv_sat_i, const uview_1d<Pack> &qv_supersat_i,
      const uview_1d<Pack> &rhofacr, const uview_1d<Pack> &rhofaci, const uview_1d<Pack> &acn,
      const uview_1d<Pack> &qv, const uview_1d<Pack> &th_atm, const uview_1d<Pack> &qc,
      const uview_1d<Pack> &nc, const uview_1d<Pack> &qr, const uview_1d<Pack> &nr,
      const uview_1d<Pack> &qi, const uview_1d<Pack> &ni, const uview_1d<Pack> &qm,
      const uview_1d<Pack> &bm, const uview_1d<Pack> &qc_incld, const uview_1d<Pack> &qr_incld,
      const uview_1d<Pack> &qi_incld, const uview_1d<Pack> &qm_incld,
      const uview_1d<Pack> &nc_incld, const uview_1d<Pack> &nr_incld,
      const uview_1d<Pack> &ni_incld, const uview_1d<Pack> &bm_incld, const uview_1d<Pack> &mu_c,
      const uview_1d<Pack> &nu, const uview_1d<Pack> &lamc, const uview_1d<Pack> &cdist,
      const uview_1d<Pack> &cdist1, const uview_1d<Pack> &cdistr, const uview_1d<Pack> &mu_r,
      const uview_1d<Pack> &lamr, const uview_1d<Pack> &logn0r,
      const uview_1d<Pack> &qv2qi_depos_tend, const uview_1d<Pack> &precip_total_tend,
      const uview_1d<Pack> &nevapr, const uview_1d<Pack> &qr_evap_tend,
      const uview_1d<Pack> &vap_liq_exchange, const uview_1d<Pack> &vap_ice_exchange,
      const uview_1d<Pack> &liq_ice_exchange, const uview_1d<Pack> &qr2qv_evap,
      const uview_1d<Pack> &qi2qv_sublim, const uview_1d<Pack> &qc2qr_accret,
      const uview_1d<Pack> &qc2qr_autoconv, const uview_1d<Pack> &qv2qi_vapdep,
      const uview_1d<Pack> &qc2qi_berg, const uview_1d<Pack> &qc2qr_ice_shed,
      const uview_1d<Pack> &qc2qi_collect, const uview_1d<Pack> &qr2qi_collect,
      const uview_1d<Pack> &qc2qi_hetero_freeze, const uview_1d<Pack> &qr2qi_immers_freeze,
      const uview_1d<Pack> &qi2qr_melt, const uview_1d<Pack> &pratot,
      const uview_1d<Pack> &prctot, bool &is_hydromet_present, const Int &nk,
      const P3Runtime &runtime_options);

#ifdef SCREAM_P3_SMALL_KERNELS
  static void p3_main_part2_disp(
      const Int &nj, const Int &nk, const Scalar &max_total_ni, const bool &do_predict_nc,
      const bool &do_prescribed_CCN, const Scalar &dt, const Scalar &inv_dt,
      const uview_2d<const Pack> &hetfrz_immersion_nucleation_tend,
      const uview_2d<const Pack> &hetfrz_contact_nucleation_tend,
      const uview_2d<const Pack> &hetfrz_deposition_nucleation_tend, const view_dnu_table &dnu,
      const view_ice_table &ice_table_vals, const view_collect_table &collect_table_vals,
      const view_2d_table &revap_table_vals, const uview_2d<const Pack> &pres,
      const uview_2d<const Pack> &dpres, const uview_2d<const Pack> &dz,
      const uview_2d<const Pack> &nc_nuceat_tend, const uview_2d<const Pack> &inv_exner,
      const uview_2d<const Pack> &exner, const uview_2d<const Pack> &inv_cld_frac_l,
      const uview_2d<const Pack> &inv_cld_frac_i, const uview_2d<const Pack> &inv_cld_frac_r,
      const uview_2d<const Pack> &ni_activated, const uview_2d<const Pack> &inv_qc_relvar,
      const uview_2d<const Pack> &cld_frac_i, const uview_2d<const Pack> &cld_frac_l,
      const uview_2d<const Pack> &cld_frac_r, const uview_2d<const Pack> &qv_prev,
      const uview_2d<const Pack> &t_prev, const uview_2d<Pack> &T_atm, const uview_2d<Pack> &rho,
      const uview_2d<Pack> &inv_rho, const uview_2d<Pack> &qv_sat_l,
      const uview_2d<Pack> &qv_sat_i, const uview_2d<Pack> &qv_supersat_i,
      const uview_2d<Pack> &rhofacr, const uview_2d<Pack> &rhofaci, const uview_2d<Pack> &acn,
      const uview_2d<Pack> &qv, const uview_2d<Pack> &th_atm, const uview_2d<Pack> &qc,
      const uview_2d<Pack> &nc, const uview_2d<Pack> &qr, const uview_2d<Pack> &nr,
      const uview_2d<Pack> &qi, const uview_2d<Pack> &ni, const uview_2d<Pack> &qm,
      const uview_2d<Pack> &bm, const uview_2d<Pack> &qc_incld, const uview_2d<Pack> &qr_incld,
      const uview_2d<Pack> &qi_incld, const uview_2d<Pack> &qm_incld,
      const uview_2d<Pack> &nc_incld, const uview_2d<Pack> &nr_incld,
      const uview_2d<Pack> &ni_incld, const uview_2d<Pack> &bm_incld, const uview_2d<Pack> &mu_c,
      const uview_2d<Pack> &nu, const uview_2d<Pack> &lamc, const uview_2d<Pack> &cdist,
      const uview_2d<Pack> &cdist1, const uview_2d<Pack> &cdistr, const uview_2d<Pack> &mu_r,
      const uview_2d<Pack> &lamr, const uview_2d<Pack> &logn0r,
      const uview_2d<Pack> &qv2qi_depos_tend, const uview_2d<Pack> &precip_total_tend,
      const uview_2d<Pack> &nevapr, const uview_2d<Pack> &qr_evap_tend,
      const uview_2d<Pack> &vap_liq_exchange, const uview_2d<Pack> &vap_ice_exchange,
      const uview_2d<Pack> &liq_ice_exchange, const uview_2d<Pack> &qr2qv_evap,
      const uview_2d<Pack> &qi2qv_sublim, const uview_2d<Pack> &qc2qr_accret,
      const uview_2d<Pack> &qc2qr_autoconv, const uview_2d<Pack> &qv2qi_vapdep,
      const uview_2d<Pack> &qc2qi_berg, const uview_2d<Pack> &qc2qr_ice_shed,
      const uview_2d<Pack> &qc2qi_collect, const uview_2d<Pack> &qr2qi_collect,
      const uview_2d<Pack> &qc2qi_hetero_freeze, const uview_2d<Pack> &qr2qi_immers_freeze,
      const uview_2d<Pack> &qi2qr_melt, const uview_2d<Pack> &pratot,
      const uview_2d<Pack> &prctot, const uview_1d<bool> &is_nucleat_possible,
      const uview_1d<bool> &is_hydromet_present, const P3Runtime &runtime_options);
#endif

  KOKKOS_FUNCTION
  static void
  p3_main_part3(const MemberType &team, const Int &nk_pack, const Scalar &max_total_ni,
                const view_dnu_table &dnu, const view_ice_table &ice_table_vals,
                const uview_1d<const Pack> &inv_exner, const uview_1d<const Pack> &cld_frac_l,
                const uview_1d<const Pack> &cld_frac_r, const uview_1d<const Pack> &cld_frac_i,
                const uview_1d<Pack> &rho, const uview_1d<Pack> &inv_rho,
                const uview_1d<Pack> &rhofaci, const uview_1d<Pack> &qv,
                const uview_1d<Pack> &th_atm, const uview_1d<Pack> &qc, const uview_1d<Pack> &nc,
                const uview_1d<Pack> &qr, const uview_1d<Pack> &nr, const uview_1d<Pack> &qi,
                const uview_1d<Pack> &ni, const uview_1d<Pack> &qm, const uview_1d<Pack> &bm,
                const uview_1d<Pack> &mu_c, const uview_1d<Pack> &nu, const uview_1d<Pack> &lamc,
                const uview_1d<Pack> &mu_r, const uview_1d<Pack> &lamr,
                const uview_1d<Pack> &vap_liq_exchange, const uview_1d<Pack> &ze_rain,
                const uview_1d<Pack> &ze_ice, const uview_1d<Pack> &diag_vm_qi,
                const uview_1d<Pack> &diag_eff_radius_qi, const uview_1d<Pack> &diag_diam_qi,
                const uview_1d<Pack> &rho_qi, const uview_1d<Pack> &diag_equiv_reflectivity,
                const uview_1d<Pack> &diag_eff_radius_qc,
                const uview_1d<Pack> &diag_eff_radius_qr, const P3Runtime &runtime_options);

#ifdef SCREAM_P3_SMALL_KERNELS
  static void p3_main_part3_disp(
      const Int &nj, const Int &nk_pack, const Scalar &max_total_ni, const view_dnu_table &dnu,
      const view_ice_table &ice_table_vals, const uview_2d<const Pack> &inv_exner,
      const uview_2d<const Pack> &cld_frac_l, const uview_2d<const Pack> &cld_frac_r,
      const uview_2d<const Pack> &cld_frac_i, const uview_2d<Pack> &rho,
      const uview_2d<Pack> &inv_rho, const uview_2d<Pack> &rhofaci, const uview_2d<Pack> &qv,
      const uview_2d<Pack> &th_atm, const uview_2d<Pack> &qc, const uview_2d<Pack> &nc,
      const uview_2d<Pack> &qr, const uview_2d<Pack> &nr, const uview_2d<Pack> &qi,
      const uview_2d<Pack> &ni, const uview_2d<Pack> &qm, const uview_2d<Pack> &bm,
      const uview_2d<Pack> &mu_c, const uview_2d<Pack> &nu, const uview_2d<Pack> &lamc,
      const uview_2d<Pack> &mu_r, const uview_2d<Pack> &lamr,
      const uview_2d<Pack> &vap_liq_exchange, const uview_2d<Pack> &ze_rain,
      const uview_2d<Pack> &ze_ice, const uview_2d<Pack> &diag_vm_qi,
      const uview_2d<Pack> &diag_eff_radius_qi, const uview_2d<Pack> &diag_diam_qi,
      const uview_2d<Pack> &rho_qi, const uview_2d<Pack> &diag_equiv_reflectivity,
      const uview_2d<Pack> &diag_eff_radius_qc, const uview_2d<Pack> &diag_eff_radius_qr,
      const uview_1d<bool> &is_nucleat_possible, const uview_1d<bool> &is_hydromet_present,
      const P3Runtime &runtime_options);
#endif

  // Return microseconds elapsed
  static Int p3_main(const P3Runtime &runtime_options, const P3PrognosticState &prognostic_state,
                     const P3DiagnosticInputs &diagnostic_inputs,
                     const P3DiagnosticOutputs &diagnostic_outputs,
                     const P3Infrastructure &infrastructure, const P3HistoryOnly &history_only,
                     const P3LookupTables &lookup_tables,
#ifdef SCREAM_P3_SMALL_KERNELS
                     const P3Temporaries &temporaries,
#endif
                     const WorkspaceManager &workspace_mgr,
                     Int nj,  // number of columns
                     Int nk); // number of vertical cells per column

  static Int
  p3_main_internal(const P3Runtime &runtime_options, const P3PrognosticState &prognostic_state,
                   const P3DiagnosticInputs &diagnostic_inputs,
                   const P3DiagnosticOutputs &diagnostic_outputs,
                   const P3Infrastructure &infrastructure, const P3HistoryOnly &history_only,
                   const P3LookupTables &lookup_tables, const WorkspaceManager &workspace_mgr,
                   Int nj,  // number of columns
                   Int nk); // number of vertical cells per column

#ifdef SCREAM_P3_SMALL_KERNELS
  static Int
  p3_main_internal_disp(const P3Runtime &runtime_options, const P3PrognosticState &prognostic_state,
                        const P3DiagnosticInputs &diagnostic_inputs,
                        const P3DiagnosticOutputs &diagnostic_outputs,
                        const P3Infrastructure &infrastructure, const P3HistoryOnly &history_only,
                        const P3LookupTables &lookup_tables, const P3Temporaries &temporaries,
                        const WorkspaceManager &workspace_mgr,
                        Int nj,  // number of columns
                        Int nk); // number of vertical cells per column
#endif

  KOKKOS_FUNCTION
  static void ice_supersat_conservation(Pack &qidep, Pack &qinuc, Pack &qinuc_cnt,
                                        const Pack &cld_frac_i, const Pack &qv,
                                        const Pack &qv_sat_i, const Pack &t_atm, const Real &dt,
                                        const Pack &qi2qv_sublim_tend,
                                        const Pack &qr2qv_evap_tend,
                                        const bool &use_hetfrz_classnuc,
                                        const Mask &context = Mask(true));

  KOKKOS_FUNCTION
  static void nc_conservation(const Pack &nc, const Pack &nc_selfcollect_tend, const Real &dt,
                              Pack &nc_collect_tend, Pack &nc2ni_immers_freeze_tend,
                              Pack &nc_accret_tend, Pack &nc2nr_autoconv_tend, Pack &ncheti_cnt,
                              Pack &nicnt, const bool &use_hetfrz_classnuc,
                              const Mask &context = Mask(true));

  KOKKOS_FUNCTION
  static void nr_conservation(const Pack &nr, const Pack &ni2nr_melt_tend,
                              const Pack &nr_ice_shed_tend, const Pack &ncshdc,
                              const Pack &nc2nr_autoconv_tend, const Real &dt,
                              const Real &nmltratio, Pack &nr_collect_tend,
//[shanyp 20260220
//			      Pack &nr2ni_immers_freeze_tend, Pack &nr_selfcollect_tend,
                              Pack &nr2ni_immers_freeze_tend, Pack &nr_selfcollect_tend, Pack &nr_breakup_tend,
//shanyp 20260220]
                              Pack &nr_evap_tend, const Mask &context = Mask(true));

  KOKKOS_FUNCTION
  static void ni_conservation(const Pack &ni, const Pack &ni_nucleat_tend,
                              const Pack &nr2ni_immers_freeze_tend,
                              const Pack &nc2ni_immers_freeze_tend, const Pack &ncheti_cnt,
                              const Pack &nicnt, const Pack &ninuc_cnt, const Real &dt,
                              Pack &ni2nr_melt_tend, Pack &ni_sublim_tend,
                              Pack &ni_selfcollect_tend, const bool &use_hetfrz_classnuc,
                              const Mask &context = Mask(true));

  KOKKOS_FUNCTION
  static void prevent_liq_supersaturation(const Pack &pres, const Pack &t_atm, const Pack &qv,
                                          const Scalar &dt, const Pack &qidep, const Pack &qinuc,
                                          Pack &qi2qv_sublim_tend, Pack &qr2qv_evap_tend,
                                          const Mask &context = Mask(true));
}; // struct Functions

template <typename ScalarT, typename DeviceT>
constexpr ScalarT Functions<ScalarT, DeviceT>::P3C::lookup_table_1a_dum1_c;

} // namespace p3
} // namespace scream

// If a GPU build, without relocatable device code enabled, make all code available
// to the translation unit; otherwise, ETI is used.
#if defined(EAMXX_ENABLE_GPU) && !defined(KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE) && \
    !defined(KOKKOS_ENABLE_HIP_RELOCATABLE_DEVICE_CODE)
#include "p3_autoconversion_impl.hpp"
#include "p3_back_to_cell_average_impl.hpp"
#include "p3_calc_liq_relaxation_timescale_impl.hpp"
#include "p3_calc_rime_density_impl.hpp"
#include "p3_check_values_impl.hpp"
#include "p3_cldliq_imm_freezing_impl.hpp"
#include "p3_cloud_rain_acc_impl.hpp"
#include "p3_cloud_sed_impl.hpp"
#include "p3_conservation_impl.hpp"
#include "p3_droplet_self_coll_impl.hpp"
#include "p3_dsd2_impl.hpp"
#include "p3_evaporate_rain_impl.hpp"
#include "p3_find_impl.hpp"
#include "p3_get_time_space_phys_variables_impl.hpp"
#include "p3_ice_classical_nucleation_impl.hpp"
#include "p3_ice_cldliq_wet_growth_impl.hpp"
#include "p3_ice_collection_impl.hpp"
#include "p3_ice_deposition_sublimation_impl.hpp"
#include "p3_ice_melting_impl.hpp"
#include "p3_ice_nucleation_impl.hpp"
#include "p3_ice_relaxation_timescale_impl.hpp"
#include "p3_ice_sed_impl.hpp"
#include "p3_ice_supersat_conservation_impl.hpp"
#include "p3_impose_max_total_ni_impl.hpp"
#include "p3_incloud_mixingratios_impl.hpp"
#include "p3_init_impl.hpp"
#include "p3_main_impl.hpp"
#include "p3_main_impl_part1.hpp"
#include "p3_main_impl_part2.hpp"
#include "p3_main_impl_part3.hpp"
#include "p3_nc_conservation_impl.hpp"
#include "p3_ni_conservation_impl.hpp"
#include "p3_nr_conservation_impl.hpp"
#include "p3_prevent_liq_supersaturation_impl.hpp"
#include "p3_rain_imm_freezing_impl.hpp"
#include "p3_rain_sed_impl.hpp"
#include "p3_rain_self_collection_impl.hpp"
#include "p3_subgrid_variance_scaling_impl.hpp"
#include "p3_table3_impl.hpp"
#include "p3_table_ice_impl.hpp"
#include "p3_update_prognostics_impl.hpp"
#include "p3_upwind_impl.hpp"
#endif // GPU && !KOKKOS_ENABLE_*_RELOCATABLE_DEVICE_CODE
#endif // P3_FUNCTIONS_HPP
