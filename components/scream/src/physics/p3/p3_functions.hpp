#ifndef P3_FUNCTIONS_HPP
#define P3_FUNCTIONS_HPP

#include "physics/share/physics_constants.hpp"

#include "share/scream_types.hpp"

#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/ekat_workspace.hpp"

namespace scream {
namespace p3 {

/*
 * Functions is a stateless struct used to encapsulate a
 * number of functions for p3. We use the ETI pattern for
 * these functions.
 *
 * P3 assumptions:
 *  - Kokkos team policies have a vector length of 1
 */

template <typename ScalarT, typename DeviceT>
struct Functions
{
  //
  // ---------- P3 constants ---------
  //
  struct P3C {
    // Constants for ice lookup tables
    enum {
      densize     = 5,
      rimsize     = 4,
      isize       = 50,
      tabsize     = 12, // number of quantities used from lookup table
      rcollsize   = 30,
      coltabsize  = 2,  // number of ice-rain collection  quantities used from lookup table

      // switch for warm-rain parameterization
      // 1 => Seifert and Beheng 2001
      // 2 => Beheng 1994
      // 3 => Khairoutdinov and Kogan 2000
      iparam      = 3,
      dnusize     = 16,
    };

    static constexpr ScalarT lookup_table_1a_dum1_c =  4.135985029041767e+00; // 1.0/(0.1*log10(261.7))
    static constexpr const char* p3_lookup_base = "./data/p3_lookup_table_1.dat-v";
    static constexpr const char* p3_version = "4"; // TODO: Change this so that the table version and table path is a runtime option.
  };

  //
  // ------- Types --------
  //

  using Scalar = ScalarT;
  using Device = DeviceT;

  template <typename S>
  using BigPack = ekat::Pack<S,SCREAM_PACK_SIZE>;
  template <typename S>
  using SmallPack = ekat::Pack<S,SCREAM_SMALL_PACK_SIZE>;

  using IntSmallPack = SmallPack<Int>;
  using Pack = BigPack<Scalar>;
  using Spack = SmallPack<Scalar>;

  using Mask = ekat::Mask<BigPack<Scalar>::n>;
  using Smask = ekat::Mask<SmallPack<Scalar>::n>;

  using KT = KokkosTypes<Device>;

  using C = scream::physics::Constants<Scalar>;

  template <typename S>
  using view_1d = typename KT::template view_1d<S>;
  template <typename S>
  using view_2d = typename KT::template view_2d<S>;

  template <typename S, int N>
  using view_1d_ptr_array = typename KT::template view_1d_ptr_carray<S, N>;

  template <typename S>
  using uview_1d = typename ekat::template Unmanaged<view_1d<S> >;
  template <typename S>
  using uview_2d = typename ekat::template Unmanaged<view_2d<S> >;

  using MemberType = typename KT::MemberType;

  using Workspace = typename ekat::WorkspaceManager<Spack, Device>::Workspace;

  // This struct stores prognostic variables evolved by P3.
  struct P3PrognosticState {
    // Cloud mass mixing ratio [kg kg-1]
    view_2d<Spack> qc;
    // Cloud number mixing ratio [# kg-1]
    view_2d<Spack> nc;
    // Rain mass mixing ratio [kg kg-1]
    view_2d<Spack> qr;
    // Rain number mixing ratio [# kg-1]
    view_2d<Spack> nr;
    // Ice total mass mixing ratio [kg kg-1]
    view_2d<Spack> qi;
    // Ice rime mass mixing ratio [kg kg-1]
    view_2d<Spack> qm;
    // Ice total number mixing ratio [# kg-1]
    view_2d<Spack> ni;
    // Ice rime volume mixing ratio [m3 kg-1]
    view_2d<Spack> bm;
    // Water vapor mixing ratio [kg kg-1]
    view_2d<Spack> qv;
    // Potential temperature [K]
    view_2d<Spack> th;
  };

  // This struct stores diagnostic variables used by P3.
  struct P3DiagnosticInputs {
    // CCN activated number tendency [kg-1 s-1]
    view_2d<const Spack> nc_nuceat_tend;
    // Activated ice nuclei concentration [1/kg]
    view_2d<const Spack> ni_activated;
    // Assumed SGS 1/(var(qc)/mean(qc)) [kg2/kg2]
    view_2d<const Spack> inv_qc_relvar;
    // Ice cloud fraction
    view_2d<const Spack> cld_frac_i;
    // Liquid cloud fraction
    view_2d<const Spack> cld_frac_l;
    // Rain cloud fraction
    view_2d<const Spack> cld_frac_r;
    // Pressure [Pa]
    view_2d<const Spack> pres;
    // Vertical grid spacing [m]
    view_2d<const Spack> dz;
    // Pressure thickness [Pa]
    view_2d<const Spack> dpres;
    // Exner expression
    view_2d<const Spack> exner;
    // qv from previous step [kg/kg]
    view_2d<const Spack> qv_prev;
    // T from previous step [K]
    view_2d<const Spack> t_prev;
  };

  // This struct stores diagnostic outputs computed by P3.
  struct P3DiagnosticOutputs {
    // Size distribution shape parameter for radiation
    view_2d<Spack> mu_c;
    // Size distribution slope parameter for radiation
    view_2d<Spack> lamc;
    // qitend due to deposition/sublimation
    view_2d<Spack> cmeiout;
    // Precipitation rate, liquid [m s-1]
    view_1d<Scalar> precip_liq_surf;
    // Precipitation rate, solid [m s-1]
    view_1d<Scalar> precip_ice_surf;
    // Effective cloud radius [m]
    view_2d<Spack> diag_effc;
    // Effective ice radius [m]
    view_2d<Spack> diag_effi;
    // Bulk density of ice [kg m-3]
    view_2d<Spack> rho_qi;
    // Total precipitation (rain + snow)
    view_2d<Spack> precip_total_tend;
    // Evaporation of total precipitation (rain + snow)
    view_2d<Spack> nevapr;
    // Evaporation of rain
    view_2d<Spack> qr_evap_tend;
    // Grid-box average rain flux [kg m^-2 s^-1] pverp
    view_2d<Spack> precip_liq_flux;
    // Grid-box average ice/snow flux [kg m^-2 s^-1] pverp
    view_2d<Spack> precip_ice_flux;
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
    // Coordinates of columns, nj x 3
    view_2d<const Scalar> col_location;
  };

  // This struct stores tendencies computed by P3 and used by other
  // parameterizations.
  struct P3HistoryOnly {
    // Sum of liq-ice phase change tendencies
    view_2d<Spack> liq_ice_exchange;
    // Sum of vap-liq phase change tendencies
    view_2d<Spack> vap_liq_exchange;
    // Sum of vap-ice phase change tendencies
    view_2d<Spack> vap_ice_exchange;
  };

  // -- Table3 --

  struct Table3 {
    IntSmallPack dumii, dumjj;
    Spack rdumii, rdumjj;
  };

  struct TableIce {
    IntSmallPack dumi, dumjj, dumii, dumzz;
    Spack dum1, dum4, dum5, dum6;
  };

  struct TableRain {
    IntSmallPack dumj;
    Spack dum3;
  };

  // lookup table values for rain shape parameter mu_r
  using view_1d_table = typename KT::template view_1d_table<Scalar, C::MU_R_TABLE_DIM>;

  // lookup table values for rain number- and mass-weighted fallspeeds and ventilation parameters
  using view_2d_table = typename KT::template view_2d_table<Scalar, C::VTABLE_DIM0, C::VTABLE_DIM1>;

  // ice lookup table values
  using view_itab_table    = typename KT::template view<const Scalar[P3C::densize][P3C::rimsize][P3C::isize][P3C::tabsize]>;

  // ice lookup table values for ice-rain collision/collection
  using view_itabcol_table = typename KT::template view<const Scalar[P3C::densize][P3C::rimsize][P3C::isize][P3C::rcollsize][P3C::coltabsize]>;

  // droplet spectral shape parameter for mass spectra, used for Seifert and Beheng (2001)
  // warm rain autoconversion/accretion option only (iparam = 1)
  using view_dnu_table = typename KT::template view_1d_table<Scalar, P3C::dnusize>;

  //
  // --------- Functions ---------
  //

  // Call from host to initialize the static table entries.
  static void init_kokkos_tables(
    view_2d_table& vn_table, view_2d_table& vm_table, view_2d_table& revap_table,
    view_1d_table& mu_r_table, view_dnu_table& dnu);

  static void init_kokkos_ice_lookup_tables(
    view_itab_table& itab, view_itabcol_table& itabcol);

  // Map (mu_r, lamr) to Table3 data.
  KOKKOS_FUNCTION
  static void lookup(const Spack& mu_r, const Spack& lamr,
                     Table3& t,
                     const Smask& context = Smask(true) );

  // Converts quantities to cell averages
  KOKKOS_FUNCTION
  static void back_to_cell_average(const Spack& cld_frac_l, const Spack& cld_frac_r,
                                   const Spack& cld_frac_i, Spack& qc2qr_accret_tend, Spack& qr2qv_evap_tend,
                                   Spack& qc2qr_autoconv_tend, Spack& nc_accret_tend, Spack& nc_selfcollect_tend,
                                   Spack& nc2nr_autoconv_tend, Spack& nr_selfcollect_tend, Spack& nr_evap_tend,
                                   Spack& ncautr, 
                                   Spack& qi2qv_sublim_tend, Spack& nr_ice_shed_tend, Spack& qc2qi_hetero_freeze_tend,
                                   Spack& qr2qi_collect_tend, Spack& qc2qr_ice_shed_tend, Spack& qi2qr_melt_tend,
                                   Spack& qc2qi_collect_tend, Spack& qr2qi_immers_freeze_tend, Spack& ni2nr_melt_tend,
                                   Spack& nc_collect_tend, Spack& ncshdc, Spack& nc2ni_immers_freeze_tend,
                                   Spack& nr_collect_tend, Spack& ni_selfcollect_tend, Spack& qv2qi_vapdep_tend,
                                   Spack& nr2ni_immers_freeze_tend, Spack& ni_sublim_tend, Spack& qv2qi_nucleat_tend,
                                   Spack& ni_nucleat_tend, Spack& qc2qi_berg_tend,
                                   const Smask& context = Smask(true) );

  // Limits ice process rates to prevent overdepletion of sources such that
  // the subsequent adjustments are done with maximum possible rates for the
  // time step.
  KOKKOS_FUNCTION
  static void prevent_ice_overdepletion(
    const Spack& pres, const Spack& t, const Spack& qv, const Spack& latent_heat_sublim, const Scalar& inv_dt,
    Spack& qv2qi_vapdep_tend, Spack& qi2qv_sublim_tend, const Smask& range_mask = Smask(true),
    const Smask& context = Smask(true) );

  //------------------------------------------------------------------------------------------!
  // Finds indices in 3D ice (only) lookup table
  // ------------------------------------------------------------------------------------------!
  KOKKOS_FUNCTION
  static void lookup_ice(const Spack& qi, const Spack& ni,
                         const Spack& qm, const Spack& rhop, TableIce& t,
                         const Smask& context = Smask(true) );

  //------------------------------------------------------------------------------------------!
  // Finds indices in 3D rain lookup table
  //------------------------------------------------------------------------------------------!
  KOKKOS_FUNCTION
  static void lookup_rain(const Spack& qr, const Spack& nr, TableRain& t,
                          const Smask& context = Smask(true) );

  // Apply Table3 data to the table to return a value. This performs bilinear
  // interpolation within the quad given by {t.dumii, t.dumjj} x {t.dumii+1,
  // t.dumjj+1}.
  KOKKOS_FUNCTION
  static Spack apply_table(const view_2d_table& table,
                           const Table3& t);

  // Apply TableIce data to the ice tables to return a value.
  KOKKOS_FUNCTION
  static Spack apply_table_ice(const int& index, const view_itab_table& itab,
                               const TableIce& t,
                               const Smask& context = Smask(true) );

  // Interpolates lookup table values for rain/ice collection processes
  KOKKOS_FUNCTION
  static Spack apply_table_coll(const int& index, const view_itabcol_table& itabcoll,
                                const TableIce& ti, const TableRain& tr,
                                const Smask& context = Smask(true) );

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
  KOKKOS_FUNCTION
  static void calc_first_order_upwind_step(
    const uview_1d<const Spack>& rho,
    const uview_1d<const Spack>& inv_rho, // 1/rho
    const uview_1d<const Spack>& inv_dz,
    const MemberType& team,
    const Int& nk, const Int& k_bot, const Int& k_top, const Int& kdir, const Scalar& dt_sub,
    const view_1d_ptr_array<Spack, nfield>& flux, // workspace
    const view_1d_ptr_array<Spack, nfield>& V,    // (behaviorally const)
    const view_1d_ptr_array<Spack, nfield>& r);   // in/out

  // Evolve 1 mixing ratio. This is a syntax-convenience version of the above.
  KOKKOS_FUNCTION
  static void calc_first_order_upwind_step(
    const uview_1d<const Spack>& rho,
    const uview_1d<const Spack>& inv_rho, // 1/rho
    const uview_1d<const Spack>& inv_dz,
    const MemberType& team,
    const Int& nk, const Int& k_bot, const Int& k_top, const Int& kdir, const Scalar& dt_sub,
    const uview_1d<Spack>& flux,
    const uview_1d<const Spack>& V,
    const uview_1d<Spack>& r);

  // This is the main routine. It can be called by the user if kdir is known at
  // compile time. So far it is not, so the above versions are called instead.
  template <Int kdir, int nfield>
  KOKKOS_FUNCTION
  static void calc_first_order_upwind_step(
    const uview_1d<const Spack>& rho,
    const uview_1d<const Spack>& inv_rho,
    const uview_1d<const Spack>& inv_dz,
    const MemberType& team,
    const Int& nk, const Int& k_bot, const Int& k_top, const Scalar& dt_sub,
    const view_1d_ptr_array<Spack, nfield>& flux,
    const view_1d_ptr_array<Spack, nfield>& V, // (behaviorally const)
    const view_1d_ptr_array<Spack, nfield>& r);

  template <int nfield>
  KOKKOS_FUNCTION
  static void generalized_sedimentation(
    const uview_1d<const Spack>& rho,
    const uview_1d<const Spack>& inv_rho,
    const uview_1d<const Spack>& inv_dz,
    const MemberType& team,
    const Int& nk, const Int& k_qxtop, Int& k_qxbot, const Int& kbot, const Int& kdir, const Scalar& Co_max, Scalar& dt_left, Scalar& prt_accum,
    const view_1d_ptr_array<Spack, nfield>& fluxes,
    const view_1d_ptr_array<Spack, nfield>& Vs, // (behaviorally const)
    const view_1d_ptr_array<Spack, nfield>& rs);

  // Cloud sedimentation
  KOKKOS_FUNCTION
  static void cloud_sedimentation(
    const uview_1d<Spack>& qc_incld,
    const uview_1d<const Spack>& rho,
    const uview_1d<const Spack>& inv_rho,
    const uview_1d<const Spack>& cld_frac_l,
    const uview_1d<const Spack>& acn,
    const uview_1d<const Spack>& inv_dz,
    const view_dnu_table& dnu,
    const MemberType& team,
    const Workspace& workspace,
    const Int& nk, const Int& ktop, const Int& kbot, const Int& kdir, const Scalar& dt, const Scalar& inv_dt,
    const bool& do_predict_nc,
    const uview_1d<Spack>& qc,
    const uview_1d<Spack>& nc,
    const uview_1d<Spack>& nc_incld,
    const uview_1d<Spack>& mu_c,
    const uview_1d<Spack>& lamc,
    const uview_1d<Spack>& qc_tend,
    const uview_1d<Spack>& nc_tend,
    Scalar& precip_liq_surf);

  // TODO: comment
  KOKKOS_FUNCTION
  static void rain_sedimentation(
    const uview_1d<const Spack>& rho,
    const uview_1d<const Spack>& inv_rho,
    const uview_1d<const Spack>& rhofacr,
    const uview_1d<const Spack>& cld_frac_r,
    const uview_1d<const Spack>& inv_dz,
    const uview_1d<Spack>& qr_incld,
    const MemberType& team,
    const Workspace& workspace,
    const view_2d_table& vn_table, const view_2d_table& vm_table,
    const Int& nk, const Int& ktop, const Int& kbot, const Int& kdir, const Scalar& dt, const Scalar& inv_dt,
    const uview_1d<Spack>& qr,
    const uview_1d<Spack>& nr,
    const uview_1d<Spack>& nr_incld,
    const uview_1d<Spack>& mu_r,
    const uview_1d<Spack>& lamr,
    const uview_1d<Spack>& precip_liq_flux,
    const uview_1d<Spack>& qr_tend,
    const uview_1d<Spack>& nr_tend,
    Scalar& precip_liq_surf);

  // TODO: comment
  KOKKOS_FUNCTION
  static void ice_sedimentation(
    const uview_1d<const Spack>& rho,
    const uview_1d<const Spack>& inv_rho,
    const uview_1d<const Spack>& rhofaci,
    const uview_1d<const Spack>& cld_frac_i,
    const uview_1d<const Spack>& inv_dz,
    const MemberType& team,
    const Workspace& workspace,
    const Int& nk, const Int& ktop, const Int& kbot, const Int& kdir, const Scalar& dt, const Scalar& inv_dt,
    const uview_1d<Spack>& qi,
    const uview_1d<Spack>& qi_incld,
    const uview_1d<Spack>& ni,
    const uview_1d<Spack>& ni_incld,
    const uview_1d<Spack>& qm,
    const uview_1d<Spack>& qm_incld,
    const uview_1d<Spack>& bm,
    const uview_1d<Spack>& bm_incld,
    const uview_1d<Spack>& qi_tend,
    const uview_1d<Spack>& ni_tend,
    const view_itab_table& itab,
    Scalar& precip_ice_surf);

  // homogeneous freezing of cloud and rain
  KOKKOS_FUNCTION
  static void homogeneous_freezing(
    const uview_1d<const Spack>& t,
    const uview_1d<const Spack>& exner,
    const uview_1d<const Spack>& latent_heat_fusion,
    const MemberType& team,
    const Int& nk, const Int& ktop, const Int& kbot, const Int& kdir,
    const uview_1d<Spack>& qc,
    const uview_1d<Spack>& nc,
    const uview_1d<Spack>& qr,
    const uview_1d<Spack>& nr,
    const uview_1d<Spack>& qi,
    const uview_1d<Spack>& ni,
    const uview_1d<Spack>& qm,
    const uview_1d<Spack>& bm,
    const uview_1d<Spack>& th);

  // -- Find layers

  // Find the bottom and top of the mixing ratio, e.g., qr. It's worth casing
  // these out in two ways: 1 thread/column vs many, and by kdir.
  KOKKOS_FUNCTION
  static Int find_bottom (
    const MemberType& team,
    const uview_1d<const Scalar>& v, const Scalar& small,
    const Int& kbot, const Int& ktop, const Int& kdir,
    bool& log_present);

  KOKKOS_FUNCTION
  static Int find_top (
    const MemberType& team,
    const uview_1d<const Scalar>& v, const Scalar& small,
    const Int& kbot, const Int& ktop, const Int& kdir,
    bool& log_present);

  KOKKOS_FUNCTION
  static void cloud_water_conservation(const Spack& qc, const Scalar dt,
    Spack& qc2qr_autoconv_tend, Spack& qc2qr_accret_tend, Spack &qc2qi_collect_tend, Spack& qc2qi_hetero_freeze_tend, 
    Spack& qc2qr_ice_shed_tend, Spack& qc2qi_berg_tend, Spack& qi2qv_sublim_tend, Spack& qv2qi_vapdep_tend,
    const Smask& context = Smask(true) );

  KOKKOS_FUNCTION
  static void rain_water_conservation(
    const Spack& qr, const Spack& qc2qr_autoconv_tend, const Spack& qc2qr_accret_tend, const Spack& qi2qr_melt_tend, const Spack& qc2qr_ice_shed_tend, const Scalar dt,
    Spack& qr2qv_evap_tend, Spack& qr2qi_collect_tend, Spack& qr2qi_immers_freeze_tend,
    const Smask& context = Smask(true) );

  KOKKOS_FUNCTION
  static void ice_water_conservation(
    const Spack& qi,const Spack& qv2qi_vapdep_tend,const Spack& qv2qi_nucleat_tend,const Spack& qc2qi_berg_tend, const Spack &qr2qi_collect_tend,
    const Spack &qc2qi_collect_tend,const Spack& qr2qi_immers_freeze_tend,const Spack& qc2qi_hetero_freeze_tend,const Scalar dt,
    Spack& qi2qv_sublim_tend, Spack& qi2qr_melt_tend,
    const Smask& context = Smask(true) );

  // TODO: comment
  KOKKOS_FUNCTION
  static void get_cloud_dsd2(
    const Spack& qc, Spack& nc, Spack& mu_c, const Spack& rho, Spack& nu,
    const view_dnu_table& dnu, Spack& lamc, Spack& cdist, Spack& cdist1, const Spack& cld_frac_l,
    const Smask& context = Smask(true) );

  // Computes and returns rain size distribution parameters
  KOKKOS_FUNCTION
  static void get_rain_dsd2 (
    const Spack& qr, Spack& nr, Spack& mu_r,
    Spack& lamr, Spack& cdistr, Spack& logn0r, const Spack& cld_frac_rm,
    const Smask& context = Smask(true) );

  // Calculates rime density
  KOKKOS_FUNCTION
  static void calc_rime_density(const Spack& t, const Spack& rhofaci,
    const Spack& table_val_qi_fallspd, const Spack& acn, const Spack& lamc,
    const Spack& mu_c, const Spack& qc_incld, const Spack& qc2qi_collect_tend,
    Spack& vtrmi1, Spack& rho_qm_cloud,
    const Smask& context = Smask(true) );

  // Computes contact and immersion freezing droplets
  KOKKOS_FUNCTION
  static void cldliq_immersion_freezing(const Spack& t, const Spack& lamc,
    const Spack& mu_c, const Spack& cdist1, const Spack& qc_incld, const Spack& inv_qc_relvar,
    Spack& qc2qi_hetero_freeze_tend, Spack& nc2ni_immers_freeze_tend,
    const Smask& context = Smask(true) );

  // Computes the immersion freezing of rain
  KOKKOS_FUNCTION
  static void rain_immersion_freezing(const Spack& t, const Spack& lamr,
    const Spack& mu_r, const Spack& cdistr, const Spack& qr_incld,
    Spack& qr2qi_immers_freeze_tend, Spack& nr2ni_immers_freeze_tend,
    const Smask& context = Smask(true) );

  // Computes droplet self collection
  KOKKOS_FUNCTION
  static void droplet_self_collection(const Spack& rho, const Spack& inv_rho,
    const Spack& qc_incld, const Spack& mu_c, const Spack& nu,
    const Spack& nc2nr_autoconv_tend, Spack& nc_selfcollect_tend,
    const Smask& context = Smask(true) );

  // Computes the accretion of clouds by rain
  KOKKOS_FUNCTION
  static void cloud_rain_accretion(const Spack& rho, const Spack& inv_rho,
    const Spack& qc_incld, const Spack& nc_incld, const Spack& qr_incld, const Spack& inv_qc_relvar,
    Spack& qc2qr_accret_tend, Spack& nc_accret_tend,
    const Smask& context = Smask(true) );

  // Computes cloud water autoconversion process rate
  KOKKOS_FUNCTION
  static void cloud_water_autoconversion(const Spack& rho,  const Spack& qc_incld,
    const Spack& nc_incld, const Spack& inv_qc_relvar,
    Spack& qc2qr_autoconv_tend, Spack& nc2nr_autoconv_tend, Spack& ncautr,
    const Smask& context = Smask(true) );

  // Computes rain self collection process rate
  KOKKOS_FUNCTION
  static void rain_self_collection(const Spack& rho, const Spack& qr_incld, const Spack& nr_incld, Spack& nr_selfcollect_tend,
                                   const Smask& context = Smask(true) );

  // Impose maximum ice number
  KOKKOS_FUNCTION
  static void impose_max_total_Ni(Spack& ni_local, const Scalar& max_total_Ni, const Spack& inv_rho_local,
                                  const Smask& context = Smask(true) );

  //--------------------------------------------------------------------------------
  //  Calculates and returns the bulk rime density from the prognostic ice variables
  //  and adjusts qm and bm appropriately.
  //--------------------------------------------------------------------------------
  KOKKOS_FUNCTION
  static Spack calc_bulk_rho_rime(
    const Spack& qi_tot, Spack& qi_rim, Spack& bi_rim,
    const Smask& context = Smask(true) );

  // TODO - comment
  KOKKOS_FUNCTION
  static void compute_rain_fall_velocity(
    const view_2d_table& vn_table, const view_2d_table& vm_table,
    const Spack& qr_incld, const Spack& cld_frac_r, const Spack& rhofacr, 
    Spack& nr_incld, Spack& mu_r, Spack& lamr, Spack& V_qr, Spack& V_nr,
    const Smask& context = Smask(true));

  //---------------------------------------------------------------------------------
  // update prognostic microphysics and thermodynamics variables
  //---------------------------------------------------------------------------------
  //-- ice-phase dependent processes:
  KOKKOS_FUNCTION
  static void update_prognostic_ice(
    const Spack& qc2qi_hetero_freeze_tend, const Spack& qc2qi_collect_tend,
    const Spack& qc2qr_ice_shed_tend,  const Spack& nc_collect_tend,  const Spack& nc2ni_immers_freeze_tend, const Spack& ncshdc,
    const Spack& qr2qi_collect_tend,  const Spack& nr_collect_tend,  const Spack& qr2qi_immers_freeze_tend, const Spack& nr2ni_immers_freeze_tend,
    const Spack& nr_ice_shed_tend, const Spack& qi2qr_melt_tend,  const Spack& ni2nr_melt_tend,  const Spack& qi2qv_sublim_tend,
    const Spack& qv2qi_vapdep_tend,  const Spack& qv2qi_nucleat_tend,  const Spack& ni_nucleat_tend,  const Spack& ni_selfcollect_tend,
    const Spack& ni_sublim_tend,  const Spack& qc2qi_berg_tend, const Spack& exner,  const Spack& latent_heat_sublim,
    const Spack& latent_heat_fusion,    const bool do_predict_nc, const Smask& log_wetgrowth, const Scalar dt,
    const Scalar& nmltratio, const Spack& rho_qm_cloud, Spack& th, Spack& qv, Spack& qi,
    Spack& ni, Spack& qm, Spack& bm, Spack& qc,  Spack& nc, Spack& qr, Spack& nr,
    const Smask& context = Smask(true));

  // TODO (comments)
  KOKKOS_FUNCTION
  static void get_time_space_phys_variables(const Spack& t, const Spack& pres, const Spack& rho,
					    const Spack& latent_heat_vapor, const Spack& latent_heat_sublim,
					    const Spack& qv_sat_l, const Spack& qv_sat_i, Spack& mu,
					    Spack& dv, Spack& sc, Spack& dqsdt, Spack& dqsidt,
					    Spack& ab, Spack& abi, Spack& kap, Spack& eii,
                                            const Smask& context = Smask(true));

  // TODO (comments)
  KOKKOS_FUNCTION
  static void ice_cldliq_collection(const Spack& rho, const Spack& temp,
                                    const Spack& rhofaci, const Spack& table_val_qc2qi_collect,
                                    const Spack& qi_incld, const Spack& qc_incld,
                                    const Spack& ni_incld, const Spack& nc_incld,
                                    Spack& qc2qi_collect_tend, Spack& nc_collect_tend, Spack& qc2qr_ice_shed_tend, Spack& ncshdc,
				    const Smask& range_mask = Smask(true),
                                    const Smask& context = Smask(true));

  // TODO (comments)
  KOKKOS_FUNCTION
  static void ice_rain_collection(const Spack& rho, const Spack& temp,
                                  const Spack& rhofaci, const Spack& logn0r,
                                  const Spack& table_val_nr_collect, const Spack& table_val_qr2qi_collect,
                                  const Spack& qi_incld, const Spack& ni_incld,
                                  const Spack& qr_incld,
                                  Spack& qr2qi_collect_tend, Spack& nr_collect_tend,
                                  const Smask& context = Smask(true));

  // TODO (comments)
  KOKKOS_FUNCTION
  static void ice_self_collection(const Spack& rho, const Spack& rhofaci,
                                  const Spack& table_val_ni_self_collect, const Spack& eii,
                                  const Spack& qm_incld, const Spack& qi_incld,
                                  const Spack& ni_incld, Spack& ni_selfcollect_tend,
                                  const Smask& context = Smask(true));

  // helper fn for evaporate_rain
  KOKKOS_FUNCTION
  static void rain_evap_tscale_weight(const Spack& dt_over_tau,
				      Spack& weight,
				      const Smask& context=Smask(true));

  // helper fn for evaporate_rain
  KOKKOS_FUNCTION
  static void rain_evap_equilib_tend(const Spack& A_c,const Spack& ab,
				     const Spack& tau_eff,const Spack& tau_r,
				     Spack& tend,const Smask& context=Smask(true));

  // helper fn for evaporate_rain
  KOKKOS_FUNCTION
  static void rain_evap_instant_tend(const Spack& ssat_r, const Spack& ab,
				     const Spack& tau_r,
				     Spack& tend, const Smask& context=Smask(true));
  
  // TODO (comments)
  KOKKOS_FUNCTION
  static void evaporate_rain(const Spack& qr_incld, const Spack& qc_incld, const Spack& nr_incld, const Spack& qi_incld,
			     const Spack& cld_frac_l, const Spack& cld_frac_r, const Spack& qv, const Spack& qv_prev,
			     const Spack& qv_sat_l, const Spack& qv_sat_i, const Spack& ab, const Spack& abi,
			     const Spack& epsr, const Spack & epsi_tot, const Spack& t, const Spack& t_prev,
			     const Spack& latent_heat_sublim, const Spack& dqsdt, const Scalar& dt,
			     Spack& qr2qv_evap_tend, Spack& nr_evap_tend,
			     const Smask& context = Smask(true));

  //get number and mass tendencies due to melting ice
  KOKKOS_FUNCTION
  static void ice_melting(const Spack& rho, const Spack& t, const Spack& pres, const Spack& rhofaci,
			  const Spack& table_val_qi2qr_melting, const Spack& table_val_qi2qr_vent_melt, const Spack& latent_heat_vapor, const Spack& latent_heat_fusion,
			  const Spack& dv, const Spack& sc, const Spack& mu, const Spack& kap,
			  const Spack& qv, const Spack& qi_incld, const Spack& ni_incld,
			  Spack& qi2qr_melt_tend, Spack& ni2nr_melt_tend, const Smask& range_mask = Smask(true),
                          const Smask& context = Smask(true));

  //liquid-phase dependent processes:
  KOKKOS_FUNCTION
  static void update_prognostic_liquid(const Spack& qc2qr_accret_tend, const Spack& nc_accret_tend,
    const Spack& qc2qr_autoconv_tend,const Spack& nc2nr_autoconv_tend, const Spack& ncautr,
    const Spack& nc_selfcollect_tend, const Spack& qr2qv_evap_tend, const Spack& nr_evap_tend, const Spack& nr_selfcollect_tend,
    const bool do_predict_nc, const Spack& inv_rho, const Spack& exner, const Spack& latent_heat_vapor,
    const Scalar dt, Spack& th, Spack& qv, Spack& qc, Spack& nc, Spack& qr, Spack& nr,
    const Smask& context = Smask(true));

  // TODO (comments)
  KOKKOS_FUNCTION
  static void ice_deposition_sublimation(const Spack& qi_incld,
    const Spack& ni_incld, const Spack& t, const Spack& qv_sat_l, const Spack& qv_sat_i,
    const Spack& epsi, const Spack& abi, const Spack& qv, Spack& qv2qi_vapdep_tend,
    Spack& qi2qv_sublim_tend, Spack& ni_sublim_tend, Spack& qc2qi_berg_tend,
    const Smask& context = Smask(true));

  KOKKOS_FUNCTION
  static void ice_relaxation_timescale(
    const Spack& rho, const Spack& temp, const Spack& rhofaci, const Spack& table_val_qi2qr_melting,
    const Spack& table_val_qi2qr_vent_melt, const Spack& dv, const Spack& mu, const Spack& sc,
    const Spack& qi_incld, const Spack& ni_incld,
    Spack& epsi, Spack& epsi_tot,
    const Smask& context = Smask(true));

  KOKKOS_FUNCTION
  static void calc_liq_relaxation_timescale(
    const view_2d_table& revap_table,
    const Spack& rho, const Scalar& f1r, const Scalar& f2r,
    const Spack& dv, const Spack& mu, const Spack& sc,
    const Spack& mu_r, const Spack& lamr, const Spack& cdistr,
    const Spack& cdist, const Spack& qr_incld, const Spack& qc_incld,
    Spack& epsr, Spack& epsc,
    const Smask& context = Smask(true));

  // ice nucleation
  KOKKOS_FUNCTION
  static void ice_nucleation(const Spack& temp, const Spack& inv_rho,
                             const Spack& ni, const Spack& ni_activated,
                             const Spack& qv_supersat_i, const Scalar& inv_dt,
                             const bool& do_predict_nc,
                             Spack& qv2qi_nucleat_tend, Spack& ni_nucleat_tend,
                             const Smask& context = Smask(true));

  KOKKOS_FUNCTION
  static Spack subgrid_variance_scaling(const Spack& relvar, const Scalar& expon);

  KOKKOS_FUNCTION
  static void ice_cldliq_wet_growth(const Spack& rho, const Spack& temp, const Spack& pres, const Spack& rhofaci, const Spack& table_val_qi2qr_melting,
                                    const Spack& table_val_qi2qr_vent_melt, const Spack& latent_heat_vapor, const Spack& latent_heat_fusion, const Spack& dv,
                                    const Spack& kap, const Spack& mu, const Spack& sc, const Spack& qv, const Spack& qc_incld,
                                    const Spack& qi_incld, const Spack& ni_incld, const Spack& qr_incld,
                                    Smask& log_wetgrowth, Spack& qr2qi_collect_tend, Spack& qc2qi_collect_tend, Spack& qc_growth_rate,
				    Spack& nr_ice_shed_tend, Spack& qc2qr_ice_shed_tend, const Smask& range_mask = Smask(true),
                                    const Smask& context = Smask(true));

  // Note: not a kernel function
  static void get_latent_heat(const Int& nj, const Int& nk, view_2d<Spack>& v, view_2d<Spack>& s, view_2d<Spack>& f);

  KOKKOS_FUNCTION
  static void check_values(const uview_1d<const Spack>& qv, const uview_1d<const Spack>& temp, const Int& ktop, const Int& kbot,
                           const Int& timestepcount, const bool& force_abort, const Int& source_ind, const MemberType& team,
                           const uview_1d<const Scalar>& col_loc);

  KOKKOS_FUNCTION
  static void calculate_incloud_mixingratios(
    const Spack& qc, const Spack& qr, const Spack& qi, const Spack& qm, const Spack& nc,
    const Spack& nr, const Spack& ni, const Spack& bm, const Spack& inv_cld_frac_l,
    const Spack& inv_cld_frac_i, const Spack& inv_cld_frac_r,
    Spack& qc_incld, Spack& qr_incld, Spack& qi_incld, Spack& qm_incld,
    Spack& nc_incld, Spack& nr_incld, Spack& ni_incld, Spack& bm_incld,
    const Smask& context = Smask(true));

  //
  // main P3 functions
  //

  KOKKOS_FUNCTION
  static void p3_main_init(
    const MemberType& team,
    const Int& nk_pack,
    const uview_1d<const Spack>& cld_frac_i,
    const uview_1d<const Spack>& cld_frac_l,
    const uview_1d<const Spack>& cld_frac_r,
    const uview_1d<const Spack>& exner,
    const uview_1d<const Spack>& th,
    const uview_1d<const Spack>& dz,
    const uview_1d<Spack>& diag_ze,
    const uview_1d<Spack>& ze_ice,
    const uview_1d<Spack>& ze_rain,
    const uview_1d<Spack>& diag_effc,
    const uview_1d<Spack>& diag_effi,
    const uview_1d<Spack>& inv_cld_frac_i,
    const uview_1d<Spack>& inv_cld_frac_l,
    const uview_1d<Spack>& inv_cld_frac_r,
    const uview_1d<Spack>& inv_exner,
    const uview_1d<Spack>& t,
    const uview_1d<Spack>& qv,
    const uview_1d<Spack>& inv_dz,
    Scalar& precip_liq_surf,
    Scalar& precip_ice_surf,
    view_1d_ptr_array<Spack, 36>& zero_init);

  KOKKOS_FUNCTION
  static void p3_main_part1(
    const MemberType& team,
    const Int& nk,
    const bool& do_predict_nc,
    const Scalar& dt,
    const uview_1d<const Spack>& pres,
    const uview_1d<const Spack>& dpres,
    const uview_1d<const Spack>& dz,
    const uview_1d<const Spack>& nc_nuceat_tend,
    const uview_1d<const Spack>& exner,
    const uview_1d<const Spack>& inv_exner,
    const uview_1d<const Spack>& inv_cld_frac_l,
    const uview_1d<const Spack>& inv_cld_frac_i,
    const uview_1d<const Spack>& inv_cld_frac_r,
    const uview_1d<const Spack>& latent_heat_vapor,
    const uview_1d<const Spack>& latent_heat_sublim,
    const uview_1d<const Spack>& latent_heat_fusion,
    const uview_1d<Spack>& t,
    const uview_1d<Spack>& rho,
    const uview_1d<Spack>& inv_rho,
    const uview_1d<Spack>& qv_sat_l,
    const uview_1d<Spack>& qv_sat_i,
    const uview_1d<Spack>& qv_supersat_i,
    const uview_1d<Spack>& rhofacr,
    const uview_1d<Spack>& rhofaci,
    const uview_1d<Spack>& acn,
    const uview_1d<Spack>& qv,
    const uview_1d<Spack>& th,
    const uview_1d<Spack>& qc,
    const uview_1d<Spack>& nc,
    const uview_1d<Spack>& qr,
    const uview_1d<Spack>& nr,
    const uview_1d<Spack>& qi,
    const uview_1d<Spack>& ni,
    const uview_1d<Spack>& qm,
    const uview_1d<Spack>& bm,
    const uview_1d<Spack>& qc_incld,
    const uview_1d<Spack>& qr_incld,
    const uview_1d<Spack>& qi_incld,
    const uview_1d<Spack>& qm_incld,
    const uview_1d<Spack>& nc_incld,
    const uview_1d<Spack>& nr_incld,
    const uview_1d<Spack>& ni_incld,
    const uview_1d<Spack>& bm_incld,
    bool& is_nucleat_possible,
    bool& is_hydromet_present);

  KOKKOS_FUNCTION
  static void p3_main_part2(
    const MemberType& team,
    const Int& nk_pack,
    const bool& do_predict_nc,
    const Scalar& dt,
    const Scalar& inv_dt,
    const view_dnu_table& dnu,
    const view_itab_table& itab,
    const view_itabcol_table& itabcol,
    const view_2d_table& revap_table,
    const uview_1d<const Spack>& pres,
    const uview_1d<const Spack>& dpres,
    const uview_1d<const Spack>& dz,
    const uview_1d<const Spack>& nc_nuceat_tend,
    const uview_1d<const Spack>& exner,
    const uview_1d<const Spack>& inv_exner,
    const uview_1d<const Spack>& inv_cld_frac_l,
    const uview_1d<const Spack>& inv_cld_frac_i,
    const uview_1d<const Spack>& inv_cld_frac_r,
    const uview_1d<const Spack>& ni_activated,
    const uview_1d<const Spack>& inv_qc_relvar,
    const uview_1d<const Spack>& cld_frac_i,
    const uview_1d<const Spack>& cld_frac_l,
    const uview_1d<const Spack>& cld_frac_r,
    const uview_1d<const Spack>& qv_prev,
    const uview_1d<const Spack>& t_prev,
    const uview_1d<Spack>& t,
    const uview_1d<Spack>& rho,
    const uview_1d<Spack>& inv_rho,
    const uview_1d<Spack>& qv_sat_l,
    const uview_1d<Spack>& qv_sat_i,
    const uview_1d<Spack>& qv_supersat_i,
    const uview_1d<Spack>& rhofacr,
    const uview_1d<Spack>& rhofaci,
    const uview_1d<Spack>& acn,
    const uview_1d<Spack>& qv,
    const uview_1d<Spack>& th,
    const uview_1d<Spack>& qc,
    const uview_1d<Spack>& nc,
    const uview_1d<Spack>& qr,
    const uview_1d<Spack>& nr,
    const uview_1d<Spack>& qi,
    const uview_1d<Spack>& ni,
    const uview_1d<Spack>& qm,
    const uview_1d<Spack>& bm,
    const uview_1d<Spack>& latent_heat_vapor,
    const uview_1d<Spack>& latent_heat_sublim,
    const uview_1d<Spack>& latent_heat_fusion,
    const uview_1d<Spack>& qc_incld,
    const uview_1d<Spack>& qr_incld,
    const uview_1d<Spack>& qi_incld,
    const uview_1d<Spack>& qm_incld,
    const uview_1d<Spack>& nc_incld,
    const uview_1d<Spack>& nr_incld,
    const uview_1d<Spack>& ni_incld,
    const uview_1d<Spack>& bm_incld,
    const uview_1d<Spack>& mu_c,
    const uview_1d<Spack>& nu,
    const uview_1d<Spack>& lamc,
    const uview_1d<Spack>& cdist,
    const uview_1d<Spack>& cdist1,
    const uview_1d<Spack>& cdistr,
    const uview_1d<Spack>& mu_r,
    const uview_1d<Spack>& lamr,
    const uview_1d<Spack>& logn0r,
    const uview_1d<Spack>& cmeiout,
    const uview_1d<Spack>& precip_total_tend,
    const uview_1d<Spack>& nevapr,
    const uview_1d<Spack>& qr_evap_tend,
    const uview_1d<Spack>& vap_liq_exchange,
    const uview_1d<Spack>& vap_ice_exchange,
    const uview_1d<Spack>& liq_ice_exchange,
    const uview_1d<Spack>& pratot,
    const uview_1d<Spack>& prctot,
    bool& is_hydromet_present,
    const Int& nk=-1);

  KOKKOS_FUNCTION
  static void p3_main_part3(
    const MemberType& team,
    const Int& nk_pack,
    const view_dnu_table& dnu,
    const view_itab_table& itab,
    const uview_1d<const Spack>& exner,
    const uview_1d<const Spack>& cld_frac_l,
    const uview_1d<const Spack>& cld_frac_r,
    const uview_1d<const Spack>& cld_frac_i,
    const uview_1d<Spack>& rho,
    const uview_1d<Spack>& inv_rho,
    const uview_1d<Spack>& rhofaci,
    const uview_1d<Spack>& qv,
    const uview_1d<Spack>& th,
    const uview_1d<Spack>& qc,
    const uview_1d<Spack>& nc,
    const uview_1d<Spack>& qr,
    const uview_1d<Spack>& nr,
    const uview_1d<Spack>& qi,
    const uview_1d<Spack>& ni,
    const uview_1d<Spack>& qm,
    const uview_1d<Spack>& bm,
    const uview_1d<Spack>& latent_heat_vapor,
    const uview_1d<Spack>& latent_heat_sublim,
    const uview_1d<Spack>& mu_c,
    const uview_1d<Spack>& nu,
    const uview_1d<Spack>& lamc,
    const uview_1d<Spack>& mu_r,
    const uview_1d<Spack>& lamr,
    const uview_1d<Spack>& vap_liq_exchange,
    const uview_1d<Spack>& ze_rain,
    const uview_1d<Spack>& ze_ice,
    const uview_1d<Spack>& diag_vmi,
    const uview_1d<Spack>& diag_effi,
    const uview_1d<Spack>& diag_di,
    const uview_1d<Spack>& rho_qi,
    const uview_1d<Spack>& diag_ze,
    const uview_1d<Spack>& diag_effc);

  static void p3_main(
    const P3PrognosticState& prognostic_state,
    const P3DiagnosticInputs& diagnostic_inputs,
    const P3DiagnosticOutputs& diagnostic_outputs,
    const P3Infrastructure& infrastructure,
    const P3HistoryOnly& history_only,
    Int nj, // number of columns
    Int nk); // number of vertical cells per column
};

template <typename ScalarT, typename DeviceT>
constexpr ScalarT Functions<ScalarT, DeviceT>::P3C::lookup_table_1a_dum1_c;

extern "C" {
// decl of fortran function for loading tables from fortran p3. This will
// continue to be a bit awkward until we have fully ported all of p3.
void init_tables_from_f90_c(Real* vn_table_data, Real* vm_table_data,
                            Real* revap_table_data, Real* mu_table_data);
}

} // namespace p3
} // namespace scream

// If a GPU build, make all code available to the translation unit; otherwise,
// ETI is used.
#ifdef KOKKOS_ENABLE_CUDA
# include "p3_table3_impl.hpp"
# include "p3_table_ice_impl.hpp"
# include "p3_back_to_cell_average_impl.hpp"
# include "p3_prevent_ice_overdepletion_impl.hpp"
# include "p3_dsd2_impl.hpp"
# include "p3_upwind_impl.hpp"
# include "p3_find_impl.hpp"
# include "p3_conservation_impl.hpp"
# include "p3_autoconversion_impl.hpp"
# include "p3_rain_self_collection_impl.hpp"
# include "p3_impose_max_total_Ni_impl.hpp"
# include "p3_calc_rime_density_impl.hpp"
# include "p3_cldliq_imm_freezing_impl.hpp"
# include "p3_droplet_self_coll_impl.hpp"
# include "p3_cloud_sed_impl.hpp"
# include "p3_cloud_rain_acc_impl.hpp"
# include "p3_ice_sed_impl.hpp"
# include "p3_rain_sed_impl.hpp"
# include "p3_rain_imm_freezing_impl.hpp"
# include "p3_get_time_space_phys_variables_impl.hpp"
# include "p3_evaporate_rain_impl.hpp"
# include "p3_update_prognostics_impl.hpp"
# include "p3_ice_collection_impl.hpp"
# include "p3_ice_deposition_sublimation_impl.hpp"
# include "p3_ice_relaxation_timescale_impl.hpp"
# include "p3_ice_nucleation_impl.hpp"
# include "p3_ice_melting_impl.hpp"
# include "p3_calc_liq_relaxation_timescale_impl.hpp"
# include "p3_ice_cldliq_wet_growth_impl.hpp"
# include "p3_get_latent_heat_impl.hpp"
# include "p3_check_values_impl.hpp"
# include "p3_incloud_mixingratios_impl.hpp"
# include "p3_subgrid_variance_scaling_impl.hpp"
# include "p3_main_impl.hpp"
#endif

#endif // P3_FUNCTIONS_HPP
