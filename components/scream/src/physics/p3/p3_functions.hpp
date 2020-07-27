#ifndef P3_FUNCTIONS_HPP
#define P3_FUNCTIONS_HPP

#include "ekat/scream_types.hpp"
#include "ekat/scream_pack_kokkos.hpp"
#include "ekat/scream_workspace.hpp"
#include "physics_constants.hpp"

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
  using BigPack = scream::pack::BigPack<S>;
  template <typename S>
  using SmallPack = scream::pack::SmallPack<S>;
  using IntSmallPack = scream::pack::IntSmallPack;

  using Pack = BigPack<Scalar>;
  using Spack = SmallPack<Scalar>;

  template <typename S>
  using Mask = scream::pack::Mask<BigPack<S>::n>;

  template <typename S>
  using SmallMask = scream::pack::Mask<SmallPack<S>::n>;

  using Smask = SmallMask<Scalar>;

  using KT = KokkosTypes<Device>;

  using C = scream::physics::Constants<Scalar>;

  template <typename S>
  using view_1d = typename KT::template view_1d<S>;
  template <typename S>
  using view_2d = typename KT::template view_2d<S>;

  template <typename S, int N>
  using view_1d_ptr_array = typename KT::template view_1d_ptr_carray<S, N>;

  template <typename S>
  using uview_1d = typename ko::template Unmanaged<view_1d<S> >;
  template <typename S>
  using uview_2d = typename ko::template Unmanaged<view_2d<S> >;

  using MemberType = typename KT::MemberType;

  using Workspace = typename WorkspaceManager<Spack, Device>::Workspace;

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
  static void back_to_cell_average(const Spack& lcldm, const Spack& rcldm,
                                   const Spack& icldm, Spack& qcacc, Spack& qrevp,
                                   Spack& qcaut, Spack& ncacc, Spack& ncslf,
                                   Spack& ncautc, Spack& nrslf, Spack& nrevp,
                                   Spack& ncautr, 
                                   Spack& qisub, Spack& nrshdr, Spack& qcheti,
                                   Spack& qrcol, Spack& qcshd, Spack& qimlt,
                                   Spack& qccol, Spack& qrheti, Spack& nimlt,
                                   Spack& nccol, Spack& ncshdc, Spack& ncheti,
                                   Spack& nrcol, Spack& nislf, Spack& qidep,
                                   Spack& nrheti, Spack& nisub, Spack& qinuc,
                                   Spack& ninuc, Spack& qiberg,
                                   const Smask& context = Smask(true) );

  // Limits ice process rates to prevent overdepletion of sources such that
  // the subsequent adjustments are done with maximum possible rates for the
  // time step.
  KOKKOS_FUNCTION
  static void prevent_ice_overdepletion(
    const Spack& pres, const Spack& t, const Spack& qv, const Spack& xxls, const Scalar& odt,
    Spack& qidep, Spack& qisub,
    const Smask& context = Smask(true) );

  //------------------------------------------------------------------------------------------!
  // Finds indices in 3D ice (only) lookup table
  // ------------------------------------------------------------------------------------------!
  KOKKOS_FUNCTION
  static void lookup_ice(const Spack& qitot, const Spack& nitot,
                         const Spack& qirim, const Spack& rhop, TableIce& t,
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
    const uview_1d<const Spack>& inv_dzq,
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
    const uview_1d<const Spack>& inv_dzq,
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
    const uview_1d<const Spack>& inv_dzq,
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
    const uview_1d<const Spack>& inv_dzq,
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
    const uview_1d<const Spack>& lcldm,
    const uview_1d<const Spack>& acn,
    const uview_1d<const Spack>& inv_dzq,
    const view_dnu_table& dnu,
    const MemberType& team,
    const Workspace& workspace,
    const Int& nk, const Int& ktop, const Int& kbot, const Int& kdir, const Scalar& dt, const Scalar& odt,
    const bool& log_predictNc,
    const uview_1d<Spack>& qc,
    const uview_1d<Spack>& nc,
    const uview_1d<Spack>& nc_incld,
    const uview_1d<Spack>& mu_c,
    const uview_1d<Spack>& lamc,
    const uview_1d<Spack>& qc_tend,
    const uview_1d<Spack>& nc_tend,
    Scalar& prt_liq);

  // TODO: comment
  KOKKOS_FUNCTION
  static void rain_sedimentation(
    const uview_1d<const Spack>& rho,
    const uview_1d<const Spack>& inv_rho,
    const uview_1d<const Spack>& rhofacr,
    const uview_1d<const Spack>& rcldm,
    const uview_1d<const Spack>& inv_dzq,
    const uview_1d<Spack>& qr_incld,
    const MemberType& team,
    const Workspace& workspace,
    const view_2d_table& vn_table, const view_2d_table& vm_table,
    const Int& nk, const Int& ktop, const Int& kbot, const Int& kdir, const Scalar& dt, const Scalar& odt,
    const uview_1d<Spack>& qr,
    const uview_1d<Spack>& nr,
    const uview_1d<Spack>& nr_incld,
    const uview_1d<Spack>& mu_r,
    const uview_1d<Spack>& lamr,
    const uview_1d<Spack>& rflx,
    const uview_1d<Spack>& qr_tend,
    const uview_1d<Spack>& nr_tend,
    Scalar& prt_liq);

  // TODO: comment
  KOKKOS_FUNCTION
  static void ice_sedimentation(
    const uview_1d<const Spack>& rho,
    const uview_1d<const Spack>& inv_rho,
    const uview_1d<const Spack>& rhofaci,
    const uview_1d<const Spack>& icldm,
    const uview_1d<const Spack>& inv_dzq,
    const MemberType& team,
    const Workspace& workspace,
    const Int& nk, const Int& ktop, const Int& kbot, const Int& kdir, const Scalar& dt, const Scalar& odt,
    const uview_1d<Spack>& qitot,
    const uview_1d<Spack>& qitot_incld,
    const uview_1d<Spack>& nitot,
    const uview_1d<Spack>& nitot_incld,
    const uview_1d<Spack>& qirim,
    const uview_1d<Spack>& qirim_incld,
    const uview_1d<Spack>& birim,
    const uview_1d<Spack>& birim_incld,
    const uview_1d<Spack>& qi_tend,
    const uview_1d<Spack>& ni_tend,
    const view_itab_table& itab,
    Scalar& prt_sol);

  // homogeneous freezing of cloud and rain
  KOKKOS_FUNCTION
  static void homogeneous_freezing(
    const uview_1d<const Spack>& t,
    const uview_1d<const Spack>& exner,
    const uview_1d<const Spack>& xlf,
    const MemberType& team,
    const Int& nk, const Int& ktop, const Int& kbot, const Int& kdir,
    const uview_1d<Spack>& qc,
    const uview_1d<Spack>& nc,
    const uview_1d<Spack>& qr,
    const uview_1d<Spack>& nr,
    const uview_1d<Spack>& qitot,
    const uview_1d<Spack>& nitot,
    const uview_1d<Spack>& qirim,
    const uview_1d<Spack>& birim,
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
    Spack& qcaut, Spack& qcacc, Spack &qccol, Spack& qcheti, Spack& qcshd, Spack& qiberg, Spack& qisub, Spack& qidep,
    const Smask& context = Smask(true) );

  KOKKOS_FUNCTION
  static void rain_water_conservation(
    const Spack& qr, const Spack& qcaut, const Spack& qcacc, const Spack& qimlt, const Spack& qcshd, const Scalar dt,
    Spack& qrevp, Spack& qrcol, Spack& qrheti,
    const Smask& context = Smask(true) );

  KOKKOS_FUNCTION
  static void ice_water_conservation(
    const Spack& qitot,const Spack& qidep,const Spack& qinuc,const Spack& qiberg, const Spack &qrcol,const Spack &qccol,const Spack& qrheti,const Spack& qcheti,const Scalar dt,
    Spack& qisub, Spack& qimlt,
    const Smask& context = Smask(true) );

  // TODO: comment
  KOKKOS_INLINE_FUNCTION
  static void get_cloud_dsd2(
    const Spack& qc, Spack& nc, Spack& mu_c, const Spack& rho, Spack& nu,
    const view_dnu_table& dnu, Spack& lamc, Spack& cdist, Spack& cdist1, const Spack& lcldm,
    const Smask& context = Smask(true) );

  // Computes and returns rain size distribution parameters
  KOKKOS_FUNCTION
  static void get_rain_dsd2 (
    const Spack& qr, Spack& nr, Spack& mu_r,
    Spack& lamr, Spack& cdistr, Spack& logn0r, const Spack& rcldmm,
    const Smask& context = Smask(true) );

  // Calculates rime density
  KOKKOS_FUNCTION
  static void calc_rime_density(const Spack& t, const Spack& rhofaci,
    const Spack& f1pr02, const Spack& acn, const Spack& lamc,
    const Spack& mu_c, const Spack& qc_incld, const Spack& qccol,
    Spack& vtrmi1, Spack& rhorime_c,
    const Smask& context = Smask(true) );

  // Computes contact and immersion freezing droplets
  KOKKOS_FUNCTION
  static void cldliq_immersion_freezing(const Spack& t, const Spack& lamc,
    const Spack& mu_c, const Spack& cdist1, const Spack& qc_incld, const Spack& qc_relvar,
    Spack& qcheti, Spack& ncheti,
    const Smask& context = Smask(true) );

  // Computes the immersion freezing of rain
  KOKKOS_FUNCTION
  static void rain_immersion_freezing(const Spack& t, const Spack& lamr,
    const Spack& mu_r, const Spack& cdistr, const Spack& qr_incld,
    Spack& qrheti, Spack& nrheti,
    const Smask& context = Smask(true) );

  // Computes droplet self collection
  KOKKOS_FUNCTION
  static void droplet_self_collection(const Spack& rho, const Spack& inv_rho,
    const Spack& qc_incld, const Spack& mu_c, const Spack& nu,
    const Spack& ncautc, Spack& ncslf,
    const Smask& context = Smask(true) );

  // Computes the accretion of clouds by rain
  KOKKOS_FUNCTION
  static void cloud_rain_accretion(const Spack& rho, const Spack& inv_rho,
    const Spack& qc_incld, const Spack& nc_incld, const Spack& qr_incld, const Spack& qc_relvar,
    Spack& qcacc, Spack& ncacc,
    const Smask& context = Smask(true) );

  // Computes cloud water autoconversion process rate
  KOKKOS_FUNCTION
  static void cloud_water_autoconversion(const Spack& rho,  const Spack& qc_incld,
    const Spack& nc_incld, const Spack& qc_relvar,
    Spack& qcaut, Spack& ncautc, Spack& ncautr,
    const Smask& context = Smask(true) );

  // Computes rain self collection process rate
  KOKKOS_FUNCTION
  static void rain_self_collection(const Spack& rho, const Spack& qr_incld, const Spack& nr_incld, Spack& nrslf,
                                   const Smask& context = Smask(true) );

  // Impose maximum ice number
  KOKKOS_FUNCTION
  static void impose_max_total_Ni(Spack& nitot_local, const Scalar& max_total_Ni, const Spack& inv_rho_local,
                                  const Smask& context = Smask(true) );

  //--------------------------------------------------------------------------------
  //  Calculates and returns the bulk rime density from the prognostic ice variables
  //  and adjusts qirim and birim appropriately.
  //--------------------------------------------------------------------------------
  KOKKOS_FUNCTION
  static Spack calc_bulk_rho_rime(
    const Spack& qi_tot, Spack& qi_rim, Spack& bi_rim,
    const Smask& context = Smask(true) );

  // TODO - comment
  KOKKOS_FUNCTION
  static void compute_rain_fall_velocity(
    const view_2d_table& vn_table, const view_2d_table& vm_table,
    const Spack& qr_incld, const Spack& rcldm, const Spack& rhofacr, 
    Spack& nr_incld, Spack& mu_r, Spack& lamr, Spack& V_qr, Spack& V_nr,
    const Smask& context = Smask(true));

  //---------------------------------------------------------------------------------
  // update prognostic microphysics and thermodynamics variables
  //---------------------------------------------------------------------------------
  //-- ice-phase dependent processes:
  KOKKOS_FUNCTION
  static void update_prognostic_ice(
    const Spack& qcheti, const Spack& qccol,
    const Spack& qcshd,  const Spack& nccol,  const Spack& ncheti, const Spack& ncshdc,
    const Spack& qrcol,  const Spack& nrcol,  const Spack& qrheti, const Spack& nrheti,
    const Spack& nrshdr, const Spack& qimlt,  const Spack& nimlt,  const Spack& qisub,
    const Spack& qidep,  const Spack& qinuc,  const Spack& ninuc,  const Spack& nislf,
    const Spack& nisub,  const Spack& qiberg, const Spack& exner,  const Spack& xxls,
    const Spack& xlf,    const bool log_predictNc, const Smask& log_wetgrowth, const Scalar dt,
    const Scalar& nmltratio, const Spack& rhorime_c, Spack& th, Spack& qv, Spack& qitot,
    Spack& nitot, Spack& qirim, Spack& birim, Spack& qc,  Spack& nc, Spack& qr, Spack& nr,
    const Smask& context = Smask(true));

  // TODO (comments)
  KOKKOS_FUNCTION
  static void get_time_space_phys_variables(const Spack& t, const Spack& pres, const Spack& rho,
					    const Spack& xxlv, const Spack& xxls,
					    const Spack& qvs, const Spack& qvi, Spack& mu,
					    Spack& dv, Spack& sc, Spack& dqsdt, Spack& dqsidt,
					    Spack& ab, Spack& abi, Spack& kap, Spack& eii,
                                            const Smask& context = Smask(true));

  // TODO (comments)
  KOKKOS_FUNCTION
  static void ice_cldliq_collection(const Spack& rho, const Spack& temp,
                                    const Spack& rhofaci, const Spack& f1pr04,
                                    const Spack& qitot_incld, const Spack& qc_incld,
                                    const Spack& nitot_incld, const Spack& nc_incld,
                                    Spack& qccol, Spack& nccol, Spack& qcshd, Spack& ncshdc,
                                    const Smask& context = Smask(true));

  // TODO (comments)
  KOKKOS_FUNCTION
  static void ice_rain_collection(const Spack& rho, const Spack& temp,
                                  const Spack& rhofaci, const Spack& logn0r,
                                  const Spack& f1pr07, const Spack& f1pr08,
                                  const Spack& qitot_incld, const Spack& nitot_incld,
                                  const Spack& qr_incld,
                                  Spack& qrcol, Spack& nrcol,
                                  const Smask& context = Smask(true));

  // TODO (comments)
  KOKKOS_FUNCTION
  static void ice_self_collection(const Spack& rho, const Spack& rhofaci,
                                  const Spack& f1pr03, const Spack& eii,
                                  const Spack& qirim_incld, const Spack& qitot_incld,
                                  const Spack& nitot_incld, Spack& nislf,
                                  const Smask& context = Smask(true));

  // TODO (comments)
  KOKKOS_FUNCTION
  static void evaporate_sublimate_precip(const Spack& qr_incld, const Spack& qc_incld,
					 const Spack& nr_incld, const Spack& qitot_incld,
					 const Spack& lcldm, const Spack& rcldm,
					 const Spack& qvs, const Spack& ab, const Spack& epsr,
					 const Spack& qv, Spack& qrevp, Spack& nrevp,
                                         const Smask& context = Smask(true));

  //get number and mass tendencies due to melting ice
  KOKKOS_FUNCTION
  static void ice_melting(const Spack& rho, const Spack& t, const Spack& pres, const Spack& rhofaci,
			  const Spack& f1pr05, const Spack& f1pr14, const Spack& xxlv, const Spack& xlf,
			  const Spack& dv, const Spack& sc, const Spack& mu, const Spack& kap,
			  const Spack& qv, const Spack& qitot_incld, const Spack& nitot_incld,
			  Spack& qimlt, Spack& nimlt,
                          const Smask& context = Smask(true));

  //liquid-phase dependent processes:
  KOKKOS_FUNCTION
  static void update_prognostic_liquid(const Spack& qcacc, const Spack& ncacc,
    const Spack& qcaut,const Spack& ncautc, const Spack& ncautr,
    const Spack& ncslf, const Spack& qrevp, const Spack& nrevp, const Spack& nrslf,
    const bool log_predictNc, const Spack& inv_rho, const Spack& exner, const Spack& xxlv,
    const Scalar dt, Spack& th, Spack& qv, Spack& qc, Spack& nc, Spack& qr, Spack& nr,
    const Smask& context = Smask(true));

  // TODO (comments)
  KOKKOS_FUNCTION
  static void ice_deposition_sublimation(const Spack& qitot_incld,
    const Spack& nitot_incld, const Spack& t, const Spack& qvs, const Spack& qvi,
    const Spack& epsi, const Spack& abi, const Spack& qv, Spack& qidep,
    Spack& qisub, Spack& nisub, Spack& qiberg,
    const Smask& context = Smask(true));

  KOKKOS_FUNCTION
  static void ice_relaxation_timescale(
    const Spack& rho, const Spack& temp, const Spack& rhofaci, const Spack& f1pr05,
    const Spack& f1pr14, const Spack& dv, const Spack& mu, const Spack& sc,
    const Spack& qitot_incld, const Spack& nitot_incld,
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
                             const Spack& nitot, const Spack& naai,
                             const Spack& supi, const Scalar& odt,
                             const bool& log_predictNc,
                             Spack& qinuc, Spack& ninuc,
                             const Smask& context = Smask(true));

  KOKKOS_FUNCTION
  static Spack subgrid_variance_scaling(const Spack& relvar, const Scalar& expon);

  KOKKOS_FUNCTION
  static void ice_cldliq_wet_growth(const Spack& rho, const Spack& temp, const Spack& pres, const Spack& rhofaci, const Spack& f1pr05,
                                    const Spack& f1pr14, const Spack& xxlv, const Spack& xlf, const Spack& dv,
                                    const Spack& kap, const Spack& mu, const Spack& sc, const Spack& qv, const Spack& qc_incld,
                                    const Spack& qitot_incld, const Spack& nitot_incld, const Spack& qr_incld,
                                    Smask& log_wetgrowth, Spack& qrcol, Spack& qccol, Spack& qwgrth, Spack& nrshdr, Spack& qcshd,
                                    const Smask& context = Smask(true));

  // Note: not a kernel function
  static void get_latent_heat(const Int& ni, const Int& nk, view_2d<Spack>& v, view_2d<Spack>& s, view_2d<Spack>& f);

  KOKKOS_FUNCTION
  static void check_values(const uview_1d<const Spack>& qv, const uview_1d<const Spack>& temp, const Int& ktop, const Int& kbot,
                           const Int& timestepcount, const bool& force_abort, const Int& source_ind, const MemberType& team,
                           const uview_1d<const Scalar>& col_loc);

  KOKKOS_FUNCTION
  static void calculate_incloud_mixingratios(
    const Spack& qc, const Spack& qr, const Spack& qitot, const Spack& qirim, const Spack& nc,
    const Spack& nr, const Spack& nitot, const Spack& birim, const Spack& inv_lcldm,
    const Spack& inv_icldm, const Spack& inv_rcldm,
    Spack& qc_incld, Spack& qr_incld, Spack& qitot_incld, Spack& qirim_incld,
    Spack& nc_incld, Spack& nr_incld, Spack& nitot_incld, Spack& birim_incld,
    const Smask& context = Smask(true));

  //
  // main P3 functions
  //

  KOKKOS_FUNCTION
  static void p3_main_init(
    const MemberType& team,
    const Int& nk_pack,
    const uview_1d<const Spack>& icldm,
    const uview_1d<const Spack>& lcldm,
    const uview_1d<const Spack>& rcldm,
    const uview_1d<const Spack>& exner,
    const uview_1d<const Spack>& th,
    const uview_1d<const Spack>& dzq,
    const uview_1d<Spack>& diag_ze,
    const uview_1d<Spack>& ze_ice,
    const uview_1d<Spack>& ze_rain,
    const uview_1d<Spack>& diag_effc,
    const uview_1d<Spack>& diag_effi,
    const uview_1d<Spack>& inv_icldm,
    const uview_1d<Spack>& inv_lcldm,
    const uview_1d<Spack>& inv_rcldm,
    const uview_1d<Spack>& inv_exner,
    const uview_1d<Spack>& t,
    const uview_1d<Spack>& qv,
    const uview_1d<Spack>& inv_dzq,
    Scalar& prt_liq,
    Scalar& prt_sol,
    view_1d_ptr_array<Spack, 40>& zero_init);

  KOKKOS_FUNCTION
  static void p3_main_part1(
    const MemberType& team,
    const Int& nk,
    const bool& log_predictNc,
    const Scalar& dt,
    const uview_1d<const Spack>& pres,
    const uview_1d<const Spack>& pdel,
    const uview_1d<const Spack>& dzq,
    const uview_1d<const Spack>& ncnuc,
    const uview_1d<const Spack>& exner,
    const uview_1d<const Spack>& inv_exner,
    const uview_1d<const Spack>& inv_lcldm,
    const uview_1d<const Spack>& inv_icldm,
    const uview_1d<const Spack>& inv_rcldm,
    const uview_1d<const Spack>& xxlv,
    const uview_1d<const Spack>& xxls,
    const uview_1d<const Spack>& xlf,
    const uview_1d<Spack>& t,
    const uview_1d<Spack>& rho,
    const uview_1d<Spack>& inv_rho,
    const uview_1d<Spack>& qvs,
    const uview_1d<Spack>& qvi,
    const uview_1d<Spack>& supi,
    const uview_1d<Spack>& rhofacr,
    const uview_1d<Spack>& rhofaci,
    const uview_1d<Spack>& acn,
    const uview_1d<Spack>& qv,
    const uview_1d<Spack>& th,
    const uview_1d<Spack>& qc,
    const uview_1d<Spack>& nc,
    const uview_1d<Spack>& qr,
    const uview_1d<Spack>& nr,
    const uview_1d<Spack>& qitot,
    const uview_1d<Spack>& nitot,
    const uview_1d<Spack>& qirim,
    const uview_1d<Spack>& birim,
    const uview_1d<Spack>& qc_incld,
    const uview_1d<Spack>& qr_incld,
    const uview_1d<Spack>& qitot_incld,
    const uview_1d<Spack>& qirim_incld,
    const uview_1d<Spack>& nc_incld,
    const uview_1d<Spack>& nr_incld,
    const uview_1d<Spack>& nitot_incld,
    const uview_1d<Spack>& birim_incld,
    bool& log_nucleationPossible,
    bool& log_hydrometeorsPresent);

  KOKKOS_FUNCTION
  static void p3_main_part2(
    const MemberType& team,
    const Int& nk_pack,
    const bool& log_predictNc,
    const Scalar& dt,
    const Scalar& odt,
    const view_dnu_table& dnu,
    const view_itab_table& itab,
    const view_itabcol_table& itabcol,
    const view_2d_table& revap_table,
    const uview_1d<const Spack>& pres,
    const uview_1d<const Spack>& pdel,
    const uview_1d<const Spack>& dzq,
    const uview_1d<const Spack>& ncnuc,
    const uview_1d<const Spack>& exner,
    const uview_1d<const Spack>& inv_exner,
    const uview_1d<const Spack>& inv_lcldm,
    const uview_1d<const Spack>& inv_icldm,
    const uview_1d<const Spack>& inv_rcldm,
    const uview_1d<const Spack>& naai,
    const uview_1d<const Spack>& qc_relvar,
    const uview_1d<const Spack>& icldm,
    const uview_1d<const Spack>& lcldm,
    const uview_1d<const Spack>& rcldm,
    const uview_1d<Spack>& t,
    const uview_1d<Spack>& rho,
    const uview_1d<Spack>& inv_rho,
    const uview_1d<Spack>& qvs,
    const uview_1d<Spack>& qvi,
    const uview_1d<Spack>& supi,
    const uview_1d<Spack>& rhofacr,
    const uview_1d<Spack>& rhofaci,
    const uview_1d<Spack>& acn,
    const uview_1d<Spack>& qv,
    const uview_1d<Spack>& th,
    const uview_1d<Spack>& qc,
    const uview_1d<Spack>& nc,
    const uview_1d<Spack>& qr,
    const uview_1d<Spack>& nr,
    const uview_1d<Spack>& qitot,
    const uview_1d<Spack>& nitot,
    const uview_1d<Spack>& qirim,
    const uview_1d<Spack>& birim,
    const uview_1d<Spack>& xxlv,
    const uview_1d<Spack>& xxls,
    const uview_1d<Spack>& xlf,
    const uview_1d<Spack>& qc_incld,
    const uview_1d<Spack>& qr_incld,
    const uview_1d<Spack>& qitot_incld,
    const uview_1d<Spack>& qirim_incld,
    const uview_1d<Spack>& nc_incld,
    const uview_1d<Spack>& nr_incld,
    const uview_1d<Spack>& nitot_incld,
    const uview_1d<Spack>& birim_incld,
    const uview_1d<Spack>& omu_c,
    const uview_1d<Spack>& nu,
    const uview_1d<Spack>& olamc,
    const uview_1d<Spack>& cdist,
    const uview_1d<Spack>& cdist1,
    const uview_1d<Spack>& cdistr,
    const uview_1d<Spack>& mu_r,
    const uview_1d<Spack>& lamr,
    const uview_1d<Spack>& logn0r,
    const uview_1d<Spack>& cmeiout,
    const uview_1d<Spack>& prain,
    const uview_1d<Spack>& nevapr,
    const uview_1d<Spack>& prer_evap,
    const uview_1d<Spack>& vap_liq_exchange,
    const uview_1d<Spack>& vap_ice_exchange,
    const uview_1d<Spack>& liq_ice_exchange,
    const uview_1d<Spack>& pratot,
    const uview_1d<Spack>& prctot,
    bool& log_hydrometeorsPresent);

  KOKKOS_FUNCTION
  static void p3_main_part3(
    const MemberType& team,
    const Int& nk_pack,
    const view_dnu_table& dnu,
    const view_itab_table& itab,
    const uview_1d<const Spack>& exner,
    const uview_1d<const Spack>& lcldm,
    const uview_1d<const Spack>& rcldm,
    const uview_1d<Spack>& rho,
    const uview_1d<Spack>& inv_rho,
    const uview_1d<Spack>& rhofaci,
    const uview_1d<Spack>& qv,
    const uview_1d<Spack>& th,
    const uview_1d<Spack>& qc,
    const uview_1d<Spack>& nc,
    const uview_1d<Spack>& qr,
    const uview_1d<Spack>& nr,
    const uview_1d<Spack>& qitot,
    const uview_1d<Spack>& nitot,
    const uview_1d<Spack>& qirim,
    const uview_1d<Spack>& birim,
    const uview_1d<Spack>& xxlv,
    const uview_1d<Spack>& xxls,
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
    const uview_1d<Spack>& diag_rhoi,
    const uview_1d<Spack>& diag_ze,
    const uview_1d<Spack>& diag_effc);

  static void p3_main(
    // inputs
    const view_2d<const Spack>& pres,          // pressure                             Pa
    const view_2d<const Spack>& dzq,           // vertical grid spacing                m
    const view_2d<const Spack>& ncnuc,         // IN ccn activated number tendency     kg-1 s-1
    const view_2d<const Spack>& naai,          // IN actived ice nuclei concentration  1/kg
    const view_2d<const Spack>& qc_relvar,     // assumed SGS 1/(var(qc)/mean(qc))     kg2/kg2
    const Real&                 dt,            // model time step                      s
    const Int&                  ni,            // num columns
    const Int&                  nk,            // column size
    const Int&                  it,            // time step counter NOTE: starts at 1 for first time step
    const bool&                 log_predictNc, // .T. (.F.) for prediction (specification) of Nc
    const view_2d<const Spack>& pdel,          // pressure thickness                   Pa
    const view_2d<const Spack>& exner,         // Exner expression

    // inputs needed for PBUF variables used by other parameterizations
    const view_2d<const Spack>& icldm,         // ice cloud fraction
    const view_2d<const Spack>& lcldm,         // liquid cloud fraction
    const view_2d<const Spack>& rcldm,         // rain cloud fraction
    const view_2d<const Scalar>& col_location, // ni x 3

    // input/output  arguments
    const view_2d<Spack>& qc,    // cloud, mass mixing ratio         kg kg-1
    const view_2d<Spack>& nc,    // cloud, number mixing ratio       #  kg-1
    const view_2d<Spack>& qr,    // rain, mass mixing ratio          kg kg-1
    const view_2d<Spack>& nr,    // rain, number mixing ratio        #  kg-1
    const view_2d<Spack>& qitot, // ice, total mass mixing ratio     kg kg-1
    const view_2d<Spack>& qirim, // ice, rime mass mixing ratio      kg kg-1
    const view_2d<Spack>& nitot, // ice, total number mixing ratio   #  kg-1
    const view_2d<Spack>& birim, // ice, rime volume mixing ratio    m3 kg-1
    const view_2d<Spack>& qv,    // water vapor mixing ratio         kg kg-1
    const view_2d<Spack>& th,    // potential temperature            K

    // output arguments
    const view_1d<Scalar>& prt_liq,  // precipitation rate, liquid       m s-1
    const view_1d<Scalar>& prt_sol,  // precipitation rate, solid        m s-1
    const view_2d<Spack>& diag_ze,   // equivalent reflectivity          dBZ
    const view_2d<Spack>& diag_effc, // effective radius, cloud          m
    const view_2d<Spack>& diag_effi, // effective radius, ice            m
    const view_2d<Spack>& diag_vmi,  // mass-weighted fall speed of ice  m s-1
    const view_2d<Spack>& diag_di,   // mean diameter of ice             m
    const view_2d<Spack>& diag_rhoi, // bulk density of ice              kg m-3
    const view_2d<Spack>& mu_c,      // Size distribution shape parameter for radiation
    const view_2d<Spack>& lamc,      // Size distribution slope parameter for radiation

    // outputs for PBUF variables used by other parameterizations
    const view_2d<Spack>& cmeiout,          // qitend due to deposition/sublimation
    const view_2d<Spack>& prain,            // Total precipitation (rain + snow)
    const view_2d<Spack>& nevapr,           // evaporation of total precipitation (rain + snow)
    const view_2d<Spack>& prer_evap,        // evaporation of rain
    const view_2d<Spack>& rflx,             // grid-box average rain flux (kg m^-2 s^-1) pverp
    const view_2d<Spack>& sflx,             // grid-box average ice/snow flux (kg m^-2 s^-1) pverp
    const view_2d<Spack>& pratot,           // accretion of cloud by rain
    const view_2d<Spack>& prctot,           // autoconversion of cloud to rain
    const view_2d<Spack>& liq_ice_exchange, // sum of liq-ice phase change tendencies
    const view_2d<Spack>& vap_liq_exchange, // sum of vap-liq phase change tendencies
    const view_2d<Spack>& vap_ice_exchange);// sum of vap-ice phase change tendencies
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
# include "p3_evaporate_sublimate_precip_impl.hpp"
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

#endif
