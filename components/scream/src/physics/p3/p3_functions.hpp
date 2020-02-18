#ifndef P3_FUNCTIONS_HPP
#define P3_FUNCTIONS_HPP

#include "share/scream_types.hpp"
#include "share/scream_pack_kokkos.hpp"
#include "share/scream_workspace.hpp"
#include "p3_constants.hpp"

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

  using C = Constants<Scalar>;

  template <typename S>
  using view_1d = typename KT::template view_1d<S>;
  template <typename S>
  using view_2d = typename KT::template view_2d<S>;

  template <typename S, int N>
  using view_1d_ptr_array = typename KT::template view_1d_ptr_carray<S, N>;

  template <typename S>
  using uview_1d = typename ko::template Unmanaged<view_1d<S> >;

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
    view_2d_table& vn_table, view_2d_table& vm_table, view_1d_table& mu_r_table, view_dnu_table& dnu);

  static void init_kokkos_ice_lookup_tables(
    view_itab_table& itab, view_itabcol_table& itabcol);

  // Map (mu_r, lamr) to Table3 data.
  KOKKOS_FUNCTION
  static void lookup(const Smask& qr_gt_small, const Spack& mu_r, const Spack& lamr,
                     Table3& t);

  //------------------------------------------------------------------------------------------!
  // Finds indices in 3D ice (only) lookup table
  // ------------------------------------------------------------------------------------------!
  KOKKOS_FUNCTION
  static void lookup_ice(const Smask& qiti_gt_small, const Spack& qitot, const Spack& nitot,
                         const Spack& qirim, const Spack& rhop, TableIce& t);

  //------------------------------------------------------------------------------------------!
  // Finds indices in 3D rain lookup table
  //------------------------------------------------------------------------------------------!
  KOKKOS_FUNCTION
  static void lookup_rain(const Smask& qiti_gt_small, const Spack& qr, const Spack& nr, TableRain& t);

  // Apply Table3 data to the table to return a value. This performs bilinear
  // interpolation within the quad given by {t.dumii, t.dumjj} x {t.dumii+1,
  // t.dumjj+1}.
  KOKKOS_FUNCTION
  static Spack apply_table(const Smask& qr_gt_small, const view_2d_table& table,
                           const Table3& t);

  // Apply TableIce data to the ice tables to return a value.
  KOKKOS_FUNCTION
  static Spack apply_table_ice(const Smask& qiti_gt_small, const int& index, const view_itab_table& itab,
                               const TableIce& t);

  // Interpolates lookup table values for rain/ice collection processes
  KOKKOS_FUNCTION
  static Spack apply_table_coll(const Smask& qiti_gt_small, const int& index, const view_itabcol_table& itabcoll,
                                const TableIce& ti, const TableRain& tr);

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
    const uview_1d<const Spack>& qc_incld,
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
    const uview_1d<const Spack>& qr_incld,
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

  //  compute saturation vapor pressure
  //  polysvp1 returned in units of pa.
  //  t is input in units of k.
  //  ice refers to saturation with respect to liquid (false) or ice (true)
  KOKKOS_FUNCTION
  static Spack polysvp1(const Spack& t, const bool ice);

  // Calls polysvp1 to obtain the saturation vapor pressure, and then computes
  // and returns the saturation mixing ratio, with respect to either liquid or ice,
  // depending on value of 'ice'
  KOKKOS_FUNCTION
  static Spack qv_sat(const Spack& t_atm, const Spack& p_atm, const bool ice);

  KOKKOS_FUNCTION
  static void cloud_water_conservation(const Spack& qc, const Spack& qcnuc,const Scalar dt,
   Spack& qcaut, Spack& qcacc, Spack &qccol, Spack& qcheti, Spack& qcshd, Spack& qiberg, Spack& qisub, Spack& qidep);

  KOKKOS_FUNCTION
  static void rain_water_conservation(const Spack& qr, const Spack& qcaut, const Spack& qcacc, const Spack& qimlt, const Spack& qcshd, const Scalar dt,
   Spack& qrevp, Spack& qrcol, Spack& qrheti);

  KOKKOS_FUNCTION
  static void ice_water_conservation(const Spack& qitot,const Spack& qidep,const Spack& qinuc,const Spack& qiberg, const Spack &qrcol,const Spack &qccol,const Spack& qrheti,const Spack& qcheti,const Scalar dt, 
   Spack& qisub, Spack& qimlt);

  // TODO: comment
  template <bool zero_out=true>
  KOKKOS_INLINE_FUNCTION
  static void get_cloud_dsd2(
    const Smask& qc_gt_small, const Spack& qc, Spack& nc, Spack& mu_c, const Spack& rho, Spack& nu,
    const view_dnu_table& dnu, Spack& lamc, Spack& cdist, Spack& cdist1, const Spack& lcldm);

  // Computes and returns rain size distribution parameters
  KOKKOS_FUNCTION
  static void get_rain_dsd2 (
    const Smask& qr_gt_small, const Spack& qr, Spack& nr, Spack& mu_r,
    Spack& lamr, Spack& cdistr, Spack& logn0r, const Spack& rcldm);

  KOKKOS_FUNCTION
  static void cloud_water_autoconversion(const Spack& rho,  const Spack& qc_incld, const Spack& nc_incld,
    Spack& qcaut, Spack& ncautc, Spack& ncautr);

  KOKKOS_FUNCTION
  static void rain_self_collection(const Spack& rho, const Spack& qr_incld, const Spack& nr_incld, Spack& nrslf);

  //--------------------------------------------------------------------------------
  //  Calculates and returns the bulk rime density from the prognostic ice variables
  //  and adjusts qirim and birim appropriately.
  //--------------------------------------------------------------------------------
  KOKKOS_FUNCTION
  static Spack calc_bulk_rho_rime(
    const Smask& qi_gt_small, const Spack& qi_tot, Spack& qi_rim, Spack& bi_rim);

  // TODO - comment
  KOKKOS_FUNCTION
  static void compute_rain_fall_velocity(
    const Smask& qr_gt_small, const view_2d_table& vn_table, const view_2d_table& vm_table,
    const Spack& qr_incld, const Spack& rcldm, const Spack& rhofacr, Spack& nr,
    Spack& nr_incld, Spack& mu_r, Spack& lamr, Spack& V_qr, Spack& V_nr);

  //---------------------------------------------------------------------------------
  // update prognostic microphysics and thermodynamics variables
  //---------------------------------------------------------------------------------
  //-- ice-phase dependent processes:
  KOKKOS_FUNCTION
  static void update_prognostic_ice(const Spack& qcheti, const Spack& qccol,
    const Spack& qcshd,  const Spack& nccol,  const Spack& ncheti, const Spack& ncshdc,
    const Spack& qrcol,  const Spack& nrcol,  const Spack& qrheti, const Spack& nrheti,
    const Spack& nrshdr, const Spack& qimlt,  const Spack& nimlt,  const Spack& qisub,
    const Spack& qidep,  const Spack& qinuc,  const Spack& ninuc,  const Spack& nislf,
    const Spack& nisub,  const Spack& qiberg, const Spack& exner,  const Spack& xxls,
    const Spack& xlf,    const bool log_predictNc, const bool log_wetgrowth, const Scalar dt,
    const Spack& nmltratio, const Spack& rhorime_c, Spack& th, Spack& qv, Spack& qitot,
    Spack& nitot, Spack& qirim, Spack& birim, Spack& qc,  Spack& nc, Spack& qr,
    Spack& nr);

  // TODO (comments)
  KOKKOS_FUNCTION
  static void ice_cldliq_collection(const Spack& rho, const Spack& temp,
                                    const Spack& rhofaci, const Spack& f1pr04,
                                    const Spack& qitot_incld, const Spack& qc_incld,
                                    const Spack& nitot_incld, const Spack& nc_incld,
                                    Spack& qccol, Spack& nccol, Spack& qcshd, Spack& ncshdc);

  // TODO (comments)
  KOKKOS_FUNCTION
  static void ice_rain_collection(const Spack& rho, const Spack& temp,
                                  const Spack& rhofaci, const Spack& logn0r,
                                  const Spack& f1pr07, const Spack& f1pr08,
                                  const Spack& qitot_incld, const Spack& nitot_incld,
                                  const Spack& qr_incld,
                                  Spack& qrcol, Spack& nrcol);
  
  // TODO (comments)
  KOKKOS_FUNCTION
  static void ice_self_collection(const Spack& rho, const Spack& rhofaci,
                                  const Spack& f1pr03, const Spack& eii,
                                  const Spack& qirim_incld, const Spack& qitot_incld,
                                  const Spack& nitot_incld, Spack& nislf);

};

template <typename ScalarT, typename DeviceT>
constexpr ScalarT Functions<ScalarT, DeviceT>::P3C::lookup_table_1a_dum1_c;

extern "C" {
// decl of fortran function for loading tables from fortran p3. This will
// continue to be a bit awkward until we have fully ported all of p3.
void init_tables_from_f90_c(Real* vn_table_data, Real* vm_table_data, Real* mu_table_data);
}

} // namespace p3
} // namespace scream

// If a GPU build, make all code available to the translation unit; otherwise,
// ETI is used.
#ifdef KOKKOS_ENABLE_CUDA
# include "p3_functions_math_impl.hpp"
# include "p3_functions_table3_impl.hpp"
# include "p3_functions_table_ice_impl.hpp"
# include "p3_functions_dsd2_impl.hpp"
# include "p3_functions_upwind_impl.hpp"
# include "p3_functions_find_impl.hpp"
# include "p3_functions_conservation_impl.hpp"
# include "p3_functions_autoconversion_impl.hpp"
# include "p3_functions_cloud_sed_impl.hpp"
# include "p3_functions_ice_sed_impl.hpp"
# include "p3_functions_rain_sed_impl.hpp"
# include "p3_functions_update_prognostics_impl.hpp"
# include "p3_functions_ice_collection_impl.hpp"
#endif

#endif
