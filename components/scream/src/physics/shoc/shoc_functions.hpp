#ifndef SHOC_FUNCTIONS_HPP
#define SHOC_FUNCTIONS_HPP

#include "physics/share/physics_constants.hpp"
#include "physics/shoc/shoc_constants.hpp"

#include "share/scream_types.hpp"

#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/ekat_workspace.hpp"

namespace scream {
namespace shoc {

/*
 * Functions is a stateless struct used to encapsulate a
 * number of functions for SHOC. We use the ETI pattern for
 * these functions.
 *
 * SHOC assumptions:
 *  - Kokkos team policies have a vector length of 1
 */

template <typename ScalarT, typename DeviceT>
struct Functions
{
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

  using Mask  = ekat::Mask<Pack::n>;
  using Smask = ekat::Mask<Spack::n>;

  using KT = ekat::KokkosTypes<Device>;

  using C  = physics::Constants<Scalar>;
  using SC = shoc::Constants<Scalar>;

  template <typename S>
  using view_1d = typename KT::template view_1d<S>;
  template <typename S>
  using view_2d = typename KT::template view_2d<S>;
  template <typename S>
  using view_3d = typename KT::template view_3d<S>;

  template <typename S, int N>
  using view_1d_ptr_array = typename KT::template view_1d_ptr_carray<S, N>;

  template <typename S>
  using uview_1d = typename ekat::template Unmanaged<view_1d<S> >;

  template <typename S>
  using uview_2d = typename ekat::template Unmanaged<view_2d<S> >;

  using MemberType = typename KT::MemberType;

  using Workspace = typename ekat::WorkspaceManager<Spack, Device>::Workspace;

  //
  // --------- Functions ---------
  //
  KOKKOS_FUNCTION
  static void calc_shoc_varorcovar(
    const MemberType&            team,
    const Int&                   nlev,
    const Scalar&                tunefac,
    const uview_1d<const Spack>& isotropy_zi,
    const uview_1d<const Spack>& tkh_zi,
    const uview_1d<const Spack>& dz_zi,
    const uview_1d<const Spack>& invar1,
    const uview_1d<const Spack>& invar2,
    const uview_1d<Spack>&       varorcovar);

  KOKKOS_FUNCTION
  static void calc_shoc_vertflux(
    const MemberType& team,
    const Int& nlev,
    const uview_1d<const Spack>& tkh_zi,
    const uview_1d<const Spack>& dz_zi,
    const uview_1d<const Spack>& invar,
    const uview_1d<Spack>& vertflux);

  KOKKOS_FUNCTION
  static void shoc_diag_second_moments_srf(
    const Scalar& wthl_sfc, const Scalar& uw_sfc, const Scalar& vw_sfc,
    Scalar& ustar2, Scalar& wstar);

  KOKKOS_FUNCTION
  static void shoc_diag_second_moments_ubycond(
    Scalar& thl_sec, Scalar& qw_sec, Scalar& wthl_sec, Scalar& wqw_sec,
    Scalar& qwthl_sec, Scalar& uw_sec, Scalar& vw_sec, Scalar& wtke_sec);

  KOKKOS_FUNCTION
  static void update_host_dse(
    const MemberType& team,
    const Int& nlev,
    const uview_1d<const Spack>& thlm,
    const uview_1d<const Spack>& shoc_ql,
    const uview_1d<const Spack>& exner,
    const uview_1d<const Spack>& zt_grid,
    const Scalar& phis,
    const uview_1d<Spack>& host_dse);

  KOKKOS_FUNCTION
  static void compute_diag_third_shoc_moment(
    const MemberType& team,
    const Int& nlev,
    const Int& nlevi,
    const uview_1d<const Spack>& w_sec,
    const uview_1d<const Spack>& thl_sec,
    const uview_1d<const Spack>& wthl_sec,
    const uview_1d<const Spack>& tke,
    const uview_1d<const Spack>& dz_zt,
    const uview_1d<const Spack>& dz_zi,
    const uview_1d<const Spack>& isotropy_zi,
    const uview_1d<const Spack>& brunt_zi,
    const uview_1d<const Spack>& w_sec_zi,
    const uview_1d<const Spack>& thetal_zi,
    const uview_1d<Spack>& w3);

  KOKKOS_FUNCTION
  static void shoc_pblintd_init_pot(
    const MemberType& team, const Int& nlev,
    const view_1d<const Spack>& thl, const view_1d<const Spack>& ql, const view_1d<const Spack>& q,
    const view_1d<Spack>& thv);

  KOKKOS_FUNCTION
  static void compute_shoc_mix_shoc_length(
    const MemberType&            team,
    const Int&                   nlev,
    const uview_1d<const Spack>& tke,
    const uview_1d<const Spack>& brunt,
    const Scalar&                tscale,
    const uview_1d<const Spack>& zt_grid,
    const Scalar&                l_inf,
    const uview_1d<Spack>&       shoc_mix);

  KOKKOS_FUNCTION
  static void check_tke(
    const MemberType& team,
    const Int& nlev,
    const uview_1d<Spack>& tke);

  KOKKOS_FUNCTION
  static void clipping_diag_third_shoc_moments(
    const MemberType& team,
    const Int& nlevi,
    const uview_1d<const Spack>& w_sec_zi,
    const uview_1d<Spack>& w3);

  KOKKOS_FUNCTION
  static void linear_interp(
    const MemberType& team,
    const uview_1d<const Spack>& x1,
    const uview_1d<const Spack>& x2,
    const uview_1d<const Spack>& y1,
    const uview_1d<Spack>& y2,
    const Int& km1,
    const Int& km2,
    const Scalar& minthresh);

  KOKKOS_FUNCTION
  static void shoc_energy_integrals(
    const MemberType&            team,
    const Int&                   nlev,
    const uview_1d<const Spack>& host_dse,
    const uview_1d<const Spack>& pdel,
    const uview_1d<const Spack>& rtm,
    const uview_1d<const Spack>& rcm,
    const uview_1d<const Spack>& u_wind,
    const uview_1d<const Spack>& v_wind,
    Scalar&                      se_int,
    Scalar&                      ke_int,
    Scalar&                      wv_int,
    Scalar&                      wl_int);

  KOKKOS_FUNCTION
  static void shoc_diag_second_moments_lbycond(
    const Scalar& wthl_sfc, const Scalar& wqw_sfc, const Scalar& uw_sfc, const Scalar& vw_sfc,
    const Scalar& ustar2, const Scalar& wstar,
    Scalar& wthl_sec, Scalar& wqw_sec, Scalar& uw_sec, Scalar& vw_sec,
    Scalar& wtke_sec, Scalar& thl_sec, Scalar& qw_sec, Scalar& qwthl_sec);

  KOKKOS_FUNCTION
  static void diag_second_moments(const MemberType& team, const Int& nlev, const Int& nlevi,
     const uview_1d<const Spack>& thetal, const uview_1d<const Spack>& qw, const uview_1d<const Spack>& u_wind,
     const uview_1d<const Spack>& v_wind, const uview_1d<const Spack>& tke, const uview_1d<const Spack>& isotropy,
     const uview_1d<const Spack>& tkh, const uview_1d<const Spack>& tk, const uview_1d<const Spack>& dz_zi,
     const uview_1d<const Spack>& zt_grid, const uview_1d<const Spack>& zi_grid, const uview_1d<const Spack>& shoc_mix,
     const uview_1d<Spack>& isotropy_zi, const uview_1d<Spack>& tkh_zi, const uview_1d<Spack>& tk_zi,
     const uview_1d<Spack>& thl_sec, const uview_1d<Spack>& qw_sec, const uview_1d<Spack>& wthl_sec,
     const uview_1d<Spack>& wqw_sec, const uview_1d<Spack>& qwthl_sec, const uview_1d<Spack>& uw_sec,
     const uview_1d<Spack>& vw_sec, const uview_1d<Spack>& wtke_sec, const uview_1d<Spack>& w_sec);

  KOKKOS_FUNCTION
  static void diag_second_shoc_moments(const MemberType& team, const Int& nlev, const Int& nlevi,
     const uview_1d<const Spack>& thetal, const uview_1d<const Spack>& qw, const uview_1d<const Spack>& u_wind,
     const uview_1d<const Spack>& v_wind, const uview_1d<const Spack>& tke, const uview_1d<const Spack>& isotropy,
     const uview_1d<const Spack>& tkh, const uview_1d<const Spack>& tk, const uview_1d<const Spack>& dz_zi,
     const uview_1d<const Spack>& zt_grid, const uview_1d<const Spack>& zi_grid, const uview_1d<const Spack>& shoc_mix,
     const Scalar& wthl_sfc, const Scalar& wqw_sfc, const Scalar& uw_sfc, const Scalar& vw_sfc, Scalar& ustar2, Scalar& wstar,
     const uview_1d<Spack>& isotropy_zi, const uview_1d<Spack>& tkh_zi, const uview_1d<Spack>& tk_zi, const uview_1d<Spack>& thl_sec,
     const uview_1d<Spack>& qw_sec, const uview_1d<Spack>& wthl_sec, const uview_1d<Spack>& wqw_sec, const uview_1d<Spack>& qwthl_sec,
     const uview_1d<Spack>& uw_sec, const uview_1d<Spack>& vw_sec, const uview_1d<Spack>& wtke_sec, const uview_1d<Spack>& w_sec);

  KOKKOS_FUNCTION
  static void compute_brunt_shoc_length(
    const MemberType&            team,
    const Int&                   nlev,
    const Int&                   nlevi,
    const uview_1d<const Spack>& dz_zt,
    const uview_1d<const Spack>& thv,
    const uview_1d<const Spack>& thv_zi,
    const uview_1d<Spack>&       brunt);

  KOKKOS_FUNCTION
  static void compute_l_inf_shoc_length(
    const MemberType&            team,
    const Int&                   nlev,
    const uview_1d<const Spack>& zt_grid,
    const uview_1d<const Spack>& dz_zt,
    const uview_1d<const Spack>& tke,
    Scalar&                      l_inf);

  KOKKOS_FUNCTION
  static void check_length_scale_shoc_length(
    const MemberType&      team,
    const Int&             nlev,
    const Scalar&          host_dx,
    const Scalar&          host_dy,
    const uview_1d<Spack>& shoc_mix);

  KOKKOS_FUNCTION
  static void compute_conv_vel_shoc_length(
    const MemberType&            team,
    const Int&                   nlev,
    const Scalar&                pblh,
    const uview_1d<const Spack>& zt_grid,
    const uview_1d<const Spack>& dz_zt,
    const uview_1d<const Spack>& thv,
    const uview_1d<const Spack>& wthv_sec,
    Scalar&                      conv_vel);

  KOKKOS_FUNCTION
  static void shoc_diag_obklen(
    const Scalar& uw_sfc,
    const Scalar& vw_sfc,
    const Scalar& wthl_sfc,
    const Scalar& wqw_sfc,
    const Scalar& thl_sfc,
    const Scalar& cldliq_sfc,
    const Scalar& qv_sfc,
    Scalar&       ustar,
    Scalar&       kbfs,
    Scalar&       obklen);

  KOKKOS_FUNCTION
  static void shoc_pblintd_cldcheck(
    const Scalar& zi, const Scalar& cldn,
    Scalar& pblh);

  KOKKOS_FUNCTION
  static void compute_conv_time_shoc_length(
    const Scalar& pblh,
    Scalar&       conv_vel,
    Scalar&       tscale);

  KOKKOS_FUNCTION
  static void shoc_length(
    const MemberType&            team,
    const Int&                   nlev,
    const Int&                   nlevi,
    const Scalar&                host_dx,
    const Scalar&                host_dy,
    const Scalar&                pblh,
    const uview_1d<const Spack>& tke,
    const uview_1d<const Spack>& zt_grid,
    const uview_1d<const Spack>& zi_grid,
    const uview_1d<const Spack>& dz_zt,
    const uview_1d<const Spack>& dz_zi,
    const uview_1d<const Spack>& wthv_sec,
    const uview_1d<const Spack>& thetal,
    const uview_1d<const Spack>& thv,
    const uview_1d<Spack>&       thv_zi,
    const uview_1d<Spack>&       brunt,
    const uview_1d<Spack>&       shoc_mix);

  KOKKOS_FUNCTION
  static void shoc_energy_fixer(
    const MemberType&            team,
    const Int&                   nlev,
    const Int&                   nlevi,
    const Scalar&                dtime,
    const Int&                   nadv,
    const uview_1d<const Spack>& zt_grid,
    const uview_1d<const Spack>& zi_grid,
    const Scalar&                se_b,
    const Scalar&                ke_b,
    const Scalar&                wv_b,
    const Scalar&                wl_b,
    const Scalar&                se_a,
    const Scalar&                ke_a,
    const Scalar&                wv_a,
    const Scalar&                wl_a,
    const Scalar&                wthl_sfc,
    const Scalar&                wqw_sfc,
    const uview_1d<const Spack>& rho_zt,
    const uview_1d<const Spack>& tke,
    const uview_1d<const Spack>& pint,
    const uview_1d<Spack>&       rho_zi,
    const uview_1d<Spack>&       host_dse);

  KOKKOS_FUNCTION
  static void compute_shoc_vapor(
    const MemberType&            team,
    const Int&                   nlev,
    const uview_1d<const Spack>& qw,
    const uview_1d<const Spack>& ql,
    const uview_1d<Spack>&       qv);

  KOKKOS_FUNCTION
  static void update_prognostics_implicit(
    const MemberType&            team,
    const Int&                   nlev,
    const Int&                   nlevi,
    const Int&                   num_tracer,
    const Scalar&                dtime,
    const uview_1d<const Spack>& dz_zt,
    const uview_1d<const Spack>& dz_zi,
    const uview_1d<const Spack>& rho_zt,
    const uview_1d<const Spack>& zt_grid,
    const uview_1d<const Spack>& zi_grid,
    const uview_1d<const Spack>& tk,
    const uview_1d<const Spack>& tkh,
    const Scalar&                uw_sfc,
    const Scalar&                vw_sfc,
    const Scalar&                wthl_sfc,
    const Scalar&                wqw_sfc,
    const uview_1d<const Spack>& wtracer_sfc,
    const uview_1d<Spack>&       rdp_zt,
    const uview_1d<Spack>&       tmpi,
    const uview_1d<Spack>&       tkh_zi,
    const uview_1d<Spack>&       tk_zi,
    const uview_1d<Spack>&       rho_zi,
    const uview_1d<Scalar>&       du,
    const uview_1d<Scalar>&       dl,
    const uview_1d<Scalar>&       d,
    const uview_2d<Spack>&       X1,
    const uview_1d<Spack>&       thetal,
    const uview_1d<Spack>&       qw,
    const uview_2d<Spack>&       tracer,
    const uview_1d<Spack>&       tke,
    const uview_1d<Spack>&       u_wind,
    const uview_1d<Spack>&       v_wind);

  KOKKOS_FUNCTION
  static void diag_third_shoc_moments(
    const MemberType&            team,
    const Int&                   nlev,
    const Int&                   nlevi,
    const uview_1d<const Spack>& w_sec,
    const uview_1d<const Spack>& thl_sec,
    const uview_1d<const Spack>& wthl_sec,
    const uview_1d<const Spack>& isotropy,
    const uview_1d<const Spack>& brunt,
    const uview_1d<const Spack>& thetal,
    const uview_1d<const Spack>& tke,
    const uview_1d<const Spack>& dz_zt,
    const uview_1d<const Spack>& dz_zi,
    const uview_1d<const Spack>& zt_grid,
    const uview_1d<const Spack>& zi_grid,
    const uview_1d<Spack>&       w_sec_zi,
    const uview_1d<Spack>&       isotropy_zi,
    const uview_1d<Spack>&       brunt_zi,
    const uview_1d<Spack>&       thetal_zi,
    const uview_1d<Spack>&       w3);

  KOKKOS_FUNCTION
  static void adv_sgs_tke(
    const MemberType&            team,
    const Int&                   nlev,
    const Real&                  dtime,
    const uview_1d<const Spack>& shoc_mix,
    const uview_1d<const Spack>& wthv_sec,
    const uview_1d<const Spack>& sterm_zt,
    const uview_1d<const Spack>& tk,
    const uview_1d<Spack>&       tke,
    const uview_1d<Spack>&       a_diss);

  KOKKOS_FUNCTION
  static void shoc_assumed_pdf(
    const MemberType&            team,
    const Int&                   nlev,
    const Int&                   nlevi,
    const uview_1d<const Spack>& thetal,
    const uview_1d<const Spack>& qw,
    const uview_1d<const Spack>& w_field,
    const uview_1d<const Spack>& thl_sec,
    const uview_1d<const Spack>& qw_sec,
    const uview_1d<const Spack>& wthl_sec,
    const uview_1d<const Spack>& w_sec,
    const uview_1d<const Spack>& wqw_sec,
    const uview_1d<const Spack>& qwthl_sec,
    const uview_1d<const Spack>& w3,
    const uview_1d<const Spack>& pres,
    const uview_1d<const Spack>& zt_grid,
    const uview_1d<const Spack>& zi_grid,
    const uview_1d<Spack>&       wthl_sec_zt,
    const uview_1d<Spack>&       wqw_sec_zt,
    const uview_1d<Spack>&       w3_zt,
    const uview_1d<Spack>&       thl_sec_zt,
    const uview_1d<Spack>&       qwthl_sec_zt,
    const uview_1d<Spack>&       qw_sec_zt,
    const uview_1d<Spack>&       shoc_cldfrac,
    const uview_1d<Spack>&       shoc_ql,
    const uview_1d<Spack>&       wqls,
    const uview_1d<Spack>&       wthv_sec,
    const uview_1d<Spack>&       shoc_ql2);

  KOKKOS_FUNCTION
  static void compute_shr_prod(
    const MemberType&            team,
    const Int&                   nlevi,
    const Int&                   nlev,
    const uview_1d<const Spack>& dz_zi,
    const uview_1d<const Spack>& u_wind,
    const uview_1d<const Spack>& v_wind,
    const uview_1d<Spack>&       sterm);

  KOKKOS_FUNCTION
  static void compute_tmpi(
    const MemberType&            team,
    const Int&                   nlevi,
    const Scalar&                dtime,
    const uview_1d<const Spack>& rho_zi,
    const uview_1d<const Spack>& dz_zi,
    const uview_1d<Spack>&       tmpi);

  KOKKOS_FUNCTION
  static void integ_column_stability(
    const MemberType&            team,
    const Int&                   nlev,
    const uview_1d<const Spack>& dz_zt,
    const uview_1d<const Spack>& pres,
    const uview_1d<const Spack>& brunt,
    Scalar&                      brunt_int);

  KOKKOS_FUNCTION
  static void isotropic_ts(
    const MemberType&            team,
    const Int&                   nlev,
    const Scalar&                brunt_int,
    const uview_1d<const Spack>& tke,
    const uview_1d<const Spack>& a_diss,
    const uview_1d<const Spack>& brunt,
    const uview_1d<Spack>&       isotropy);

  KOKKOS_FUNCTION
  static void dp_inverse(
    const MemberType&            team,
    const Int&                   nlev,
    const uview_1d<const Spack>& rho_zt,
    const uview_1d<const Spack>& dz_zt,
    const uview_1d<Spack>&       rdp_zt);

  KOKKOS_FUNCTION
  static void shoc_main(const Int& shcol, const Int& nlev, const Int& nlevi, const Spack& dtime, const Int& nadv, const uview_1d<const Spack>& host_dx, const uview_1d<const Spack>& host_dy, const uview_1d<const Spack>& thv, const uview_1d<const Spack>& zt_grid, const uview_1d<const Spack>& zi_grid, const uview_1d<const Spack>& pres, const uview_1d<const Spack>& presi, const uview_1d<const Spack>& pdel, const uview_1d<const Spack>& wthl_sfc, const uview_1d<const Spack>& wqw_sfc, const uview_1d<const Spack>& uw_sfc, const uview_1d<const Spack>& vw_sfc, const uview_1d<const Spack>& wtracer_sfc, const Int& num_qtracers, const uview_1d<const Spack>& w_field, const uview_1d<const Spack>& exner, const uview_1d<const Spack>& phis, const uview_1d<Spack>& host_dse, const uview_1d<Spack>& tke, const uview_1d<Spack>& thetal, const uview_1d<Spack>& qw, const uview_1d<Spack>& u_wind, const uview_1d<Spack>& v_wind, const uview_1d<Spack>& qtracers, const uview_1d<Spack>& wthv_sec, const uview_1d<Spack>& tkh, const uview_1d<Spack>& tk, const uview_1d<Spack>& shoc_ql, const uview_1d<Spack>& shoc_cldfrac, const uview_1d<Spack>& pblh, const uview_1d<Spack>& shoc_mix, const uview_1d<Spack>& isotropy, const uview_1d<Spack>& w_sec, const uview_1d<Spack>& thl_sec, const uview_1d<Spack>& qw_sec, const uview_1d<Spack>& qwthl_sec, const uview_1d<Spack>& wthl_sec, const uview_1d<Spack>& wqw_sec, const uview_1d<Spack>& wtke_sec, const uview_1d<Spack>& uw_sec, const uview_1d<Spack>& vw_sec, const uview_1d<Spack>& w3, const uview_1d<Spack>& wqls_sec, const uview_1d<Spack>& brunt, const uview_1d<Spack>& shoc_ql2);

  KOKKOS_FUNCTION
  static void pblintd_height(
    const MemberType& team,
    const Int& nlev,
    const Int& npbl,
    const uview_1d<const Spack>& z,
    const uview_1d<const Spack>& u,
    const uview_1d<const Spack>& v,
    const Scalar& ustar,
    const uview_1d<const Spack>& thv,
    const Scalar& thv_ref,
    Scalar& pblh,
    const uview_1d<Spack>& rino,
    bool& check);
  
  KOKKOS_FUNCTION
  static void vd_shoc_decomp(
    const MemberType&            team,
    const Int&                   nlev,
    const uview_1d<const Spack>& kv_term,
    const uview_1d<const Spack>& tmpi,
    const uview_1d<const Spack>& rdp_zt,
    const Scalar&                dtime,
    const Scalar&                flux,
    const uview_1d<Scalar>&      du,
    const uview_1d<Scalar>&      dl,
    const uview_1d<Scalar>&      d);

  KOKKOS_FUNCTION
  static void vd_shoc_solve(
    const MemberType&       team,
    const uview_1d<Scalar>& du,
    const uview_1d<Scalar>& dl,
    const uview_1d<Scalar>& d,
    const uview_2d<Spack>&  var);

  KOKKOS_FUNCTION
  static void pblintd_surf_temp(const Int& nlev, const Int& nlevi, const Int& npbl,
      const uview_1d<const Spack>& z, const Scalar& ustar,
      const Scalar& obklen, const Scalar& kbfs,
      const uview_1d<const Spack>& thv, Scalar& tlv,
      Scalar& pblh, bool& check, const uview_1d<Spack>& rino);

  KOKKOS_FUNCTION
  static void pblintd_check_pblh(const Int& nlevi, const Int& npbl,
              const uview_1d<const Spack>& z, const Scalar& ustar, const bool& check, Scalar& pblh);

  KOKKOS_FUNCTION
  static void pblintd(
    const MemberType&            team,
    const Int&                   nlev,
    const Int&                   nlevi,
    const Int&                   npbl,
    const uview_1d<const Spack>& z,
    const uview_1d<const Spack>& zi,
    const uview_1d<const Spack>& thl,
    const uview_1d<const Spack>& ql,
    const uview_1d<const Spack>& q,
    const uview_1d<const Spack>& u,
    const uview_1d<const Spack>& v,
    const Scalar&                ustar,
    const Scalar&                obklen,
    const Scalar&                kbfs,
    const uview_1d<const Spack>& cldn,
    const uview_1d<Spack>&       rino,
    const uview_1d<Spack>&       thv,
    Scalar&                      pblh);

  KOKKOS_FUNCTION
  static void shoc_grid(
    const MemberType&            team,
    const Int&                   nlev,
    const Int&                   nlevi,
    const uview_1d<const Spack>& zt_grid,
    const uview_1d<const Spack>& zi_grid,
    const uview_1d<const Spack>& pdel,
    const uview_1d<Spack>&       dz_zt,
    const uview_1d<Spack>&       dz_zi,
    const uview_1d<Spack>&       rho_zt);

  KOKKOS_FUNCTION
  static void eddy_diffusivities(
    const MemberType&            team,
    const Int&                   nlev,
    const Scalar&                obklen,
    const Scalar&                pblh,
    const uview_1d<const Spack>& zt_grid,
    const uview_1d<const Spack>& shoc_mix,
    const uview_1d<const Spack>& sterm_zt,
    const uview_1d<const Spack>& isotropy,
    const uview_1d<const Spack>& tke,
    const uview_1d<Spack>&       tkh,
    const uview_1d<Spack>&       tk);

  KOKKOS_FUNCTION
  static void shoc_tke(
    const MemberType&            team,
    const Int&                   nlev,
    const Int&                   nlevi,
    const Scalar&                dtime,
    const uview_1d<const Spack>& wthv_sec,
    const uview_1d<const Spack>& shoc_mix,
    const uview_1d<const Spack>& dz_zi,
    const uview_1d<const Spack>& dz_zt,
    const uview_1d<const Spack>& pres,
    const uview_1d<const Spack>& u_wind,
    const uview_1d<const Spack>& v_wind,
    const uview_1d<const Spack>& brunt,
    const Scalar&                obklen,
    const uview_1d<const Spack>& zt_grid,
    const uview_1d<const Spack>& zi_grid,
    const Scalar&                pblh,
    const uview_1d<Spack>&       sterm,
    const uview_1d<Spack>&       sterm_zt,
    const uview_1d<Spack>&       a_diss,
    const uview_1d<Spack>&       tke,
    const uview_1d<Spack>&       tk,
    const uview_1d<Spack>&       tkh,
    const uview_1d<Spack>&       isotropy);
}; // struct Functions

} // namespace shoc
} // namespace scream

// If a GPU build, make all code available to the translation unit; otherwise,
// ETI is used.
#ifdef KOKKOS_ENABLE_CUDA
# include "shoc_calc_shoc_varorcovar_impl.hpp"
# include "shoc_calc_shoc_vertflux_impl.hpp"
# include "shoc_diag_second_moments_srf_impl.hpp"
# include "shoc_diag_second_moments_ubycond_impl.hpp"
# include "shoc_update_host_dse_impl.hpp"
# include "shoc_compute_diag_third_shoc_moment_impl.hpp"
# include "shoc_pblintd_init_pot_impl.hpp"
# include "shoc_compute_shoc_mix_shoc_length_impl.hpp"
# include "shoc_check_tke_impl.hpp"
# include "shoc_linear_interp_impl.hpp"
# include "shoc_clipping_diag_third_shoc_moments_impl.hpp"
# include "shoc_energy_integrals_impl.hpp"
# include "shoc_diag_second_moments_lbycond_impl.hpp"
# include "shoc_diag_second_moments_impl.hpp"
# include "shoc_diag_second_shoc_moments_impl.hpp"
# include "shoc_compute_shr_prod_impl.hpp"
# include "shoc_compute_brunt_shoc_length_impl.hpp"
# include "shoc_compute_l_inf_shoc_length_impl.hpp"
# include "shoc_check_length_scale_shoc_length_impl.hpp"
# include "shoc_compute_conv_vel_shoc_length_impl.hpp"
# include "shoc_diag_obklen_impl.hpp"
# include "shoc_pblintd_cldcheck_impl.hpp"
# include "shoc_compute_conv_time_shoc_length_impl.hpp"
# include "shoc_length_impl.hpp"
# include "shoc_energy_fixer_impl.hpp"
# include "shoc_compute_shoc_vapor_impl.hpp"
# include "shoc_update_prognostics_implicit_impl.hpp"
# include "shoc_diag_third_shoc_moments_impl.hpp"
# include "shoc_assumed_pdf_impl.hpp"
# include "shoc_adv_sgs_tke_impl.hpp"
# include "shoc_compute_tmpi_impl.hpp"
# include "shoc_integ_column_stability_impl.hpp"
# include "shoc_isotropic_ts_impl.hpp"
# include "shoc_dp_inverse_impl.hpp"
# include "shoc_shoc_main_impl.hpp"
# include "shoc_pblintd_height_impl.hpp"
# include "shoc_tridiag_solver_impl.hpp"
# include "shoc_pblintd_surf_temp_impl.hpp"
# include "shoc_pblintd_check_pblh_impl.hpp"
# include "shoc_pblintd_impl.hpp"
# include "shoc_grid_impl.hpp"
# include "shoc_eddy_diffusivities_impl.hpp"
# include "shoc_tke_impl.hpp"
#endif // KOKKOS_ENABLE_CUDA

#endif // SHOC_FUNCTIONS_HPP
