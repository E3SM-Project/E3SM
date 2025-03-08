#ifndef SHOC_FUNCTIONS_HPP
#define SHOC_FUNCTIONS_HPP

#include "physics/share/physics_constants.hpp"
#include "physics/shoc/shoc_constants.hpp"

#include "share/eamxx_types.hpp"

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

  using WorkspaceMgr = typename ekat::WorkspaceManager<Spack,  Device>;
  using Workspace    = typename WorkspaceMgr::Workspace;

  // This struct stores runtime options for shoc_main
 struct SHOCRuntime {
   SHOCRuntime() = default;
   // Runtime options for isotropic_ts
   Scalar lambda_low;
   Scalar lambda_high;
   Scalar lambda_slope;
   Scalar lambda_thresh;
   Scalar thl2tune;
   Scalar qw2tune;
   Scalar qwthl2tune;
   Scalar w2tune;
   Scalar length_fac;
   Scalar c_diag_3rd_mom;
   Scalar Ckh;
   Scalar Ckm;
 };

  // This struct stores input views for shoc_main.
  struct SHOCInput {
    SHOCInput() = default;

    // Grid spacing of host model in x direction [m]
    view_1d<const Scalar> dx;
    // grid spacing of host model in y direction [m]
    view_1d<const Scalar> dy;
    // heights, for thermo grid [m]
    view_2d<const Spack>  zt_grid;
    // heights, for interface grid [m]
    view_2d<const Spack>  zi_grid;
    // pressure levels on thermo grid [Pa]
    view_2d<const Spack>  pres;
    // pressure levels on interface grid [Pa]
    view_2d<const Spack>  presi;
    // Differences in pressure levels [Pa]
    view_2d<const Spack>  pdel;
    // virtual potential temperature [K]
    view_2d<const Spack>  thv;
    // large scale vertical velocity [m/s]
    view_2d<const Spack>  w_field;
    // Surface sensible heat flux [K m/s]
    view_1d<const Scalar> wthl_sfc;
    // Surface latent heat flux [kg/kg m/s]
    view_1d<const Scalar> wqw_sfc;
    // Surface momentum flux (u-direction) [m2/s2]
    view_1d<const Scalar> uw_sfc;
    // Surface momentum flux (v-direction) [m2/s2]
    view_1d<const Scalar> vw_sfc;
    // Surface flux for tracers [varies]
    view_2d<const Spack>  wtracer_sfc;
    // Inverse of the exner function [-]
    view_2d<const Spack>  inv_exner;
    // Host model surface geopotential height
    view_1d<const Scalar> phis;
  };

  // This struct stores input/outputs views for shoc_main.
  struct SHOCInputOutput {
    SHOCInputOutput() = default;

    // prognostic temp variable of host model
    // dry static energy [J/kg]
    // dse = Cp*T + g*z + phis
    view_2d<Spack>  host_dse;
    // turbulent kinetic energy [m2/s2]
    view_2d<Spack>  tke;
    // liquid water potential temperature [K]
    view_2d<Spack>  thetal;
    // total water mixing ratio [kg/kg]
    view_2d<Spack>  qw;
    // Vector-valued wind (u,v) [m/s]
    view_3d<Spack>  horiz_wind;
    // buoyancy flux [K m/s]
    view_2d<Spack>  wthv_sec;
    // tracers [varies]
    view_3d<Spack>  qtracers;
    // eddy coefficient for momentum [m2/s]
    view_2d<Spack>  tk;
    // Cloud fraction [-]
    view_2d<Spack>  shoc_cldfrac;
    // cloud liquid mixing ratio [kg/kg]
    view_2d<Spack>  shoc_ql;
  };

  // This struct stores output only views for shoc_main.
  struct SHOCOutput {
    SHOCOutput() = default;

    // planetary boundary layer depth [m]
    view_1d<Scalar> pblh;
    // surface friction velocity [m/s]
    view_1d<Scalar> ustar;
    // Monin Obukhov length [m]
    view_1d<Scalar> obklen;
    // cloud liquid mixing ratio variance [kg^2/kg^2]
    view_2d<Spack>  shoc_ql2;
    // eddy coefficient for heat [m2/s]
    view_2d<Spack>  tkh;
  };

  // This struct stores output views for SHOC diagnostics for shoc_main.
  struct SHOCHistoryOutput {
    SHOCHistoryOutput() = default;

    // Turbulent length scale [m]
    view_2d<Spack>  shoc_mix;
    // vertical velocity variance [m2/s2]
    view_2d<Spack>  w_sec;
    // temperature variance [K^2]
    view_2d<Spack>  thl_sec;
    // moisture variance [kg2/kg2]
    view_2d<Spack>  qw_sec;
    // temp moisture covariance [K kg/kg]
    view_2d<Spack>  qwthl_sec;
    // vertical heat flux [K m/s]
    view_2d<Spack>  wthl_sec;
    // vertical moisture flux [K m/s]
    view_2d<Spack>  wqw_sec;
    // vertical tke flux [m3/s3]
    view_2d<Spack>  wtke_sec;
    // vertical zonal momentum flux [m2/s2]
    view_2d<Spack>  uw_sec;
    // vertical meridional momentum flux [m2/s2]
    view_2d<Spack>  vw_sec;
    // third moment vertical velocity [m3/s3]
    view_2d<Spack>  w3;
    // liquid water flux [kg/kg m/s]
    view_2d<Spack>  wqls_sec;
    // brunt vaisala frequency [s-1]
    view_2d<Spack>  brunt;
    // return to isotropic timescale [s]
    view_2d<Spack>  isotropy;
  };

#ifdef SCREAM_SHOC_SMALL_KERNELS
  struct SHOCTemporaries {
    SHOCTemporaries() = default;

    view_1d<Scalar> se_b;
    view_1d<Scalar> ke_b;
    view_1d<Scalar> wv_b;
    view_1d<Scalar> wl_b;
    view_1d<Scalar> se_a;
    view_1d<Scalar> ke_a;
    view_1d<Scalar> wv_a;
    view_1d<Scalar> wl_a;
    view_1d<Scalar> kbfs;
    view_1d<Scalar> ustar2;
    view_1d<Scalar> wstar;

    view_2d<Spack> rho_zt;
    view_2d<Spack> shoc_qv;
    view_2d<Spack> tabs;
    view_2d<Spack> dz_zt;
    view_2d<Spack> dz_zi;
    view_2d<Spack> tkh;
  };
#endif

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
    const uview_1d<const Spack>& inv_exner,
    const uview_1d<const Spack>& zt_grid,
    const Scalar& phis,
    const uview_1d<Spack>& host_dse);
#ifdef SCREAM_SHOC_SMALL_KERNELS
  static void update_host_dse_disp(
    const Int& shcol,
    const Int& nlev,
    const view_2d<const Spack>& thlm,
    const view_2d<const Spack>& shoc_ql,
    const view_2d<const Spack>& inv_exner,
    const view_2d<const Spack>& zt_grid,
    const view_1d<const Scalar>& phis,
    const view_2d<Spack>& host_dse);
#endif

  KOKKOS_FUNCTION
  static void compute_diag_third_shoc_moment(
    const MemberType& team,
    const Int& nlev,
    const Int& nlevi,
    const Scalar& c_diag_3rd_mom,
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
    const Scalar&                length_fac,
    const uview_1d<const Spack>& tke,
    const uview_1d<const Spack>& brunt,
    const uview_1d<const Spack>& zt_grid,
    const Scalar&                l_inf,
    const uview_1d<Spack>&       shoc_mix);

  KOKKOS_FUNCTION
  static void check_tke(
    const MemberType& team,
    const Int& nlev,
    const uview_1d<Spack>& tke);
#ifdef SCREAM_SHOC_SMALL_KERNELS
  static void check_tke_disp(
    const Int& schol,
    const Int& nlev,
    const view_2d<Spack>& tke);
#endif

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
#ifdef SCREAM_SHOC_SMALL_KERNELS
  static void shoc_energy_integrals_disp(
    const Int&                   shcol,
    const Int&                   nlev,
    const view_2d<const Spack>& host_dse,
    const view_2d<const Spack>& pdel,
    const view_2d<const Spack>& rtm,
    const view_2d<const Spack>& rcm,
    const uview_2d<const Spack>& u_wind,
    const uview_2d<const Spack>& v_wind,
    const view_1d<Scalar>& se_b_slot,
    const view_1d<Scalar>& ke_b_slot,
    const view_1d<Scalar>& wv_b_slot,
    const view_1d<Scalar>& wl_b_slot);
#endif

  KOKKOS_FUNCTION
  static void shoc_diag_second_moments_lbycond(
    const Scalar& wthl_sfc, const Scalar& wqw_sfc, const Scalar& uw_sfc, const Scalar& vw_sfc,
    const Scalar& ustar2, const Scalar& wstar,
    Scalar& wthl_sec, Scalar& wqw_sec, Scalar& uw_sec, Scalar& vw_sec,
    Scalar& wtke_sec, Scalar& thl_sec, Scalar& qw_sec, Scalar& qwthl_sec);

  KOKKOS_FUNCTION
  static void diag_second_moments(const MemberType& team, const Int& nlev, const Int& nlevi,
     const Real& thl2tune, const Real& qw2tune, const Real& qwthl2tune, const Real& w2tune,
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
     const Scalar& thl2tune, const Scalar& qw2tune, const Scalar& qwthl2tune, const Scalar& w2tune,
     const uview_1d<const Spack>& thetal, const uview_1d<const Spack>& qw, const uview_1d<const Spack>& u_wind,
     const uview_1d<const Spack>& v_wind, const uview_1d<const Spack>& tke, const uview_1d<const Spack>& isotropy,
     const uview_1d<const Spack>& tkh, const uview_1d<const Spack>& tk, const uview_1d<const Spack>& dz_zi,
     const uview_1d<const Spack>& zt_grid, const uview_1d<const Spack>& zi_grid, const uview_1d<const Spack>& shoc_mix,
     const Scalar& wthl_sfc, const Scalar& wqw_sfc, const Scalar& uw_sfc, const Scalar& vw_sfc, Scalar& ustar2, Scalar& wstar,
     const Workspace& workspace, const uview_1d<Spack>& thl_sec,
     const uview_1d<Spack>& qw_sec, const uview_1d<Spack>& wthl_sec, const uview_1d<Spack>& wqw_sec, const uview_1d<Spack>& qwthl_sec,
     const uview_1d<Spack>& uw_sec, const uview_1d<Spack>& vw_sec, const uview_1d<Spack>& wtke_sec, const uview_1d<Spack>& w_sec);
#ifdef SCREAM_SHOC_SMALL_KERNELS
  static void diag_second_shoc_moments_disp(
    const Int& shcol, const Int& nlev, const Int& nlevi,
    const Scalar& thl2tune,
    const Scalar& qw2tune,
    const Scalar& qwthl2tune,
    const Scalar& w2tune,
    const view_2d<const Spack>& thetal,
    const view_2d<const Spack>& qw,
    const view_2d<const Spack>& u_wind,
    const view_2d<const Spack>& v_wind,
    const view_2d<const Spack>& tke,
    const view_2d<const Spack>& isotropy,
    const view_2d<const Spack>& tkh,
    const view_2d<const Spack>& tk,
    const view_2d<const Spack>& dz_zi,
    const view_2d<const Spack>& zt_grid,
    const view_2d<const Spack>& zi_grid,
    const view_2d<const Spack>& shoc_mix,
    const view_1d<const Scalar>& wthl_sfc,
    const view_1d<const Scalar>& wqw_sfc,
    const view_1d<const Scalar>& uw_sfc,
    const view_1d<const Scalar>& vw_sfc,
    const view_1d<Scalar>& ustar2,
    const view_1d<Scalar>& wstar,
    const WorkspaceMgr& workspace_mgr,
    const view_2d<Spack>& thl_sec,
    const view_2d<Spack>& qw_sec,
    const view_2d<Spack>& wthl_sec,
    const view_2d<Spack>& wqw_sec,
    const view_2d<Spack>& qwthl_sec,
    const view_2d<Spack>& uw_sec,
    const view_2d<Spack>& vw_sec,
    const view_2d<Spack>& wtke_sec,
    const view_2d<Spack>& w_sec);
#endif

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
    const Scalar&          dx,
    const Scalar&          dy,
    const uview_1d<Spack>& shoc_mix);

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
#ifdef SCREAM_SHOC_SMALL_KERNELS
  static void shoc_diag_obklen_disp(
    const Int&                   shcol,
    const Int&                   nlev,
    const view_1d<const Scalar>& uw_sfc,
    const view_1d<const Scalar>& vw_sfc,
    const view_1d<const Scalar>& wthl_sfc,
    const view_1d<const Scalar>& wqw_sfc,
    const view_2d<const Scalar>& thl_sfc,
    const view_2d<const Scalar>& cldliq_sfc,
    const view_2d<const Scalar>& qv_sfc,
    const view_1d<Scalar>&       ustar,
    const view_1d<Scalar>&       kbfs,
    const view_1d<Scalar>&       obklen);
#endif

  KOKKOS_FUNCTION
  static void shoc_pblintd_cldcheck(
    const Scalar& zi, const Scalar& cldn,
    Scalar& pblh);

  KOKKOS_FUNCTION
  static void shoc_length(
    const MemberType&            team,
    const Int&                   nlev,
    const Int&                   nlevi,
    const Scalar&                length_fac,
    const Scalar&                dx,
    const Scalar&                dy,
    const uview_1d<const Spack>& zt_grid,
    const uview_1d<const Spack>& zi_grid,
    const uview_1d<const Spack>& dz_zt,
    const uview_1d<const Spack>& tke,
    const uview_1d<const Spack>& thv,
    const Workspace&             workspace,
    const uview_1d<Spack>&       brunt,
    const uview_1d<Spack>&       shoc_mix);
#ifdef SCREAM_SHOC_SMALL_KERNELS
  static void shoc_length_disp(
    const Int&                   shcol,
    const Int&                   nlev,
    const Int&                   nlevi,
    const Scalar&                length_fac,
    const view_1d<const Scalar>& dx,
    const view_1d<const Scalar>& dy,
    const view_2d<const Spack>&  zt_grid,
    const view_2d<const Spack>&  zi_grid,
    const view_2d<const Spack>&  dz_zt,
    const view_2d<const Spack>&  tke,
    const view_2d<const Spack>&  thv,
    const WorkspaceMgr&          workspace_mgr,
    const view_2d<Spack>&        brunt,
    const view_2d<Spack>&        shoc_mix);
#endif

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
    const Workspace&             workspace,
    const uview_1d<Spack>&       host_dse);
#ifdef SCREAM_SHOC_SMALL_KERNELS
  static void shoc_energy_fixer_disp(
    const Int&                   shcol,
    const Int&                   nlev,
    const Int&                   nlevi,
    const Scalar&                dtime,
    const Int&                   nadv,
    const view_2d<const Spack>&  zt_grid,
    const view_2d<const Spack>&  zi_grid,
    const view_1d<const Scalar>& se_b,
    const view_1d<const Scalar>& ke_b,
    const view_1d<const Scalar>& wv_b,
    const view_1d<const Scalar>& wl_b,
    const view_1d<const Scalar>& se_a,
    const view_1d<const Scalar>& ke_a,
    const view_1d<const Scalar>& wv_a,
    const view_1d<const Scalar>& wl_a,
    const view_1d<const Scalar>& wthl_sfc,
    const view_1d<const Scalar>& wqw_sfc,
    const view_2d<const Spack>&  rho_zt,
    const view_2d<const Spack>&  tke,
    const view_2d<const Spack>&  pint,
    const WorkspaceMgr&          workspace_mgr,
    const view_2d<Spack>&        host_dse);
#endif

  KOKKOS_FUNCTION
  static void compute_shoc_vapor(
    const MemberType&            team,
    const Int&                   nlev,
    const uview_1d<const Spack>& qw,
    const uview_1d<const Spack>& ql,
    const uview_1d<Spack>&       qv);
#ifdef SCREAM_SHOC_SMALL_KERNELS
  static void compute_shoc_vapor_disp(
    const Int&                  shcol,
    const Int&                  nlev,
    const view_2d<const Spack>& qw,
    const view_2d<const Spack>& ql,
    const view_2d<Spack>&       qv);
#endif

  KOKKOS_FUNCTION
  static void compute_shoc_temperature(
    const MemberType&            team,
    const Int&                   nlev,
    const uview_1d<const Spack>& thetal,
    const uview_1d<const Spack>& ql,
    const uview_1d<const Spack>& inv_exner,
    const uview_1d<Spack>&       tabs);
#ifdef SCREAM_SHOC_SMALL_KERNELS
  static void compute_shoc_temperature_disp(
    const Int&                  shcol,
    const Int&                  nlev,
    const view_2d<const Spack>& thetal,
    const view_2d<const Spack>& ql,
    const view_2d<const Spack>& inv_exner,
    const view_2d<Spack>&       tabs);
#endif

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
    const Workspace&             workspace,
    const uview_1d<Spack>&       thetal,
    const uview_1d<Spack>&       qw,
    const uview_2d<Spack>&       tracer,
    const uview_1d<Spack>&       tke,
    const uview_1d<Spack>&       u_wind,
    const uview_1d<Spack>&       v_wind);
#ifdef SCREAM_SHOC_SMALL_KERNELS
  static void update_prognostics_implicit_disp(
    const Int&                   shcol,
    const Int&                   nlev,
    const Int&                   nlevi,
    const Int&                   num_tracer,
    const Scalar&                dtime,
    const view_2d<const Spack>&  dz_zt,
    const view_2d<const Spack>&  dz_zi,
    const view_2d<const Spack>&  rho_zt,
    const view_2d<const Spack>&  zt_grid,
    const view_2d<const Spack>&  zi_grid,
    const view_2d<const Spack>&  tk,
    const view_2d<const Spack>&  tkh,
    const view_1d<const Scalar>& uw_sfc,
    const view_1d<const Scalar>& vw_sfc,
    const view_1d<const Scalar>& wthl_sfc,
    const view_1d<const Scalar>& wqw_sfc,
    const view_2d<const Spack>&  wtracer_sfc,
    const WorkspaceMgr&          workspace_mgr,
    const view_2d<Spack>&        thetal,
    const view_2d<Spack>&        qw,
    const view_3d<Spack>&        tracer,
    const view_2d<Spack>&        tke,
    const view_2d<Spack>&        u_wind,
    const view_2d<Spack>&        v_wind);
#endif

  KOKKOS_FUNCTION
  static void diag_third_shoc_moments(
    const MemberType&            team,
    const Int&                   nlev,
    const Int&                   nlevi,
    const Scalar&                c_diag_3rd_mom,
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
    const Workspace&             workspace,
    const uview_1d<Spack>&       w3);
#ifdef SCREAM_SHOC_SMALL_KERNELS
  static void diag_third_shoc_moments_disp(
    const Int&                  shcol,
    const Int&                  nlev,
    const Int&                  nlevi,
    const Scalar&               c_diag_3rd_mom,
    const view_2d<const Spack>& w_sec,
    const view_2d<const Spack>& thl_sec,
    const view_2d<const Spack>& wthl_sec,
    const view_2d<const Spack>& isotropy,
    const view_2d<const Spack>& brunt,
    const view_2d<const Spack>& thetal,
    const view_2d<const Spack>& tke,
    const view_2d<const Spack>& dz_zt,
    const view_2d<const Spack>& dz_zi,
    const view_2d<const Spack>& zt_grid,
    const view_2d<const Spack>& zi_grid,
    const WorkspaceMgr&         workspace_mgr,
    const view_2d<Spack>&       w3);
#endif

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
    const Workspace&             workspace,
    const uview_1d<Spack>&       shoc_cldfrac,
    const uview_1d<Spack>&       shoc_ql,
    const uview_1d<Spack>&       wqls,
    const uview_1d<Spack>&       wthv_sec,
    const uview_1d<Spack>&       shoc_ql2);
#ifdef SCREAM_SHOC_SMALL_KERNELS
  static void shoc_assumed_pdf_disp(
    const Int&                  shcol,
    const Int&                  nlev,
    const Int&                  nlevi,
    const view_2d<const Spack>& thetal,
    const view_2d<const Spack>& qw,
    const view_2d<const Spack>& w_field,
    const view_2d<const Spack>& thl_sec,
    const view_2d<const Spack>& qw_sec,
    const view_2d<const Spack>& wthl_sec,
    const view_2d<const Spack>& w_sec,
    const view_2d<const Spack>& wqw_sec,
    const view_2d<const Spack>& qwthl_sec,
    const view_2d<const Spack>& w3,
    const view_2d<const Spack>& pres,
    const view_2d<const Spack>& zt_grid,
    const view_2d<const Spack>& zi_grid,
    const WorkspaceMgr&         workspace_mgr,
    const view_2d<Spack>&       shoc_cldfrac,
    const view_2d<Spack>&       shoc_ql,
    const view_2d<Spack>&       wqls,
    const view_2d<Spack>&       wthv_sec,
    const view_2d<Spack>&       shoc_ql2);
#endif

  KOKKOS_INLINE_FUNCTION
  static void shoc_assumed_pdf_compute_buoyancy_flux(
    const Spack& wthlsec,
    const Spack& wqwsec,
    const Spack& pval,
    const Spack& wqls,
    Spack&       wthv_sec);

  KOKKOS_INLINE_FUNCTION
  static void shoc_assumed_pdf_compute_cloud_liquid_variance(
    const Spack& a,
    const Spack& s1,
    const Spack& ql1,
    const Spack& C1,
    const Spack& std_s1,
    const Spack& s2,
    const Spack& ql2,
    const Spack& C2,
    const Spack& std_s2,
    const Spack& shoc_ql,
    Spack&       shoc_ql2);


  KOKKOS_INLINE_FUNCTION
  static void shoc_assumed_pdf_compute_liquid_water_flux(
    const Spack& a,
    const Spack& w1_1,
    const Spack& w_first,
    const Spack& ql1,
    const Spack& w1_2,
    const Spack& ql2,
    Spack&       wqls);

  KOKKOS_INLINE_FUNCTION
  static void shoc_assumed_pdf_compute_qs(
    const Spack& Tl1_1,
    const Spack& Tl1_2,
    const Spack& pval,
    const Smask& active_entries,
    Spack&       qs1,
    Spack&       beta1,
    Spack&       qs2,
    Spack&       beta2);

  KOKKOS_INLINE_FUNCTION
  static void shoc_assumed_pdf_compute_s(
    const Spack& qw1,
    const Spack& qs,
    const Spack& beta,
    const Spack& pval,
    const Spack& thl2,
    const Spack& qw2,
    const Spack& sqrtthl2,
    const Spack& sqrtqw2,
    const Spack& r_qwthl,
    Spack&       s,
    Spack&       std_s,
    Spack&       qn,
    Spack&       C);

  KOKKOS_INLINE_FUNCTION
  static void shoc_assumed_pdf_compute_sgs_liquid(
    const Spack& a,
    const Spack& ql1,
    const Spack& ql2,
    Spack&       shoc_ql);

  KOKKOS_INLINE_FUNCTION
  static void shoc_assumed_pdf_compute_temperature(
    const Spack& thl1,
    const Spack& pval,
    Spack&       Tl1);

  KOKKOS_INLINE_FUNCTION
  static void shoc_assumed_pdf_inplume_correlations(
    const Spack& sqrtqw2_1,
    const Spack& sqrtthl2_1,
    const Spack& a,
    const Spack& sqrtqw2_2,
    const Spack& sqrtthl2_2,
    const Spack& qwthlsec,
    const Spack& qw1_1,
    const Spack& qw_first,
    const Spack& thl1_1,
    const Spack& thl_first,
    const Spack& qw1_2,
    const Spack& thl1_2,
    Spack&       r_qwthl_1);

  KOKKOS_INLINE_FUNCTION
  static void shoc_assumed_pdf_qw_parameters(
    const Spack& wqwsec,
    const Spack& sqrtw2,
    const Spack& Skew_w,
    const Spack& sqrtqt,
    const Spack& qwsec,
    const Spack& w1_2,
    const Spack& w1_1,
    const Spack& qw_first,
    const Spack& a,
    const Scalar rt_tol,
    const Scalar w_thresh,
    Spack&       qw1_1,
    Spack&       qw1_2,
    Spack&       qw2_1,
    Spack&       qw2_2,
    Spack&       sqrtqw2_1,
    Spack&       sqrtqw2_2);

  KOKKOS_INLINE_FUNCTION
  static void shoc_assumed_pdf_thl_parameters(
    const Spack& wthlsec,
    const Spack& sqrtw2,
    const Spack& sqrtthl,
    const Spack& thlsec,
    const Spack& thl_first,
    const Spack& w1_1,
    const Spack& w1_2,
    const Spack& Skew_w,
    const Spack& a,
    const Scalar thl_tol,
    const Scalar w_thresh,
    Spack&       thl1_1,
    Spack&       thl1_2,
    Spack&       thl2_1,
    Spack&       thl2_2,
    Spack&       sqrtthl2_1,
    Spack&       sqrtthl2_2);

  KOKKOS_INLINE_FUNCTION
  static void shoc_assumed_pdf_tilde_to_real(
    const Spack& w_first,
    const Spack& sqrtw2,
    Spack&       w1);

  KOKKOS_INLINE_FUNCTION
  static void shoc_assumed_pdf_vv_parameters(
    const Spack& w_first,
    const Spack& w_sec,
    const Spack& w3var,
    const Scalar w_tol_sqd,
    Spack&       Skew_w,
    Spack&       w1_1,
    Spack&       w1_2,
    Spack&       w2_1,
    Spack&       w2_2,
    Spack&       a);

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
    const Scalar&                lambda_low,
    const Scalar&                lambda_high,
    const Scalar&                lambda_slope,
    const Scalar&                lambda_thresh,
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

  static Int shoc_init(
    const Int&                  nbot_shoc,
    const Int&                  ntop_shoc,
    const view_1d<const Spack>& pref_mid);

#ifndef SCREAM_SHOC_SMALL_KERNELS
  KOKKOS_FUNCTION
  static void shoc_main_internal(
    const MemberType&            team,
    const Int&                   nlev,         // Number of levels
    const Int&                   nlevi,        // Number of levels on interface grid
    const Int&                   npbl,         // Maximum number of levels in pbl from surface
    const Int&                   nadv,         // Number of times to loop SHOC
    const Int&                   num_qtracers, // Number of tracers
    const Scalar&                dtime,        // SHOC timestep [s]
    // Runtime Parameters
    const Scalar&                lambda_low,
    const Scalar&                lambda_high,
    const Scalar&                lambda_slope,
    const Scalar&                lambda_thresh,
    const Scalar&                thl2tune,
    const Scalar&                qw2tune,
    const Scalar&                qwthl2tune,
    const Scalar&                w2tune,
    const Scalar&                length_fac,
    const Scalar&                c_diag_3rd_mom,
    const Scalar&                Ckh,
    const Scalar&                Ckm,
    // Input Variables
    const Scalar&                host_dx,
    const Scalar&                host_dy,
    const uview_1d<const Spack>& zt_grid,
    const uview_1d<const Spack>& zi_grid,
    const uview_1d<const Spack>& pres,
    const uview_1d<const Spack>& presi,
    const uview_1d<const Spack>& pdel,
    const uview_1d<const Spack>& thv,
    const uview_1d<const Spack>& w_field,
    const Scalar&                wthl_sfc,
    const Scalar&                wqw_sfc,
    const Scalar&                uw_sfc,
    const Scalar&                vw_sfc,
    const uview_1d<const Spack>& wtracer_sfc,
    const uview_1d<const Spack>& inv_exner,
    const Scalar&                phis,
    // Local Workspace
    const Workspace&             workspace,
    // Input/Output Variables
    const uview_1d<Spack>&       host_dse,
    const uview_1d<Spack>&       tke,
    const uview_1d<Spack>&       thetal,
    const uview_1d<Spack>&       qw,
    const uview_1d<Spack>&       u_wind,
    const uview_1d<Spack>&       v_wind,
    const uview_1d<Spack>&       wthv_sec,
    const uview_2d<Spack>&       qtracers,
    const uview_1d<Spack>&       tk,
    const uview_1d<Spack>&       shoc_cldfrac,
    const uview_1d<Spack>&       shoc_ql,
    // Output Variables
    Scalar&                      pblh,
    Scalar&                      ustar,
    Scalar&                      obklen,
    const uview_1d<Spack>&       shoc_ql2,
    const uview_1d<Spack>&       tkh,
    // Diagnostic Output Variables
    const uview_1d<Spack>&       shoc_mix,
    const uview_1d<Spack>&       w_sec,
    const uview_1d<Spack>&       thl_sec,
    const uview_1d<Spack>&       qw_sec,
    const uview_1d<Spack>&       qwthl_sec,
    const uview_1d<Spack>&       wthl_sec,
    const uview_1d<Spack>&       wqw_sec,
    const uview_1d<Spack>&       wtke_sec,
    const uview_1d<Spack>&       uw_sec,
    const uview_1d<Spack>&       vw_sec,
    const uview_1d<Spack>&       w3,
    const uview_1d<Spack>&       wqls_sec,
    const uview_1d<Spack>&       brunt,
    const uview_1d<Spack>&       isotropy);
#else
  static void shoc_main_internal(
    const Int&                   shcol,        // Number of columns
    const Int&                   nlev,         // Number of levels
    const Int&                   nlevi,        // Number of levels on interface grid
    const Int&                   npbl,         // Maximum number of levels in pbl from surface
    const Int&                   nadv,         // Number of times to loop SHOC
    const Int&                   num_qtracers, // Number of tracers
    const Scalar&                dtime,        // SHOC timestep [s]
    // Runtime Parameters
    const Scalar&                lambda_low,
    const Scalar&                lambda_high,
    const Scalar&                lambda_slope,
    const Scalar&                lambda_thresh,
    const Scalar&                thl2tune,
    const Scalar&                qw2tune,
    const Scalar&                qwthl2tune,
    const Scalar&                w2tune,
    const Scalar&                length_fac,
    const Scalar&                c_diag_3rd_mom,
    const Scalar&                Ckh,
    const Scalar&                Ckm,
    // Input Variables
    const view_1d<const Scalar>& host_dx,
    const view_1d<const Scalar>& host_dy,
    const view_2d<const Spack>& zt_grid,
    const view_2d<const Spack>& zi_grid,
    const view_2d<const Spack>& pres,
    const view_2d<const Spack>& presi,
    const view_2d<const Spack>& pdel,
    const view_2d<const Spack>& thv,
    const view_2d<const Spack>& w_field,
    const view_1d<const Scalar>& wthl_sfc,
    const view_1d<const Scalar>& wqw_sfc,
    const view_1d<const Scalar>& uw_sfc,
    const view_1d<const Scalar>& vw_sfc,
    const view_2d<const Spack>& wtracer_sfc,
    const view_2d<const Spack>& inv_exner,
    const view_1d<const Scalar>& phis,
    // Workspace Manager
    WorkspaceMgr&               workspace_mgr,
    // Input/Output Variables
    const view_2d<Spack>&       host_dse,
    const view_2d<Spack>&       tke,
    const view_2d<Spack>&       thetal,
    const view_2d<Spack>&       qw,
    const uview_2d<Spack>&      u_wind,
    const uview_2d<Spack>&      v_wind,
    const view_2d<Spack>&       wthv_sec,
    const view_3d<Spack>&       qtracers,
    const view_2d<Spack>&       tk,
    const view_2d<Spack>&       shoc_cldfrac,
    const view_2d<Spack>&       shoc_ql,
    // Output Variables
    const view_1d<Scalar>&      pblh,
    const view_1d<Scalar>&      ustar,
    const view_1d<Scalar>&      obklen,
    const view_2d<Spack>&       shoc_ql2,
    const view_2d<Spack>&       tkh,
    // Diagnostic Output Variables
    const view_2d<Spack>&       shoc_mix,
    const view_2d<Spack>&       w_sec,
    const view_2d<Spack>&       thl_sec,
    const view_2d<Spack>&       qw_sec,
    const view_2d<Spack>&       qwthl_sec,
    const view_2d<Spack>&       wthl_sec,
    const view_2d<Spack>&       wqw_sec,
    const view_2d<Spack>&       wtke_sec,
    const view_2d<Spack>&       uw_sec,
    const view_2d<Spack>&       vw_sec,
    const view_2d<Spack>&       w3,
    const view_2d<Spack>&       wqls_sec,
    const view_2d<Spack>&       brunt,
    const view_2d<Spack>&       isotropy,
    // Temporaries
    const view_1d<Scalar>& se_b,
    const view_1d<Scalar>& ke_b,
    const view_1d<Scalar>& wv_b,
    const view_1d<Scalar>& wl_b,
    const view_1d<Scalar>& se_a,
    const view_1d<Scalar>& ke_a,
    const view_1d<Scalar>& wv_a,
    const view_1d<Scalar>& wl_a,
    const view_1d<Scalar>& kbfs,
    const view_1d<Scalar>& ustar2,
    const view_1d<Scalar>& wstar,
    const view_2d<Spack>& rho_zt,
    const view_2d<Spack>& shoc_qv,
    const view_2d<Spack>& tabs,
    const view_2d<Spack>& dz_zt,
    const view_2d<Spack>& dz_zi);
#endif

  // Return microseconds elapsed
  static Int shoc_main(
    const Int&               shcol,                // Number of SHOC columns in the array
    const Int&               nlev,                 // Number of levels
    const Int&               nlevi,                // Number of levels on interface grid
    const Int&               npbl,                 // Maximum number of levels in pbl from surface
    const Int&               nadv,                 // Number of times to loop SHOC
    const Int&               num_q_tracers,        // Number of tracers
    const Scalar&            dtime,                // SHOC timestep [s]
    WorkspaceMgr&            workspace_mgr,        // WorkspaceManager for local variables
    const SHOCRuntime&       shoc_runtime,         // Runtime options
    const SHOCInput&         shoc_input,           // Input
    const SHOCInputOutput&   shoc_input_output,    // Input/Output
    const SHOCOutput&        shoc_output,          // Output
    const SHOCHistoryOutput& shoc_history_output   // Output (diagnostic)
#ifdef SCREAM_SHOC_SMALL_KERNELS
    , const SHOCTemporaries& shoc_temporaries      // Temporaries for small kernels
#endif
                       );

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
    const Workspace&             workspace,
    Scalar&                      pblh);
#ifdef SCREAM_SHOC_SMALL_KERNELS
  static void pblintd_disp(
    const Int&                   shcol,
    const Int&                   nlev,
    const Int&                   nlevi,
    const Int&                   npbl,
    const view_2d<const Spack>&  z,
    const view_2d<const Spack>&  zi,
    const view_2d<const Spack>&  thl,
    const view_2d<const Spack>&  ql,
    const view_2d<const Spack>&  q,
    const view_2d<const Spack>&  u,
    const view_2d<const Spack>&  v,
    const view_1d<const Scalar>& ustar,
    const view_1d<const Scalar>& obklen,
    const view_1d<const Scalar>& kbfs,
    const view_2d<const Spack>&  cldn,
    const WorkspaceMgr&          workspace_mgr,
    const view_1d<Scalar>&       pblh);
#endif

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
#ifdef SCREAM_SHOC_SMALL_KERNELS
  static void shoc_grid_disp(
    const Int&                  shcol,
    const Int&                  nlev,
    const Int&                  nlevi,
    const view_2d<const Spack>& zt_grid,
    const view_2d<const Spack>& zi_grid,
    const view_2d<const Spack>& pdel,
    const view_2d<Spack>&       dz_zt,
    const view_2d<Spack>&       dz_zi,
    const view_2d<Spack>&       rho_zt);
#endif

  KOKKOS_FUNCTION
  static void eddy_diffusivities(
    const MemberType&            team,
    const Int&                   nlev,
    const Scalar&                Ckh,
    const Scalar&                Ckm,
    const Scalar&                pblh,
    const uview_1d<const Spack>& zt_grid,
    const uview_1d<const Spack>& tabs,
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
    const Scalar&                lambda_low,
    const Scalar&                lambda_high,
    const Scalar&                lambda_slope,
    const Scalar&                lambda_thresh,
    const Scalar&                Ckh,
    const Scalar&                Ckm,
    const uview_1d<const Spack>& wthv_sec,
    const uview_1d<const Spack>& shoc_mix,
    const uview_1d<const Spack>& dz_zi,
    const uview_1d<const Spack>& dz_zt,
    const uview_1d<const Spack>& pres,
    const uview_1d<const Spack>& tabs,
    const uview_1d<const Spack>& u_wind,
    const uview_1d<const Spack>& v_wind,
    const uview_1d<const Spack>& brunt,
    const uview_1d<const Spack>& zt_grid,
    const uview_1d<const Spack>& zi_grid,
    const Scalar&                pblh,
    const Workspace&             workspace,
    const uview_1d<Spack>&       tke,
    const uview_1d<Spack>&       tk,
    const uview_1d<Spack>&       tkh,
    const uview_1d<Spack>&       isotropy);
#ifdef SCREAM_SHOC_SMALL_KERNELS
  static void shoc_tke_disp(
    const Int&                   shcol,
    const Int&                   nlev,
    const Int&                   nlevi,
    const Scalar&                dtime,
    const Scalar&                lambda_low,
    const Scalar&                lambda_high,
    const Scalar&                lambda_slope,
    const Scalar&                lambda_thresh,
    const Scalar&                Ckh,
    const Scalar&                Ckm,
    const view_2d<const Spack>&  wthv_sec,
    const view_2d<const Spack>&  shoc_mix,
    const view_2d<const Spack>&  dz_zi,
    const view_2d<const Spack>&  dz_zt,
    const view_2d<const Spack>&  pres,
    const view_2d<const Spack>&  tabs,
    const view_2d<const Spack>&  u_wind,
    const view_2d<const Spack>&  v_wind,
    const view_2d<const Spack>&  brunt,
    const view_2d<const Spack>&  zt_grid,
    const view_2d<const Spack>&  zi_grid,
    const view_1d<const Scalar>& pblh,
    const WorkspaceMgr&          workspace_mgr,
    const view_2d<Spack>&        tke,
    const view_2d<Spack>&        tk,
    const view_2d<Spack>&        tkh,
    const view_2d<Spack>&        isotropy);
#endif
}; // struct Functions

} // namespace shoc
} // namespace scream

// If a GPU build, without relocatable device code enabled, make all code available
// to the translation unit; otherwise, ETI is used.
#if defined(EAMXX_ENABLE_GPU) && !defined(KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE)  \
                                && !defined(KOKKOS_ENABLE_HIP_RELOCATABLE_DEVICE_CODE)

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
# include "shoc_diag_obklen_impl.hpp"
# include "shoc_pblintd_cldcheck_impl.hpp"
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
# include "shoc_main_impl.hpp"
# include "shoc_pblintd_height_impl.hpp"
# include "shoc_tridiag_solver_impl.hpp"
# include "shoc_pblintd_surf_temp_impl.hpp"
# include "shoc_pblintd_check_pblh_impl.hpp"
# include "shoc_pblintd_impl.hpp"
# include "shoc_grid_impl.hpp"
# include "shoc_eddy_diffusivities_impl.hpp"
# include "shoc_tke_impl.hpp"
# include "shoc_compute_shoc_temperature_impl.hpp"

#endif // GPU && !KOKKOS_ENABLE_*_RELOCATABLE_DEVICE_CODE

// Some functions should be inlined, thus do not use ETI
# include "shoc_assumed_pdf_compute_buoyancy_flux_impl.hpp"
# include "shoc_assumed_pdf_compute_cloud_liquid_variance_impl.hpp"
# include "shoc_assumed_pdf_compute_liquid_water_flux_impl.hpp"
# include "shoc_assumed_pdf_compute_qs_impl.hpp"
# include "shoc_assumed_pdf_compute_s_impl.hpp"
# include "shoc_assumed_pdf_compute_sgs_liquid_impl.hpp"
# include "shoc_assumed_pdf_compute_temperature_impl.hpp"
# include "shoc_assumed_pdf_inplume_correlations_impl.hpp"
# include "shoc_assumed_pdf_qw_parameters_impl.hpp"
# include "shoc_assumed_pdf_compute_s_impl.hpp"
# include "shoc_assumed_pdf_thl_parameters_impl.hpp"
# include "shoc_assumed_pdf_tilde_to_real_impl.hpp"
# include "shoc_assumed_pdf_vv_parameters_impl.hpp"

#endif // SHOC_FUNCTIONS_HPP
