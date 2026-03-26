#ifndef SHOC_FUNCTIONS_HPP
#define SHOC_FUNCTIONS_HPP

#include "share/physics/physics_constants.hpp"
#include "physics/shoc/shoc_constants.hpp"

#include "share/core/eamxx_types.hpp"

#include <ekat_pack_kokkos.hpp>
#include <ekat_workspace.hpp>

namespace scream
{
namespace shoc
{

/*
 * Functions is a stateless struct used to encapsulate a
 * number of functions for SHOC. We use the ETI pattern for
 * these functions.
 *
 * SHOC assumptions:
 *  - Kokkos team policies have a vector length of 1
 */

template <typename ScalarT, typename DeviceT> struct Functions {
  //
  // ------- Types --------
  //

  using Scalar = ScalarT;
  using Device = DeviceT;

  using Pack         = ekat::Pack<Scalar, SCREAM_PACK_SIZE>;
  using IntPack = ekat::Pack<Int, SCREAM_PACK_SIZE>;

  using Mask = ekat::Mask<Pack::n>;

  using KT = ekat::KokkosTypes<Device>;

  using C  = physics::Constants<Scalar>;
  using SC = shoc::Constants<Scalar>;

  template <typename S> using view_1d = typename KT::template view_1d<S>;
  template <typename S> using view_2d = typename KT::template view_2d<S>;
  template <typename S> using view_3d = typename KT::template view_3d<S>;

  template <typename S> using view_2d_strided = typename KT::template sview<S **>;
  template <typename S> using view_3d_strided = typename KT::template sview<S ***>;

  template <typename S, int N>
  using view_1d_ptr_array = typename KT::template view_1d_ptr_carray<S, N>;

  template <typename S> using uview_1d = typename ekat::template Unmanaged<view_1d<S>>;

  template <typename S> using uview_2d = typename ekat::template Unmanaged<view_2d<S>>;

  template <typename S>
  using uview_2d_strided = typename ekat::template Unmanaged<view_2d_strided<S>>;

  using MemberType = typename KT::MemberType;

  using WorkspaceMgr = typename ekat::WorkspaceManager<Pack, Device>;
  using Workspace    = typename WorkspaceMgr::Workspace;

  // This struct stores runtime options for shoc_main
  struct SHOCRuntime {
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
    bool shoc_1p5tke;
    bool extra_diags;
    bool do_3d_turb;
  };

  // This struct stores input views for shoc_main.
  struct SHOCInput {
    // Grid spacing of host model in x direction [m]
    view_1d<const Scalar> dx;
    // grid spacing of host model in y direction [m]
    view_1d<const Scalar> dy;
    // heights, for thermo grid [m]
    view_2d<const Pack> zt_grid;
    // heights, for interface grid [m]
    view_2d<const Pack> zi_grid;
    // pressure levels on thermo grid [Pa]
    view_2d<const Pack> pres;
    // pressure levels on interface grid [Pa]
    view_2d<const Pack> presi;
    // Differences in pressure levels [Pa]
    view_2d<const Pack> pdel;
    // virtual potential temperature [K]
    view_2d<const Pack> thv;
    // large scale vertical velocity [m/s]
    view_2d<const Pack> w_field;
    // Surface sensible heat flux [K m/s]
    view_1d<const Scalar> wthl_sfc;
    // Surface latent heat flux [kg/kg m/s]
    view_1d<const Scalar> wqw_sfc;
    // Surface momentum flux (u-direction) [m2/s2]
    view_1d<const Scalar> uw_sfc;
    // Surface momentum flux (v-direction) [m2/s2]
    view_1d<const Scalar> vw_sfc;
    // Surface flux for tracers [varies]
    view_2d<const Pack> wtracer_sfc;
    // Inverse of the exner function [-]
    view_2d<const Pack> inv_exner;
    // Host model surface geopotential height
    view_1d<const Scalar> phis;
    // 3D strain term for shear production of TKE [/s2]
    view_2d<const Pack> strain2;
  };

  // This struct stores input/outputs views for shoc_main.
  struct SHOCInputOutput {
    // prognostic temp variable of host model
    // dry static energy [J/kg]
    // dse = Cp*T + g*z + phis
    view_2d<Pack> host_dse;
    // turbulent kinetic energy [m2/s2]
    view_2d<Pack> tke;
    // liquid water potential temperature [K]
    view_2d<Pack> thetal;
    // total water mixing ratio [kg/kg]
    view_2d<Pack> qw;
    // Vector-valued wind (u,v) [m/s]
    view_3d<Pack> horiz_wind;
    // buoyancy flux [K m/s]
    view_2d<Pack> wthv_sec;
    // tracers [varies]
    view_3d_strided<Pack> qtracers;
    // eddy coefficient for momentum [m2/s]
    view_2d<Pack> tk;
    // Cloud fraction [-]
    view_2d<Pack> shoc_cldfrac;
    // cloud liquid mixing ratio [kg/kg]
    view_2d<Pack> shoc_ql;
  };

  // This struct stores output only views for shoc_main.
  struct SHOCOutput {
    // planetary boundary layer depth [m]
    view_1d<Scalar> pblh;
    // surface friction velocity [m/s]
    view_1d<Scalar> ustar;
    // Monin Obukhov length [m]
    view_1d<Scalar> obklen;
    // cloud liquid mixing ratio variance [kg^2/kg^2]
    view_2d<Pack> shoc_ql2;
    // eddy coefficient for heat [m2/s]
    view_2d<Pack> tkh;
  };

  // This struct stores output views for SHOC diagnostics for shoc_main.
  struct SHOCHistoryOutput {
    // Turbulent length scale [m]
    view_2d<Pack> shoc_mix;
    // vertical velocity variance [m2/s2]
    view_2d<Pack> w_sec;
    // temperature variance [K^2]
    view_2d<Pack> thl_sec;
    // moisture variance [kg2/kg2]
    view_2d<Pack> qw_sec;
    // temp moisture covariance [K kg/kg]
    view_2d<Pack> qwthl_sec;
    // vertical heat flux [K m/s]
    view_2d<Pack> wthl_sec;
    // vertical moisture flux [K m/s]
    view_2d<Pack> wqw_sec;
    // vertical tke flux [m3/s3]
    view_2d<Pack> wtke_sec;
    // vertical zonal momentum flux [m2/s2]
    view_2d<Pack> uw_sec;
    // vertical meridional momentum flux [m2/s2]
    view_2d<Pack> vw_sec;
    // third moment vertical velocity [m3/s3]
    view_2d<Pack> w3;
    // liquid water flux [kg/kg m/s]
    view_2d<Pack> wqls_sec;
    // brunt vaisala frequency [s-1]
    view_2d<Pack> brunt;
    // return to isotropic timescale [s]
    view_2d<Pack> isotropy;
    // shoc condensation kg/kg/s
    view_2d<Pack> shoc_cond;
    // shoc evaporation kg/kg/s
    view_2d<Pack> shoc_evap;
  };

#ifdef SCREAM_SHOC_SMALL_KERNELS
  struct SHOCTemporaries {
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

    view_2d<Pack> rho_zt;
    view_2d<Pack> shoc_qv;
    view_2d<Pack> tabs;
    view_2d<Pack> dz_zt;
    view_2d<Pack> dz_zi;
    view_2d<Pack> tkh;
  };
#endif

  //
  // --------- Functions ---------
  //
  KOKKOS_FUNCTION
  static void calc_shoc_varorcovar(const MemberType &team, const Int &nlev, const Scalar &tunefac,
                                   const uview_1d<const Pack> &isotropy_zi,
                                   const uview_1d<const Pack> &tkh_zi,
                                   const uview_1d<const Pack> &dz_zi,
                                   const uview_1d<const Pack> &invar1,
                                   const uview_1d<const Pack> &invar2,
                                   const uview_1d<Pack> &varorcovar);

  KOKKOS_FUNCTION
  static void calc_shoc_vertflux(const MemberType &team, const Int &nlev,
                                 const uview_1d<const Pack> &tkh_zi,
                                 const uview_1d<const Pack> &dz_zi,
                                 const uview_1d<const Pack> &invar,
                                 const uview_1d<Pack> &vertflux);

  KOKKOS_FUNCTION
  static void shoc_diag_second_moments_srf(const Scalar &wthl_sfc, const Scalar &uw_sfc,
                                           const Scalar &vw_sfc, Scalar &ustar2, Scalar &wstar);

  KOKKOS_FUNCTION
  static void shoc_diag_second_moments_ubycond(Scalar &thl_sec, Scalar &qw_sec, Scalar &wthl_sec,
                                               Scalar &wqw_sec, Scalar &qwthl_sec, Scalar &uw_sec,
                                               Scalar &vw_sec, Scalar &wtke_sec);

  KOKKOS_FUNCTION
  static void update_host_dse(const MemberType &team, const Int &nlev,
                              const uview_1d<const Pack> &thlm,
                              const uview_1d<const Pack> &shoc_ql,
                              const uview_1d<const Pack> &inv_exner,
                              const uview_1d<const Pack> &zt_grid, const Scalar &phis,
                              const uview_1d<Pack> &host_dse);
#ifdef SCREAM_SHOC_SMALL_KERNELS
  static void
  update_host_dse_disp(const Int &shcol, const Int &nlev, const view_2d<const Pack> &thlm,
                       const view_2d<const Pack> &shoc_ql, const view_2d<const Pack> &inv_exner,
                       const view_2d<const Pack> &zt_grid, const view_1d<const Scalar> &phis,
                       const view_2d<Pack> &host_dse);
#endif

  KOKKOS_FUNCTION
  static void compute_diag_third_shoc_moment(
      const MemberType &team, const Int &nlev, const Int &nlevi, const Scalar &c_diag_3rd_mom,
      const bool &shoc_1p5tke, const uview_1d<const Pack> &w_sec,
      const uview_1d<const Pack> &thl_sec, const uview_1d<const Pack> &wthl_sec,
      const uview_1d<const Pack> &tke, const uview_1d<const Pack> &dz_zt,
      const uview_1d<const Pack> &dz_zi, const uview_1d<const Pack> &isotropy_zi,
      const uview_1d<const Pack> &brunt_zi, const uview_1d<const Pack> &w_sec_zi,
      const uview_1d<const Pack> &thetal_zi, const uview_1d<Pack> &w3);

  KOKKOS_FUNCTION
  static void shoc_pblintd_init_pot(const MemberType &team, const Int &nlev,
                                    const view_1d<const Pack> &thl, const view_1d<const Pack> &ql,
                                    const view_1d<const Pack> &q, const view_1d<Pack> &thv);

  KOKKOS_FUNCTION
  static void compute_shoc_mix_shoc_length(
      const MemberType &team, const Int &nlev, const Scalar &length_fac, const bool &shoc_1p5tke,
      const uview_1d<const Pack> &tke, const uview_1d<const Pack> &brunt,
      const uview_1d<const Pack> &zt_grid, const uview_1d<const Pack> &dz_zt,
      const uview_1d<const Pack> &tk, const Scalar &l_inf, const uview_1d<Pack> &shoc_mix);

  KOKKOS_FUNCTION
  static void check_tke(const MemberType &team, const Int &nlev, const uview_1d<Pack> &tke);
#ifdef SCREAM_SHOC_SMALL_KERNELS
  static void check_tke_disp(const Int &schol, const Int &nlev, const view_2d<Pack> &tke);
#endif

  KOKKOS_FUNCTION
  static void clipping_diag_third_shoc_moments(const MemberType &team, const Int &nlevi,
                                               const uview_1d<const Pack> &w_sec_zi,
                                               const uview_1d<Pack> &w3);

  KOKKOS_FUNCTION
  static void linear_interp(const MemberType &team, const uview_1d<const Pack> &x1,
                            const uview_1d<const Pack> &x2, const uview_1d<const Pack> &y1,
                            const uview_1d<Pack> &y2, const Int &km1, const Int &km2,
                            const Scalar &minthresh);

  KOKKOS_FUNCTION
  static void
  shoc_energy_integrals(const MemberType &team, const Int &nlev,
                        const uview_1d<const Pack> &host_dse, const uview_1d<const Pack> &pdel,
                        const uview_1d<const Pack> &rtm, const uview_1d<const Pack> &rcm,
                        const uview_1d<const Pack> &u_wind, const uview_1d<const Pack> &v_wind,
                        Scalar &se_int, Scalar &ke_int, Scalar &wv_int, Scalar &wl_int);
#ifdef SCREAM_SHOC_SMALL_KERNELS
  static void
  shoc_energy_integrals_disp(const Int &shcol, const Int &nlev,
                             const view_2d<const Pack> &host_dse, const view_2d<const Pack> &pdel,
                             const view_2d<const Pack> &rtm, const view_2d<const Pack> &rcm,
                             const uview_2d<const Pack> &u_wind,
                             const uview_2d<const Pack> &v_wind, const view_1d<Scalar> &se_b_slot,
                             const view_1d<Scalar> &ke_b_slot, const view_1d<Scalar> &wv_b_slot,
                             const view_1d<Scalar> &wl_b_slot);
#endif

  KOKKOS_FUNCTION
  static void shoc_diag_second_moments_lbycond(const Scalar &wthl_sfc, const Scalar &wqw_sfc,
                                               const Scalar &uw_sfc, const Scalar &vw_sfc,
                                               const Scalar &ustar2, const Scalar &wstar,
                                               Scalar &wthl_sec, Scalar &wqw_sec, Scalar &uw_sec,
                                               Scalar &vw_sec, Scalar &wtke_sec, Scalar &thl_sec,
                                               Scalar &qw_sec, Scalar &qwthl_sec);

  KOKKOS_FUNCTION
  static void diag_second_moments(
      const MemberType &team, const Int &nlev, const Int &nlevi, const Real &thl2tune,
      const Real &qw2tune, const Real &qwthl2tune, const Real &w2tune, const bool &shoc_1p5tke,
      const uview_1d<const Pack> &thetal, const uview_1d<const Pack> &qw,
      const uview_1d<const Pack> &u_wind, const uview_1d<const Pack> &v_wind,
      const uview_1d<const Pack> &tke, const uview_1d<const Pack> &isotropy,
      const uview_1d<const Pack> &tkh, const uview_1d<const Pack> &tk,
      const uview_1d<const Pack> &dz_zi, const uview_1d<const Pack> &zt_grid,
      const uview_1d<const Pack> &zi_grid, const uview_1d<const Pack> &shoc_mix,
      const uview_1d<Pack> &isotropy_zi, const uview_1d<Pack> &tkh_zi,
      const uview_1d<Pack> &tk_zi, const uview_1d<Pack> &thl_sec, const uview_1d<Pack> &qw_sec,
      const uview_1d<Pack> &wthl_sec, const uview_1d<Pack> &wqw_sec,
      const uview_1d<Pack> &qwthl_sec, const uview_1d<Pack> &uw_sec,
      const uview_1d<Pack> &vw_sec, const uview_1d<Pack> &wtke_sec, const uview_1d<Pack> &w_sec);

  KOKKOS_FUNCTION
  static void diag_second_shoc_moments(
      const MemberType &team, const Int &nlev, const Int &nlevi, const Scalar &thl2tune,
      const Scalar &qw2tune, const Scalar &qwthl2tune, const Scalar &w2tune,
      const bool &shoc_1p5tke, const uview_1d<const Pack> &thetal, const uview_1d<const Pack> &qw,
      const uview_1d<const Pack> &u_wind, const uview_1d<const Pack> &v_wind,
      const uview_1d<const Pack> &tke, const uview_1d<const Pack> &isotropy,
      const uview_1d<const Pack> &tkh, const uview_1d<const Pack> &tk,
      const uview_1d<const Pack> &dz_zi, const uview_1d<const Pack> &zt_grid,
      const uview_1d<const Pack> &zi_grid, const uview_1d<const Pack> &shoc_mix,
      const Scalar &wthl_sfc, const Scalar &wqw_sfc, const Scalar &uw_sfc, const Scalar &vw_sfc,
      Scalar &ustar2, Scalar &wstar, const Workspace &workspace, const uview_1d<Pack> &thl_sec,
      const uview_1d<Pack> &qw_sec, const uview_1d<Pack> &wthl_sec,
      const uview_1d<Pack> &wqw_sec, const uview_1d<Pack> &qwthl_sec,
      const uview_1d<Pack> &uw_sec, const uview_1d<Pack> &vw_sec, const uview_1d<Pack> &wtke_sec,
      const uview_1d<Pack> &w_sec);
#ifdef SCREAM_SHOC_SMALL_KERNELS
  static void diag_second_shoc_moments_disp(
      const Int &shcol, const Int &nlev, const Int &nlevi, const Scalar &thl2tune,
      const Scalar &qw2tune, const Scalar &qwthl2tune, const Scalar &w2tune,
      const bool &shoc_1p5tke, const view_2d<const Pack> &thetal, const view_2d<const Pack> &qw,
      const view_2d<const Pack> &u_wind, const view_2d<const Pack> &v_wind,
      const view_2d<const Pack> &tke, const view_2d<const Pack> &isotropy,
      const view_2d<const Pack> &tkh, const view_2d<const Pack> &tk,
      const view_2d<const Pack> &dz_zi, const view_2d<const Pack> &zt_grid,
      const view_2d<const Pack> &zi_grid, const view_2d<const Pack> &shoc_mix,
      const view_1d<const Scalar> &wthl_sfc, const view_1d<const Scalar> &wqw_sfc,
      const view_1d<const Scalar> &uw_sfc, const view_1d<const Scalar> &vw_sfc,
      const view_1d<Scalar> &ustar2, const view_1d<Scalar> &wstar,
      const WorkspaceMgr &workspace_mgr, const view_2d<Pack> &thl_sec,
      const view_2d<Pack> &qw_sec, const view_2d<Pack> &wthl_sec, const view_2d<Pack> &wqw_sec,
      const view_2d<Pack> &qwthl_sec, const view_2d<Pack> &uw_sec, const view_2d<Pack> &vw_sec,
      const view_2d<Pack> &wtke_sec, const view_2d<Pack> &w_sec);
#endif

  KOKKOS_FUNCTION
  static void compute_brunt_shoc_length(const MemberType &team, const Int &nlev, const Int &nlevi,
                                        const uview_1d<const Pack> &dz_zt,
                                        const uview_1d<const Pack> &thv,
                                        const uview_1d<const Pack> &thv_zi,
                                        const uview_1d<Pack> &brunt);

  KOKKOS_FUNCTION
  static void compute_l_inf_shoc_length(const MemberType &team, const Int &nlev,
                                        const uview_1d<const Pack> &zt_grid,
                                        const uview_1d<const Pack> &dz_zt,
                                        const uview_1d<const Pack> &tke, Scalar &l_inf);

  KOKKOS_FUNCTION
  static void check_length_scale_shoc_length(const MemberType &team, const Int &nlev,
                                             const Scalar &dx, const Scalar &dy,
                                             const uview_1d<Pack> &shoc_mix);

  KOKKOS_FUNCTION
  static void shoc_diag_obklen(const Scalar &uw_sfc, const Scalar &vw_sfc, const Scalar &wthl_sfc,
                               const Scalar &wqw_sfc, const Scalar &thl_sfc,
                               const Scalar &cldliq_sfc, const Scalar &qv_sfc, Scalar &ustar,
                               Scalar &kbfs, Scalar &obklen);
#ifdef SCREAM_SHOC_SMALL_KERNELS
  static void
  shoc_diag_obklen_disp(const Int &shcol, const Int &nlev, const view_1d<const Scalar> &uw_sfc,
                        const view_1d<const Scalar> &vw_sfc, const view_1d<const Scalar> &wthl_sfc,
                        const view_1d<const Scalar> &wqw_sfc, const view_2d<const Scalar> &thl_sfc,
                        const view_2d<const Scalar> &cldliq_sfc,
                        const view_2d<const Scalar> &qv_sfc, const view_1d<Scalar> &ustar,
                        const view_1d<Scalar> &kbfs, const view_1d<Scalar> &obklen);
#endif

  KOKKOS_FUNCTION
  static void shoc_pblintd_cldcheck(const Scalar &zi, const Scalar &cldn, Scalar &pblh);

  KOKKOS_FUNCTION
  static void shoc_length(const MemberType &team, const Int &nlev, const Int &nlevi,
                          const Scalar &length_fac, const bool &shoc_1p5tke, const Scalar &dx,
                          const Scalar &dy, const uview_1d<const Pack> &zt_grid,
                          const uview_1d<const Pack> &zi_grid, const uview_1d<const Pack> &dz_zt,
                          const uview_1d<const Pack> &tke, const uview_1d<const Pack> &thv,
                          const uview_1d<const Pack> &tk, const Workspace &workspace,
                          const uview_1d<Pack> &brunt, const uview_1d<Pack> &shoc_mix);
#ifdef SCREAM_SHOC_SMALL_KERNELS
  static void shoc_length_disp(const Int &shcol, const Int &nlev, const Int &nlevi,
                               const Scalar &length_fac, const bool &tke_1p5_closure,
                               const view_1d<const Scalar> &dx, const view_1d<const Scalar> &dy,
                               const view_2d<const Pack> &zt_grid,
                               const view_2d<const Pack> &zi_grid,
                               const view_2d<const Pack> &dz_zt, const view_2d<const Pack> &tke,
                               const view_2d<const Pack> &thv, const view_2d<const Pack> &tk,
                               const WorkspaceMgr &workspace_mgr, const view_2d<Pack> &brunt,
                               const view_2d<Pack> &shoc_mix);
#endif

  KOKKOS_FUNCTION
  static void shoc_energy_fixer(const MemberType &team, const Int &nlev, const Int &nlevi,
                                const Scalar &dtime, const Int &nadv,
                                const uview_1d<const Pack> &zt_grid,
                                const uview_1d<const Pack> &zi_grid, const Scalar &se_b,
                                const Scalar &ke_b, const Scalar &wv_b, const Scalar &wl_b,
                                const Scalar &se_a, const Scalar &ke_a, const Scalar &wv_a,
                                const Scalar &wl_a, const Scalar &wthl_sfc, const Scalar &wqw_sfc,
                                const uview_1d<const Pack> &rho_zt,
                                const uview_1d<const Pack> &tke, const uview_1d<const Pack> &pint,
                                const Workspace &workspace, const uview_1d<Pack> &host_dse);
#ifdef SCREAM_SHOC_SMALL_KERNELS
  static void
  shoc_energy_fixer_disp(const Int &shcol, const Int &nlev, const Int &nlevi, const Scalar &dtime,
                         const Int &nadv, const view_2d<const Pack> &zt_grid,
                         const view_2d<const Pack> &zi_grid, const view_1d<const Scalar> &se_b,
                         const view_1d<const Scalar> &ke_b, const view_1d<const Scalar> &wv_b,
                         const view_1d<const Scalar> &wl_b, const view_1d<const Scalar> &se_a,
                         const view_1d<const Scalar> &ke_a, const view_1d<const Scalar> &wv_a,
                         const view_1d<const Scalar> &wl_a, const view_1d<const Scalar> &wthl_sfc,
                         const view_1d<const Scalar> &wqw_sfc, const view_2d<const Pack> &rho_zt,
                         const view_2d<const Pack> &tke, const view_2d<const Pack> &pint,
                         const WorkspaceMgr &workspace_mgr, const view_2d<Pack> &host_dse);
#endif

  KOKKOS_FUNCTION
  static void compute_shoc_vapor(const MemberType &team, const Int &nlev,
                                 const uview_1d<const Pack> &qw, const uview_1d<const Pack> &ql,
                                 const uview_1d<Pack> &qv);
#ifdef SCREAM_SHOC_SMALL_KERNELS
  static void compute_shoc_vapor_disp(const Int &shcol, const Int &nlev,
                                      const view_2d<const Pack> &qw,
                                      const view_2d<const Pack> &ql, const view_2d<Pack> &qv);
#endif

  KOKKOS_FUNCTION
  static void compute_shoc_temperature(const MemberType &team, const Int &nlev,
                                       const uview_1d<const Pack> &thetal,
                                       const uview_1d<const Pack> &ql,
                                       const uview_1d<const Pack> &inv_exner,
                                       const uview_1d<Pack> &tabs);
#ifdef SCREAM_SHOC_SMALL_KERNELS
  static void compute_shoc_temperature_disp(const Int &shcol, const Int &nlev,
                                            const view_2d<const Pack> &thetal,
                                            const view_2d<const Pack> &ql,
                                            const view_2d<const Pack> &inv_exner,
                                            const view_2d<Pack> &tabs);
#endif

  KOKKOS_FUNCTION
  static void update_prognostics_implicit(
      const MemberType &team, const Int &nlev, const Int &nlevi, const Int &num_tracer,
      const Scalar &dtime, const uview_1d<const Pack> &dz_zt, const uview_1d<const Pack> &dz_zi,
      const uview_1d<const Pack> &rho_zt, const uview_1d<const Pack> &zt_grid,
      const uview_1d<const Pack> &zi_grid, const uview_1d<const Pack> &tk,
      const uview_1d<const Pack> &tkh, const Scalar &uw_sfc, const Scalar &vw_sfc,
      const Scalar &wthl_sfc, const Scalar &wqw_sfc, const uview_1d<const Pack> &wtracer_sfc,
      const Workspace &workspace, const uview_1d<Pack> &thetal, const uview_1d<Pack> &qw,
      const uview_2d_strided<Pack> &tracer, const uview_1d<Pack> &tke,
      const uview_1d<Pack> &u_wind, const uview_1d<Pack> &v_wind);
#ifdef SCREAM_SHOC_SMALL_KERNELS
  static void update_prognostics_implicit_disp(
      const Int &shcol, const Int &nlev, const Int &nlevi, const Int &num_tracer,
      const Scalar &dtime, const view_2d<const Pack> &dz_zt, const view_2d<const Pack> &dz_zi,
      const view_2d<const Pack> &rho_zt, const view_2d<const Pack> &zt_grid,
      const view_2d<const Pack> &zi_grid, const view_2d<const Pack> &tk,
      const view_2d<const Pack> &tkh, const view_1d<const Scalar> &uw_sfc,
      const view_1d<const Scalar> &vw_sfc, const view_1d<const Scalar> &wthl_sfc,
      const view_1d<const Scalar> &wqw_sfc, const view_2d<const Pack> &wtracer_sfc,
      const WorkspaceMgr &workspace_mgr, const view_2d<Pack> &thetal, const view_2d<Pack> &qw,
      const view_3d_strided<Pack> &tracer, const view_2d<Pack> &tke, const view_2d<Pack> &u_wind,
      const view_2d<Pack> &v_wind);
#endif

  KOKKOS_FUNCTION
  static void diag_third_shoc_moments(
      const MemberType &team, const Int &nlev, const Int &nlevi, const Scalar &c_diag_3rd_mom,
      const bool &shoc_1p5tke, const uview_1d<const Pack> &w_sec,
      const uview_1d<const Pack> &thl_sec, const uview_1d<const Pack> &wthl_sec,
      const uview_1d<const Pack> &isotropy, const uview_1d<const Pack> &brunt,
      const uview_1d<const Pack> &thetal, const uview_1d<const Pack> &tke,
      const uview_1d<const Pack> &dz_zt, const uview_1d<const Pack> &dz_zi,
      const uview_1d<const Pack> &zt_grid, const uview_1d<const Pack> &zi_grid,
      const Workspace &workspace, const uview_1d<Pack> &w3);
#ifdef SCREAM_SHOC_SMALL_KERNELS
  static void diag_third_shoc_moments_disp(
      const Int &shcol, const Int &nlev, const Int &nlevi, const Scalar &c_diag_3rd_mom,
      const bool &shoc_1p5tke, const view_2d<const Pack> &w_sec,
      const view_2d<const Pack> &thl_sec, const view_2d<const Pack> &wthl_sec,
      const view_2d<const Pack> &isotropy, const view_2d<const Pack> &brunt,
      const view_2d<const Pack> &thetal, const view_2d<const Pack> &tke,
      const view_2d<const Pack> &dz_zt, const view_2d<const Pack> &dz_zi,
      const view_2d<const Pack> &zt_grid, const view_2d<const Pack> &zi_grid,
      const WorkspaceMgr &workspace_mgr, const view_2d<Pack> &w3);
#endif

  KOKKOS_FUNCTION
  static void adv_sgs_tke(const MemberType &team, const Int &nlev, const Real &dtime,
                          const bool &shoc_1p5tke, const bool &do_3d_turb,
                          const uview_1d<const Pack> &shoc_mix, const uview_1d<const Pack> &wthv_sec,
                          const uview_1d<const Pack> &sterm_zt, const uview_1d<const Pack> &tk,
                          const uview_1d<const Pack> &brunt, const uview_1d<const Pack> &strain2,
                          const uview_1d<Pack> &tke, const uview_1d<Pack> &a_diss);

  KOKKOS_FUNCTION
  static void
  shoc_assumed_pdf(const MemberType &team, const Int &nlev, const Int &nlevi,
                   const uview_1d<const Pack> &thetal, const uview_1d<const Pack> &qw,
                   const uview_1d<const Pack> &w_field, const uview_1d<const Pack> &thl_sec,
                   const uview_1d<const Pack> &qw_sec, const Scalar &dtime,
                   const bool &extra_diags, const uview_1d<const Pack> &wthl_sec,
                   const uview_1d<const Pack> &w_sec, const uview_1d<const Pack> &wqw_sec,
                   const uview_1d<const Pack> &qwthl_sec, const uview_1d<const Pack> &w3,
                   const uview_1d<const Pack> &pres, const uview_1d<const Pack> &zt_grid,
                   const uview_1d<const Pack> &zi_grid, const Workspace &workspace,
                   const uview_1d<Pack> &shoc_cond, const uview_1d<Pack> &shoc_evap,
                   const uview_1d<Pack> &shoc_cldfrac, const uview_1d<Pack> &shoc_ql,
                   const uview_1d<Pack> &wqls, const uview_1d<Pack> &wthv_sec,
                   const uview_1d<Pack> &shoc_ql2);
#ifdef SCREAM_SHOC_SMALL_KERNELS
  static void shoc_assumed_pdf_disp(
      const Int &shcol, const Int &nlev, const Int &nlevi, const view_2d<const Pack> &thetal,
      const view_2d<const Pack> &qw, const view_2d<const Pack> &w_field,
      const view_2d<const Pack> &thl_sec, const view_2d<const Pack> &qw_sec, const Scalar &dtime,
      const bool &extra_diags, const view_2d<const Pack> &wthl_sec,
      const view_2d<const Pack> &w_sec, const view_2d<const Pack> &wqw_sec,
      const view_2d<const Pack> &qwthl_sec, const view_2d<const Pack> &w3,
      const view_2d<const Pack> &pres, const view_2d<const Pack> &zt_grid,
      const view_2d<const Pack> &zi_grid, const WorkspaceMgr &workspace_mgr,
      const view_2d<Pack> &shoc_cond, const view_2d<Pack> &shoc_evap,
      const view_2d<Pack> &shoc_cldfrac, const view_2d<Pack> &shoc_ql, const view_2d<Pack> &wqls,
      const view_2d<Pack> &wthv_sec, const view_2d<Pack> &shoc_ql2);
#endif

  KOKKOS_INLINE_FUNCTION
  static void shoc_assumed_pdf_compute_buoyancy_flux(const Pack &wthlsec, const Pack &wqwsec,
                                                     const Pack &pval, const Pack &wqls,
                                                     Pack &wthv_sec);

  KOKKOS_INLINE_FUNCTION
  static void shoc_assumed_pdf_compute_cloud_liquid_variance(const Pack &a, const Pack &s1,
                                                             const Pack &ql1, const Pack &C1,
                                                             const Pack &std_s1, const Pack &s2,
                                                             const Pack &ql2, const Pack &C2,
                                                             const Pack &std_s2,
                                                             const Pack &shoc_ql, Pack &shoc_ql2);

  KOKKOS_INLINE_FUNCTION
  static void shoc_assumed_pdf_compute_liquid_water_flux(const Pack &a, const Pack &w1_1,
                                                         const Pack &w_first, const Pack &ql1,
                                                         const Pack &w1_2, const Pack &ql2,
                                                         Pack &wqls);

  KOKKOS_INLINE_FUNCTION
  static void shoc_assumed_pdf_compute_qs(const Pack &Tl1_1, const Pack &Tl1_2, const Pack &pval,
                                          const Mask &active_entries, Pack &qs1, Pack &beta1,
                                          Pack &qs2, Pack &beta2);

  KOKKOS_INLINE_FUNCTION
  static void shoc_assumed_pdf_compute_s(const Pack &qw1, const Pack &qs, const Pack &beta,
                                         const Pack &pval, const Pack &thl2, const Pack &qw2,
                                         const Pack &sqrtthl2, const Pack &sqrtqw2,
                                         const Pack &r_qwthl, Pack &s, Pack &std_s, Pack &qn,
                                         Pack &C);

  KOKKOS_INLINE_FUNCTION
  static void shoc_assumed_pdf_compute_sgs_liquid(const Pack &a, const Pack &ql1,
                                                  const Pack &ql2, Pack &shoc_ql);

  KOKKOS_INLINE_FUNCTION
  static void shoc_assumed_pdf_compute_temperature(const Pack &thl1, const Pack &pval,
                                                   Pack &Tl1);

  KOKKOS_INLINE_FUNCTION
  static void shoc_assumed_pdf_inplume_correlations(const Pack &sqrtqw2_1, const Pack &sqrtthl2_1,
                                                    const Pack &a, const Pack &sqrtqw2_2,
                                                    const Pack &sqrtthl2_2, const Pack &qwthlsec,
                                                    const Pack &qw1_1, const Pack &qw_first,
                                                    const Pack &thl1_1, const Pack &thl_first,
                                                    const Pack &qw1_2, const Pack &thl1_2,
                                                    Pack &r_qwthl_1);

  KOKKOS_INLINE_FUNCTION
  static void shoc_assumed_pdf_qw_parameters(
      const Pack &wqwsec, const Pack &sqrtw2, const Pack &Skew_w, const Pack &sqrtqt,
      const Pack &qwsec, const Pack &w1_2, const Pack &w1_1, const Pack &qw_first,
      const Pack &a, const Scalar rt_tol, const Scalar w_thresh, Pack &qw1_1, Pack &qw1_2,
      Pack &qw2_1, Pack &qw2_2, Pack &sqrtqw2_1, Pack &sqrtqw2_2);

  KOKKOS_INLINE_FUNCTION
  static void shoc_assumed_pdf_thl_parameters(
      const Pack &wthlsec, const Pack &sqrtw2, const Pack &sqrtthl, const Pack &thlsec,
      const Pack &thl_first, const Pack &w1_1, const Pack &w1_2, const Pack &Skew_w,
      const Pack &a, const Scalar thl_tol, const Scalar w_thresh, Pack &thl1_1, Pack &thl1_2,
      Pack &thl2_1, Pack &thl2_2, Pack &sqrtthl2_1, Pack &sqrtthl2_2);

  KOKKOS_INLINE_FUNCTION
  static void shoc_assumed_pdf_tilde_to_real(const Pack &w_first, const Pack &sqrtw2, Pack &w1);

  KOKKOS_INLINE_FUNCTION
  static void shoc_assumed_pdf_vv_parameters(const Pack &w_first, const Pack &w_sec,
                                             const Pack &w3var, const Scalar w_tol_sqd,
                                             Pack &Skew_w, Pack &w1_1, Pack &w1_2, Pack &w2_1,
                                             Pack &w2_2, Pack &a);

  KOKKOS_FUNCTION
  static void compute_shr_prod(const MemberType &team, const Int &nlevi, const Int &nlev,
                               const uview_1d<const Pack> &dz_zi,
                               const uview_1d<const Pack> &u_wind,
                               const uview_1d<const Pack> &v_wind, const uview_1d<Pack> &sterm);

  KOKKOS_FUNCTION
  static void compute_tmpi(const MemberType &team, const Int &nlevi, const Scalar &dtime,
                           const uview_1d<const Pack> &rho_zi, const uview_1d<const Pack> &dz_zi,
                           const uview_1d<Pack> &tmpi);

  KOKKOS_FUNCTION
  static void integ_column_stability(const MemberType &team, const Int &nlev,
                                     const uview_1d<const Pack> &dz_zt,
                                     const uview_1d<const Pack> &pres,
                                     const uview_1d<const Pack> &brunt, Scalar &brunt_int);

  KOKKOS_FUNCTION
  static void isotropic_ts(const MemberType &team, const Int &nlev, const Scalar &lambda_low,
                           const Scalar &lambda_high, const Scalar &lambda_slope,
                           const Scalar &lambda_thresh, const Scalar &brunt_int,
                           const uview_1d<const Pack> &tke, const uview_1d<const Pack> &a_diss,
                           const uview_1d<const Pack> &brunt, const uview_1d<Pack> &isotropy);

  KOKKOS_FUNCTION
  static void dp_inverse(const MemberType &team, const Int &nlev,
                         const uview_1d<const Pack> &rho_zt, const uview_1d<const Pack> &dz_zt,
                         const uview_1d<Pack> &rdp_zt);

  static Int shoc_init(const Int &nbot_shoc, const Int &ntop_shoc,
                       const view_1d<const Pack> &pref_mid);

#ifndef SCREAM_SHOC_SMALL_KERNELS
  KOKKOS_FUNCTION
  static void shoc_main_internal(
      const MemberType &team,
      const Int &nlev,         // Number of levels
      const Int &nlevi,        // Number of levels on interface grid
      const Int &npbl,         // Maximum number of levels in pbl from surface
      const Int &nadv,         // Number of times to loop SHOC
      const Int &num_qtracers, // Number of tracers
      const Scalar &dtime,     // SHOC timestep [s]
      // Runtime Parameters
      const Scalar &lambda_low, const Scalar &lambda_high, const Scalar &lambda_slope,
      const Scalar &lambda_thresh, const Scalar &thl2tune, const Scalar &qw2tune,
      const Scalar &qwthl2tune, const Scalar &w2tune, const Scalar &length_fac,
      const Scalar &c_diag_3rd_mom, const Scalar &Ckh, const Scalar &Ckm, const bool &shoc_1p5tke,
      const bool &do_3d_turb, const bool &extra_diags,
      // Input Variables
      const Scalar &host_dx, const Scalar &host_dy, const uview_1d<const Pack> &zt_grid,
      const uview_1d<const Pack> &zi_grid, const uview_1d<const Pack> &pres,
      const uview_1d<const Pack> &presi, const uview_1d<const Pack> &pdel,
      const uview_1d<const Pack> &thv, const uview_1d<const Pack> &w_field,
      const Scalar &wthl_sfc, const Scalar &wqw_sfc, const Scalar &uw_sfc, const Scalar &vw_sfc,
      const uview_1d<const Pack> &wtracer_sfc, const uview_1d<const Pack> &inv_exner,
      const Scalar &phis,
      const uview_1d<const Pack> &strain2,
      // Local Workspace
      const Workspace &workspace,
      // Input/Output Variables
      const uview_1d<Pack> &host_dse, const uview_1d<Pack> &tke, const uview_1d<Pack> &thetal,
      const uview_1d<Pack> &qw, const uview_1d<Pack> &u_wind, const uview_1d<Pack> &v_wind,
      const uview_1d<Pack> &wthv_sec, const uview_2d_strided<Pack> &qtracers,
      const uview_1d<Pack> &tk, const uview_1d<Pack> &shoc_cldfrac,
      const uview_1d<Pack> &shoc_ql,
      // Output Variables
      Scalar &pblh, Scalar &ustar, Scalar &obklen, const uview_1d<Pack> &shoc_ql2,
      const uview_1d<Pack> &tkh,
      // Diagnostic Output Variables
      const uview_1d<Pack> &shoc_cond, const uview_1d<Pack> &shoc_evap,
      const uview_1d<Pack> &shoc_mix, const uview_1d<Pack> &w_sec, const uview_1d<Pack> &thl_sec,
      const uview_1d<Pack> &qw_sec, const uview_1d<Pack> &qwthl_sec,
      const uview_1d<Pack> &wthl_sec, const uview_1d<Pack> &wqw_sec,
      const uview_1d<Pack> &wtke_sec, const uview_1d<Pack> &uw_sec, const uview_1d<Pack> &vw_sec,
      const uview_1d<Pack> &w3, const uview_1d<Pack> &wqls_sec, const uview_1d<Pack> &brunt,
      const uview_1d<Pack> &isotropy);
#else
  static void shoc_main_internal(
      const Int &shcol,        // Number of columns
      const Int &nlev,         // Number of levels
      const Int &nlevi,        // Number of levels on interface grid
      const Int &npbl,         // Maximum number of levels in pbl from surface
      const Int &nadv,         // Number of times to loop SHOC
      const Int &num_qtracers, // Number of tracers
      const Scalar &dtime,     // SHOC timestep [s]
      // Runtime Parameters
      const Scalar &lambda_low, const Scalar &lambda_high, const Scalar &lambda_slope,
      const Scalar &lambda_thresh, const Scalar &thl2tune, const Scalar &qw2tune,
      const Scalar &qwthl2tune, const Scalar &w2tune, const Scalar &length_fac,
      const Scalar &c_diag_3rd_mom, const Scalar &Ckh, const Scalar &Ckm, const bool &shoc_1p5tke,
      const bool &do_3d_turb, const bool &extra_diags,
      // Input Variables
      const view_1d<const Scalar> &host_dx, const view_1d<const Scalar> &host_dy,
      const view_2d<const Pack> &zt_grid, const view_2d<const Pack> &zi_grid,
      const view_2d<const Pack> &pres, const view_2d<const Pack> &presi,
      const view_2d<const Pack> &pdel, const view_2d<const Pack> &thv,
      const view_2d<const Pack> &w_field, const view_1d<const Scalar> &wthl_sfc,
      const view_1d<const Scalar> &wqw_sfc, const view_1d<const Scalar> &uw_sfc,
      const view_1d<const Scalar> &vw_sfc, const view_2d<const Pack> &wtracer_sfc,
      const view_2d<const Pack> &inv_exner, const view_1d<const Scalar> &phis,
      const view_2d<const Pack> &strain2,
      // Workspace Manager
      WorkspaceMgr &workspace_mgr,
      // Input/Output Variables
      const view_2d<Pack> &host_dse, const view_2d<Pack> &tke, const view_2d<Pack> &thetal,
      const view_2d<Pack> &qw, const uview_2d<Pack> &u_wind, const uview_2d<Pack> &v_wind,
      const view_2d<Pack> &wthv_sec, const view_3d_strided<Pack> &qtracers,
      const view_2d<Pack> &tk, const view_2d<Pack> &shoc_cldfrac, const view_2d<Pack> &shoc_ql,
      // Output Variables
      const view_1d<Scalar> &pblh, const view_1d<Scalar> &ustar, const view_1d<Scalar> &obklen,
      const view_2d<Pack> &shoc_ql2, const view_2d<Pack> &tkh,
      // Diagnostic Output Variables
      const view_2d<Pack> &shoc_evap, const view_2d<Pack> &shoc_cond,
      const view_2d<Pack> &shoc_mix, const view_2d<Pack> &w_sec, const view_2d<Pack> &thl_sec,
      const view_2d<Pack> &qw_sec, const view_2d<Pack> &qwthl_sec, const view_2d<Pack> &wthl_sec,
      const view_2d<Pack> &wqw_sec, const view_2d<Pack> &wtke_sec, const view_2d<Pack> &uw_sec,
      const view_2d<Pack> &vw_sec, const view_2d<Pack> &w3, const view_2d<Pack> &wqls_sec,
      const view_2d<Pack> &brunt, const view_2d<Pack> &isotropy,
      // Temporaries
      const view_1d<Scalar> &se_b, const view_1d<Scalar> &ke_b, const view_1d<Scalar> &wv_b,
      const view_1d<Scalar> &wl_b, const view_1d<Scalar> &se_a, const view_1d<Scalar> &ke_a,
      const view_1d<Scalar> &wv_a, const view_1d<Scalar> &wl_a, const view_1d<Scalar> &kbfs,
      const view_1d<Scalar> &ustar2, const view_1d<Scalar> &wstar, const view_2d<Pack> &rho_zt,
      const view_2d<Pack> &shoc_qv, const view_2d<Pack> &tabs, const view_2d<Pack> &dz_zt,
      const view_2d<Pack> &dz_zi);
#endif

  // Return microseconds elapsed
  static Int shoc_main(const Int &shcol,            // Number of SHOC columns in the array
                       const Int &nlev,             // Number of levels
                       const Int &nlevi,            // Number of levels on interface grid
                       const Int &npbl,             // Maximum number of levels in pbl from surface
                       const Int &nadv,             // Number of times to loop SHOC
                       const Int &num_q_tracers,    // Number of tracers
                       const Scalar &dtime,         // SHOC timestep [s]
                       WorkspaceMgr &workspace_mgr, // WorkspaceManager for local variables
                       const SHOCRuntime &shoc_runtime,             // Runtime options
                       const SHOCInput &shoc_input,                 // Input
                       const SHOCInputOutput &shoc_input_output,    // Input/Output
                       const SHOCOutput &shoc_output,               // Output
                       const SHOCHistoryOutput &shoc_history_output // Output (diagnostic)
#ifdef SCREAM_SHOC_SMALL_KERNELS
                       ,
                       const SHOCTemporaries &shoc_temporaries // Temporaries for small kernels
#endif
  );

  KOKKOS_FUNCTION
  static void pblintd_height(const MemberType &team, const Int &nlev, const Int &npbl,
                             const uview_1d<const Pack> &z, const uview_1d<const Pack> &u,
                             const uview_1d<const Pack> &v, const Scalar &ustar,
                             const uview_1d<const Pack> &thv, const Scalar &thv_ref, Scalar &pblh,
                             const uview_1d<Pack> &rino, bool &check);

  KOKKOS_FUNCTION
  static void vd_shoc_decomp(const MemberType &team, const Int &nlev,
                             const uview_1d<const Pack> &kv_term,
                             const uview_1d<const Pack> &tmpi, const uview_1d<const Pack> &rdp_zt,
                             const Scalar &dtime, const Scalar &flux, const uview_1d<Scalar> &du,
                             const uview_1d<Scalar> &dl, const uview_1d<Scalar> &d);

  KOKKOS_FUNCTION
  static void vd_shoc_solve(const MemberType &team, const uview_1d<Scalar> &du,
                            const uview_1d<Scalar> &dl, const uview_1d<Scalar> &d,
                            const uview_2d<Pack> &var);

  KOKKOS_FUNCTION
  static void pblintd_surf_temp(const Int &nlev, const Int &nlevi, const Int &npbl,
                                const uview_1d<const Pack> &z, const Scalar &ustar,
                                const Scalar &obklen, const Scalar &kbfs,
                                const uview_1d<const Pack> &thv, Scalar &tlv, Scalar &pblh,
                                bool &check, const uview_1d<Pack> &rino);

  KOKKOS_FUNCTION
  static void pblintd_check_pblh(const Int &nlevi, const Int &npbl, const uview_1d<const Pack> &z,
                                 const Scalar &ustar, const bool &check, Scalar &pblh);

  KOKKOS_FUNCTION
  static void pblintd(const MemberType &team, const Int &nlev, const Int &nlevi, const Int &npbl,
                      const uview_1d<const Pack> &z, const uview_1d<const Pack> &zi,
                      const uview_1d<const Pack> &thl, const uview_1d<const Pack> &ql,
                      const uview_1d<const Pack> &q, const uview_1d<const Pack> &u,
                      const uview_1d<const Pack> &v, const Scalar &ustar, const Scalar &obklen,
                      const Scalar &kbfs, const uview_1d<const Pack> &cldn,
                      const Workspace &workspace, Scalar &pblh);
#ifdef SCREAM_SHOC_SMALL_KERNELS
  static void pblintd_disp(const Int &shcol, const Int &nlev, const Int &nlevi, const Int &npbl,
                           const view_2d<const Pack> &z, const view_2d<const Pack> &zi,
                           const view_2d<const Pack> &thl, const view_2d<const Pack> &ql,
                           const view_2d<const Pack> &q, const view_2d<const Pack> &u,
                           const view_2d<const Pack> &v, const view_1d<const Scalar> &ustar,
                           const view_1d<const Scalar> &obklen, const view_1d<const Scalar> &kbfs,
                           const view_2d<const Pack> &cldn, const WorkspaceMgr &workspace_mgr,
                           const view_1d<Scalar> &pblh);
#endif

  KOKKOS_FUNCTION
  static void shoc_grid(const MemberType &team, const Int &nlev, const Int &nlevi,
                        const uview_1d<const Pack> &zt_grid, const uview_1d<const Pack> &zi_grid,
                        const uview_1d<const Pack> &pdel, const uview_1d<Pack> &dz_zt,
                        const uview_1d<Pack> &dz_zi, const uview_1d<Pack> &rho_zt);
#ifdef SCREAM_SHOC_SMALL_KERNELS
  static void shoc_grid_disp(const Int &shcol, const Int &nlev, const Int &nlevi,
                             const view_2d<const Pack> &zt_grid,
                             const view_2d<const Pack> &zi_grid, const view_2d<const Pack> &pdel,
                             const view_2d<Pack> &dz_zt, const view_2d<Pack> &dz_zi,
                             const view_2d<Pack> &rho_zt);
#endif

  KOKKOS_FUNCTION
  static void
  eddy_diffusivities(const MemberType &team, const Int &nlev, const bool &shoc_1p5tke,
                     const Scalar &Ckh, const Scalar &Ckm, const Scalar &pblh,
                     const uview_1d<const Pack> &zt_grid, const uview_1d<const Pack> &tabs,
                     const uview_1d<const Pack> &shoc_mix, const uview_1d<const Pack> &sterm_zt,
                     const uview_1d<const Pack> &isotropy, const uview_1d<const Pack> &tke,
                     const uview_1d<Pack> &tkh, const uview_1d<Pack> &tk);

  KOKKOS_FUNCTION
  static void shoc_tke(const MemberType &team, const Int &nlev, const Int &nlevi,
                       const Scalar &dtime, const Scalar &lambda_low, const Scalar &lambda_high,
                       const Scalar &lambda_slope, const Scalar &lambda_thresh, const Scalar &Ckh,
                       const Scalar &Ckm, const bool &shoc_1p5tke, const bool &do_3d_turb,
                       const uview_1d<const Pack> &wthv_sec, const uview_1d<const Pack> &strain2,
                       const uview_1d<const Pack> &shoc_mix,
                       const uview_1d<const Pack> &dz_zi, const uview_1d<const Pack> &dz_zt,
                       const uview_1d<const Pack> &pres, const uview_1d<const Pack> &tabs,
                       const uview_1d<const Pack> &u_wind, const uview_1d<const Pack> &v_wind,
                       const uview_1d<const Pack> &brunt, const uview_1d<const Pack> &zt_grid,
                       const uview_1d<const Pack> &zi_grid, const Scalar &pblh,
                       const Workspace &workspace, const uview_1d<Pack> &tke,
                       const uview_1d<Pack> &tk, const uview_1d<Pack> &tkh,
                       const uview_1d<Pack> &isotropy);
#ifdef SCREAM_SHOC_SMALL_KERNELS
  static void shoc_tke_disp(const Int &shcol, const Int &nlev, const Int &nlevi,
                            const Scalar &dtime, const Scalar &lambda_low,
                            const Scalar &lambda_high, const Scalar &lambda_slope,
                            const Scalar &lambda_thresh, const Scalar &Ckh, const Scalar &Ckm,
                            const bool &shoc_1p5tke, const view_2d<const Pack> &wthv_sec,
                            const view_2d<const Pack> &shoc_mix, const view_2d<const Pack> &dz_zi,
                            const view_2d<const Pack> &dz_zt, const view_2d<const Pack> &pres,
                            const view_2d<const Pack> &tabs, const view_2d<const Pack> &u_wind,
                            const view_2d<const Pack> &v_wind, const view_2d<const Pack> &brunt,
                            const view_2d<const Pack> &zt_grid,
                            const view_2d<const Pack> &zi_grid, const view_1d<const Scalar> &pblh,
                            const WorkspaceMgr &workspace_mgr, const view_2d<Pack> &tke,
                            const view_2d<Pack> &tk, const view_2d<Pack> &tkh,
                            const view_2d<Pack> &isotropy);
#endif
}; // struct Functions

} // namespace shoc
} // namespace scream

// If a GPU build, without relocatable device code enabled, make all code available
// to the translation unit; otherwise, ETI is used.
#if defined(EAMXX_ENABLE_GPU) && !defined(KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE) && \
    !defined(KOKKOS_ENABLE_HIP_RELOCATABLE_DEVICE_CODE)

#include "shoc_adv_sgs_tke_impl.hpp"
#include "shoc_assumed_pdf_impl.hpp"
#include "shoc_calc_shoc_varorcovar_impl.hpp"
#include "shoc_calc_shoc_vertflux_impl.hpp"
#include "shoc_check_length_scale_shoc_length_impl.hpp"
#include "shoc_check_tke_impl.hpp"
#include "shoc_clipping_diag_third_shoc_moments_impl.hpp"
#include "shoc_compute_brunt_shoc_length_impl.hpp"
#include "shoc_compute_diag_third_shoc_moment_impl.hpp"
#include "shoc_compute_l_inf_shoc_length_impl.hpp"
#include "shoc_compute_shoc_mix_shoc_length_impl.hpp"
#include "shoc_compute_shoc_temperature_impl.hpp"
#include "shoc_compute_shoc_vapor_impl.hpp"
#include "shoc_compute_shr_prod_impl.hpp"
#include "shoc_compute_tmpi_impl.hpp"
#include "shoc_diag_obklen_impl.hpp"
#include "shoc_diag_second_moments_impl.hpp"
#include "shoc_diag_second_moments_lbycond_impl.hpp"
#include "shoc_diag_second_moments_srf_impl.hpp"
#include "shoc_diag_second_moments_ubycond_impl.hpp"
#include "shoc_diag_second_shoc_moments_impl.hpp"
#include "shoc_diag_third_shoc_moments_impl.hpp"
#include "shoc_dp_inverse_impl.hpp"
#include "shoc_eddy_diffusivities_impl.hpp"
#include "shoc_energy_fixer_impl.hpp"
#include "shoc_energy_integrals_impl.hpp"
#include "shoc_grid_impl.hpp"
#include "shoc_integ_column_stability_impl.hpp"
#include "shoc_isotropic_ts_impl.hpp"
#include "shoc_length_impl.hpp"
#include "shoc_linear_interp_impl.hpp"
#include "shoc_main_impl.hpp"
#include "shoc_pblintd_check_pblh_impl.hpp"
#include "shoc_pblintd_cldcheck_impl.hpp"
#include "shoc_pblintd_height_impl.hpp"
#include "shoc_pblintd_impl.hpp"
#include "shoc_pblintd_init_pot_impl.hpp"
#include "shoc_pblintd_surf_temp_impl.hpp"
#include "shoc_tke_impl.hpp"
#include "shoc_tridiag_solver_impl.hpp"
#include "shoc_update_host_dse_impl.hpp"
#include "shoc_update_prognostics_implicit_impl.hpp"

#endif // GPU && !KOKKOS_ENABLE_*_RELOCATABLE_DEVICE_CODE

// Some functions should be inlined, thus do not use ETI
#include "shoc_assumed_pdf_compute_buoyancy_flux_impl.hpp"
#include "shoc_assumed_pdf_compute_cloud_liquid_variance_impl.hpp"
#include "shoc_assumed_pdf_compute_liquid_water_flux_impl.hpp"
#include "shoc_assumed_pdf_compute_qs_impl.hpp"
#include "shoc_assumed_pdf_compute_s_impl.hpp"
#include "shoc_assumed_pdf_compute_sgs_liquid_impl.hpp"
#include "shoc_assumed_pdf_compute_temperature_impl.hpp"
#include "shoc_assumed_pdf_inplume_correlations_impl.hpp"
#include "shoc_assumed_pdf_qw_parameters_impl.hpp"
#include "shoc_assumed_pdf_thl_parameters_impl.hpp"
#include "shoc_assumed_pdf_tilde_to_real_impl.hpp"
#include "shoc_assumed_pdf_vv_parameters_impl.hpp"

#endif // SHOC_FUNCTIONS_HPP
