#include "p3_test_data.hpp"
#include "ekat/kokkos/ekat_kokkos_types.hpp"
#include "p3_data.hpp"

#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/ekat_assert.hpp"

#include <random>

using scream::Real;
using scream::Int;


namespace scream {
namespace p3 {

void BackToCellAverageData::randomize(std::mt19937_64& engine)
{
  // Populate the struct with numbers between 0 and 1.
  std::uniform_real_distribution<Real> data_dist(0.0, 1.0);
  cld_frac_l = data_dist(engine);
  cld_frac_r = data_dist(engine);
  cld_frac_i = data_dist(engine);
  qc2qr_accret_tend = data_dist(engine);
  qr2qv_evap_tend = data_dist(engine);
  qc2qr_autoconv_tend = data_dist(engine);
  nc_accret_tend = data_dist(engine);
  nc_selfcollect_tend = data_dist(engine);
  nc2nr_autoconv_tend = data_dist(engine);
  nr_selfcollect_tend = data_dist(engine);
  nr_evap_tend = data_dist(engine);
  ncautr = data_dist(engine);
  qcnuc = data_dist(engine);
  nc_nuceat_tend = data_dist(engine);
  qi2qv_sublim_tend = data_dist(engine);
  nr_ice_shed_tend = data_dist(engine);
  qc2qi_hetero_freeze_tend = data_dist(engine);
  qr2qi_collect_tend = data_dist(engine);
  qc2qr_ice_shed_tend = data_dist(engine);
  qi2qr_melt_tend = data_dist(engine);
  qc2qi_collect_tend = data_dist(engine);
  qr2qi_immers_freeze_tend = data_dist(engine);
  ni2nr_melt_tend = data_dist(engine);
  nc_collect_tend = data_dist(engine);
  ncshdc = data_dist(engine);
  nc2ni_immers_freeze_tend = data_dist(engine);
  nr_collect_tend = data_dist(engine);
  ni_selfcollect_tend = data_dist(engine);
  qv2qi_vapdep_tend = data_dist(engine);
  nr2ni_immers_freeze_tend = data_dist(engine);
  ni_sublim_tend = data_dist(engine);
  qv2qi_nucleat_tend = data_dist(engine);
  ni_nucleat_tend = data_dist(engine);
  qc2qi_berg_tend = data_dist(engine);
}

void CalcLiqRelaxationData::randomize(std::mt19937_64& engine)
{
  // Populate the struct's input fields with numbers between 0 and 1.
  std::uniform_real_distribution<Real> data_dist(0.0, 1.0);
  rho = data_dist(engine);
  f1r = data_dist(engine);
  f2r = data_dist(engine);
  dv = data_dist(engine);
  mu = data_dist(engine);
  sc = data_dist(engine);
  mu_r = data_dist(engine);
  lamr = data_dist(engine);
  cdistr = data_dist(engine);
  cdist = data_dist(engine);
  qr_incld = data_dist(engine);
  qc_incld = data_dist(engine);
}

CheckValuesData::CheckValuesData(
  Int kts_, Int kte_, Int timestepcount_, Int source_ind_, bool force_abort_) :
  PhysicsTestData( { {(kte_-kts_)+1} },
                   { {&qv, &temp, &col_loc} }),
  kts(kts_), kte(kte_), timestepcount(timestepcount_), source_ind(source_ind_), force_abort(force_abort_)
{
  EKAT_REQUIRE_MSG(nk() >= 3 || (kte == 1 && kts == 1), "nk too small to use for col_loc");
}

CalcUpwindData::CalcUpwindData(
  Int kts_, Int kte_, Int kdir_, Int kbot_, Int k_qxtop_, Int num_arrays_, Real dt_sub_) :
  PhysicsTestData({ {(kte_ - kts_)+1, num_arrays_}, {(kte_ - kts_)+1} },
                  { {&vs, &qnx, &fluxes},           {&rho, &inv_rho, &inv_dz} }),
  kts(kts_), kte(kte_), kdir(kdir_), kbot(kbot_), k_qxtop(k_qxtop_), num_arrays(num_arrays_), dt_sub(dt_sub_)
{}

void CalcUpwindData::convert_to_ptr_arr(std::vector<Real*>& mem_space, Real**& fluxes_, Real**& vs_, Real**& qnx_)
{
  mem_space.resize(num_arrays*3);
  for (Int i = 0; i < num_arrays; ++i) {
    mem_space[i]              = fluxes + (i*nk());
    mem_space[i+num_arrays]   = vs     + (i*nk());
    mem_space[i+num_arrays*2] = qnx    + (i*nk());
  }
  fluxes_ = mem_space.data();
  vs_     = mem_space.data() + num_arrays;
  qnx_    = mem_space.data() + num_arrays*2;
}

GenSedData::GenSedData(
  Int kts_, Int kte_, Int kdir_, Int k_qxtop_, Int k_qxbot_, Int kbot_, Real Co_max_, Real dt_left_,
  Real prt_accum_, Int num_arrays_) :
  CalcUpwindData(kts_, kte_, kdir_, kbot_, k_qxtop_, num_arrays_, 0.0),
  Co_max(Co_max_), k_qxbot(k_qxbot_), dt_left(dt_left_), prt_accum(prt_accum_)
{ }

CloudSedData::CloudSedData(
  Int kts_, Int kte_, Int ktop_, Int kbot_, Int kdir_,
  Real dt_, Real inv_dt_, bool do_predict_nc_, Real precip_liq_surf_) :
  PhysicsTestData( { {(kte_ - kts_) + 1} },
                   { {&qc_incld, &rho, &inv_rho, &cld_frac_l, &acn, &inv_dz, &qc, &nc, &nc_incld, &mu_c, &lamc, &qc_tend, &nc_tend} }),
  kts(kts_), kte(kte_), ktop(ktop_), kbot(kbot_), kdir(kdir_),
  dt(dt_), inv_dt(inv_dt_), do_predict_nc(do_predict_nc_), precip_liq_surf(precip_liq_surf_)
{}

IceSedData::IceSedData(
  Int kts_, Int kte_, Int ktop_, Int kbot_, Int kdir_,
  Real dt_, Real inv_dt_, Real precip_ice_surf_) :
  PhysicsTestData( { {(kte_ - kts_) + 1} },
                   { {&rho, &inv_rho, &rhofaci, &cld_frac_i, &inv_dz, &qi, &qi_incld, &ni, &ni_incld, &qm, &qm_incld, &bm, &bm_incld, &qi_tend, &ni_tend} }),
  kts(kts_), kte(kte_), ktop(ktop_), kbot(kbot_), kdir(kdir_),
  dt(dt_), inv_dt(inv_dt_), precip_ice_surf(precip_ice_surf_)
{}

RainSedData::RainSedData(
  Int kts_, Int kte_, Int ktop_, Int kbot_, Int kdir_,
  Real dt_, Real inv_dt_, Real precip_liq_surf_) :
  PhysicsTestData({ {(kte_ - kts_) + 2} }, // extra real at end for precip_liq_flux, so just add 1 to all
                  { {&rho, &inv_rho, &rhofacr, &cld_frac_r, &inv_dz, &qr_incld, &qr, &nr, &nr_incld, &mu_r, &lamr, &qr_tend, &nr_tend, &precip_liq_flux} }),
  kts(kts_), kte(kte_), ktop(ktop_), kbot(kbot_), kdir(kdir_),
  dt(dt_), inv_dt(inv_dt_), precip_liq_surf(precip_liq_surf_)
{}

HomogeneousFreezingData::HomogeneousFreezingData(
  Int kts_, Int kte_, Int ktop_, Int kbot_, Int kdir_) :
  PhysicsTestData( { {(kte_ - kts_) + 1} },
                   { {&T_atm, &inv_exner, &latent_heat_fusion, &qc, &nc, &qr, &nr, &qi, &ni, &qm, &bm, &th_atm} }),
  kts(kts_), kte(kte_), ktop(ktop_), kbot(kbot_), kdir(kdir_)
{}

P3MainPart1Data::P3MainPart1Data(
  Int kts_, Int kte_, Int kbot_, Int ktop_, Int kdir_,
  bool do_predict_nc_, bool do_prescribed_CCN_, Real dt_, bool, bool) :
  PhysicsTestData( { {(kte_ - kts_) + 1} }, { {
    &pres, &dpres, &dz, &nc_nuceat_tend, &inv_exner, &exner, &inv_cld_frac_l, &inv_cld_frac_i, &inv_cld_frac_r, &latent_heat_vapor, &latent_heat_sublim, &latent_heat_fusion, &nccn_prescribed,
    &T_atm, &rho, &inv_rho, &qv_sat_l, &qv_sat_i, &qv_supersat_i, &rhofacr, &rhofaci,
    &acn, &qv, &th_atm, &qc, &nc, &qr, &nr, &qi, &ni, &qm, &bm, &qc_incld, &qr_incld, &qi_incld,
    &qm_incld, &nc_incld, &nr_incld, &ni_incld, &bm_incld} }),
  kts(kts_), kte(kte_), kbot(kbot_), ktop(ktop_), kdir(kdir_),
  do_predict_nc(do_predict_nc_), do_prescribed_CCN(do_prescribed_CCN_), dt(dt_)
{}

///////////////////////////////////////////////////////////////////////////////

P3MainPart2Data::P3MainPart2Data(
  Int kts_, Int kte_, Int kbot_, Int ktop_, Int kdir_,
  bool do_predict_nc_, bool do_prescribed_CCN_, Real dt_, Real, bool) :
  PhysicsTestData( { {(kte_ - kts_) + 1} }, { {
    &pres, &dpres, &dz, &nc_nuceat_tend, &inv_exner, &exner, &inv_cld_frac_l, &inv_cld_frac_i, &inv_cld_frac_r, &ni_activated, &inv_qc_relvar, &cld_frac_i, &cld_frac_l, &cld_frac_r, &qv_prev, &t_prev,
    &T_atm, &rho, &inv_rho, &qv_sat_l, &qv_sat_i, &qv_supersat_i, &rhofacr, &rhofaci, &acn,
    &qv, &th_atm, &qc, &nc, &qr, &nr, &qi, &ni, &qm, &bm, &latent_heat_vapor, &latent_heat_sublim, &latent_heat_fusion, &qc_incld, &qr_incld,
    &qi_incld, &qm_incld, &nc_incld, &nr_incld, &ni_incld, &bm_incld, &mu_c, &nu, &lamc, &cdist, &cdist1,
    &cdistr, &mu_r, &lamr, &logn0r, &qv2qi_depos_tend, &precip_total_tend, &nevapr, &qr_evap_tend, &vap_liq_exchange,
    &vap_ice_exchange, &liq_ice_exchange,
    &P3_qr2qv_evap, &P3_qi2qv_sublim, &P3_qc2qr_accret, &P3_qc2qr_autoconv, &P3_qv2qi_vapdep, &P3_qc2qi_berg, &P3_qc2qr_ice_shed, &P3_qc2qi_collect, &P3_qr2qi_collect, &P3_qc2qi_hetero_freeze, &P3_qr2qi_immers_freeze, &P3_qi2qr_melt,
    &pratot, &prctot} }),
  kts(kts_), kte(kte_), kbot(kbot_), ktop(ktop_), kdir(kdir_),
  do_predict_nc(do_predict_nc_), do_prescribed_CCN(do_prescribed_CCN_), dt(dt_), inv_dt(1 / dt)
{}

///////////////////////////////////////////////////////////////////////////////

P3MainPart3Data::P3MainPart3Data(
  Int kts_, Int kte_, Int kbot_, Int ktop_, Int kdir_) :
  PhysicsTestData( { {(kte_ - kts_) + 1} }, { {
    &inv_exner, &cld_frac_l, &cld_frac_r, &cld_frac_i,
    &rho, &inv_rho, &rhofaci,
    &qv, &th_atm, &qc, &nc, &qr, &nr, &qi, &ni, &qm, &bm, &latent_heat_vapor, &latent_heat_sublim,
    &mu_c, &nu, &lamc, &mu_r,
    &lamr, &vap_liq_exchange,
    &ze_rain, &ze_ice, &diag_vm_qi, &diag_eff_radius_qi, &diag_diam_qi, &rho_qi, &diag_equiv_reflectivity,
    &diag_eff_radius_qc, &diag_eff_radius_qr} }),
  kts(kts_), kte(kte_), kbot(kbot_), ktop(ktop_), kdir(kdir_)
{}

///////////////////////////////////////////////////////////////////////////////

P3MainData::P3MainData(
  Int its_, Int ite_, Int kts_, Int kte_, Int it_, Real dt_, bool do_predict_nc_, bool do_prescribed_CCN_, Real) :
  PhysicsTestData( { {(ite_ - its_) + 1, (kte_ - kts_) + 1}, {(ite_ - its_) + 1, (kte_ - kts_) + 2} }, { {
    &pres, &dz, &nc_nuceat_tend, &nccn_prescribed, &ni_activated, &dpres, &inv_exner, &cld_frac_i, &cld_frac_l, &cld_frac_r,
    &inv_qc_relvar, &qc, &nc, &qr, &nr, &qi, &qm, &ni, &bm, &qv, &th_atm, &qv_prev, &t_prev, 
    &diag_eff_radius_qc, &diag_eff_radius_qi, &diag_eff_radius_qr, &rho_qi, &mu_c, &lamc, &qv2qi_depos_tend, &precip_total_tend, &nevapr,
    &qr_evap_tend, &liq_ice_exchange,
    &P3_qr2qv_evap, &P3_qi2qv_sublim, &P3_qc2qr_accret, &P3_qc2qr_autoconv, &P3_qv2qi_vapdep, &P3_qc2qi_berg, &P3_qc2qr_ice_shed, &P3_qc2qi_collect, &P3_qr2qi_collect, &P3_qc2qi_hetero_freeze, &P3_qr2qi_immers_freeze, &P3_qi2qr_melt, &P3_qr_sed, &P3_qc_sed, &P3_qi_sed,
    &vap_liq_exchange, &vap_ice_exchange, &precip_liq_flux,
    &precip_ice_flux},
    {&precip_liq_surf, &precip_ice_surf} }), // these two are (ni, nk+1)
  its(its_), ite(ite_), kts(kts_), kte(kte_), it(it_), dt(dt_), do_predict_nc(do_predict_nc_), do_prescribed_CCN(do_prescribed_CCN_), use_hetfrz_classnuc(false),
  hetfrz_immersion_nucleation_tend(nullptr), hetfrz_contact_nucleation_tend(nullptr), hetfrz_deposition_nucleation_tend(nullptr) 
{}

void IceSupersatConservationData::randomize(std::mt19937_64& engine)
{
  std::uniform_real_distribution<Real> data_dist(0.0, 1.0);

  cld_frac_i         = data_dist(engine);
  qv                 = data_dist(engine);
  qv_sat_i           = data_dist(engine);
  latent_heat_sublim = data_dist(engine);
  t_atm              = data_dist(engine);
  dt                 = data_dist(engine);
  qi2qv_sublim_tend  = data_dist(engine);
  qr2qv_evap_tend    = data_dist(engine);
  qidep              = data_dist(engine);
  qinuc              = data_dist(engine);
}

void NcConservationData::randomize(std::mt19937_64& engine)
{
  std::uniform_real_distribution<Real> data_dist(0.0, 1.0);

  nc                       = data_dist(engine);
  nc_selfcollect_tend      = data_dist(engine);
  dt                       = data_dist(engine);
  nc_collect_tend          = data_dist(engine);
  nc2ni_immers_freeze_tend = data_dist(engine);
  nc_accret_tend           = data_dist(engine);
  nc2nr_autoconv_tend      = data_dist(engine);
}

void NrConservationData::randomize(std::mt19937_64& engine)
{
  std::uniform_real_distribution<Real> data_dist(0.0, 1.0);

  nr                       = data_dist(engine);
  ni2nr_melt_tend          = data_dist(engine);
  nr_ice_shed_tend         = data_dist(engine);
  ncshdc                   = data_dist(engine);
  nc2nr_autoconv_tend      = data_dist(engine);
  dt                       = data_dist(engine);
  nmltratio                = data_dist(engine);
  nr_collect_tend          = data_dist(engine);
  nr2ni_immers_freeze_tend = data_dist(engine);
  nr_selfcollect_tend      = data_dist(engine);
  nr_evap_tend             = data_dist(engine);
}

void NiConservationData::randomize(std::mt19937_64& engine)
{
  std::uniform_real_distribution<Real> data_dist(0.0, 1.0);

  ni                       = data_dist(engine);
  ni_nucleat_tend          = data_dist(engine);
  nr2ni_immers_freeze_tend = data_dist(engine);
  nc2ni_immers_freeze_tend = data_dist(engine);
  dt                       = data_dist(engine);
  ni2nr_melt_tend          = data_dist(engine);
  ni_sublim_tend           = data_dist(engine);
  ni_selfcollect_tend      = data_dist(engine);
}

void PreventLiqSupersaturationData::randomize(std::mt19937_64& engine)
{
  /*Create random test data which changes each invocation, yet is
    physically reasonable. Follows examples in common_physics_functions_tests.cpp.
    Note that rates must be chosen carefully to prevent qv and T_atm from going negative.
    Note also that I'm using crude approx of latent heats here because these tests shouldn't
    care about slight numerical inaccuracies. I'm hardcoding cp here because I can't figure out
    how to make physics_constants.hpp available here. It may also be nice/needed to add eps to
    "make T and q exactly zero" calculations to prevent neg values due to roundoff error.
  */

  // Construct random input data
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF pdf_qv(1e-5,1e-3),
    pdf_pres(0.1,102000),
    pdf_temp(200.0,300.0),
    //pdf_dt(0.1,300.), //since dt is Scalar, always gets set to 1st index of Pack... so can't use rand val here.
    pdf_rate(0.,1e-3);

  Real cp=1004; //approx cp is good enough for testing.

  pres                     = pdf_pres(engine);
  t_atm                    = pdf_temp(engine);
  qv                       = pdf_qv(engine);
  latent_heat_vapor        = 2.5e6; //approx val is good enough for testing
  latent_heat_sublim       = 2.838e6; //approx val is good enough for testing
  dt                       = 60; //pdf_dt(engine);

  //qv sinks: don't let qv go neg.
  qidep                    = std::min(pdf_rate(engine), qv/dt ); //don't let dep make qv neg by itself
  qinuc                    = std::min(pdf_rate(engine), qv/dt - qidep); //don't let qidep+qinuc make qv neg

  //qv sources: don't let T go neg.
  qi2qv_sublim_tend        = std::min(pdf_rate(engine), cp/latent_heat_sublim * t_atm/dt + qidep+qinuc ); //don't let sublim make T neg.

  qr2qv_evap_tend          = std::min(pdf_rate(engine),
				      cp/latent_heat_vapor*t_atm/dt
    				      +(qidep+qinuc-qi2qv_sublim_tend)*latent_heat_sublim/latent_heat_vapor );

  /*
  pres                     = data_dist(engine);
  t_atm                    = data_dist(engine);
  qv                       = data_dist(engine);
  latent_heat_vapor        = data_dist(engine);
  latent_heat_sublim       = data_dist(engine);
  dt                       = data_dist(engine);
  qidep                    = data_dist(engine);
  qinuc                    = data_dist(engine);
  qi2qv_sublim_tend        = data_dist(engine);
  qr2qv_evap_tend          = data_dist(engine);
  */
}

///////////////////////////////////////////////////////////////////////////////

//
// _host function definitions
//

template <typename T>
std::vector<T*> ptr_to_arr(T** data, int n)
{
  std::vector<T*> result(n);
  for (int i = 0; i < n; ++i) result[i] = data[i];

  return result;
}

template <int N>
void calc_first_order_upwind_step_host_impl(
  Int kts, Int kte, Int kdir, Int kbot, Int k_qxtop, Real dt_sub,
  Real* rho, Real* inv_rho, Real* inv_dz,
  Real** fluxes, Real** vs, Real** qnx)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack = typename P3F::Spack;
  using view_1d = typename P3F::view_1d<Spack>;
  using KT = typename P3F::KT;
  using ExeSpace = typename KT::ExeSpace;
  using MemberType = typename P3F::MemberType;
  using view_1d_ptr_array = typename P3F::view_1d_ptr_array<Spack, N>;
  using uview_1d = typename P3F::uview_1d<Spack>;

  EKAT_REQUIRE_MSG(kts == 1, "kts must be 1, got " << kts);

  // Adjust for 0-based indexing
  kts -= 1;
  kte -= 1;
  kbot -= 1;
  k_qxtop -= 1;

  const Int nk = (kte - kts) + 1;
  const Int nk_pack = ekat::npack<Spack>(nk);

  // Setup views
  std::vector<view_1d> temp_d(3);
  std::vector<view_1d> fluxes_d(N), vs_d(N), qnx_d(N);

  ekat::host_to_device({rho, inv_rho, inv_dz}, nk, temp_d);

  view_1d rho_d(temp_d[0]), inv_rho_d(temp_d[1]), inv_dz_d(temp_d[2]);

  ekat::host_to_device(ptr_to_arr((const Real**)fluxes, N), nk, fluxes_d);
  ekat::host_to_device(ptr_to_arr((const Real**)vs, N)    , nk, vs_d);
  ekat::host_to_device(ptr_to_arr((const Real**)qnx, N)   , nk, qnx_d);

  Kokkos::Array<view_1d, N> fluxes_a, vs_a, qnx_a;
  for (Int i = 0; i < N; ++i) {
    fluxes_a[i] = fluxes_d[i];
    vs_a[i]     = vs_d[i];
    qnx_a[i]    = qnx_d[i];
  }

  // Call core function from kernel
  auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    view_1d_ptr_array fluxes_ptr, vs_ptr, qnx_ptr;
    for (int i = 0; i < N; ++i) {
      fluxes_ptr[i] = (uview_1d*)(&fluxes_a[i]);
      vs_ptr[i]     = (uview_1d*)(&vs_a[i]);
      qnx_ptr[i]    = (uview_1d*)(&qnx_a[i]);
    }
    uview_1d urho_d(rho_d), uinv_rho_d(inv_rho_d), uinv_dz_d(inv_dz_d);
    P3F::calc_first_order_upwind_step<N>(urho_d, uinv_rho_d, uinv_dz_d, team, nk, kbot, k_qxtop, kdir, dt_sub, fluxes_ptr, vs_ptr, qnx_ptr);
  });

  // Sync back to host
  ekat::device_to_host(ptr_to_arr(fluxes, N), nk, fluxes_d);
  ekat::device_to_host(ptr_to_arr(qnx, N), nk, qnx_d);
}

template <int N>
void generalized_sedimentation_host_impl(
  Int kts, Int kte, Int kdir, Int k_qxtop, Int* k_qxbot, Int kbot, Real Co_max,
  Real* dt_left, Real* prt_accum, Real* inv_dz, Real* inv_rho, Real* rho,
  Real** vs, Real** fluxes, Real** qnx)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack = typename P3F::Spack;
  using Singlep = typename ekat::Pack<Real, 1>;
  using view_1d = typename P3F::view_1d<Spack>;
  using view_1ds = typename P3F::view_1d<Singlep>;
  using KT = typename P3F::KT;
  using ExeSpace = typename KT::ExeSpace;
  using MemberType = typename P3F::MemberType;
  using view_1d_ptr_array = typename P3F::view_1d_ptr_array<Spack, N>;
  using uview_1d = typename P3F::uview_1d<Spack>;
  using ekat::host_to_device;
  using ekat::device_to_host;

  EKAT_REQUIRE_MSG(kts == 1, "kts must be 1, got " << kts);

  // Adjust for 0-based indexing
  kts -= 1;
  kte -= 1;
  kbot -= 1;
  k_qxtop -= 1;
  *k_qxbot -= 1;

  const Int nk = (kte - kts) + 1;
  const Int nk_pack = ekat::npack<Spack>(nk);

  // Set up views
  std::vector<view_1d> temp_d(3);
  std::vector<view_1d> fluxes_d(N), vs_d(N), qnx_d(N);
  std::vector<view_1ds> scalar_temp(1);
  std::vector<Real> scalars = {*prt_accum, *dt_left, static_cast<Real>(*k_qxbot)};

  host_to_device({rho, inv_rho, inv_dz}, nk, temp_d);
  host_to_device({scalars.data()}, scalars.size(), scalar_temp);

  view_1d rho_d(temp_d[0]), inv_rho_d(temp_d[1]), inv_dz_d(temp_d[2]);
  view_1ds scalars_d(scalar_temp[0]);

  host_to_device(ptr_to_arr((const Real**)fluxes, N), nk, fluxes_d);
  host_to_device(ptr_to_arr((const Real**)vs, N)    , nk, vs_d);
  host_to_device(ptr_to_arr((const Real**)qnx, N)   , nk, qnx_d);

  Kokkos::Array<view_1d, N> fluxes_a, vs_a, qnx_a;
  for (Int i = 0; i < N; ++i) {
    fluxes_a[i] = fluxes_d[i];
    vs_a[i]     = vs_d[i];
    qnx_a[i]    = qnx_d[i];
  }

  // Call core function from kernel
  auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    view_1d_ptr_array fluxes_ptr, vs_ptr, qnx_ptr;
    for (int i = 0; i < N; ++i) {
      fluxes_ptr[i] = (uview_1d*)(&fluxes_a[i]);
      vs_ptr[i]     = (uview_1d*)(&vs_a[i]);
      qnx_ptr[i]    = (uview_1d*)(&qnx_a[i]);
    }
    uview_1d urho_d(rho_d), uinv_rho_d(inv_rho_d), uinv_dz_d(inv_dz_d);

    // Each thread needs their own copy, like we expect in the main program, or else we will hit
    // data race issues
    Real prt_accum_k = scalars_d(0)[0];
    Real dt_left_k   = scalars_d(1)[0];
    Int k_qxbot_k    = static_cast<int>(scalars_d(2)[0]);

    P3F::generalized_sedimentation<N>(urho_d, uinv_rho_d, uinv_dz_d, team, nk, k_qxtop, k_qxbot_k, kbot, kdir, Co_max, dt_left_k, prt_accum_k, fluxes_ptr, vs_ptr, qnx_ptr);

    scalars_d(0)[0] = prt_accum_k;
    scalars_d(1)[0] = dt_left_k;
    scalars_d(2)[0] = k_qxbot_k;
  });

  // Sync back to host
  device_to_host(ptr_to_arr(fluxes, N), nk, fluxes_d);
  device_to_host(ptr_to_arr(qnx, N), nk, qnx_d);
  device_to_host({scalars.data()}, scalars.size(), scalar_temp);

  // Set scalars
  *prt_accum = scalars[0];
  *dt_left   = scalars[1];
  *k_qxbot   = scalars[2] + 1;
}

void calc_first_order_upwind_step_host(
  Int kts, Int kte, Int kdir, Int kbot, Int k_qxtop, Real dt_sub,
  Real* rho, Real* inv_rho, Real* inv_dz,
  Int num_arrays, Real** fluxes, Real** vs, Real** qnx)
{
  if (num_arrays == 1) {
    calc_first_order_upwind_step_host_impl<1>(kts, kte, kdir, kbot, k_qxtop, dt_sub, rho, inv_rho, inv_dz, fluxes, vs, qnx);
  }
  else if (num_arrays == 2) {
    calc_first_order_upwind_step_host_impl<2>(kts, kte, kdir, kbot, k_qxtop, dt_sub, rho, inv_rho, inv_dz, fluxes, vs, qnx);
  }
  else if (num_arrays == 4) {
    calc_first_order_upwind_step_host_impl<4>(kts, kte, kdir, kbot, k_qxtop, dt_sub, rho, inv_rho, inv_dz, fluxes, vs, qnx);
  }
  else {
    EKAT_REQUIRE_MSG(false, "Unsupported num arrays in bridge calc_first_order_upwind_step_f: " << num_arrays);
  }
}

void generalized_sedimentation_host(
  Int kts, Int kte, Int kdir, Int k_qxtop, Int* k_qxbot, Int kbot, Real Co_max,
  Real* dt_left, Real* prt_accum, Real* inv_dz, Real* inv_rho, Real* rho,
  Int num_arrays, Real** vs, Real** fluxes, Real** qnx)
{
  if (num_arrays == 1) {
    generalized_sedimentation_host_impl<1>(kts, kte, kdir, k_qxtop, k_qxbot, kbot, Co_max, dt_left, prt_accum,
                                        inv_dz, inv_rho, rho, vs, fluxes, qnx);
  }
  else if (num_arrays == 2) {
    generalized_sedimentation_host_impl<2>(kts, kte, kdir, k_qxtop, k_qxbot, kbot, Co_max, dt_left, prt_accum,
                                        inv_dz, inv_rho, rho, vs, fluxes, qnx);
  }
  else if (num_arrays == 4) {
    generalized_sedimentation_host_impl<4>(kts, kte, kdir, k_qxtop, k_qxbot, kbot, Co_max, dt_left, prt_accum,
                                        inv_dz, inv_rho, rho, vs, fluxes, qnx);
  }
  else {
    EKAT_REQUIRE_MSG(false, "Unsupported num arrays in bridge calc_first_order_upwind_step_f: " << num_arrays);
  }
}

void cloud_sedimentation_host(
  Int kts, Int kte, Int ktop, Int kbot, Int kdir,
  Real* qc_incld, Real* rho, Real* inv_rho, Real* cld_frac_l, Real* acn, Real* inv_dz,
  Real dt, Real inv_dt, bool do_predict_nc,
  Real* qc, Real* nc, Real* nc_incld, Real* mu_c, Real* lamc, Real* precip_liq_surf, Real* qc_tend, Real* nc_tend)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack = typename P3F::Spack;
  using view_1d = typename P3F::view_1d<Spack>;
  using KT = typename P3F::KT;
  using ExeSpace = typename KT::ExeSpace;
  using MemberType = typename P3F::MemberType;

  EKAT_REQUIRE_MSG(kts == 1, "kts must be 1, got " << kts);

  // Adjust for 0-based indexing
  kts -= 1;
  kte -= 1;
  ktop -= 1;
  kbot -= 1;

  const Int nk = (kte - kts) + 1;
  const Int nk_pack = ekat::npack<Spack>(nk);

  // Set up views
  const auto dnu = P3F::p3_init().dnu_table_vals;

  std::vector<view_1d> temp_d(CloudSedData::NUM_ARRAYS);

  ekat::host_to_device({qc_incld, rho, inv_rho, cld_frac_l, acn, inv_dz, qc, nc, nc_incld, mu_c, lamc, qc_tend, nc_tend},
                       nk, temp_d);

  view_1d
    qc_incld_d(temp_d[0]),
    rho_d     (temp_d[1]),
    inv_rho_d (temp_d[2]),
    cld_frac_l_d   (temp_d[3]),
    acn_d     (temp_d[4]),
    inv_dz_d (temp_d[5]),
    qc_d      (temp_d[6]),
    nc_d      (temp_d[7]),
    nc_incld_d(temp_d[8]),
    mu_c_d    (temp_d[9]),
    lamc_d    (temp_d[10]),
    qc_tend_d (temp_d[11]),
    nc_tend_d (temp_d[12]);

  // Call core function from kernel
  auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, nk_pack);
  ekat::WorkspaceManager<Spack> wsm(rho_d.extent(0), 4, policy);
  Kokkos::parallel_reduce(policy, KOKKOS_LAMBDA(const MemberType& team, Real& precip_liq_surf_k) {

    P3F::cloud_sedimentation(
      qc_incld_d, rho_d, inv_rho_d, cld_frac_l_d, acn_d, inv_dz_d, dnu,
      team, wsm.get_workspace(team),
      nk, ktop, kbot, kdir, dt, inv_dt, do_predict_nc,
      qc_d, nc_d, nc_incld_d, mu_c_d, lamc_d, qc_tend_d, nc_tend_d,
      precip_liq_surf_k);

  }, *precip_liq_surf);

  // Sync back to host
  std::vector<view_1d> inout_views = {qc_d, nc_d, nc_incld_d, mu_c_d, lamc_d, qc_tend_d, nc_tend_d};
  ekat::device_to_host({qc, nc, nc_incld, mu_c, lamc, qc_tend, nc_tend}, nk, inout_views);
}

void ice_sedimentation_host(
  Int kts, Int kte, Int ktop, Int kbot, Int kdir,
  Real* rho, Real* inv_rho, Real* rhofaci, Real* cld_frac_i, Real* inv_dz,
  Real dt, Real inv_dt,
  Real* qi, Real* qi_incld, Real* ni, Real* qm, Real* qm_incld, Real* bm, Real* bm_incld,
  Real* ni_incld, Real* precip_ice_surf, Real* qi_tend, Real* ni_tend)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack      = typename P3F::Spack;
  using view_1d    = typename P3F::view_1d<Spack>;
  using KT         = typename P3F::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename P3F::MemberType;

  EKAT_REQUIRE_MSG(kts == 1, "kts must be 1, got " << kts);

  // Adjust for 0-based indexing
  kts  -= 1;
  kte  -= 1;
  ktop -= 1;
  kbot -= 1;

  const Int nk = (kte - kts) + 1;
  const Int nk_pack = ekat::npack<Spack>(nk);

  // Set up views
  std::vector<view_1d> temp_d(IceSedData::NUM_ARRAYS);

  ekat::host_to_device({rho, inv_rho, rhofaci, cld_frac_i, inv_dz, qi, qi_incld, ni, qm, qm_incld, bm, bm_incld, ni_incld, qi_tend, ni_tend},
                       nk, temp_d);

  view_1d
    rho_d        (temp_d[0]),
    inv_rho_d    (temp_d[1]),
    rhofaci_d    (temp_d[2]),
    cld_frac_i_d (temp_d[3]),
    inv_dz_d     (temp_d[4]),
    qi_d         (temp_d[5]),
    qi_incld_d   (temp_d[6]),
    ni_d         (temp_d[7]),
    qm_d         (temp_d[8]),
    qm_incld_d   (temp_d[9]),
    bm_d         (temp_d[10]),
    bm_incld_d   (temp_d[11]),
    ni_incld_d   (temp_d[12]),
    qi_tend_d    (temp_d[13]),
    ni_tend_d    (temp_d[14]);

  // Call core function from kernel
  auto ice_table_vals = P3F::p3_init().ice_table_vals;
  auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, nk_pack);
  ekat::WorkspaceManager<Spack> wsm(rho_d.extent(0), 6, policy);
  Real my_precip_ice_surf = 0;
  Kokkos::parallel_reduce(policy, KOKKOS_LAMBDA(const MemberType& team, Real& precip_ice_surf_k) {

    P3F::ice_sedimentation(
      rho_d, inv_rho_d, rhofaci_d, cld_frac_i_d, inv_dz_d,
      team, wsm.get_workspace(team),
      nk, ktop, kbot, kdir, dt, inv_dt,
      qi_d, qi_incld_d, ni_d, ni_incld_d, qm_d, qm_incld_d, bm_d, bm_incld_d,
      qi_tend_d, ni_tend_d, ice_table_vals,
      precip_ice_surf_k, P3F::P3Runtime());

  }, my_precip_ice_surf);
  *precip_ice_surf += my_precip_ice_surf;

  // Sync back to host
  std::vector<view_1d> inout_views = {qi_d, qi_incld_d, ni_d, ni_incld_d, qm_d, qm_incld_d,
                                      bm_d, bm_incld_d, qi_tend_d, ni_tend_d};
  ekat::device_to_host({qi, qi_incld, ni, ni_incld, qm, qm_incld, bm, bm_incld, qi_tend, ni_tend}, nk, inout_views);
}

void rain_sedimentation_host(
  Int kts, Int kte, Int ktop, Int kbot, Int kdir,
  Real* qr_incld, Real* rho, Real* inv_rho, Real* rhofacr, Real* cld_frac_r, Real* inv_dz,
  Real dt, Real inv_dt,
  Real* qr, Real* nr, Real* nr_incld, Real* mu_r, Real* lamr, Real* precip_liq_surf, Real* precip_liq_flux, Real* qr_tend, Real* nr_tend)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack      = typename P3F::Spack;
  using view_1d    = typename P3F::view_1d<Spack>;
  using KT         = typename P3F::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename P3F::MemberType;

  EKAT_REQUIRE_MSG(kts == 1, "kts must be 1, got " << kts);

  // Adjust for 0-based indexing
  kts  -= 1;
  kte  -= 1;
  ktop -= 1;
  kbot -= 1;

  const Int nk = (kte - kts) + 1;
  const Int nk_pack = ekat::npack<Spack>(nk);

  // Set up views
  std::vector<view_1d> temp_d(RainSedData::NUM_ARRAYS);
  std::vector<size_t> sizes(RainSedData::NUM_ARRAYS, nk);
  sizes[RainSedData::NUM_ARRAYS - 1] = nk+1;

  ekat::host_to_device({qr_incld, rho, inv_rho, rhofacr, cld_frac_r, inv_dz, qr, nr, nr_incld, mu_r, lamr, qr_tend, nr_tend, precip_liq_flux},
                       sizes, temp_d);

  view_1d
    qr_incld_d       (temp_d[0]),
    rho_d            (temp_d[1]),
    inv_rho_d        (temp_d[2]),
    rhofacr_d        (temp_d[3]),
    cld_frac_r_d     (temp_d[4]),
    inv_dz_d         (temp_d[5]),
    qr_d             (temp_d[6]),
    nr_d             (temp_d[7]),
    nr_incld_d       (temp_d[8]),
    mu_r_d           (temp_d[9]),
    lamr_d           (temp_d[10]),
    qr_tend_d        (temp_d[11]),
    nr_tend_d        (temp_d[12]),
    precip_liq_flux_d(temp_d[13]);

  // Call core function from kernel
  auto tables = P3F::p3_init();
  auto vn_table_vals = tables.vn_table_vals;
  auto vm_table_vals = tables.vm_table_vals;
  auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, nk_pack);
  ekat::WorkspaceManager<Spack> wsm(rho_d.extent(0), 4, policy);
  Real my_precip_liq_surf = 0;
  Kokkos::parallel_reduce(policy, KOKKOS_LAMBDA(const MemberType& team, Real& precip_liq_surf_k) {

    P3F::rain_sedimentation(
      rho_d, inv_rho_d, rhofacr_d, cld_frac_r_d, inv_dz_d, qr_incld_d,
      team, wsm.get_workspace(team), vn_table_vals, vm_table_vals,
      nk, ktop, kbot, kdir, dt, inv_dt,
      qr_d, nr_d, nr_incld_d, mu_r_d, lamr_d, precip_liq_flux_d, qr_tend_d, nr_tend_d,
      precip_liq_surf_k, P3F::P3Runtime());

  }, my_precip_liq_surf);
  *precip_liq_surf += my_precip_liq_surf;

  // Sync back to host
  std::vector<size_t> sizes_out(8, nk);
  sizes_out[7] = nk+1;

  std::vector<view_1d> inout_views = {qr_d, nr_d, nr_incld_d, mu_r_d, lamr_d, qr_tend_d, nr_tend_d, precip_liq_flux_d};
  ekat::device_to_host({qr, nr, nr_incld, mu_r, lamr, qr_tend, nr_tend, precip_liq_flux}, sizes_out, inout_views);
}

void homogeneous_freezing_host(
  Int kts, Int kte, Int ktop, Int kbot, Int kdir,
  Real* T_atm, Real* inv_exner,
  Real* qc, Real* nc, Real* qr, Real* nr, Real* qi, Real* ni, Real* qm, Real* bm, Real* th_atm)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack      = typename P3F::Spack;
  using view_1d    = typename P3F::view_1d<Spack>;
  using KT         = typename P3F::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename P3F::MemberType;

  EKAT_REQUIRE_MSG(kts == 1, "kts must be 1, got " << kts);

  // Adjust for 0-based indexing
  kts  -= 1;
  kte  -= 1;
  ktop -= 1;
  kbot -= 1;

  const Int nk = (kte - kts) + 1;
  const Int nk_pack = ekat::npack<Spack>(nk);

  // Set up views
  std::vector<view_1d> temp_d(HomogeneousFreezingData::NUM_ARRAYS-1);

  ekat::host_to_device({T_atm, inv_exner, qc, nc, qr, nr, qi, ni, qm, bm, th_atm},
                       nk, temp_d);

  int current_index = 0;
  view_1d
    t_d                   (temp_d[current_index++]),
    inv_exner_d           (temp_d[current_index++]),
    qc_d                  (temp_d[current_index++]),
    nc_d                  (temp_d[current_index++]),
    qr_d                  (temp_d[current_index++]),
    nr_d                  (temp_d[current_index++]),
    qi_d                  (temp_d[current_index++]),
    ni_d                  (temp_d[current_index++]),
    qm_d                  (temp_d[current_index++]),
    bm_d                  (temp_d[current_index++]),
    th_atm_d              (temp_d[current_index++]);

  // Call core function from kernel
  auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {

    P3F::homogeneous_freezing(
      t_d, inv_exner_d,
      team,
      nk, ktop, kbot, kdir,
      qc_d, nc_d, qr_d, nr_d, qi_d, ni_d, qm_d, bm_d, th_atm_d);
  });

  // Sync back to host
  std::vector<view_1d> inout_views = {qc_d, nc_d, qr_d, nr_d, qi_d, ni_d, qm_d, bm_d, th_atm_d};

  ekat::device_to_host({qc, nc, qr, nr, qi, ni, qm, bm, th_atm}, nk, inout_views);
}

void check_values_host(Real* qv, Real* temp, Int kstart, Int kend,
                    Int timestepcount, bool force_abort, Int source_ind, Real* col_loc)
{
  using P3F        = Functions<Real, DefaultDevice>;
  using Spack      = typename P3F::Spack;
  using view_1d    = typename P3F::view_1d<Spack>;
  using suview_1d  = typename P3F::uview_1d<Real>;
  using KT         = typename P3F::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename P3F::MemberType;

  EKAT_REQUIRE_MSG(kend > kstart,
                    "ktop must be larger than kstart, kstart, kend " << kend << kstart);

  kstart -= 1;
  kend -= 1;
  const auto nk = (kend - kstart) + 1;
  const Int nk_pack = ekat::npack<Spack>(nk);
  std::vector<view_1d> cvd_d(CheckValuesData::NUM_ARRAYS+1);

  ekat::host_to_device<Int>({qv, temp, col_loc}, {nk, nk, 3}, cvd_d);

  view_1d qv_d(cvd_d[0]), temp_d(cvd_d[1]), col_loc_d(cvd_d[2]);
  suview_1d ucol_loc_d(reinterpret_cast<Real*>(col_loc_d.data()), 3);

  auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {

    P3F::check_values(qv_d, temp_d, kstart, kend, timestepcount, force_abort, source_ind, team,
                      ucol_loc_d);
  });
}

void p3_main_part1_host(
  Int kts, Int kte, Int kbot, Int ktop, Int kdir,
  bool do_predict_nc, bool do_prescribed_CCN,
  Real dt,
  Real* pres, Real* dpres, Real* dz, Real* nc_nuceat_tend, Real* nccn_prescribed, Real* inv_exner, Real* exner, Real* inv_cld_frac_l, Real* inv_cld_frac_i,
  Real* inv_cld_frac_r,
  Real* T_atm, Real* rho, Real* inv_rho, Real* qv_sat_l, Real* qv_sat_i, Real* qv_supersat_i, Real* rhofacr, Real* rhofaci,
  Real* acn, Real* qv, Real* th_atm, Real* qc, Real* nc, Real* qr, Real* nr, Real* qi, Real* ni, Real* qm, Real* bm, Real* qc_incld, Real* qr_incld, Real* qi_incld,
  Real* qm_incld, Real* nc_incld, Real* nr_incld, Real* ni_incld, Real* bm_incld,
  bool* is_nucleat_possible, bool* is_hydromet_present)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack      = typename P3F::Spack;
  using view_1d    = typename P3F::view_1d<Spack>;
  using bview_1d   = typename P3F::view_1d<bool>;
  using KT         = typename P3F::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename P3F::MemberType;

  EKAT_REQUIRE_MSG(kts == 1, "kts must be 1, got " << kts);

  // Adjust for 0-based indexing
  kts  -= 1;
  kte  -= 1;
  ktop -= 1;
  kbot -= 1;

  const Int nk = (kte - kts) + 1;
  const Int nk_pack = ekat::npack<Spack>(nk);

  // Set up views
  std::vector<view_1d> temp_d(P3MainPart1Data::NUM_ARRAYS);

  ekat::host_to_device({pres, dpres, dz, nc_nuceat_tend, inv_exner, exner, inv_cld_frac_l, inv_cld_frac_i, inv_cld_frac_r,
        T_atm, rho, inv_rho, qv_sat_l, qv_sat_i, qv_supersat_i, rhofacr, rhofaci,
        acn, qv, th_atm, qc, nc, qr, nr, qi, ni, qm, bm, qc_incld, qr_incld, qi_incld,
        qm_incld, nc_incld, nr_incld, ni_incld, bm_incld, nccn_prescribed},
    nk, temp_d);

  view_1d
    pres_d               (temp_d[0]),
    dpres_d              (temp_d[1]),
    dz_d                 (temp_d[2]),
    nc_nuceat_tend_d     (temp_d[3]),
    inv_exner_d          (temp_d[4]),
    exner_d              (temp_d[5]),
    inv_cld_frac_l_d     (temp_d[6]),
    inv_cld_frac_i_d     (temp_d[7]),
    inv_cld_frac_r_d     (temp_d[8]),
    t_d                  (temp_d[9]),
    rho_d                (temp_d[10]),
    inv_rho_d            (temp_d[11]),
    qv_sat_l_d           (temp_d[12]),
    qv_sat_i_d           (temp_d[13]),
    qv_supersat_i_d      (temp_d[14]),
    rhofacr_d            (temp_d[15]),
    rhofaci_d            (temp_d[16]),
    acn_d                (temp_d[17]),
    qv_d                 (temp_d[18]),
    th_atm_d             (temp_d[19]),
    qc_d                 (temp_d[20]),
    nc_d                 (temp_d[21]),
    qr_d                 (temp_d[22]),
    nr_d                 (temp_d[23]),
    qi_d                 (temp_d[24]),
    ni_d                 (temp_d[25]),
    qm_d                 (temp_d[26]),
    bm_d                 (temp_d[27]),
    qc_incld_d           (temp_d[28]),
    qr_incld_d           (temp_d[29]),
    qi_incld_d           (temp_d[30]),
    qm_incld_d           (temp_d[31]),
    nc_incld_d           (temp_d[32]),
    nr_incld_d           (temp_d[33]),
    ni_incld_d           (temp_d[34]),
    bm_incld_d           (temp_d[35]),
    nccn_prescribed_d    (temp_d[36]);

  // Call core function from kernel
  bview_1d bools_d("bools", 2);
  auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {

    P3F::p3_main_part1(
      team, nk, do_predict_nc, do_prescribed_CCN, dt,
      pres_d, dpres_d, dz_d, nc_nuceat_tend_d, nccn_prescribed_d, inv_exner_d, exner_d, inv_cld_frac_l_d, inv_cld_frac_i_d,
      inv_cld_frac_r_d,
      t_d, rho_d, inv_rho_d, qv_sat_l_d, qv_sat_i_d, qv_supersat_i_d, rhofacr_d, rhofaci_d,
      acn_d, qv_d, th_atm_d, qc_d, nc_d, qr_d, nr_d, qi_d, ni_d, qm_d, bm_d, qc_incld_d, qr_incld_d, qi_incld_d,
      qm_incld_d, nc_incld_d, nr_incld_d, ni_incld_d, bm_incld_d,
      bools_d(0), bools_d(1), P3F::P3Runtime());
  });

  // Sync back to host
  std::vector<view_1d> inout_views = {
    t_d, rho_d, inv_rho_d, qv_sat_l_d, qv_sat_i_d, qv_supersat_i_d, rhofacr_d, rhofaci_d,
    acn_d, qv_d, th_atm_d, qc_d, nc_d, qr_d, nr_d, qi_d, ni_d, qm_d, bm_d, qc_incld_d, qr_incld_d, qi_incld_d,
    qm_incld_d, nc_incld_d, nr_incld_d, ni_incld_d, bm_incld_d};

  ekat::device_to_host({T_atm, rho, inv_rho, qv_sat_l, qv_sat_i, qv_supersat_i, rhofacr, rhofaci,
        acn, qv, th_atm, qc, nc, qr, nr, qi, ni, qm, bm, qc_incld, qr_incld, qi_incld,
        qm_incld, nc_incld, nr_incld, ni_incld, bm_incld},
    nk, inout_views);

  const auto bools_h = Kokkos::create_mirror_view(bools_d);
  Kokkos::deep_copy(bools_h, bools_d);

  *is_nucleat_possible  = bools_h(0);
  *is_hydromet_present = bools_h(1);
}

void p3_main_part2_host(
  Int kts, Int kte, Int kbot, Int ktop, Int kdir, bool do_predict_nc, bool do_prescribed_CCN, Real dt, Real inv_dt,
  const Real *hetfrz_immersion_nucleation_tend, const Real *hetfrz_contact_nucleation_tend, const Real *hetfrz_deposition_nucleation_tend,
  Real* pres, Real* dpres, Real* dz, Real* nc_nuceat_tend, Real* inv_exner, Real* exner, Real* inv_cld_frac_l, Real* inv_cld_frac_i,
  Real* inv_cld_frac_r, Real* ni_activated, Real* inv_qc_relvar, Real* cld_frac_i, Real* cld_frac_l, Real* cld_frac_r, Real* qv_prev, Real* t_prev,
  Real* T_atm, Real* rho, Real* inv_rho, Real* qv_sat_l, Real* qv_sat_i, Real* qv_supersat_i, Real* rhofacr, Real* rhofaci,
  Real* acn, Real* qv, Real* th_atm, Real* qc, Real* nc, Real* qr, Real* nr, Real* qi, Real* ni,
  Real* qm, Real* bm, Real* qc_incld, Real* qr_incld, Real* qi_incld, Real* qm_incld, Real* nc_incld, Real* nr_incld,
  Real* ni_incld, Real* bm_incld, Real* mu_c, Real* nu, Real* lamc, Real* cdist, Real* cdist1, Real* cdistr, Real* mu_r, Real* lamr, Real* logn0r, Real* qv2qi_depos_tend, Real* precip_total_tend,
  Real* nevapr, Real* qr_evap_tend, Real* vap_liq_exchange, Real* vap_ice_exchange, Real* liq_ice_exchange,
  Real* P3_qr2qv_evap, Real* P3_qi2qv_sublim, Real* P3_qc2qr_accret, Real* P3_qc2qr_autoconv, Real* P3_qv2qi_vapdep, Real* P3_qc2qi_berg, Real* P3_qc2qr_ice_shed, Real* P3_qc2qi_collect, Real* P3_qr2qi_collect, Real* P3_qc2qi_hetero_freeze, Real* P3_qr2qi_immers_freeze, Real* P3_qi2qr_melt,
  Real* pratot,
  Real* prctot, bool* is_hydromet_present)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack      = typename P3F::Spack;
  using view_1d    = typename P3F::view_1d<Spack>;
  using bview_1d   = typename P3F::view_1d<bool>;
  using KT         = typename P3F::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename P3F::MemberType;

  EKAT_REQUIRE_MSG(kts == 1, "kts must be 1, got " << kts);

  // Adjust for 0-based indexing
  kts  -= 1;
  kte  -= 1;
  ktop -= 1;
  kbot -= 1;

  const Int nk = (kte - kts) + 1;
  const Int nk_pack = ekat::npack<Spack>(nk);
  const Real max_total_ni = 740.0e3;  // Hard-code this value for F90 comparison

  // Set up views
  std::vector<view_1d> temp_d(P3MainPart2Data::NUM_ARRAYS);
  std::vector<Real> hetfrz_0(nk,0), hetfrz_1(nk,0), hetfrz_2(nk,0);
  hetfrz_immersion_nucleation_tend = hetfrz_0.data();
  hetfrz_contact_nucleation_tend = hetfrz_1.data();
  hetfrz_deposition_nucleation_tend = hetfrz_2.data();

  ekat::host_to_device({hetfrz_immersion_nucleation_tend, hetfrz_contact_nucleation_tend, hetfrz_deposition_nucleation_tend,
        pres, dpres, dz, nc_nuceat_tend, inv_exner, exner, inv_cld_frac_l, inv_cld_frac_i, inv_cld_frac_r, ni_activated, inv_qc_relvar, cld_frac_i, cld_frac_l, cld_frac_r,
        T_atm, rho, inv_rho, qv_sat_l, qv_sat_i, qv_supersat_i, rhofacr, rhofaci, acn,
        qv, th_atm, qc, nc, qr, nr, qi, ni, qm, bm, qc_incld, qr_incld,
        qi_incld, qm_incld, nc_incld, nr_incld, ni_incld, bm_incld, mu_c, nu, lamc, cdist, cdist1,
        cdistr, mu_r, lamr, logn0r, qv2qi_depos_tend, precip_total_tend, nevapr, qr_evap_tend, vap_liq_exchange,
        vap_ice_exchange, liq_ice_exchange,
        P3_qr2qv_evap, P3_qi2qv_sublim, P3_qc2qr_accret, P3_qc2qr_autoconv, P3_qv2qi_vapdep, P3_qc2qi_berg, P3_qc2qr_ice_shed, P3_qc2qi_collect, P3_qr2qi_collect, P3_qc2qi_hetero_freeze, P3_qr2qi_immers_freeze, P3_qi2qr_melt,
        pratot, prctot, qv_prev, t_prev
        },
    nk, temp_d);

  int current_index = 0;
  view_1d
    hetfrz_immersion_nucleation_tend_d(temp_d[current_index++]), 
    hetfrz_contact_nucleation_tend_d(temp_d[current_index++]), 
    hetfrz_deposition_nucleation_tend_d(temp_d[current_index++]),
    pres_d              (temp_d[current_index++]),
    dpres_d             (temp_d[current_index++]),
    dz_d                (temp_d[current_index++]),
    nc_nuceat_tend_d    (temp_d[current_index++]),
    inv_exner_d         (temp_d[current_index++]),
    exner_d             (temp_d[current_index++]),
    inv_cld_frac_l_d    (temp_d[current_index++]),
    inv_cld_frac_i_d    (temp_d[current_index++]),
    inv_cld_frac_r_d    (temp_d[current_index++]),
    ni_activated_d      (temp_d[current_index++]),
    inv_qc_relvar_d     (temp_d[current_index++]),
    cld_frac_i_d        (temp_d[current_index++]),
    cld_frac_l_d        (temp_d[current_index++]),
    cld_frac_r_d        (temp_d[current_index++]),
    t_d                 (temp_d[current_index++]),
    rho_d               (temp_d[current_index++]),
    inv_rho_d           (temp_d[current_index++]),
    qv_sat_l_d          (temp_d[current_index++]),
    qv_sat_i_d          (temp_d[current_index++]),
    qv_supersat_i_d     (temp_d[current_index++]),
    rhofacr_d           (temp_d[current_index++]),
    rhofaci_d           (temp_d[current_index++]),
    acn_d               (temp_d[current_index++]),
    qv_d                (temp_d[current_index++]),
    th_atm_d            (temp_d[current_index++]),
    qc_d                (temp_d[current_index++]),
    nc_d                (temp_d[current_index++]),
    qr_d                (temp_d[current_index++]),
    nr_d                (temp_d[current_index++]),
    qi_d                (temp_d[current_index++]),
    ni_d                (temp_d[current_index++]),
    qm_d                (temp_d[current_index++]),
    bm_d                (temp_d[current_index++]),
    qc_incld_d          (temp_d[current_index++]),
    qr_incld_d          (temp_d[current_index++]),
    qi_incld_d          (temp_d[current_index++]),
    qm_incld_d          (temp_d[current_index++]),
    nc_incld_d          (temp_d[current_index++]),
    nr_incld_d          (temp_d[current_index++]),
    ni_incld_d          (temp_d[current_index++]),
    bm_incld_d          (temp_d[current_index++]),
    mu_c_d              (temp_d[current_index++]),
    nu_d                (temp_d[current_index++]),
    lamc_d              (temp_d[current_index++]),
    cdist_d             (temp_d[current_index++]),
    cdist1_d            (temp_d[current_index++]),
    cdistr_d            (temp_d[current_index++]),
    mu_r_d              (temp_d[current_index++]),
    lamr_d              (temp_d[current_index++]),
    logn0r_d            (temp_d[current_index++]),
    qv2qi_depos_tend_d  (temp_d[current_index++]),
    precip_total_tend_d (temp_d[current_index++]),
    nevapr_d            (temp_d[current_index++]),
    qr_evap_tend_d      (temp_d[current_index++]),
    vap_liq_exchange_d  (temp_d[current_index++]),
    vap_ice_exchange_d  (temp_d[current_index++]),
    liq_ice_exchange_d  (temp_d[current_index++]),
    P3_qr2qv_evap_d     (temp_d[current_index++]),
    P3_qi2qv_sublim_d   (temp_d[current_index++]),
    P3_qc2qr_accret_d   (temp_d[current_index++]),
    P3_qc2qr_autoconv_d (temp_d[current_index++]),
    P3_qv2qi_vapdep_d   (temp_d[current_index++]),
    P3_qc2qi_berg_d     (temp_d[current_index++]),
    P3_qc2qr_ice_shed_d (temp_d[current_index++]),
    P3_qc2qi_collect_d  (temp_d[current_index++]),
    P3_qr2qi_collect_d  (temp_d[current_index++]),
    P3_qc2qi_hetero_freeze_d (temp_d[current_index++]),
    P3_qr2qi_immers_freeze_d (temp_d[current_index++]),
    P3_qi2qr_melt_d     (temp_d[current_index++]),
    pratot_d            (temp_d[current_index++]),
    prctot_d            (temp_d[current_index++]),
    qv_prev_d           (temp_d[current_index++]),
    t_prev_d            (temp_d[current_index++]);

  // Call core function from kernel
  auto tables = P3F::p3_init();
  const auto dnu         = tables.dnu_table_vals;
  const auto ice_table_vals        = tables.ice_table_vals;
  const auto collect_table_vals     = tables.collect_table_vals;
  const auto revap_table_vals = tables.revap_table_vals;
  bview_1d bools_d("bools", 1);
  auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {

    P3F::p3_main_part2(
      team, nk_pack, max_total_ni, do_predict_nc, do_prescribed_CCN, dt, inv_dt, 
      hetfrz_immersion_nucleation_tend_d, hetfrz_contact_nucleation_tend_d, hetfrz_deposition_nucleation_tend_d, 
      dnu, ice_table_vals, collect_table_vals, revap_table_vals,
      pres_d, dpres_d, dz_d, nc_nuceat_tend_d, inv_exner_d, exner_d, inv_cld_frac_l_d,
      inv_cld_frac_i_d, inv_cld_frac_r_d, ni_activated_d, inv_qc_relvar_d, cld_frac_i_d, cld_frac_l_d, cld_frac_r_d,
      qv_prev_d, t_prev_d, t_d, rho_d, inv_rho_d, qv_sat_l_d, qv_sat_i_d, qv_supersat_i_d, rhofacr_d, rhofaci_d, acn_d,
      qv_d, th_atm_d, qc_d, nc_d, qr_d, nr_d, qi_d, ni_d, qm_d, bm_d,
      qc_incld_d, qr_incld_d, qi_incld_d,
      qm_incld_d, nc_incld_d, nr_incld_d, ni_incld_d, bm_incld_d,
      mu_c_d, nu_d, lamc_d, cdist_d, cdist1_d, cdistr_d, mu_r_d, lamr_d,
      logn0r_d, qv2qi_depos_tend_d, precip_total_tend_d, nevapr_d, qr_evap_tend_d, vap_liq_exchange_d,
      vap_ice_exchange_d, liq_ice_exchange_d,
      P3_qr2qv_evap_d, P3_qi2qv_sublim_d, P3_qc2qr_accret_d, P3_qc2qr_autoconv_d, P3_qv2qi_vapdep_d,
      P3_qc2qi_berg_d, P3_qc2qr_ice_shed_d, P3_qc2qi_collect_d, P3_qr2qi_collect_d, P3_qc2qi_hetero_freeze_d,
      P3_qr2qi_immers_freeze_d, P3_qi2qr_melt_d,
      pratot_d, prctot_d, bools_d(0),nk, P3F::P3Runtime());
  });

  // Sync back to host. Skip intent in variables.
  std::vector<view_1d> inout_views = {
    t_d, rho_d, inv_rho_d, qv_sat_l_d, qv_sat_i_d, qv_supersat_i_d, rhofacr_d, rhofaci_d, acn_d,
    qv_d, th_atm_d, qc_d, nc_d, qr_d, nr_d, qi_d, ni_d, qm_d, bm_d,
    qc_incld_d, qr_incld_d, qi_incld_d, qm_incld_d,
    nc_incld_d, nr_incld_d, ni_incld_d, bm_incld_d, mu_c_d, nu_d, lamc_d,
    cdist_d, cdist1_d, cdistr_d, mu_r_d, lamr_d, logn0r_d, qv2qi_depos_tend_d, precip_total_tend_d,
    nevapr_d, qr_evap_tend_d, vap_liq_exchange_d, vap_ice_exchange_d,
    liq_ice_exchange_d, 
    P3_qr2qv_evap_d, P3_qi2qv_sublim_d, P3_qc2qr_accret_d, P3_qc2qr_autoconv_d, P3_qv2qi_vapdep_d,
    P3_qc2qi_berg_d, P3_qc2qr_ice_shed_d, P3_qc2qi_collect_d, P3_qr2qi_collect_d, P3_qc2qi_hetero_freeze_d,
    P3_qr2qi_immers_freeze_d, P3_qi2qr_melt_d,
    pratot_d, prctot_d
  };

  ekat::device_to_host({
      T_atm, rho, inv_rho, qv_sat_l, qv_sat_i, qv_supersat_i, rhofacr, rhofaci, acn, qv, th_atm, qc, nc,
      qr, nr, qi, ni, qm, bm, qc_incld, qr_incld,
      qi_incld, qm_incld, nc_incld, nr_incld, ni_incld, bm_incld,
      mu_c, nu, lamc, cdist, cdist1, cdistr, mu_r, lamr, logn0r, qv2qi_depos_tend, precip_total_tend,
      nevapr, qr_evap_tend, vap_liq_exchange, vap_ice_exchange, liq_ice_exchange,
      P3_qr2qv_evap, P3_qi2qv_sublim, P3_qc2qr_accret, P3_qc2qr_autoconv, P3_qv2qi_vapdep,
      P3_qc2qi_berg, P3_qc2qr_ice_shed, P3_qc2qi_collect, P3_qr2qi_collect, P3_qc2qi_hetero_freeze,
      P3_qr2qi_immers_freeze, P3_qi2qr_melt,
      pratot, prctot},
    nk, inout_views);

  const auto bools_h = Kokkos::create_mirror_view(bools_d);
  Kokkos::deep_copy(bools_h, bools_d);

  *is_hydromet_present = bools_h(0);
}

void p3_main_part3_host(
  Int kts, Int kte, Int kbot, Int ktop, Int kdir,
  Real* inv_exner, Real* cld_frac_l, Real* cld_frac_r, Real* cld_frac_i,
  Real* rho, Real* inv_rho, Real* rhofaci, Real* qv, Real* th_atm, Real* qc,
  Real* nc, Real* qr, Real* nr, Real* qi, Real* ni, Real* qm,
  Real* bm, Real* mu_c, Real* nu, Real* lamc,
  Real* mu_r, Real* lamr, Real* vap_liq_exchange, Real* ze_rain, Real* ze_ice,
  Real* diag_vm_qi, Real* diag_eff_radius_qi, Real* diag_diam_qi, Real* rho_qi,
  Real* diag_equiv_reflectivity, Real* diag_eff_radius_qc, Real* diag_eff_radius_qr)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack      = typename P3F::Spack;
  using view_1d    = typename P3F::view_1d<Spack>;
  using KT         = typename P3F::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename P3F::MemberType;

  EKAT_REQUIRE_MSG(kts == 1, "kts must be 1, got " << kts);

  // Adjust for 0-based indexing
  kts  -= 1;
  kte  -= 1;
  ktop -= 1;
  kbot -= 1;

  const Int nk = (kte - kts) + 1;
  const Int nk_pack = ekat::npack<Spack>(nk);
  const Real max_total_ni = 740.0e3;  // Hard-code this value for F90 comparison

  // Set up views
  std::vector<view_1d> temp_d(P3MainPart3Data::NUM_ARRAYS-2);

  ekat::host_to_device({
      inv_exner, cld_frac_l, cld_frac_r, cld_frac_i, rho, inv_rho, rhofaci, qv, th_atm, qc,
      nc, qr, nr, qi, ni, qm, bm, mu_c, nu, lamc, mu_r,
      lamr, vap_liq_exchange, ze_rain, ze_ice, diag_vm_qi, diag_eff_radius_qi, diag_diam_qi,
      rho_qi, diag_equiv_reflectivity, diag_eff_radius_qc, diag_eff_radius_qr},
    nk, temp_d);

  int current_index = 0;
  view_1d
    inv_exner_d                (temp_d[current_index++]),
    cld_frac_l_d               (temp_d[current_index++]),
    cld_frac_r_d               (temp_d[current_index++]),
    cld_frac_i_d               (temp_d[current_index++]),
    rho_d                      (temp_d[current_index++]),
    inv_rho_d                  (temp_d[current_index++]),
    rhofaci_d                  (temp_d[current_index++]),
    qv_d                       (temp_d[current_index++]),
    th_atm_d                   (temp_d[current_index++]),
    qc_d                       (temp_d[current_index++]),
    nc_d                       (temp_d[current_index++]),
    qr_d                       (temp_d[current_index++]),
    nr_d                       (temp_d[current_index++]),
    qi_d                       (temp_d[current_index++]),
    ni_d                       (temp_d[current_index++]),
    qm_d                       (temp_d[current_index++]),
    bm_d                       (temp_d[current_index++]),
    mu_c_d                     (temp_d[current_index++]),
    nu_d                       (temp_d[current_index++]),
    lamc_d                     (temp_d[current_index++]),
    mu_r_d                     (temp_d[current_index++]),
    lamr_d                     (temp_d[current_index++]),
    vap_liq_exchange_d         (temp_d[current_index++]),
    ze_rain_d                  (temp_d[current_index++]),
    ze_ice_d                   (temp_d[current_index++]),
    diag_vm_qi_d               (temp_d[current_index++]),
    diag_eff_radius_qi_d       (temp_d[current_index++]),
    diag_diam_qi_d             (temp_d[current_index++]),
    rho_qi_d                   (temp_d[current_index++]),
    diag_equiv_reflectivity_d  (temp_d[current_index++]),
    diag_eff_radius_qc_d       (temp_d[current_index++]),
    diag_eff_radius_qr_d       (temp_d[current_index++]);

  // Call core function from kernel
  auto tables = P3F::p3_init();
  const auto dnu            = tables.dnu_table_vals;
  const auto ice_table_vals = tables.ice_table_vals;
  auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, nk_pack);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {

    P3F::p3_main_part3(team, nk_pack, max_total_ni, dnu, ice_table_vals,
                       inv_exner_d, cld_frac_l_d, cld_frac_r_d, cld_frac_i_d, rho_d, inv_rho_d,
                       rhofaci_d, qv_d, th_atm_d, qc_d, nc_d, qr_d, nr_d,
                       qi_d, ni_d, qm_d, bm_d,
                       mu_c_d, nu_d, lamc_d, mu_r_d, lamr_d,
                       vap_liq_exchange_d, ze_rain_d, ze_ice_d,
                       diag_vm_qi_d, diag_eff_radius_qi_d, diag_diam_qi_d, rho_qi_d,
                       diag_equiv_reflectivity_d, diag_eff_radius_qc_d, diag_eff_radius_qr_d, P3F::P3Runtime());
  });

  // Sync back to host
  std::vector<view_1d> inout_views = {
    rho_d, inv_rho_d, rhofaci_d, qv_d, th_atm_d, qc_d, nc_d, qr_d, nr_d, qi_d,
    ni_d, qm_d, bm_d, mu_c_d, nu_d, lamc_d, mu_r_d,
    lamr_d, vap_liq_exchange_d, ze_rain_d, ze_ice_d, diag_vm_qi_d, diag_eff_radius_qi_d,
    diag_diam_qi_d, rho_qi_d, diag_equiv_reflectivity_d, diag_eff_radius_qc_d, diag_eff_radius_qr_d
  };

  ekat::device_to_host({
      rho, inv_rho, rhofaci, qv, th_atm, qc, nc, qr, nr, qi, ni, qm, bm,
      mu_c, nu, lamc, mu_r, lamr, vap_liq_exchange, ze_rain, ze_ice,
      diag_vm_qi, diag_eff_radius_qi, diag_diam_qi, rho_qi, diag_equiv_reflectivity, diag_eff_radius_qc,
      diag_eff_radius_qr
    },
    nk, inout_views);
}

Int p3_main_host(
  Real* qc, Real* nc, Real* qr, Real* nr, Real* th_atm, Real* qv, Real dt,
  Real* qi, Real* qm, Real* ni, Real* bm, Real* pres, Real* dz,
  Real* nc_nuceat_tend, Real* nccn_prescribed, Real* ni_activated, Real* inv_qc_relvar, Int it, Real* precip_liq_surf,
  Real* precip_ice_surf, Int its, Int ite, Int kts, Int kte, Real* diag_eff_radius_qc,
  Real* diag_eff_radius_qi, Real* diag_eff_radius_qr, Real* rho_qi, bool do_predict_nc, bool do_prescribed_CCN, bool use_hetfrz_classnuc, Real* dpres, Real* inv_exner,
  Real* qv2qi_depos_tend, Real* precip_liq_flux, Real* precip_ice_flux, Real* cld_frac_r, Real* cld_frac_l, Real* cld_frac_i,
  Real* liq_ice_exchange, Real* vap_liq_exchange, Real* vap_ice_exchange, 
  Real* P3_qr2qv_evap, Real* P3_qi2qv_sublim, Real* P3_qc2qr_accret, Real* P3_qc2qr_autoconv, Real* P3_qv2qi_vapdep, Real* P3_qc2qi_berg, Real* P3_qc2qr_ice_shed, Real* P3_qc2qi_collect, Real* P3_qr2qi_collect, Real* P3_qc2qi_hetero_freeze, Real* P3_qr2qi_immers_freeze, Real* P3_qi2qr_melt, Real* P3_qr_sed, Real* P3_qc_sed, Real* P3_qi_sed,
  Real* qv_prev, Real* t_prev)
{
  using P3F  = Functions<Real, DefaultDevice>;

  using Spack      = typename P3F::Spack;
  using KT         = typename P3F::KT;
  using view_2d    = typename P3F::view_2d<Spack>;
  using sview_1d   = typename P3F::view_1d<Real>;
  using sview_2d   = typename P3F::view_2d<Real>;

  EKAT_REQUIRE_MSG(its == 1, "its must be 1, got " << its);
  EKAT_REQUIRE_MSG(kts == 1, "kts must be 1, got " << kts);

  // Adjust for 0-based indexing
  its  -= 1;
  ite  -= 1;
  kts  -= 1;
  kte  -= 1;

  const Int nj    = (ite - its) + 1;
  const Int nk    = (kte - kts) + 1;

  // Set up views, pretend all views are input views for the sake of initializing kokkos views
  std::vector<view_2d> temp_d(P3MainData::NUM_ARRAYS);
  std::vector<size_t> dim1_sizes(P3MainData::NUM_ARRAYS, nj);
  std::vector<size_t> dim2_sizes(P3MainData::NUM_ARRAYS, nk);
  std::vector<const Real*> ptr_array = {
    pres, dz, nc_nuceat_tend, nccn_prescribed, ni_activated, dpres, inv_exner, cld_frac_i, cld_frac_l, cld_frac_r, inv_qc_relvar,
    qc, nc, qr, nr, qi, qm, ni, bm, qv, th_atm, qv_prev, t_prev, diag_eff_radius_qc, diag_eff_radius_qi, diag_eff_radius_qr,
    rho_qi, qv2qi_depos_tend,
    liq_ice_exchange, vap_liq_exchange, vap_ice_exchange, 
    
    P3_qr2qv_evap, P3_qi2qv_sublim, P3_qc2qr_accret, P3_qc2qr_autoconv, P3_qv2qi_vapdep, P3_qc2qi_berg, P3_qc2qr_ice_shed, P3_qc2qi_collect, P3_qr2qi_collect, P3_qc2qi_hetero_freeze, P3_qr2qi_immers_freeze, P3_qi2qr_melt, P3_qr_sed, P3_qc_sed, P3_qi_sed,
    precip_liq_flux, precip_ice_flux, 
    precip_liq_surf, precip_ice_surf
  };

  int dim_sizes_len = dim1_sizes.size();
  dim2_sizes[dim_sizes_len-4] = nk+1; // precip_liq_flux
  dim2_sizes[dim_sizes_len-3] = nk+1; // precip_ice_flux
  dim1_sizes[dim_sizes_len-2] = 1; dim2_sizes[dim_sizes_len-2] = nj; // precip_liq_surf
  dim1_sizes[dim_sizes_len-1] = 1; dim2_sizes[dim_sizes_len-1] = nj; // precip_ice_surf

  // Initialize outputs to avoid uninitialized read warnings in memory checkers
  for (size_t i = P3MainData::NUM_INPUT_ARRAYS; i < P3MainData::NUM_ARRAYS; ++i) {
    for (size_t j = 0; j < dim1_sizes[i]*dim2_sizes[i]; ++j) {
      const_cast<Real*>(ptr_array[i])[j] = 0;
    }
  }

  ekat::host_to_device(ptr_array, dim1_sizes, dim2_sizes, temp_d);

  int counter = 0;
  view_2d
    pres_d                 (temp_d[counter++]), //0
    dz_d                   (temp_d[counter++]),
    nc_nuceat_tend_d       (temp_d[counter++]),
    nccn_prescribed_d      (temp_d[counter++]),
    ni_activated_d         (temp_d[counter++]),
    dpres_d                (temp_d[counter++]), //5
    inv_exner_d            (temp_d[counter++]),
    cld_frac_i_d           (temp_d[counter++]),
    cld_frac_l_d           (temp_d[counter++]),
    cld_frac_r_d           (temp_d[counter++]),
    inv_qc_relvar_d        (temp_d[counter++]), //10
    qc_d                   (temp_d[counter++]),
    nc_d                   (temp_d[counter++]),
    qr_d                   (temp_d[counter++]),
    nr_d                   (temp_d[counter++]),
    qi_d                   (temp_d[counter++]), //15
    qm_d                   (temp_d[counter++]),
    ni_d                   (temp_d[counter++]),
    bm_d                   (temp_d[counter++]),
    qv_d                   (temp_d[counter++]),
    th_atm_d               (temp_d[counter++]), //20
    qv_prev_d              (temp_d[counter++]),
    t_prev_d               (temp_d[counter++]),
    diag_eff_radius_qc_d   (temp_d[counter++]),
    diag_eff_radius_qi_d   (temp_d[counter++]),
    diag_eff_radius_qr_d   (temp_d[counter++]), //25
    rho_qi_d               (temp_d[counter++]),
    qv2qi_depos_tend_d     (temp_d[counter++]),
    liq_ice_exchange_d     (temp_d[counter++]),
    vap_liq_exchange_d     (temp_d[counter++]),
    vap_ice_exchange_d     (temp_d[counter++]), //30
    P3_qr2qv_evap_d           (temp_d[counter++]), //35
    P3_qi2qv_sublim_d         (temp_d[counter++]),
    P3_qc2qr_accret_d         (temp_d[counter++]),
    P3_qc2qr_autoconv_d       (temp_d[counter++]),
    P3_qv2qi_vapdep_d         (temp_d[counter++]),
    P3_qc2qi_berg_d           (temp_d[counter++]), //40
    P3_qc2qr_ice_shed_d       (temp_d[counter++]),
    P3_qc2qi_collect_d        (temp_d[counter++]),
    P3_qr2qi_collect_d        (temp_d[counter++]),
    P3_qc2qi_hetero_freeze_d  (temp_d[counter++]),
    P3_qr2qi_immers_freeze_d  (temp_d[counter++]), //45
    P3_qi2qr_melt_d           (temp_d[counter++]), //46
    P3_qr_sed_d               (temp_d[counter++]),
    P3_qc_sed_d               (temp_d[counter++]),
    P3_qi_sed_d               (temp_d[counter++]), //49
    precip_liq_flux_d      (temp_d[counter++]),
    precip_ice_flux_d      (temp_d[counter++]), // 35
    precip_liq_surf_temp_d (temp_d[counter++]),
    precip_ice_surf_temp_d (temp_d[counter++]); 

  std::vector<Real> hetfrz_immersion_nucleation_tend(nj*nk, 0.0);
  std::vector<Real> hetfrz_contact_nucleation_tend(nj*nk, 0.0);
  std::vector<Real> hetfrz_deposition_nucleation_tend(nj*nk, 0.0);
  std::vector<const Real*> pointers = {
    hetfrz_immersion_nucleation_tend.data(), 
    hetfrz_contact_nucleation_tend.data(), 
    hetfrz_deposition_nucleation_tend.data()};
  std::vector<size_t> dim1(3, nj);
  std::vector<size_t> dim2(3, nk);
  std::vector<view_2d> view(3);
  ekat::host_to_device(pointers, dim1, dim2, view);
  view_2d hetfrz_immersion_nucleation_tend_d(view[0]);
  view_2d hetfrz_contact_nucleation_tend_d(view[1]);
  view_2d hetfrz_deposition_nucleation_tend_d(view[2]);

  // Special cases: precip_liq_surf=1d<scalar>(ni), precip_ice_surf=1d<scalar>(ni), col_location=2d<scalar>(nj, 3)
  sview_1d precip_liq_surf_d("precip_liq_surf_d", nj), precip_ice_surf_d("precip_ice_surf_d", nj);
  sview_2d col_location_d("col_location_d", nj, 3);

  view_2d mu_c_d("mu_c_d",nj,nk);
  view_2d lamc_d("lamc_d",nj,nk);
  view_2d precip_total_tend_d("precip_total_tend_d",nj,nk);
  view_2d nevapr_d("nevapr_d",nj,nk);
  view_2d diag_equiv_reflectivity_d("diag_equiv_reflectivity_d",nj,nk);
  view_2d qr_evap_tend_d("qr_evap_tend_d",nj,nk);

  Kokkos::parallel_for(nj, KOKKOS_LAMBDA(const Int& i) {
    precip_liq_surf_d(i) = precip_liq_surf_temp_d(0, i / Spack::n)[i % Spack::n];
    precip_ice_surf_d(i) = precip_ice_surf_temp_d(0, i / Spack::n)[i % Spack::n];

    for (int j = 0; j < 3; ++j) {
      col_location_d(i, j) = i+1;
    }
  });

  // Pack our data into structs and ship it off to p3_main.
  P3F::P3PrognosticState prog_state{qc_d, nc_d, qr_d, nr_d, qi_d, qm_d,
                                    ni_d, bm_d, qv_d, th_atm_d};
  P3F::P3DiagnosticInputs diag_inputs{nc_nuceat_tend_d, nccn_prescribed_d, ni_activated_d, inv_qc_relvar_d, cld_frac_i_d,
                                      cld_frac_l_d, cld_frac_r_d, pres_d, dz_d, dpres_d,
                                      inv_exner_d, qv_prev_d, t_prev_d, hetfrz_immersion_nucleation_tend_d, hetfrz_contact_nucleation_tend_d, hetfrz_deposition_nucleation_tend_d};
  P3F::P3DiagnosticOutputs diag_outputs{qv2qi_depos_tend_d, precip_liq_surf_d,
                                        precip_ice_surf_d, diag_eff_radius_qc_d, diag_eff_radius_qi_d, diag_eff_radius_qr_d,
                                        rho_qi_d,precip_liq_flux_d, precip_ice_flux_d, precip_total_tend_d, nevapr_d, diag_equiv_reflectivity_d};
  P3F::P3Infrastructure infrastructure{dt, it, its, ite, kts, kte,
                                       do_predict_nc, do_prescribed_CCN, col_location_d};
  P3F::P3HistoryOnly history_only{liq_ice_exchange_d, vap_liq_exchange_d,
                                  vap_ice_exchange_d,
                                  P3_qr2qv_evap_d, P3_qi2qv_sublim_d, P3_qc2qr_accret_d, P3_qc2qr_autoconv_d, P3_qv2qi_vapdep_d,
                                  P3_qc2qi_berg_d, P3_qc2qr_ice_shed_d, P3_qc2qi_collect_d, P3_qr2qi_collect_d, P3_qc2qi_hetero_freeze_d,
                                  P3_qr2qi_immers_freeze_d, P3_qi2qr_melt_d, P3_qr_sed_d, P3_qc_sed_d, P3_qi_sed_d
                                  };

  const Int nk_pack = ekat::npack<Spack>(nk);
#ifdef SCREAM_P3_SMALL_KERNELS
  view_2d
    mu_r("mu_r", nj, nk_pack), T_atm("T_atm", nj, nk_pack), lamr("lamr", nj, nk_pack), logn0r("logn0r", nj, nk_pack), nu("nu", nj, nk_pack),
    cdist("cdist", nj, nk_pack), cdist1("cdist1", nj, nk_pack), cdistr("cdistr", nj, nk_pack), inv_cld_frac_i("inv_cld_frac_i", nj, nk_pack),
    inv_cld_frac_l("inv_cld_frac_l", nj, nk_pack), inv_cld_frac_r("inv_cld_frac_r", nj, nk_pack), qc_incld("qc_incld", nj, nk_pack),
    qr_incld("qr_incld", nj, nk_pack), qi_incld("qi_incld", nj, nk_pack), qm_incld("qm_incld", nj, nk_pack), nc_incld("nc_incld", nj, nk_pack),
    nr_incld("nr_incld", nj, nk_pack), ni_incld("ni_incld", nj, nk_pack), bm_incld("bm_incld", nj, nk_pack), inv_dz("inv_dz", nj, nk_pack),
    inv_rho("inv_rho", nj, nk_pack), ze_ice("ze_ice", nj, nk_pack), ze_rain("ze_rain", nj, nk_pack), prec("prec", nj, nk_pack),
    rho("rho", nj, nk_pack), rhofacr("rhofacr", nj, nk_pack), rhofaci("rhofaci", nj, nk_pack),  acn("acn", nj, nk_pack), qv_sat_l("qv_sat", nj, nk_pack),
    qv_sat_i("qv_sat_i", nj, nk_pack), sup("sup", nj, nk_pack), qv_supersat_i("qv_supersat", nj, nk_pack), tmparr2("tmparr2", nj, nk_pack),
    exner("exner", nj, nk_pack), diag_vm_qi("diag_vm_qi", nj, nk_pack),
    diag_diam_qi("diag_diam_qi", nj, nk_pack), pratot("pratot", nj, nk_pack), prctot("prctot", nj, nk_pack), qtend_ignore("qtend_ignore", nj, nk_pack),
    ntend_ignore("ntend_ignore", nj, nk_pack), mu_c("mu_c", nj, nk_pack), lamc("lamc", nj, nk_pack), qr_evap_tend("qr_evap_tend", nj, nk_pack),
    v_qc("v_qc", nj, nk_pack), v_nc("v_nc", nj, nk_pack), flux_qx("flux_qx", nj, nk_pack), flux_nx("flux_nx", nj, nk_pack), v_qit("v_qit", nj, nk_pack),
    v_nit("v_nit", nj, nk_pack), flux_nit("flux_nit", nj, nk_pack), flux_bir("flux_bir", nj, nk_pack), flux_qir("flux_qir", nj, nk_pack),
    flux_qit("flux_qit", nj, nk_pack), v_qr("v_qr", nj, nk_pack), v_nr("v_nr", nj, nk_pack);

  P3F::P3Temporaries temporaries{
    mu_r, T_atm, lamr, logn0r, nu, cdist, cdist1, cdistr, inv_cld_frac_i,
    inv_cld_frac_l, inv_cld_frac_r, qc_incld, qr_incld, qi_incld, qm_incld,
    nc_incld, nr_incld, ni_incld, bm_incld, inv_dz, inv_rho, ze_ice, ze_rain,
    prec, rho, rhofacr, rhofaci, acn, qv_sat_l, qv_sat_i, sup, qv_supersat_i,
    tmparr2, exner, diag_vm_qi, diag_diam_qi,
    pratot, prctot, qtend_ignore, ntend_ignore, mu_c, lamc, qr_evap_tend,
    v_qc, v_nc, flux_qx, flux_nx, v_qit, v_nit, flux_nit, flux_bir, flux_qir,
    flux_qit, v_qr, v_nr
  };
#endif

  // load tables
  auto lookup_tables = P3F::p3_init();
  P3F::P3Runtime runtime_options{740.0e3};

  // Create local workspace
  const auto policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(nj, nk_pack);
  ekat::WorkspaceManager<Spack, KT::Device> workspace_mgr(nk_pack, 52, policy);

  auto elapsed_microsec = P3F::p3_main(runtime_options, prog_state, diag_inputs, diag_outputs, infrastructure,
                                       history_only, lookup_tables,
#ifdef SCREAM_P3_SMALL_KERNELS
                                       temporaries,
#endif
                                       workspace_mgr, nj, nk);

  Kokkos::parallel_for(nj, KOKKOS_LAMBDA(const Int& i) {
    precip_liq_surf_temp_d(0, i / Spack::n)[i % Spack::n] = precip_liq_surf_d(i);
    precip_ice_surf_temp_d(0, i / Spack::n)[i % Spack::n] = precip_ice_surf_d(i);
  });

  // Sync back to host
  std::vector<view_2d> inout_views = {
    qc_d, nc_d, qr_d, nr_d, qi_d, qm_d, ni_d, bm_d, qv_d, th_atm_d,
    diag_eff_radius_qc_d, diag_eff_radius_qi_d, diag_eff_radius_qr_d, rho_qi_d,
    qv2qi_depos_tend_d,
    liq_ice_exchange_d, vap_liq_exchange_d, vap_ice_exchange_d,
    P3_qr2qv_evap_d, P3_qi2qv_sublim_d, P3_qc2qr_accret_d, P3_qc2qr_autoconv_d, P3_qv2qi_vapdep_d,
    P3_qc2qi_berg_d, P3_qc2qr_ice_shed_d, P3_qc2qi_collect_d, P3_qr2qi_collect_d, P3_qc2qi_hetero_freeze_d,
    P3_qr2qi_immers_freeze_d, P3_qi2qr_melt_d, P3_qr_sed_d, P3_qc_sed_d, P3_qi_sed_d,
    precip_liq_flux_d, precip_ice_flux_d, precip_liq_surf_temp_d, precip_ice_surf_temp_d
  };
  std::vector<size_t> dim1_sizes_out(P3MainData::NUM_ARRAYS - 13, nj);
  std::vector<size_t> dim2_sizes_out(P3MainData::NUM_ARRAYS - 13, nk);
  int dim_sizes_out_len = dim1_sizes_out.size();
  dim2_sizes_out[dim_sizes_out_len-4] = nk+1; // precip_liq_flux
  dim2_sizes_out[dim_sizes_out_len-3] = nk+1; // precip_ice_flux
  dim1_sizes_out[dim_sizes_out_len-2] = 1; dim2_sizes_out[dim_sizes_out_len-2] = nj; // precip_liq_surf
  dim1_sizes_out[dim_sizes_out_len-1] = 1; dim2_sizes_out[dim_sizes_out_len-1] = nj; // precip_ice_surf

  ekat::device_to_host({
      qc, nc, qr, nr, qi, qm, ni, bm, qv, th_atm, diag_eff_radius_qc, diag_eff_radius_qi, diag_eff_radius_qr,
      rho_qi, qv2qi_depos_tend,
      liq_ice_exchange, vap_liq_exchange, vap_ice_exchange,
      P3_qr2qv_evap, P3_qi2qv_sublim, P3_qc2qr_accret, P3_qc2qr_autoconv, P3_qv2qi_vapdep,
      P3_qc2qi_berg, P3_qc2qr_ice_shed, P3_qc2qi_collect, P3_qr2qi_collect, P3_qc2qi_hetero_freeze,
      P3_qr2qi_immers_freeze, P3_qi2qr_melt, P3_qr_sed, P3_qc_sed, P3_qi_sed,
      precip_liq_flux, precip_ice_flux, precip_liq_surf, precip_ice_surf
    },
    dim1_sizes_out, dim2_sizes_out, inout_views);

  return elapsed_microsec;
}

} // namespace p3
} // namespace scream
