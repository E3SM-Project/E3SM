#include "zm_test_data.hpp"
#include "zm_functions.hpp"

#include <ekat_pack_kokkos.hpp>
#include <ekat_subview_utils.hpp>

#include <random>

using scream::Real;
using scream::Int;

//
// A C++ interface to ZM fortran calls and vice versa
//

namespace scream {
namespace zm {

using ZMF = Functions<Real, DefaultDevice>;
using ZMC = typename ZMF::ZMC;

using ExeSpace   = typename ZMF::KT::ExeSpace;
using MemberType = typename ZMF::KT::MemberType;

using view0dr_d = ZMF::view_0d<Real>;
using view0di_d = ZMF::view_0d<Int>;
using view1di_d = ZMF::view_1d<Int>;
using view1db_d = ZMF::view_1d<bool>;
using view1dr_d = ZMF::view_1d<Real>;
using view2dr_d = ZMF::view_2d<Real>;
using view3dr_d = ZMF::view_3d<Real>;

using WSM = typename ZMF::WorkspaceManager;

extern "C" {

void zm_common_init_bridge_f();

void zm_common_finalize_bridge_f();

void ientropy_bridge_f(Real s, Real p, Real qt, Real* t, Real* qst, Real tfg);

void entropy_bridge_f(Real tk, Real p, Real qtot, Real* entropy);

void zm_transport_tracer_bridge_f(Int pcols, Int pver, bool* doconvtran, Real* q, Int ncnst, Real* mu, Real* md, Real* du, Real* eu, Real* ed, Real* dp, Int* jt, Int* mx, Int* ideep, Int il1g, Int il2g, Real* fracis, Real* dqdt, Real* dpdry, Real dt);

void zm_transport_momentum_bridge_f(Int pcols, Int ncol, Int pver, Int pverp, Real* wind_in, Int nwind, Real* mu, Real* md, Real* du, Real* eu, Real* ed, Real* dp, Int* jt, Int* mx, Int* ideep, Int il1g, Int il2g, Real* wind_tend, Real* pguall, Real* pgdall, Real* icwu, Real* icwd, Real dt, Real* seten);

void compute_dilute_cape_bridge_f(Int pcols, Int ncol, Int pver, Int pverp, Int num_cin, Int num_msg, Real* sp_humidity_in, Real* temperature_in, Real* zmid, Real* pmid, Real* pint, Int* pblt, Real* tpert, Real* parcel_temp, Real* parcel_qsat, Int* msemax_klev, Real* lcl_temperature, Int* lcl_klev, Int* eql_klev, Real* cape, bool calc_msemax_klev, Int* prev_msemax_klev, bool use_input_tq_mx, Real* q_mx, Real* t_mx);

void find_mse_max_bridge_f(Int pcols, Int ncol, Int pver, Int num_msg, Int* msemax_top_k, bool pergro_active, Real* temperature, Real* zmid, Real* sp_humidity, Int* msemax_klev, Real* mse_max_val);

void compute_dilute_parcel_bridge_f(Int pcols, Int ncol, Int pver, Int num_msg, Int* klaunch, Real* pmid, Real* temperature, Real* sp_humidity, Real* tpert, Int* pblt, Real* parcel_temp, Real* parcel_vtemp, Real* parcel_qsat, Real* lcl_pmid, Real* lcl_temperature, Int* lcl_klev);

void compute_cape_from_parcel_bridge_f(Int pcols, Int ncol, Int pver, Int pverp, Int num_cin, Int num_msg, Real* temperature, Real* tv, Real* sp_humidity, Real* pint, Int* msemax_klev, Real* lcl_pmid, Int* lcl_klev, Real* parcel_qsat, Real* parcel_temp, Real* parcel_vtemp, Int* eql_klev, Real* cape);

void zm_conv_mcsp_calculate_shear_bridge_f(Int pcols, Int ncol, Int pver, Real* state_pmid, Real* state_u, Real* state_v, Real* mcsp_shear);

void zm_conv_mcsp_tend_bridge_f(Int pcols, Int ncol, Int pver, Int pverp, Real ztodt, Int* jctop, Real* state_pmid, Real* state_pint, Real* state_pdel, Real* state_s, Real* state_q, Real* state_u, Real* state_v, Real* ptend_zm_s, Real* ptend_zm_q, Real* ptend_s, Real* ptend_q, Real* ptend_u, Real* ptend_v, Real* mcsp_dt_out, Real* mcsp_dq_out, Real* mcsp_du_out, Real* mcsp_dv_out, Real* mcsp_freq, Real* mcsp_shear, Real* zm_depth);

void zm_conv_main_bridge_f(Int pcols, Int ncol, Int pver, Int pverp, bool is_first_step, Real time_step, Real* t_mid, Real* q_mid_in, Real* omega, Real* p_mid_in, Real* p_int_in, Real* p_del_in, Real* geos, Real* z_mid_in, Real* z_int_in, Real* pbl_hgt, Real* tpert, Real* landfrac, Real* t_star, Real* q_star, Int* lengath, Int* gather_index, Int* msemax_klev_g, Int* jctop, Int* jcbot, Int* jt, Real* prec, Real* heat, Real* qtnd, Real* cape, Real* dcape, Real* mcon, Real* pflx, Real* zdu, Real* mflx_up, Real* entr_up, Real* detr_up, Real* mflx_dn, Real* entr_dn, Real* p_del, Real* dsubcld, Real* ql, Real* rliq, Real* rprd, Real* dlf);

void zm_conv_evap_bridge_f(Int pcols, Int ncol, Int pver, Int pverp, Real time_step, Real* p_mid, Real* p_del, Real* t_mid, Real* q_mid, Real* prdprec, Real* cldfrc, Real* tend_s, Real* tend_q, Real* tend_s_snwprd, Real* tend_s_snwevmlt, Real* prec, Real* snow, Real* ntprprd, Real* ntsnprd, Real* flxprec, Real* flxsnow);

void zm_calc_fractional_entrainment_bridge_f(Int pcols, Int ncol, Int pver, Int pverp, Int msg, Int* jb, Int* jt, Int* j0, Real* z_mid, Real* z_int, Real* dz, Real* h_env, Real* h_env_sat, Real* h_env_min, Real* lambda, Real* lambda_max);

void zm_downdraft_properties_bridge_f(Int pcols, Int ncol, Int pver, Int pverp, Int msg, Int* jb, Int* jt, Int* j0, Int* jd, Real* z_int, Real* dz, Real* s_mid, Real* q_mid, Real* h_env, Real* lambda, Real* lambda_max, Real* qsthat, Real* hsthat, Real* gamhat, Real* rprd, Real* mflx_up, Real* mflx_dn, Real* entr_dn, Real* s_dnd, Real* q_dnd, Real* h_dnd, Real* q_dnd_sat, Real* evp, Real* totevp);

void zm_cloud_properties_bridge_f(Int pcols, Int ncol, Int pver, Int pverp, Int msg, Int limcnv, Real* p_mid, Real* z_mid, Real* z_int, Real* t_mid, Real* s_mid, Real* s_int, Real* q_mid, Real* landfrac, Real* tpert_g, Int* jb, Int* lel, Int* jt, Int* jlcl, Int* j0, Int* jd, Real* mflx_up, Real* entr_up, Real* detr_up, Real* mflx_dn, Real* entr_dn, Real* mflx_net, Real* s_upd, Real* q_upd, Real* ql, Real* s_dnd, Real* q_dnd, Real* qst, Real* cu, Real* evp, Real* pflx, Real* rprd);

void zm_closure_bridge_f(Int pcols, Int ncol, Int pver, Int pverp, Int msg, Real cape_threshold_in, Int* lcl, Int* lel, Int* jt, Int* mx, Real* dsubcld, Real* z_mid, Real* z_int, Real* p_mid, Real* p_del, Real* t_mid, Real* s_mid, Real* q_mid, Real* qs, Real* ql, Real* s_int, Real* q_int, Real* t_pcl_lcl, Real* t_pcl, Real* q_pcl_sat, Real* s_upd, Real* q_upd, Real* mflx_net, Real* detr_up, Real* mflx_up, Real* mflx_dn, Real* q_dnd, Real* s_dnd, Real* cape, Real* cld_base_mass_flux);

void zm_calc_output_tend_bridge_f(Int pcols, Int ncol, Int pver, Int pverp, Int msg, Int* jt, Int* mx, Real* dsubcld, Real* p_del, Real* s_int, Real* q_int, Real* s_upd, Real* q_upd, Real* mflx_up, Real* detr_up, Real* mflx_dn, Real* s_dnd, Real* q_dnd, Real* ql, Real* evp, Real* cu, Real* dsdt, Real* dqdt, Real* dl);
} // extern "C" : end _f decls

// Inits and finalizes are not intended to be called outside this comp unit
namespace {

void zm_common_init_f()
{
  zm_common_init_bridge_f();
}

void zm_common_finalize_f()
{
  zm_common_finalize_bridge_f();
}


// Wrapper around zm_common_init for cxx
void zm_common_init()
{
  ZMF::zm_common_init();
}

void zm_finalize_cxx()
{
  ZMF::zm_finalize();
}

}

void ientropy_f(IentropyData& d)
{
  d.transition<ekat::TransposeDirection::c2f>();
  zm_common_init_f();
  ientropy_bridge_f(d.s, d.p, d.qt, &d.t, &d.qst, d.tfg);
  zm_common_finalize_f();
  d.transition<ekat::TransposeDirection::f2c>();
}

void ientropy(IentropyData& d)
{
  zm_common_init();

  // create device views and copy
  const auto policy = ekat::TeamPolicyFactory<ExeSpace>::get_default_team_policy(1, 1);

  // unpack data scalars because we do not want the lambda to capture d
  const Real p = d.p;
  const Real qt = d.qt;
  const Real s = d.s;
  const Real tfg = d.tfg;
  view0dr_d qst_d("qst_d");
  view0dr_d t_d("t_d");
  auto qst_h = Kokkos::create_mirror_view(qst_d);
  auto t_h   = Kokkos::create_mirror_view(t_d);

  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    ZMF::invert_entropy(
      team,
      s,
      p,
      qt,
      tfg,
      t_d(),
      qst_d());
  });

  // Get outputs back, start with scalars
  Kokkos::deep_copy(qst_h, qst_d);
  d.qst = qst_h();
  Kokkos::deep_copy(t_h, t_d);
  d.t = t_h();

  zm_finalize_cxx();
}

void entropy_f(EntropyData& d)
{
  d.transition<ekat::TransposeDirection::c2f>();
  zm_common_init_f();
  entropy_bridge_f(d.tk, d.p, d.qtot, &d.entropy);
  zm_common_finalize_f();
  d.transition<ekat::TransposeDirection::f2c>();
}

void entropy(EntropyData& d)
{
  zm_common_init();

  // create device views and copy
  const auto policy = ekat::TeamPolicyFactory<ExeSpace>::get_default_team_policy(1, 1);

  // unpack data scalars because we do not want the lambda to capture d
  const Real p = d.p;
  const Real qtot = d.qtot;
  const Real tk = d.tk;
  view0dr_d entropy_d("entropy_d");
  auto entropy_h = Kokkos::create_mirror_view(entropy_d);

  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType&) {
    entropy_d() = ZMF::entropy(
      tk,
      p,
      qtot);
  });

  // Get outputs back, start with scalars
  Kokkos::deep_copy(entropy_h, entropy_d);
  d.entropy = entropy_h();

  zm_finalize_cxx();
}

void zm_transport_tracer_f(ZmTransportTracerData& d)
{
  d.transition<ekat::TransposeDirection::c2f>();
  zm_common_init_f();
  zm_transport_tracer_bridge_f(d.pcols, d.pver, d.doconvtran, d.q, d.ncnst, d.mu, d.md, d.du, d.eu, d.ed, d.dp, d.jt, d.mx, d.ideep, d.il1g, d.il2g, d.fracis, d.dqdt, d.dpdry, d.dt);
  zm_common_finalize_f();
  d.transition<ekat::TransposeDirection::f2c>();
}

void zm_transport_tracer(ZmTransportTracerData& d)
{
  zm_common_init();

  // create device views and copy
  std::vector<view2dr_d> vec2dr_in(7);
  ekat::host_to_device({d.dp, d.dpdry, d.du, d.ed, d.eu, d.md, d.mu}, d.pcols, d.pver, vec2dr_in);

  std::vector<view3dr_d> vec3dr_in(3);
  ekat::host_to_device({d.dqdt, d.fracis, d.q}, d.pcols, d.pver, d.ncnst, vec3dr_in);

  std::vector<view1di_d> vec1di_in(3);
  ekat::host_to_device({d.ideep, d.jt, d.mx}, d.pcols, vec1di_in);

  std::vector<view1db_d> vec1db_in(1);
  ekat::host_to_device({d.doconvtran}, d.ncnst, vec1db_in);

  view2dr_d
    dp_d(vec2dr_in[0]),
    dpdry_d(vec2dr_in[1]),
    du_d(vec2dr_in[2]),
    ed_d(vec2dr_in[3]),
    eu_d(vec2dr_in[4]),
    md_d(vec2dr_in[5]),
    mu_d(vec2dr_in[6]);

  view3dr_d
    dqdt_d(vec3dr_in[0]),
    fracis_d(vec3dr_in[1]),
    q_d(vec3dr_in[2]);

  view1di_d
    ideep_d(vec1di_in[0]),
    jt_d(vec1di_in[1]),
    mx_d(vec1di_in[2]);

  view1db_d
    doconvtran_d(vec1db_in[0]);

  const auto policy = ekat::TeamPolicyFactory<ExeSpace>::get_default_team_policy(d.pcols, d.pver);

  WSM wsm(d.pver * d.ncnst, 10, policy);
  ZMF::ZmRuntimeOpt init_cp = ZMF::s_common_init;

  // unpack data scalars because we do not want the lambda to capture d
  const Real dt = d.dt;
  const Int ncnst = d.ncnst;
  const Int pver = d.pver;

  // Find min top and bottom
  assert(d.il1g == 0 && d.il2g == d.pcols-1);
  Int ktm, kbm;
  Kokkos::RangePolicy<ExeSpace> rpolicy(0, d.pcols);
  Kokkos::parallel_reduce("FindMinJt", rpolicy, KOKKOS_LAMBDA(const int i, Int& update) {
    if (jt_d(i) < update) {
      update = jt_d(i);
    }
  }, Kokkos::Min<Int>(ktm));

  Kokkos::parallel_reduce("FindMinMx", rpolicy, KOKKOS_LAMBDA(const int i, Int& update) {
    if (mx_d(i) < update) {
      update = mx_d(i);
    }
  }, Kokkos::Min<Int>(kbm));

  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    // Get single-column subviews of all inputs, shouldn't need any i-indexing
    // after this.
    const auto q_c = ekat::subview(q_d, i);
    const auto mu_c = ekat::subview(mu_d, i);
    const auto md_c = ekat::subview(md_d, i);
    const auto du_c = ekat::subview(du_d, i);
    const auto eu_c = ekat::subview(eu_d, i);
    const auto ed_c = ekat::subview(ed_d, i);
    const auto dp_c = ekat::subview(dp_d, i);
    const auto fracis_c = ekat::subview(fracis_d, i);
    const auto dqdt_c = ekat::subview(dqdt_d, i);
    const auto dpdry_c = ekat::subview(dpdry_d, i);

    ZMF::zm_transport_tracer(
      team,
      wsm.get_workspace(team),
      init_cp,
      pver,
      doconvtran_d,
      q_c,
      ncnst,
      mu_c,
      md_c,
      du_c,
      eu_c,
      ed_c,
      dp_c,
      jt_d(i),
      mx_d(i),
      ktm,
      kbm,
      fracis_c,
      dpdry_c,
      dt,
      dqdt_c);
  });

  // Now get arrays
  std::vector<view3dr_d> vec3dr_out = {dqdt_d};
  ekat::device_to_host({d.dqdt}, d.pcols, d.pver, d.ncnst, vec3dr_out);

  zm_finalize_cxx();
}

void zm_transport_momentum_f(ZmTransportMomentumData& d)
{
  d.transition<ekat::TransposeDirection::c2f>();
  zm_common_init_f();
  zm_transport_momentum_bridge_f(d.pcols, d.ncol, d.pver, d.pverp, d.wind_in, d.nwind, d.mu, d.md, d.du, d.eu, d.ed, d.dp, d.jt, d.mx, d.ideep, d.il1g, d.il2g, d.wind_tend, d.pguall, d.pgdall, d.icwu, d.icwd, d.dt, d.seten);
  zm_common_finalize_f();
  d.transition<ekat::TransposeDirection::f2c>();
}

void zm_transport_momentum(ZmTransportMomentumData& d)
{
  zm_common_init();

  // create device views and copy
  std::vector<view2dr_d> vec2dr_in(7);
  ekat::host_to_device({d.dp, d.du, d.ed, d.eu, d.md, d.mu, d.seten}, d.pcols, d.pver, vec2dr_in);

  std::vector<view3dr_d> vec3dr_in(6);
  ekat::host_to_device({d.icwd, d.icwu, d.pgdall, d.pguall, d.wind_in, d.wind_tend}, d.pcols, d.pver, d.nwind, vec3dr_in);

  std::vector<view1di_d> vec1di_in(3);
  ekat::host_to_device({d.ideep, d.jt, d.mx}, d.pcols, vec1di_in);

  view2dr_d
    dp_d(vec2dr_in[0]),
    du_d(vec2dr_in[1]),
    ed_d(vec2dr_in[2]),
    eu_d(vec2dr_in[3]),
    md_d(vec2dr_in[4]),
    mu_d(vec2dr_in[5]),
    seten_d(vec2dr_in[6]);

  view3dr_d
    icwd_d(vec3dr_in[0]),
    icwu_d(vec3dr_in[1]),
    pgdall_d(vec3dr_in[2]),
    pguall_d(vec3dr_in[3]),
    wind_in_d(vec3dr_in[4]),
    wind_tend_d(vec3dr_in[5]);

  view1di_d
    ideep_d(vec1di_in[0]),
    jt_d(vec1di_in[1]),
    mx_d(vec1di_in[2]);

  const auto policy = ekat::TeamPolicyFactory<ExeSpace>::get_default_team_policy(d.pcols, d.pver);

  WSM wsm(d.pverp * d.nwind, 12, policy);

  // unpack data scalars because we do not want the lambda to capture d
  const Real dt = d.dt;
  const Int il1g = d.il1g;
  const Int il2g = d.il2g;
  const Int nwind = d.nwind;
  const Int pver = d.pver;
  const Int pverp = d.pverp;

  // Find min top and bottom
  assert(d.pverp == d.pver + 1);
  assert(d.il1g == 0 && d.il2g == d.pcols-1);
  Int ktm, kbm;
  Kokkos::RangePolicy<ExeSpace> rpolicy(0, d.pcols);
  Kokkos::parallel_reduce("FindMinJt", rpolicy, KOKKOS_LAMBDA(const int i, Int& update) {
    if (jt_d(i) < update) {
      update = jt_d(i);
    }
  }, Kokkos::Min<Int>(ktm));

  Kokkos::parallel_reduce("FindMinMx", rpolicy, KOKKOS_LAMBDA(const int i, Int& update) {
    if (mx_d(i) < update) {
      update = mx_d(i);
    }
  }, Kokkos::Min<Int>(kbm));


  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    // Get single-column subviews of all inputs, shouldn't need any i-indexing
    // after this.
    const auto wind_in_c = ekat::subview(wind_in_d, i);
    const auto mu_c = ekat::subview(mu_d, i);
    const auto md_c = ekat::subview(md_d, i);
    const auto du_c = ekat::subview(du_d, i);
    const auto eu_c = ekat::subview(eu_d, i);
    const auto ed_c = ekat::subview(ed_d, i);
    const auto dp_c = ekat::subview(dp_d, i);
    const auto wind_tend_c = ekat::subview(wind_tend_d, i);
    const auto pguall_c = ekat::subview(pguall_d, i);
    const auto pgdall_c = ekat::subview(pgdall_d, i);
    const auto icwu_c = ekat::subview(icwu_d, i);
    const auto icwd_c = ekat::subview(icwd_d, i);
    const auto seten_c = ekat::subview(seten_d, i);

    ZMF::zm_transport_momentum(
      team,
      wsm.get_workspace(team),
      pver,
      pverp,
      wind_in_c,
      nwind,
      mu_c,
      md_c,
      du_c,
      eu_c,
      ed_c,
      dp_c,
      jt_d(i),
      mx_d(i),
      ideep_d(i),
      il1g,
      il2g,
      dt,
      ktm,
      kbm,
      wind_tend_c,
      pguall_c,
      pgdall_c,
      icwu_c,
      icwd_c,
      seten_c);
  });

  // Now get arrays
  std::vector<view2dr_d> vec2dr_out = {seten_d};
  ekat::device_to_host({d.seten}, d.pcols, d.pver, vec2dr_out);

  std::vector<view3dr_d> vec3dr_out = {icwd_d, icwu_d, pgdall_d, pguall_d, wind_tend_d};
  ekat::device_to_host({d.icwd, d.icwu, d.pgdall, d.pguall, d.wind_tend}, d.pcols, d.pver, d.nwind, vec3dr_out);

  zm_finalize_cxx();
}

void compute_dilute_cape_f(ComputeDiluteCapeData& d)
{
  d.transition<ekat::TransposeDirection::c2f>();
  zm_common_init_f();
  compute_dilute_cape_bridge_f(d.pcols, d.ncol, d.pver, d.pverp, d.num_cin, d.num_msg, d.sp_humidity_in, d.temperature_in, d.zmid, d.pmid, d.pint, d.pblt, d.tpert, d.parcel_temp, d.parcel_qsat, d.msemax_klev, d.lcl_temperature, d.lcl_klev, d.eql_klev, d.cape, d.calc_msemax_klev, d.prev_msemax_klev, d.use_input_tq_mx, d.q_mx, d.t_mx);
  zm_common_finalize_f();
  d.transition<ekat::TransposeDirection::f2c>();
}

void compute_dilute_cape(ComputeDiluteCapeData& d)
{
  zm_common_init();

  // create device views and copy
  std::vector<view1dr_d> vec1dr_in(5);
  ekat::host_to_device({d.cape, d.lcl_temperature, d.q_mx, d.t_mx, d.tpert}, d.pcols, vec1dr_in);

  std::vector<view2dr_d> vec2dr_in(7);
  std::vector<int> vec2dr_in_0_sizes = {d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols};
  std::vector<int> vec2dr_in_1_sizes = {d.pver, d.pver, d.pverp, d.pver, d.pver, d.pver, d.pver};
  ekat::host_to_device({d.parcel_qsat, d.parcel_temp, d.pint, d.pmid, d.sp_humidity_in, d.temperature_in, d.zmid}, vec2dr_in_0_sizes, vec2dr_in_1_sizes, vec2dr_in);

  std::vector<view1di_d> vec1di_in(5);
  ekat::host_to_device({d.eql_klev, d.lcl_klev, d.msemax_klev, d.pblt, d.prev_msemax_klev}, d.pcols, vec1di_in);

  view1dr_d
    cape_d(vec1dr_in[0]),
    lcl_temperature_d(vec1dr_in[1]),
    q_mx_d(vec1dr_in[2]),
    t_mx_d(vec1dr_in[3]),
    tpert_d(vec1dr_in[4]);

  view2dr_d
    parcel_qsat_d(vec2dr_in[0]),
    parcel_temp_d(vec2dr_in[1]),
    pint_d(vec2dr_in[2]),
    pmid_d(vec2dr_in[3]),
    sp_humidity_in_d(vec2dr_in[4]),
    temperature_in_d(vec2dr_in[5]),
    zmid_d(vec2dr_in[6]);

  view1di_d
    eql_klev_d(vec1di_in[0]),
    lcl_klev_d(vec1di_in[1]),
    msemax_klev_d(vec1di_in[2]),
    pblt_d(vec1di_in[3]),
    prev_msemax_klev_d(vec1di_in[4]);

  const auto policy = ekat::TeamPolicyFactory<ExeSpace>::get_default_team_policy(d.pcols, d.pver);

  WSM wsm(d.pver, 11, policy);
  ZMF::ZmRuntimeOpt init_cp = ZMF::s_common_init;

  // unpack data scalars because we do not want the lambda to capture d
  const Int num_cin = d.num_cin;
  const Int num_msg = d.num_msg;
  const Int pver = d.pver;
  const Int pverp = d.pverp;
  const bool calc_msemax_klev = d.calc_msemax_klev;
  const bool use_input_tq_mx = d.use_input_tq_mx;

  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    // Get single-column subviews of all inputs, shouldn't need any i-indexing
    // after this.
    const auto sp_humidity_in_c = ekat::subview(sp_humidity_in_d, i);
    const auto temperature_in_c = ekat::subview(temperature_in_d, i);
    const auto zmid_c = ekat::subview(zmid_d, i);
    const auto pmid_c = ekat::subview(pmid_d, i);
    const auto pint_c = ekat::subview(pint_d, i);
    const auto parcel_temp_c = ekat::subview(parcel_temp_d, i);
    const auto parcel_qsat_c = ekat::subview(parcel_qsat_d, i);

    ZMF::compute_dilute_cape(
      team,
      wsm.get_workspace(team),
      init_cp,
      pver,
      pverp,
      num_cin,
      num_msg,
      sp_humidity_in_c,
      temperature_in_c,
      zmid_c,
      pmid_c,
      pint_c,
      pblt_d(i),
      tpert_d(i),
      calc_msemax_klev,
      prev_msemax_klev_d(i),
      use_input_tq_mx,
      parcel_qsat_c,
      msemax_klev_d(i),
      lcl_klev_d(i),
      eql_klev_d(i),
      cape_d(i),
      q_mx_d(i),
      t_mx_d(i),
      parcel_temp_c,
      lcl_temperature_d(i));
  });

  // Now get arrays
  std::vector<view1dr_d> vec1dr_out = {cape_d, lcl_temperature_d, q_mx_d, t_mx_d};
  ekat::device_to_host({d.cape, d.lcl_temperature, d.q_mx, d.t_mx}, d.pcols, vec1dr_out);

  std::vector<view2dr_d> vec2dr_out = {parcel_qsat_d, parcel_temp_d};
  ekat::device_to_host({d.parcel_qsat, d.parcel_temp}, d.pcols, d.pver, vec2dr_out);

  std::vector<view1di_d> vec1di_out = {eql_klev_d, lcl_klev_d, msemax_klev_d};
  ekat::device_to_host({d.eql_klev, d.lcl_klev, d.msemax_klev}, d.pcols, vec1di_out);

  zm_finalize_cxx();
}

void find_mse_max_f(FindMseMaxData& d)
{
  d.transition<ekat::TransposeDirection::c2f>();
  zm_common_init_f();
  find_mse_max_bridge_f(d.pcols, d.ncol, d.pver, d.num_msg, d.msemax_top_k, d.pergro_active, d.temperature, d.zmid, d.sp_humidity, d.msemax_klev, d.mse_max_val);
  zm_common_finalize_f();
  d.transition<ekat::TransposeDirection::f2c>();
}

void find_mse_max(FindMseMaxData& d)
{
  zm_common_init();

  // create device views and copy
  std::vector<view1dr_d> vec1dr_in(1);
  ekat::host_to_device({d.mse_max_val}, d.pcols, vec1dr_in);

  std::vector<view2dr_d> vec2dr_in(3);
  ekat::host_to_device({d.sp_humidity, d.temperature, d.zmid}, d.pcols, d.pver, vec2dr_in);

  std::vector<view1di_d> vec1di_in(2);
  ekat::host_to_device({d.msemax_klev, d.msemax_top_k}, d.pcols, vec1di_in);

  view1dr_d
    mse_max_val_d(vec1dr_in[0]);

  view2dr_d
    sp_humidity_d(vec2dr_in[0]),
    temperature_d(vec2dr_in[1]),
    zmid_d(vec2dr_in[2]);

  view1di_d
    msemax_klev_d(vec1di_in[0]),
    msemax_top_k_d(vec1di_in[1]);

  const auto policy = ekat::TeamPolicyFactory<ExeSpace>::get_default_team_policy(d.pcols, d.pver);

  ZMF::ZmRuntimeOpt init_cp = ZMF::s_common_init;

  // unpack data scalars because we do not want the lambda to capture d
  const Int num_msg = d.num_msg;
  const Int pver = d.pver;
  const bool pergro_active = d.pergro_active;

  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    // Get single-column subviews of all inputs, shouldn't need any i-indexing
    // after this.
    const auto temperature_c = ekat::subview(temperature_d, i);
    const auto zmid_c = ekat::subview(zmid_d, i);
    const auto sp_humidity_c = ekat::subview(sp_humidity_d, i);

    ZMF::find_mse_max(
      team,
      init_cp,
      pver,
      num_msg,
      msemax_top_k_d(i),
      pergro_active,
      temperature_c,
      zmid_c,
      sp_humidity_c,
      msemax_klev_d(i),
      mse_max_val_d(i));
  });

  // Now get arrays
  std::vector<view1dr_d> vec1dr_out = {mse_max_val_d};
  ekat::device_to_host({d.mse_max_val}, d.pcols, vec1dr_out);

  std::vector<view1di_d> vec1di_out = {msemax_klev_d};
  ekat::device_to_host({d.msemax_klev}, d.pcols, vec1di_out);

  zm_finalize_cxx();
}

void compute_dilute_parcel_f(ComputeDiluteParcelData& d)
{
  d.transition<ekat::TransposeDirection::c2f>();
  zm_common_init_f();
  compute_dilute_parcel_bridge_f(d.pcols, d.ncol, d.pver, d.num_msg, d.klaunch, d.pmid, d.temperature, d.sp_humidity, d.tpert, d.pblt, d.parcel_temp, d.parcel_vtemp, d.parcel_qsat, d.lcl_pmid, d.lcl_temperature, d.lcl_klev);
  zm_common_finalize_f();
  d.transition<ekat::TransposeDirection::f2c>();
}

void compute_dilute_parcel(ComputeDiluteParcelData& d)
{
  zm_common_init();

  // create device views and copy
  std::vector<view1dr_d> vec1dr_in(3);
  ekat::host_to_device({d.lcl_pmid, d.lcl_temperature, d.tpert}, d.pcols, vec1dr_in);

  std::vector<view2dr_d> vec2dr_in(6);
  ekat::host_to_device({d.parcel_qsat, d.parcel_temp, d.parcel_vtemp, d.pmid, d.sp_humidity, d.temperature}, d.pcols, d.pver, vec2dr_in);

  std::vector<view1di_d> vec1di_in(3);
  ekat::host_to_device({d.klaunch, d.lcl_klev, d.pblt}, d.pcols, vec1di_in);

  view1dr_d
    lcl_pmid_d(vec1dr_in[0]),
    lcl_temperature_d(vec1dr_in[1]),
    tpert_d(vec1dr_in[2]);

  view2dr_d
    parcel_qsat_d(vec2dr_in[0]),
    parcel_temp_d(vec2dr_in[1]),
    parcel_vtemp_d(vec2dr_in[2]),
    pmid_d(vec2dr_in[3]),
    sp_humidity_d(vec2dr_in[4]),
    temperature_d(vec2dr_in[5]);

  view1di_d
    klaunch_d(vec1di_in[0]),
    lcl_klev_d(vec1di_in[1]),
    pblt_d(vec1di_in[2]);

  const auto policy = ekat::TeamPolicyFactory<ExeSpace>::get_default_team_policy(d.pcols, d.pver);

  WSM wsm(d.pver, 7, policy);
  ZMF::ZmRuntimeOpt init_cp = ZMF::s_common_init;

  // unpack data scalars because we do not want the lambda to capture d
  const Int num_msg = d.num_msg;
  const Int pver = d.pver;

  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    // Get single-column subviews of all inputs, shouldn't need any i-indexing
    // after this.
    const auto pmid_c = ekat::subview(pmid_d, i);
    const auto temperature_c = ekat::subview(temperature_d, i);
    const auto sp_humidity_c = ekat::subview(sp_humidity_d, i);
    const auto parcel_temp_c = ekat::subview(parcel_temp_d, i);
    const auto parcel_vtemp_c = ekat::subview(parcel_vtemp_d, i);
    const auto parcel_qsat_c = ekat::subview(parcel_qsat_d, i);

    ZMF::compute_dilute_parcel(
      team,
      wsm.get_workspace(team),
      init_cp,
      pver,
      num_msg,
      klaunch_d(i),
      pmid_c,
      temperature_c,
      sp_humidity_c,
      tpert_d(i),
      pblt_d(i),
      parcel_temp_c,
      parcel_vtemp_c,
      parcel_qsat_c,
      lcl_pmid_d(i),
      lcl_temperature_d(i),
      lcl_klev_d(i));
  });

  // Now get arrays
  std::vector<view1dr_d> vec1dr_out = {lcl_pmid_d, lcl_temperature_d};
  ekat::device_to_host({d.lcl_pmid, d.lcl_temperature}, d.pcols, vec1dr_out);

  std::vector<view2dr_d> vec2dr_out = {parcel_qsat_d, parcel_temp_d, parcel_vtemp_d};
  ekat::device_to_host({d.parcel_qsat, d.parcel_temp, d.parcel_vtemp}, d.pcols, d.pver, vec2dr_out);

  std::vector<view1di_d> vec1di_out = {lcl_klev_d};
  ekat::device_to_host({d.lcl_klev}, d.pcols, vec1di_out);

  zm_finalize_cxx();
}

void compute_cape_from_parcel_f(ComputeCapeFromParcelData& d)
{
  d.transition<ekat::TransposeDirection::c2f>();
  zm_common_init_f();
  compute_cape_from_parcel_bridge_f(d.pcols, d.ncol, d.pver, d.pverp, d.num_cin, d.num_msg, d.temperature, d.tv, d.sp_humidity, d.pint, d.msemax_klev, d.lcl_pmid, d.lcl_klev, d.parcel_qsat, d.parcel_temp, d.parcel_vtemp, d.eql_klev, d.cape);
  zm_common_finalize_f();
  d.transition<ekat::TransposeDirection::f2c>();
}

void compute_cape_from_parcel(ComputeCapeFromParcelData& d)
{
  zm_common_init();

  // create device views and copy
  std::vector<view1dr_d> vec1dr_in(2);
  ekat::host_to_device({d.cape, d.lcl_pmid}, d.pcols, vec1dr_in);

  std::vector<view2dr_d> vec2dr_in(7);
  std::vector<int> vec2dr_in_0_sizes = {d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols};
  std::vector<int> vec2dr_in_1_sizes = {d.pver, d.pver, d.pver, d.pverp, d.pver, d.pver, d.pver};
  ekat::host_to_device({d.parcel_qsat, d.parcel_temp, d.parcel_vtemp, d.pint, d.sp_humidity, d.temperature, d.tv}, vec2dr_in_0_sizes, vec2dr_in_1_sizes, vec2dr_in);

  std::vector<view1di_d> vec1di_in(3);
  ekat::host_to_device({d.eql_klev, d.lcl_klev, d.msemax_klev}, d.pcols, vec1di_in);

  view1dr_d
    cape_d(vec1dr_in[0]),
    lcl_pmid_d(vec1dr_in[1]);

  view2dr_d
    parcel_qsat_d(vec2dr_in[0]),
    parcel_temp_d(vec2dr_in[1]),
    parcel_vtemp_d(vec2dr_in[2]),
    pint_d(vec2dr_in[3]),
    sp_humidity_d(vec2dr_in[4]),
    temperature_d(vec2dr_in[5]),
    tv_d(vec2dr_in[6]);

  view1di_d
    eql_klev_d(vec1di_in[0]),
    lcl_klev_d(vec1di_in[1]),
    msemax_klev_d(vec1di_in[2]);

  const auto policy = ekat::TeamPolicyFactory<ExeSpace>::get_default_team_policy(d.pcols, d.pver);

  WSM wsm(d.pver, 3, policy);
  ZMF::ZmRuntimeOpt init_cp = ZMF::s_common_init;

  // unpack data scalars because we do not want the lambda to capture d
  const Int num_cin = d.num_cin;
  const Int num_msg = d.num_msg;
  const Int pver = d.pver;
  const Int pverp = d.pverp;

  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    // Get single-column subviews of all inputs, shouldn't need any i-indexing
    // after this.
    const auto temperature_c = ekat::subview(temperature_d, i);
    const auto tv_c = ekat::subview(tv_d, i);
    const auto sp_humidity_c = ekat::subview(sp_humidity_d, i);
    const auto pint_c = ekat::subview(pint_d, i);
    const auto parcel_qsat_c = ekat::subview(parcel_qsat_d, i);
    const auto parcel_temp_c = ekat::subview(parcel_temp_d, i);
    const auto parcel_vtemp_c = ekat::subview(parcel_vtemp_d, i);

    ZMF::compute_cape_from_parcel(
      team,
      wsm.get_workspace(team),
      init_cp,
      pver,
      pverp,
      num_cin,
      num_msg,
      temperature_c,
      tv_c,
      sp_humidity_c,
      pint_c,
      msemax_klev_d(i),
      lcl_pmid_d(i),
      lcl_klev_d(i),
      parcel_qsat_c,
      parcel_temp_c,
      parcel_vtemp_c,
      eql_klev_d(i),
      cape_d(i));
  });

  // Now get arrays
  std::vector<view1dr_d> vec1dr_out = {cape_d};
  ekat::device_to_host({d.cape}, d.pcols, vec1dr_out);

  std::vector<view2dr_d> vec2dr_out = {parcel_qsat_d, parcel_temp_d, parcel_vtemp_d};
  ekat::device_to_host({d.parcel_qsat, d.parcel_temp, d.parcel_vtemp}, d.pcols, d.pver, vec2dr_out);

  std::vector<view1di_d> vec1di_out = {eql_klev_d};
  ekat::device_to_host({d.eql_klev}, d.pcols, vec1di_out);

  zm_finalize_cxx();
}

void zm_conv_mcsp_calculate_shear_f(ZmConvMcspCalculateShearData& d)
{
  d.transition<ekat::TransposeDirection::c2f>();
  zm_common_init_f();
  zm_conv_mcsp_calculate_shear_bridge_f(d.pcols, d.ncol, d.pver, d.state_pmid, d.state_u, d.state_v, d.mcsp_shear);
  zm_common_finalize_f();
  d.transition<ekat::TransposeDirection::f2c>();
}

void zm_conv_mcsp_calculate_shear(ZmConvMcspCalculateShearData& d)
{
  zm_common_init();

  // create device views and copy
  std::vector<view1dr_d> vec1dr_in(1);
  ekat::host_to_device({d.mcsp_shear}, d.pcols, vec1dr_in);

  std::vector<view2dr_d> vec2dr_in(3);
  ekat::host_to_device({d.state_pmid, d.state_u, d.state_v}, d.pcols, d.pver, vec2dr_in);

  view1dr_d
    mcsp_shear_d(vec1dr_in[0]);

  view2dr_d
    state_pmid_d(vec2dr_in[0]),
    state_u_d(vec2dr_in[1]),
    state_v_d(vec2dr_in[2]);

  const auto policy = ekat::TeamPolicyFactory<ExeSpace>::get_default_team_policy(d.pcols, d.pver);

  // unpack data scalars because we do not want the lambda to capture d
  const Int pver = d.pver;

  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    // Get single-column subviews of all inputs, shouldn't need any i-indexing
    // after this.
    const auto state_pmid_c = ekat::subview(state_pmid_d, i);
    const auto state_u_c = ekat::subview(state_u_d, i);

    ZMF::zm_conv_mcsp_calculate_shear(
      team,
      pver,
      state_pmid_c,
      state_u_c,
      mcsp_shear_d(i));
  });

  // Now get arrays
  std::vector<view1dr_d> vec1dr_out = {mcsp_shear_d};
  ekat::device_to_host({d.mcsp_shear}, d.pcols, vec1dr_out);

  zm_finalize_cxx();
}

void zm_conv_mcsp_tend_f(ZmConvMcspTendData& d)
{
  d.transition<ekat::TransposeDirection::c2f>();
  zm_common_init_f();
  zm_conv_mcsp_tend_bridge_f(d.pcols, d.ncol, d.pver, d.pverp, d.ztodt, d.jctop, d.state_pmid, d.state_pint, d.state_pdel, d.state_s, d.state_q, d.state_u, d.state_v, d.ptend_zm_s, d.ptend_zm_q, d.ptend_s, d.ptend_q, d.ptend_u, d.ptend_v, d.mcsp_dt_out, d.mcsp_dq_out, d.mcsp_du_out, d.mcsp_dv_out, d.mcsp_freq, d.mcsp_shear, d.zm_depth);
  zm_common_finalize_f();
  d.transition<ekat::TransposeDirection::f2c>();
}

void zm_conv_mcsp_tend(ZmConvMcspTendData& d)
{
  zm_common_init();

  // create device views and copy
  std::vector<view1dr_d> vec1dr_in(3);
  ekat::host_to_device({d.mcsp_freq, d.mcsp_shear, d.zm_depth}, d.pcols, vec1dr_in);

  std::vector<view2dr_d> vec2dr_in(17);
  std::vector<int> vec2dr_in_0_sizes = {d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols};
  std::vector<int> vec2dr_in_1_sizes = {d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pverp, d.pver, d.pver, d.pver, d.pver, d.pver};
  ekat::host_to_device({d.mcsp_dq_out, d.mcsp_dt_out, d.mcsp_du_out, d.mcsp_dv_out, d.ptend_q, d.ptend_s, d.ptend_u, d.ptend_v, d.ptend_zm_q, d.ptend_zm_s, d.state_pdel, d.state_pint, d.state_pmid, d.state_q, d.state_s, d.state_u, d.state_v}, vec2dr_in_0_sizes, vec2dr_in_1_sizes, vec2dr_in);

  std::vector<view1di_d> vec1di_in(1);
  ekat::host_to_device({d.jctop}, d.pcols, vec1di_in);

  view1dr_d
    mcsp_freq_d(vec1dr_in[0]),
    mcsp_shear_d(vec1dr_in[1]),
    zm_depth_d(vec1dr_in[2]);

  view2dr_d
    mcsp_dq_out_d(vec2dr_in[0]),
    mcsp_dt_out_d(vec2dr_in[1]),
    mcsp_du_out_d(vec2dr_in[2]),
    mcsp_dv_out_d(vec2dr_in[3]),
    ptend_q_d(vec2dr_in[4]),
    ptend_s_d(vec2dr_in[5]),
    ptend_u_d(vec2dr_in[6]),
    ptend_v_d(vec2dr_in[7]),
    ptend_zm_q_d(vec2dr_in[8]),
    ptend_zm_s_d(vec2dr_in[9]),
    state_pdel_d(vec2dr_in[10]),
    state_pint_d(vec2dr_in[11]),
    state_pmid_d(vec2dr_in[12]),
    state_q_d(vec2dr_in[13]),
    state_s_d(vec2dr_in[14]),
    state_u_d(vec2dr_in[15]),
    state_v_d(vec2dr_in[16]);

  view1di_d
    jctop_d(vec1di_in[0]);

  const auto policy = ekat::TeamPolicyFactory<ExeSpace>::get_default_team_policy(d.pcols, d.pver);

  WSM wsm(d.pver, 4, policy);
  ZMF::ZmRuntimeOpt init_cp = ZMF::s_common_init;

  // unpack data scalars because we do not want the lambda to capture d
  const Real ztodt = d.ztodt;
  const Int pver = d.pver;
  const Int pverp = d.pverp;

  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    // Get single-column subviews of all inputs, shouldn't need any i-indexing
    // after this.
    const auto state_pmid_c = ekat::subview(state_pmid_d, i);
    const auto state_pint_c = ekat::subview(state_pint_d, i);
    const auto state_pdel_c = ekat::subview(state_pdel_d, i);
    const auto state_s_c = ekat::subview(state_s_d, i);
    const auto state_q_c = ekat::subview(state_q_d, i);
    const auto state_u_c = ekat::subview(state_u_d, i);
    const auto state_v_c = ekat::subview(state_v_d, i);
    const auto ptend_zm_s_c = ekat::subview(ptend_zm_s_d, i);
    const auto ptend_zm_q_c = ekat::subview(ptend_zm_q_d, i);
    const auto ptend_s_c = ekat::subview(ptend_s_d, i);
    const auto ptend_q_c = ekat::subview(ptend_q_d, i);
    const auto ptend_u_c = ekat::subview(ptend_u_d, i);
    const auto ptend_v_c = ekat::subview(ptend_v_d, i);
    const auto mcsp_dt_out_c = ekat::subview(mcsp_dt_out_d, i);
    const auto mcsp_dq_out_c = ekat::subview(mcsp_dq_out_d, i);
    const auto mcsp_du_out_c = ekat::subview(mcsp_du_out_d, i);
    const auto mcsp_dv_out_c = ekat::subview(mcsp_dv_out_d, i);

    ZMF::zm_conv_mcsp_tend(
      team,
      wsm.get_workspace(team),
      init_cp,
      pver,
      pverp,
      ztodt,
      jctop_d(i),
      state_pmid_c,
      state_pint_c,
      state_pdel_c,
      state_s_c,
      state_q_c,
      state_u_c,
      state_v_c,
      ptend_zm_s_c,
      ptend_zm_q_c,
      ptend_s_c,
      ptend_q_c,
      ptend_u_c,
      ptend_v_c,
      mcsp_dt_out_c,
      mcsp_dq_out_c,
      mcsp_du_out_c,
      mcsp_dv_out_c,
      mcsp_freq_d(i),
      mcsp_shear_d(i),
      zm_depth_d(i));
  });

  // Now get arrays
  std::vector<view1dr_d> vec1dr_out = {mcsp_freq_d, mcsp_shear_d, zm_depth_d};
  ekat::device_to_host({d.mcsp_freq, d.mcsp_shear, d.zm_depth}, d.pcols, vec1dr_out);

  std::vector<view2dr_d> vec2dr_out = {mcsp_dq_out_d, mcsp_dt_out_d, mcsp_du_out_d, mcsp_dv_out_d, ptend_q_d, ptend_s_d, ptend_u_d, ptend_v_d};
  ekat::device_to_host({d.mcsp_dq_out, d.mcsp_dt_out, d.mcsp_du_out, d.mcsp_dv_out, d.ptend_q, d.ptend_s, d.ptend_u, d.ptend_v}, d.pcols, d.pver, vec2dr_out);

  zm_finalize_cxx();
}

void zm_conv_main_f(ZmConvMainData& d)
{
  d.transition<ekat::TransposeDirection::c2f>();
  zm_common_init_f();
  zm_conv_main_bridge_f(d.pcols, d.ncol, d.pver, d.pverp, d.is_first_step, d.time_step, d.t_mid, d.q_mid_in, d.omega, d.p_mid_in, d.p_int_in, d.p_del_in, d.geos, d.z_mid_in, d.z_int_in, d.pbl_hgt, d.tpert, d.landfrac, d.t_star, d.q_star, &d.lengath, d.gather_index, d.msemax_klev, d.jctop, d.jcbot, d.jt, d.prec, d.heat, d.qtnd, d.cape, d.dcape, d.mcon, d.pflx, d.zdu, d.mflx_up, d.entr_up, d.detr_up, d.mflx_dn, d.entr_dn, d.p_del, d.dsubcld, d.ql, d.rliq, d.rprd, d.dlf);
  zm_common_finalize_f();
  d.transition<ekat::TransposeDirection::f2c>();
}

std::vector<bool> zm_conv_main(ZmConvMainData& d)
{
  zm_common_init();

  // Upload all 2D input/output arrays [pcols, pver] or [pcols, pverp]
  std::vector<view2dr_d> vec2dr(24);
  std::vector<int> dim0(24, d.pcols);
  std::vector<int> dim1 = {d.pver, d.pver, d.pver, d.pver, d.pver, d.pverp,
                             d.pver, d.pver, d.pver, d.pver, d.pver, d.pverp,
                             d.pver, d.pverp, d.pver, d.pver, d.pver, d.pver,
                             d.pver, d.pver, d.pver, d.pverp, d.pver, d.pver};
  ekat::host_to_device({d.detr_up, d.dlf, d.entr_dn, d.entr_up, d.heat, d.mcon,
                         d.mflx_dn, d.mflx_up, d.omega, d.p_del, d.p_del_in,
                         d.p_int_in, d.p_mid_in, d.pflx, d.q_mid_in, d.q_star,
                         d.ql, d.qtnd, d.rprd, d.t_mid, d.t_star,
                         d.z_int_in, d.z_mid_in, d.zdu},
                        dim0, dim1, vec2dr);
  view2dr_d
    detr_up_d  (vec2dr[0]),  dlf_d      (vec2dr[1]),
    entr_dn_d  (vec2dr[2]),  entr_up_d  (vec2dr[3]),
    heat_d     (vec2dr[4]),  mcon_d     (vec2dr[5]),
    mflx_dn_d  (vec2dr[6]),  mflx_up_d  (vec2dr[7]),
    omega_d    (vec2dr[8]),  p_del_d    (vec2dr[9]),
    p_del_in_d (vec2dr[10]), p_int_in_d (vec2dr[11]),
    p_mid_in_d (vec2dr[12]), pflx_d     (vec2dr[13]),
    q_mid_in_d (vec2dr[14]), q_star_d   (vec2dr[15]),
    ql_d       (vec2dr[16]), qtnd_d     (vec2dr[17]),
    rprd_d     (vec2dr[18]), t_mid_d    (vec2dr[19]),
    t_star_d   (vec2dr[20]), z_int_in_d (vec2dr[21]),
    z_mid_in_d (vec2dr[22]), zdu_d      (vec2dr[23]);

  // Upload 1D scalar arrays [pcols]
  std::vector<view1dr_d> vec1dr(9);
  ekat::host_to_device({d.cape, d.dcape, d.dsubcld, d.geos, d.landfrac,
                         d.pbl_hgt, d.prec, d.rliq, d.tpert},
                        d.pcols, vec1dr);
  view1dr_d cape_d(vec1dr[0]), dcape_d(vec1dr[1]), dsubcld_d(vec1dr[2]),
            geos_d(vec1dr[3]), landfrac_d(vec1dr[4]), pbl_hgt_d(vec1dr[5]),
            prec_d(vec1dr[6]), rliq_d(vec1dr[7]), tpert_d(vec1dr[8]);

  std::vector<view1di_d> vec1di(4);
  ekat::host_to_device({d.jcbot, d.jctop, d.jt, d.msemax_klev}, d.pcols, vec1di);
  view1di_d jcbot_d(vec1di[0]), jctop_d(vec1di[1]), jt_d(vec1di[2]),
            msemax_klev_d(vec1di[3]);

  ZMF::ZmRuntimeOpt init_cp = ZMF::s_common_init;

  auto active = ZMF::zm_conv_main(
    init_cp,
    d.ncol, d.pver, d.pverp, d.is_first_step, d.time_step,
    t_mid_d, q_mid_in_d, omega_d, p_mid_in_d, p_int_in_d, p_del_in_d,
    geos_d, z_mid_in_d, z_int_in_d, pbl_hgt_d, tpert_d, landfrac_d,
    t_star_d, q_star_d,
    msemax_klev_d, jctop_d, jcbot_d, jt_d, prec_d,
    heat_d, qtnd_d, cape_d, dcape_d,
    mcon_d, pflx_d, zdu_d, mflx_up_d, entr_up_d, detr_up_d,
    mflx_dn_d, entr_dn_d, p_del_d, dsubcld_d, ql_d, rliq_d, rprd_d, dlf_d);

  // Determine active columns for gather_index (matches is_conv_active logic)
  view1di_d gather_index_d("gather_index", d.pcols);
  Kokkos::deep_copy(gather_index_d, 0);
  Int num_active = 0;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<ExeSpace>(0, d.ncol),
    KOKKOS_LAMBDA(const Int i, Int& cnt) {
      const bool col_active = active(i);
      if (col_active) {
        gather_index_d(i) = i;
        cnt++;
      }
    }, num_active);
  d.lengath = num_active;

  // Copy results back to host
  auto active_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), active);
  std::vector<bool> active_v(active_h.data(), active_h.data() + active_h.size());

  std::vector<view1dr_d> vec1dr_out = {cape_d, dcape_d, dsubcld_d, prec_d, rliq_d};
  ekat::device_to_host({d.cape, d.dcape, d.dsubcld, d.prec, d.rliq}, d.pcols, vec1dr_out);

  std::vector<view2dr_d> vec2dr_out = {detr_up_d, dlf_d, entr_dn_d, entr_up_d, heat_d,
                                        mcon_d, mflx_dn_d, mflx_up_d, p_del_d, pflx_d,
                                        ql_d, qtnd_d, rprd_d, zdu_d};
  std::vector<int> out0(14, d.pcols);
  std::vector<int> out1 = {d.pver, d.pver, d.pver, d.pver, d.pver, d.pverp,
                             d.pver, d.pver, d.pver, d.pverp, d.pver, d.pver,
                             d.pver, d.pver};
  ekat::device_to_host({d.detr_up, d.dlf, d.entr_dn, d.entr_up, d.heat, d.mcon,
                         d.mflx_dn, d.mflx_up, d.p_del, d.pflx, d.ql, d.qtnd, d.rprd, d.zdu},
                        out0, out1, vec2dr_out);

  std::vector<view1di_d> vec1di_out = {gather_index_d, jcbot_d, jctop_d, jt_d, msemax_klev_d};
  ekat::device_to_host({d.gather_index, d.jcbot, d.jctop, d.jt, d.msemax_klev},
                        d.pcols, vec1di_out);

  zm_finalize_cxx();

  return active_v;
}

void zm_conv_evap_f(ZmConvEvapData& d)
{
  d.transition<ekat::TransposeDirection::c2f>();
  zm_common_init_f();
  zm_conv_evap_bridge_f(d.pcols, d.ncol, d.pver, d.pverp, d.time_step, d.p_mid, d.p_del, d.t_mid, d.q_mid, d.prdprec, d.cldfrc, d.tend_s, d.tend_q, d.tend_s_snwprd, d.tend_s_snwevmlt, d.prec, d.snow, d.ntprprd, d.ntsnprd, d.flxprec, d.flxsnow);
  zm_common_finalize_f();
  d.transition<ekat::TransposeDirection::f2c>();
}

void zm_conv_evap(ZmConvEvapData& d)
{
  zm_common_init();

  // create device views and copy
  std::vector<view1dr_d> vec1dr_in(2);
  ekat::host_to_device({d.prec, d.snow}, d.pcols, vec1dr_in);

  std::vector<view2dr_d> vec2dr_in(14);
  std::vector<int> vec2dr_in_0_sizes = {d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols};
  std::vector<int> vec2dr_in_1_sizes = {d.pver, d.pverp, d.pverp, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver};
  ekat::host_to_device({d.cldfrc, d.flxprec, d.flxsnow, d.ntprprd, d.ntsnprd, d.p_del, d.p_mid, d.prdprec, d.q_mid, d.t_mid, d.tend_q, d.tend_s, d.tend_s_snwevmlt, d.tend_s_snwprd}, vec2dr_in_0_sizes, vec2dr_in_1_sizes, vec2dr_in);

  view1dr_d
    prec_d(vec1dr_in[0]),
    snow_d(vec1dr_in[1]);

  view2dr_d
    cldfrc_d(vec2dr_in[0]),
    flxprec_d(vec2dr_in[1]),
    flxsnow_d(vec2dr_in[2]),
    ntprprd_d(vec2dr_in[3]),
    ntsnprd_d(vec2dr_in[4]),
    p_del_d(vec2dr_in[5]),
    p_mid_d(vec2dr_in[6]),
    prdprec_d(vec2dr_in[7]),
    q_mid_d(vec2dr_in[8]),
    t_mid_d(vec2dr_in[9]),
    tend_q_d(vec2dr_in[10]),
    tend_s_d(vec2dr_in[11]),
    tend_s_snwevmlt_d(vec2dr_in[12]),
    tend_s_snwprd_d(vec2dr_in[13]);

  const auto policy = ekat::TeamPolicyFactory<ExeSpace>::get_default_team_policy(d.pcols, d.pver);

  // unpack data scalars because we do not want the lambda to capture d
  const Real time_step = d.time_step;
  const Int pver = d.pver;
  const Int pverp = d.pverp;

  ZMF::ZmRuntimeOpt init_cp = ZMF::s_common_init;

  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    // Get single-column subviews of all inputs, shouldn't need any i-indexing
    // after this.
    const auto p_mid_c = ekat::subview(p_mid_d, i);
    const auto p_del_c = ekat::subview(p_del_d, i);
    const auto t_mid_c = ekat::subview(t_mid_d, i);
    const auto q_mid_c = ekat::subview(q_mid_d, i);
    const auto prdprec_c = ekat::subview(prdprec_d, i);
    const auto cldfrc_c = ekat::subview(cldfrc_d, i);
    const auto tend_s_c = ekat::subview(tend_s_d, i);
    const auto tend_q_c = ekat::subview(tend_q_d, i);
    const auto tend_s_snwprd_c = ekat::subview(tend_s_snwprd_d, i);
    const auto tend_s_snwevmlt_c = ekat::subview(tend_s_snwevmlt_d, i);
    const auto ntprprd_c = ekat::subview(ntprprd_d, i);
    const auto ntsnprd_c = ekat::subview(ntsnprd_d, i);
    const auto flxprec_c = ekat::subview(flxprec_d, i);
    const auto flxsnow_c = ekat::subview(flxsnow_d, i);

    ZMF::zm_conv_evap(
      team,
      init_cp,
      pver,
      pverp,
      time_step,
      p_mid_c,
      p_del_c,
      t_mid_c,
      q_mid_c,
      prdprec_c,
      cldfrc_c,
      tend_s_c,
      tend_q_c,
      tend_s_snwprd_c,
      tend_s_snwevmlt_c,
      prec_d(i),
      snow_d(i),
      ntprprd_c,
      ntsnprd_c,
      flxprec_c,
      flxsnow_c);
  });

  // Now get arrays
  std::vector<view1dr_d> vec1dr_out = {prec_d, snow_d};
  ekat::device_to_host({d.prec, d.snow}, d.pcols, vec1dr_out);

  std::vector<view2dr_d> vec2dr_out = {flxprec_d, flxsnow_d, ntprprd_d, ntsnprd_d, tend_q_d, tend_s_d, tend_s_snwevmlt_d, tend_s_snwprd_d};
  std::vector<int> vec2dr_out_0_sizes = {d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols};
  std::vector<int> vec2dr_out_1_sizes = {d.pverp, d.pverp, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver};
  ekat::device_to_host({d.flxprec, d.flxsnow, d.ntprprd, d.ntsnprd, d.tend_q, d.tend_s, d.tend_s_snwevmlt, d.tend_s_snwprd}, vec2dr_out_0_sizes, vec2dr_out_1_sizes, vec2dr_out);

  zm_finalize_cxx();
}

void zm_calc_fractional_entrainment_f(ZmCalcFractionalEntrainmentData& d)
{
  d.transition<ekat::TransposeDirection::c2f>();
  zm_common_init_f();
  zm_calc_fractional_entrainment_bridge_f(d.pcols, d.ncol, d.pver, d.pverp, d.msg, d.jb, d.jt, d.j0, d.z_mid, d.z_int, d.dz, d.h_env, d.h_env_sat, d.h_env_min, d.lambda, d.lambda_max);
  zm_common_finalize_f();
  d.transition<ekat::TransposeDirection::f2c>();
}

void zm_calc_fractional_entrainment(ZmCalcFractionalEntrainmentData& d)
{
  zm_common_init();

  // create device views and copy
  std::vector<view1dr_d> vec1dr_in(2);
  ekat::host_to_device({d.h_env_min, d.lambda_max}, d.pcols, vec1dr_in);

  std::vector<view2dr_d> vec2dr_in(6);
  std::vector<int> vec2dr_in_0_sizes = {d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols};
  std::vector<int> vec2dr_in_1_sizes = {d.pver, d.pver, d.pver, d.pver, d.pverp, d.pver};
  ekat::host_to_device({d.dz, d.h_env, d.h_env_sat, d.lambda, d.z_int, d.z_mid}, vec2dr_in_0_sizes, vec2dr_in_1_sizes, vec2dr_in);

  std::vector<view1di_d> vec1di_in(3);
  ekat::host_to_device({d.j0, d.jb, d.jt}, d.pcols, vec1di_in);

  view1dr_d
    h_env_min_d(vec1dr_in[0]),
    lambda_max_d(vec1dr_in[1]);

  view2dr_d
    dz_d(vec2dr_in[0]),
    h_env_d(vec2dr_in[1]),
    h_env_sat_d(vec2dr_in[2]),
    lambda_d(vec2dr_in[3]),
    z_int_d(vec2dr_in[4]),
    z_mid_d(vec2dr_in[5]);

  view1di_d
    j0_d(vec1di_in[0]),
    jb_d(vec1di_in[1]),
    jt_d(vec1di_in[2]);

  const auto policy = ekat::TeamPolicyFactory<ExeSpace>::get_default_team_policy(d.pcols, d.pver);

  WSM wsm(d.pver, 5, policy);

  // unpack data scalars because we do not want the lambda to capture d
  const Int msg = d.msg;
  const Int pver = d.pver;
  const Int pverp = d.pverp;

  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    // Get single-column subviews of all inputs, shouldn't need any i-indexing
    // after this.
    const auto z_mid_c = ekat::subview(z_mid_d, i);
    const auto z_int_c = ekat::subview(z_int_d, i);
    const auto dz_c = ekat::subview(dz_d, i);
    const auto h_env_c = ekat::subview(h_env_d, i);
    const auto h_env_sat_c = ekat::subview(h_env_sat_d, i);
    const auto lambda_c = ekat::subview(lambda_d, i);

    ZMF::zm_calc_fractional_entrainment(
      team,
      wsm.get_workspace(team),
      pver,
      pverp,
      msg,
      jb_d(i),
      jt_d(i),
      j0_d(i),
      z_mid_c,
      z_int_c,
      dz_c,
      h_env_c,
      h_env_sat_c,
      h_env_min_d(i),
      lambda_c,
      lambda_max_d(i));
  });

  // Now get arrays
  std::vector<view1dr_d> vec1dr_out = {h_env_min_d, lambda_max_d};
  ekat::device_to_host({d.h_env_min, d.lambda_max}, d.pcols, vec1dr_out);

  std::vector<view2dr_d> vec2dr_out = {lambda_d};
  ekat::device_to_host({d.lambda}, d.pcols, d.pver, vec2dr_out);

  std::vector<view1di_d> vec1di_out = {j0_d};
  ekat::device_to_host({d.j0}, d.pcols, vec1di_out);

  zm_finalize_cxx();
}

void zm_downdraft_properties_f(ZmDowndraftPropertiesData& d)
{
  d.transition<ekat::TransposeDirection::c2f>();
  zm_common_init_f();
  zm_downdraft_properties_bridge_f(d.pcols, d.ncol, d.pver, d.pverp, d.msg, d.jb, d.jt, d.j0, d.jd, d.z_int, d.dz, d.s_mid, d.q_mid, d.h_env, d.lambda, d.lambda_max, d.qsthat, d.hsthat, d.gamhat, d.rprd, d.mflx_up, d.mflx_dn, d.entr_dn, d.s_dnd, d.q_dnd, d.h_dnd, d.q_dnd_sat, d.evp, d.totevp);
  zm_common_finalize_f();
  d.transition<ekat::TransposeDirection::f2c>();
}

void zm_downdraft_properties(ZmDowndraftPropertiesData& d)
{
  zm_common_init();

  // create device views and copy
  std::vector<view1dr_d> vec1dr_in(2);
  ekat::host_to_device({d.lambda_max, d.totevp}, d.pcols, vec1dr_in);

  std::vector<view2dr_d> vec2dr_in(18);
  std::vector<int> vec2dr_in_0_sizes = {d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols};
  std::vector<int> vec2dr_in_1_sizes = {d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pverp};
  ekat::host_to_device({d.dz, d.entr_dn, d.evp, d.gamhat, d.h_dnd, d.h_env, d.hsthat, d.lambda, d.mflx_dn, d.mflx_up, d.q_dnd, d.q_dnd_sat, d.q_mid, d.qsthat, d.rprd, d.s_dnd, d.s_mid, d.z_int}, vec2dr_in_0_sizes, vec2dr_in_1_sizes, vec2dr_in);

  std::vector<view1di_d> vec1di_in(4);
  ekat::host_to_device({d.j0, d.jb, d.jd, d.jt}, d.pcols, vec1di_in);

  view1dr_d
    lambda_max_d(vec1dr_in[0]),
    totevp_d(vec1dr_in[1]);

  view2dr_d
    dz_d(vec2dr_in[0]),
    entr_dn_d(vec2dr_in[1]),
    evp_d(vec2dr_in[2]),
    gamhat_d(vec2dr_in[3]),
    h_dnd_d(vec2dr_in[4]),
    h_env_d(vec2dr_in[5]),
    hsthat_d(vec2dr_in[6]),
    lambda_d(vec2dr_in[7]),
    mflx_dn_d(vec2dr_in[8]),
    mflx_up_d(vec2dr_in[9]),
    q_dnd_d(vec2dr_in[10]),
    q_dnd_sat_d(vec2dr_in[11]),
    q_mid_d(vec2dr_in[12]),
    qsthat_d(vec2dr_in[13]),
    rprd_d(vec2dr_in[14]),
    s_dnd_d(vec2dr_in[15]),
    s_mid_d(vec2dr_in[16]),
    z_int_d(vec2dr_in[17]);

  view1di_d
    j0_d(vec1di_in[0]),
    jb_d(vec1di_in[1]),
    jd_d(vec1di_in[2]),
    jt_d(vec1di_in[3]);

  const auto policy = ekat::TeamPolicyFactory<ExeSpace>::get_default_team_policy(d.pcols, d.pver);
  ZMF::ZmRuntimeOpt init_cp = ZMF::s_common_init;

  // unpack data scalars because we do not want the lambda to capture d
  const Int msg = d.msg;
  const Int pver = d.pver;
  const Int pverp = d.pverp;

  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    // Get single-column subviews of all inputs, shouldn't need any i-indexing
    // after this.
    const auto z_int_c = ekat::subview(z_int_d, i);
    const auto dz_c = ekat::subview(dz_d, i);
    const auto s_mid_c = ekat::subview(s_mid_d, i);
    const auto q_mid_c = ekat::subview(q_mid_d, i);
    const auto h_env_c = ekat::subview(h_env_d, i);
    const auto lambda_c = ekat::subview(lambda_d, i);
    const auto qsthat_c = ekat::subview(qsthat_d, i);
    const auto hsthat_c = ekat::subview(hsthat_d, i);
    const auto gamhat_c = ekat::subview(gamhat_d, i);
    const auto rprd_c = ekat::subview(rprd_d, i);
    const auto mflx_up_c = ekat::subview(mflx_up_d, i);
    const auto mflx_dn_c = ekat::subview(mflx_dn_d, i);
    const auto entr_dn_c = ekat::subview(entr_dn_d, i);
    const auto s_dnd_c = ekat::subview(s_dnd_d, i);
    const auto q_dnd_c = ekat::subview(q_dnd_d, i);
    const auto h_dnd_c = ekat::subview(h_dnd_d, i);
    const auto q_dnd_sat_c = ekat::subview(q_dnd_sat_d, i);
    const auto evp_c = ekat::subview(evp_d, i);

    ZMF::zm_downdraft_properties(
      team,
      init_cp,
      pver,
      pverp,
      msg,
      jb_d(i),
      jt_d(i),
      j0_d(i),
      jd_d(i),
      z_int_c,
      dz_c,
      s_mid_c,
      q_mid_c,
      h_env_c,
      lambda_c,
      lambda_max_d(i),
      qsthat_c,
      hsthat_c,
      gamhat_c,
      rprd_c,
      mflx_up_c,
      mflx_dn_c,
      entr_dn_c,
      s_dnd_c,
      q_dnd_c,
      h_dnd_c,
      q_dnd_sat_c,
      evp_c,
      totevp_d(i));
  });

  // Now get arrays
  std::vector<view1dr_d> vec1dr_out = {totevp_d};
  ekat::device_to_host({d.totevp}, d.pcols, vec1dr_out);

  std::vector<view2dr_d> vec2dr_out = {entr_dn_d, evp_d, h_dnd_d, mflx_dn_d, q_dnd_d, q_dnd_sat_d, s_dnd_d};
  ekat::device_to_host({d.entr_dn, d.evp, d.h_dnd, d.mflx_dn, d.q_dnd, d.q_dnd_sat, d.s_dnd}, d.pcols, d.pver, vec2dr_out);

  std::vector<view1di_d> vec1di_out = {jd_d, jt_d};
  ekat::device_to_host({d.jd, d.jt}, d.pcols, vec1di_out);

  zm_finalize_cxx();
}

void zm_cloud_properties_f(ZmCloudPropertiesData& d)
{
  d.transition<ekat::TransposeDirection::c2f>();
  zm_common_init_f();
  zm_cloud_properties_bridge_f(d.pcols, d.ncol, d.pver, d.pverp, d.msg, d.limcnv, d.p_mid, d.z_mid, d.z_int, d.t_mid, d.s_mid, d.s_int, d.q_mid, d.landfrac, d.tpert_g, d.jb, d.lel, d.jt, d.jlcl, d.j0, d.jd, d.mflx_up, d.entr_up, d.detr_up, d.mflx_dn, d.entr_dn, d.mflx_net, d.s_upd, d.q_upd, d.ql, d.s_dnd, d.q_dnd, d.qst, d.cu, d.evp, d.pflx, d.rprd);
  zm_common_finalize_f();
  d.transition<ekat::TransposeDirection::f2c>();
}

void zm_cloud_properties(ZmCloudPropertiesData& d)
{
  zm_common_init();

  // create device views and copy
  std::vector<view1dr_d> vec1dr_in(2);
  ekat::host_to_device({d.landfrac, d.tpert_g}, d.pcols, vec1dr_in);

  std::vector<view2dr_d> vec2dr_in(23);
  std::vector<int> vec2dr_in_0_sizes = {d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols};
  std::vector<int> vec2dr_in_1_sizes = {d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pverp, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pverp, d.pver};
  ekat::host_to_device({d.cu, d.detr_up, d.entr_dn, d.entr_up, d.evp, d.mflx_dn, d.mflx_net, d.mflx_up, d.p_mid, d.pflx, d.q_dnd, d.q_mid, d.q_upd, d.ql, d.qst, d.rprd, d.s_dnd, d.s_int, d.s_mid, d.s_upd, d.t_mid, d.z_int, d.z_mid}, vec2dr_in_0_sizes, vec2dr_in_1_sizes, vec2dr_in);

  std::vector<view1di_d> vec1di_in(6);
  ekat::host_to_device({d.j0, d.jb, d.jd, d.jlcl, d.jt, d.lel}, d.pcols, vec1di_in);

  view1dr_d
    landfrac_d(vec1dr_in[0]),
    tpert_g_d(vec1dr_in[1]);

  view2dr_d
    cu_d(vec2dr_in[0]),
    detr_up_d(vec2dr_in[1]),
    entr_dn_d(vec2dr_in[2]),
    entr_up_d(vec2dr_in[3]),
    evp_d(vec2dr_in[4]),
    mflx_dn_d(vec2dr_in[5]),
    mflx_net_d(vec2dr_in[6]),
    mflx_up_d(vec2dr_in[7]),
    p_mid_d(vec2dr_in[8]),
    pflx_d(vec2dr_in[9]),
    q_dnd_d(vec2dr_in[10]),
    q_mid_d(vec2dr_in[11]),
    q_upd_d(vec2dr_in[12]),
    ql_d(vec2dr_in[13]),
    qst_d(vec2dr_in[14]),
    rprd_d(vec2dr_in[15]),
    s_dnd_d(vec2dr_in[16]),
    s_int_d(vec2dr_in[17]),
    s_mid_d(vec2dr_in[18]),
    s_upd_d(vec2dr_in[19]),
    t_mid_d(vec2dr_in[20]),
    z_int_d(vec2dr_in[21]),
    z_mid_d(vec2dr_in[22]);

  view1di_d
    j0_d(vec1di_in[0]),
    jb_d(vec1di_in[1]),
    jd_d(vec1di_in[2]),
    jlcl_d(vec1di_in[3]),
    jt_d(vec1di_in[4]),
    lel_d(vec1di_in[5]);

  const auto policy = ekat::TeamPolicyFactory<ExeSpace>::get_default_team_policy(d.pcols, d.pver);
  ZMF::ZmRuntimeOpt init_cp = ZMF::s_common_init;
  WSM wsm(d.pverp, 20, policy);

  // unpack data scalars because we do not want the lambda to capture d
  const Int limcnv = d.limcnv;
  const Int msg = d.msg;
  const Int pver = d.pver;
  const Int pverp = d.pverp;

  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    // Get single-column subviews of all inputs, shouldn't need any i-indexing
    // after this.
    const auto p_mid_c = ekat::subview(p_mid_d, i);
    const auto z_mid_c = ekat::subview(z_mid_d, i);
    const auto z_int_c = ekat::subview(z_int_d, i);
    const auto t_mid_c = ekat::subview(t_mid_d, i);
    const auto s_mid_c = ekat::subview(s_mid_d, i);
    const auto s_int_c = ekat::subview(s_int_d, i);
    const auto q_mid_c = ekat::subview(q_mid_d, i);
    const auto mflx_up_c = ekat::subview(mflx_up_d, i);
    const auto entr_up_c = ekat::subview(entr_up_d, i);
    const auto detr_up_c = ekat::subview(detr_up_d, i);
    const auto mflx_dn_c = ekat::subview(mflx_dn_d, i);
    const auto entr_dn_c = ekat::subview(entr_dn_d, i);
    const auto mflx_net_c = ekat::subview(mflx_net_d, i);
    const auto s_upd_c = ekat::subview(s_upd_d, i);
    const auto q_upd_c = ekat::subview(q_upd_d, i);
    const auto ql_c = ekat::subview(ql_d, i);
    const auto s_dnd_c = ekat::subview(s_dnd_d, i);
    const auto q_dnd_c = ekat::subview(q_dnd_d, i);
    const auto qst_c = ekat::subview(qst_d, i);
    const auto cu_c = ekat::subview(cu_d, i);
    const auto evp_c = ekat::subview(evp_d, i);
    const auto pflx_c = ekat::subview(pflx_d, i);
    const auto rprd_c = ekat::subview(rprd_d, i);

    ZMF::zm_cloud_properties(
      team,
      wsm.get_workspace(team),
      init_cp,
      pver,
      pverp,
      msg,
      limcnv,
      p_mid_c,
      z_mid_c,
      z_int_c,
      t_mid_c,
      s_mid_c,
      s_int_c,
      q_mid_c,
      landfrac_d(i),
      tpert_g_d(i),
      jb_d(i),
      lel_d(i),
      jt_d(i),
      jlcl_d(i),
      j0_d(i),
      jd_d(i),
      mflx_up_c,
      entr_up_c,
      detr_up_c,
      mflx_dn_c,
      entr_dn_c,
      mflx_net_c,
      s_upd_c,
      q_upd_c,
      ql_c,
      s_dnd_c,
      q_dnd_c,
      qst_c,
      cu_c,
      evp_c,
      pflx_c,
      rprd_c);
  });

  // Now get arrays
  std::vector<view2dr_d> vec2dr_out = {cu_d, detr_up_d, entr_dn_d, entr_up_d, evp_d, mflx_dn_d, mflx_net_d, mflx_up_d, pflx_d, q_dnd_d, q_upd_d, ql_d, qst_d, rprd_d, s_dnd_d, s_upd_d};
  std::vector<int> vec2dr_out_0_sizes = {d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols};
  std::vector<int> vec2dr_out_1_sizes = {d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pverp, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver};
  ekat::device_to_host({d.cu, d.detr_up, d.entr_dn, d.entr_up, d.evp, d.mflx_dn, d.mflx_net, d.mflx_up, d.pflx, d.q_dnd, d.q_upd, d.ql, d.qst, d.rprd, d.s_dnd, d.s_upd}, vec2dr_out_0_sizes, vec2dr_out_1_sizes, vec2dr_out);

  std::vector<view1di_d> vec1di_out = {j0_d, jd_d, jlcl_d, jt_d};
  ekat::device_to_host({d.j0, d.jd, d.jlcl, d.jt}, d.pcols, vec1di_out);

  zm_finalize_cxx();
}

void zm_closure_f(ZmClosureData& d)
{
  d.transition<ekat::TransposeDirection::c2f>();
  zm_common_init_f();
  zm_closure_bridge_f(d.pcols, d.ncol, d.pver, d.pverp, d.msg, d.cape_threshold_in, d.lcl, d.lel, d.jt, d.mx, d.dsubcld, d.z_mid, d.z_int, d.p_mid, d.p_del, d.t_mid, d.s_mid, d.q_mid, d.qs, d.ql, d.s_int, d.q_int, d.t_pcl_lcl, d.t_pcl, d.q_pcl_sat, d.s_upd, d.q_upd, d.mflx_net, d.detr_up, d.mflx_up, d.mflx_dn, d.q_dnd, d.s_dnd, d.cape, d.cld_base_mass_flux);
  zm_common_finalize_f();
  d.transition<ekat::TransposeDirection::f2c>();
}

void zm_closure(ZmClosureData& d)
{
  zm_common_init();

  // create device views and copy
  std::vector<view1dr_d> vec1dr_in(4);
  ekat::host_to_device({d.cape, d.cld_base_mass_flux, d.dsubcld, d.t_pcl_lcl}, d.pcols, vec1dr_in);

  std::vector<view2dr_d> vec2dr_in(21);
  std::vector<int> vec2dr_in_0_sizes = {d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols};
  std::vector<int> vec2dr_in_1_sizes = {d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pverp, d.pver};
  ekat::host_to_device({d.detr_up, d.mflx_dn, d.mflx_net, d.mflx_up, d.p_del, d.p_mid, d.q_dnd, d.q_int, d.q_mid, d.q_pcl_sat, d.q_upd, d.ql, d.qs, d.s_dnd, d.s_int, d.s_mid, d.s_upd, d.t_mid, d.t_pcl, d.z_int, d.z_mid}, vec2dr_in_0_sizes, vec2dr_in_1_sizes, vec2dr_in);

  std::vector<view1di_d> vec1di_in(4);
  ekat::host_to_device({d.jt, d.lcl, d.lel, d.mx}, d.pcols, vec1di_in);

  view1dr_d
    cape_d(vec1dr_in[0]),
    cld_base_mass_flux_d(vec1dr_in[1]),
    dsubcld_d(vec1dr_in[2]),
    t_pcl_lcl_d(vec1dr_in[3]);

  view2dr_d
    detr_up_d(vec2dr_in[0]),
    mflx_dn_d(vec2dr_in[1]),
    mflx_net_d(vec2dr_in[2]),
    mflx_up_d(vec2dr_in[3]),
    p_del_d(vec2dr_in[4]),
    p_mid_d(vec2dr_in[5]),
    q_dnd_d(vec2dr_in[6]),
    q_int_d(vec2dr_in[7]),
    q_mid_d(vec2dr_in[8]),
    q_pcl_sat_d(vec2dr_in[9]),
    q_upd_d(vec2dr_in[10]),
    ql_d(vec2dr_in[11]),
    qs_d(vec2dr_in[12]),
    s_dnd_d(vec2dr_in[13]),
    s_int_d(vec2dr_in[14]),
    s_mid_d(vec2dr_in[15]),
    s_upd_d(vec2dr_in[16]),
    t_mid_d(vec2dr_in[17]),
    t_pcl_d(vec2dr_in[18]),
    z_int_d(vec2dr_in[19]),
    z_mid_d(vec2dr_in[20]);

  view1di_d
    jt_d(vec1di_in[0]),
    lcl_d(vec1di_in[1]),
    lel_d(vec1di_in[2]),
    mx_d(vec1di_in[3]);

  const auto policy = ekat::TeamPolicyFactory<ExeSpace>::get_default_team_policy(d.pcols, d.pver);

  WSM wsm(d.pver, 3, policy);
  ZMF::ZmRuntimeOpt init_cp = ZMF::s_common_init;

  // unpack data scalars because we do not want the lambda to capture d
  const Real cape_threshold_in = d.cape_threshold_in;
  const Int msg = d.msg;
  const Int pver = d.pver;
  const Int pverp = d.pverp;

  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    // Get single-column subviews of all inputs, shouldn't need any i-indexing
    // after this.
    const auto z_mid_c = ekat::subview(z_mid_d, i);
    const auto z_int_c = ekat::subview(z_int_d, i);
    const auto p_mid_c = ekat::subview(p_mid_d, i);
    const auto p_del_c = ekat::subview(p_del_d, i);
    const auto t_mid_c = ekat::subview(t_mid_d, i);
    const auto s_mid_c = ekat::subview(s_mid_d, i);
    const auto q_mid_c = ekat::subview(q_mid_d, i);
    const auto qs_c = ekat::subview(qs_d, i);
    const auto ql_c = ekat::subview(ql_d, i);
    const auto s_int_c = ekat::subview(s_int_d, i);
    const auto q_int_c = ekat::subview(q_int_d, i);
    const auto t_pcl_c = ekat::subview(t_pcl_d, i);
    const auto q_pcl_sat_c = ekat::subview(q_pcl_sat_d, i);
    const auto s_upd_c = ekat::subview(s_upd_d, i);
    const auto q_upd_c = ekat::subview(q_upd_d, i);
    const auto mflx_net_c = ekat::subview(mflx_net_d, i);
    const auto detr_up_c = ekat::subview(detr_up_d, i);
    const auto mflx_up_c = ekat::subview(mflx_up_d, i);
    const auto mflx_dn_c = ekat::subview(mflx_dn_d, i);
    const auto q_dnd_c = ekat::subview(q_dnd_d, i);
    const auto s_dnd_c = ekat::subview(s_dnd_d, i);

    ZMF::zm_closure(
      team,
      wsm.get_workspace(team),
      init_cp,
      pver,
      pverp,
      msg,
      cape_threshold_in,
      lcl_d(i),
      lel_d(i),
      jt_d(i),
      mx_d(i),
      dsubcld_d(i),
      z_mid_c,
      z_int_c,
      p_mid_c,
      p_del_c,
      t_mid_c,
      s_mid_c,
      q_mid_c,
      qs_c,
      ql_c,
      s_int_c,
      q_int_c,
      t_pcl_lcl_d(i),
      t_pcl_c,
      q_pcl_sat_c,
      s_upd_c,
      q_upd_c,
      mflx_net_c,
      detr_up_c,
      mflx_up_c,
      mflx_dn_c,
      q_dnd_c,
      s_dnd_c,
      cape_d(i),
      cld_base_mass_flux_d(i));
  });

  // Now get arrays
  std::vector<view1dr_d> vec1dr_out = {cld_base_mass_flux_d};
  ekat::device_to_host({d.cld_base_mass_flux}, d.pcols, vec1dr_out);

  zm_finalize_cxx();
}

void zm_calc_output_tend_f(ZmCalcOutputTendData& d)
{
  d.transition<ekat::TransposeDirection::c2f>();
  zm_common_init_f();
  zm_calc_output_tend_bridge_f(d.pcols, d.ncol, d.pver, d.pverp, d.msg, d.jt, d.mx, d.dsubcld, d.p_del, d.s_int, d.q_int, d.s_upd, d.q_upd, d.mflx_up, d.detr_up, d.mflx_dn, d.s_dnd, d.q_dnd, d.ql, d.evp, d.cu, d.dsdt, d.dqdt, d.dl);
  zm_common_finalize_f();
  d.transition<ekat::TransposeDirection::f2c>();
}

void zm_calc_output_tend(ZmCalcOutputTendData& d)
{
  zm_common_init();

  // create device views and copy
  std::vector<view1dr_d> vec1dr_in(1);
  ekat::host_to_device({d.dsubcld}, d.pcols, vec1dr_in);

  std::vector<view2dr_d> vec2dr_in(16);
  ekat::host_to_device({d.cu, d.detr_up, d.dl, d.dqdt, d.dsdt, d.evp, d.mflx_dn, d.mflx_up, d.p_del, d.q_dnd, d.q_int, d.q_upd, d.ql, d.s_dnd, d.s_int, d.s_upd}, d.pcols, d.pver, vec2dr_in);

  std::vector<view1di_d> vec1di_in(2);
  ekat::host_to_device({d.jt, d.mx}, d.pcols, vec1di_in);

  view1dr_d
    dsubcld_d(vec1dr_in[0]);

  view2dr_d
    cu_d(vec2dr_in[0]),
    detr_up_d(vec2dr_in[1]),
    dl_d(vec2dr_in[2]),
    dqdt_d(vec2dr_in[3]),
    dsdt_d(vec2dr_in[4]),
    evp_d(vec2dr_in[5]),
    mflx_dn_d(vec2dr_in[6]),
    mflx_up_d(vec2dr_in[7]),
    p_del_d(vec2dr_in[8]),
    q_dnd_d(vec2dr_in[9]),
    q_int_d(vec2dr_in[10]),
    q_upd_d(vec2dr_in[11]),
    ql_d(vec2dr_in[12]),
    s_dnd_d(vec2dr_in[13]),
    s_int_d(vec2dr_in[14]),
    s_upd_d(vec2dr_in[15]);

  view1di_d
    jt_d(vec1di_in[0]),
    mx_d(vec1di_in[1]);

  const auto policy = ekat::TeamPolicyFactory<ExeSpace>::get_default_team_policy(d.pcols, d.pver);

  Int ktm, kbm;
  Kokkos::RangePolicy<ExeSpace> rpolicy(0, d.pcols);
  Kokkos::parallel_reduce("FindMinJt", rpolicy, KOKKOS_LAMBDA(const int i, Int& update) {
    if (jt_d(i) < update) {
      update = jt_d(i);
    }
  }, Kokkos::Min<Int>(ktm));

  Kokkos::parallel_reduce("FindMinMx", rpolicy, KOKKOS_LAMBDA(const int i, Int& update) {
    if (mx_d(i) < update) {
      update = mx_d(i);
    }
  }, Kokkos::Min<Int>(kbm));

  // unpack data scalars because we do not want the lambda to capture d
  const Int msg = d.msg;
  const Int pver = d.pver;
  const Int pverp = d.pverp;

  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    // Get single-column subviews of all inputs, shouldn't need any i-indexing
    // after this.
    const auto p_del_c = ekat::subview(p_del_d, i);
    const auto s_int_c = ekat::subview(s_int_d, i);
    const auto q_int_c = ekat::subview(q_int_d, i);
    const auto s_upd_c = ekat::subview(s_upd_d, i);
    const auto q_upd_c = ekat::subview(q_upd_d, i);
    const auto mflx_up_c = ekat::subview(mflx_up_d, i);
    const auto detr_up_c = ekat::subview(detr_up_d, i);
    const auto mflx_dn_c = ekat::subview(mflx_dn_d, i);
    const auto s_dnd_c = ekat::subview(s_dnd_d, i);
    const auto q_dnd_c = ekat::subview(q_dnd_d, i);
    const auto ql_c = ekat::subview(ql_d, i);
    const auto evp_c = ekat::subview(evp_d, i);
    const auto cu_c = ekat::subview(cu_d, i);
    const auto dsdt_c = ekat::subview(dsdt_d, i);
    const auto dqdt_c = ekat::subview(dqdt_d, i);
    const auto dl_c = ekat::subview(dl_d, i);

    ZMF::zm_calc_output_tend(
      team,
      pver,
      pverp,
      msg,
      jt_d(i),
      mx_d(i),
      ktm,
      kbm,
      dsubcld_d(i),
      p_del_c,
      s_int_c,
      q_int_c,
      s_upd_c,
      q_upd_c,
      mflx_up_c,
      detr_up_c,
      mflx_dn_c,
      s_dnd_c,
      q_dnd_c,
      ql_c,
      evp_c,
      cu_c,
      dsdt_c,
      dqdt_c,
      dl_c);
  });

  // Now get arrays
  std::vector<view2dr_d> vec2dr_out = {dl_d, dqdt_d, dsdt_d};
  ekat::device_to_host({d.dl, d.dqdt, d.dsdt}, d.pcols, d.pver, vec2dr_out);

  zm_finalize_cxx();
}
// end glue impls

} // namespace zm
} // namespace scream
