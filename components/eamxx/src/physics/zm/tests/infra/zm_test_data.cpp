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
  zm_conv_main_bridge_f(d.pcols, d.ncol, d.pver, d.pverp, d.is_first_step, d.time_step, d.t_mid, d.q_mid_in, d.omega, d.p_mid_in, d.p_int_in, d.p_del_in, d.geos, d.z_mid_in, d.z_int_in, d.pbl_hgt, d.tpert, d.landfrac, d.t_star, d.q_star, &d.lengath, d.gather_index, d.msemax_klev_g, d.jctop, d.jcbot, d.jt, d.prec, d.heat, d.qtnd, d.cape, d.dcape, d.mcon, d.pflx, d.zdu, d.mflx_up, d.entr_up, d.detr_up, d.mflx_dn, d.entr_dn, d.p_del, d.dsubcld, d.ql, d.rliq, d.rprd, d.dlf);
  zm_common_finalize_f();
  d.transition<ekat::TransposeDirection::f2c>();
}

void zm_conv_main(ZmConvMainData& d)
{
  zm_common_init(); // Might need more specific init

  // create device views and copy
  std::vector<view1dr_d> vec1dr_in(9);
  ekat::host_to_device({d.cape, d.dcape, d.dsubcld, d.geos, d.landfrac, d.pbl_hgt, d.prec, d.rliq, d.tpert}, d.pcols, vec1dr_in);

  std::vector<view2dr_d> vec2dr_in(24);
  std::vector<int> vec2dr_in_0_sizes = {d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols};
  std::vector<int> vec2dr_in_1_sizes = {d.pver, d.pver, d.pver, d.pver, d.pver, d.pverp, d.pver, d.pver, d.pver, d.pver, d.pver, d.pverp, d.pver, d.pverp, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pver, d.pverp, d.pver, d.pver};
  ekat::host_to_device({d.detr_up, d.dlf, d.entr_dn, d.entr_up, d.heat, d.mcon, d.mflx_dn, d.mflx_up, d.omega, d.p_del, d.p_del_in, d.p_int_in, d.p_mid_in, d.pflx, d.q_mid_in, d.q_star, d.ql, d.qtnd, d.rprd, d.t_mid, d.t_star, d.z_int_in, d.z_mid_in, d.zdu}, vec2dr_in_0_sizes, vec2dr_in_1_sizes, vec2dr_in);

  std::vector<view1di_d> vec1di_in(5);
  ekat::host_to_device({d.gather_index, d.jcbot, d.jctop, d.jt, d.msemax_klev_g}, d.pcols, vec1di_in);

  view1dr_d
    cape_d(vec1dr_in[0]),
    dcape_d(vec1dr_in[1]),
    dsubcld_d(vec1dr_in[2]),
    geos_d(vec1dr_in[3]),
    landfrac_d(vec1dr_in[4]),
    pbl_hgt_d(vec1dr_in[5]),
    prec_d(vec1dr_in[6]),
    rliq_d(vec1dr_in[7]),
    tpert_d(vec1dr_in[8]);

  view2dr_d
    detr_up_d(vec2dr_in[0]),
    dlf_d(vec2dr_in[1]),
    entr_dn_d(vec2dr_in[2]),
    entr_up_d(vec2dr_in[3]),
    heat_d(vec2dr_in[4]),
    mcon_d(vec2dr_in[5]),
    mflx_dn_d(vec2dr_in[6]),
    mflx_up_d(vec2dr_in[7]),
    omega_d(vec2dr_in[8]),
    p_del_d(vec2dr_in[9]),
    p_del_in_d(vec2dr_in[10]),
    p_int_in_d(vec2dr_in[11]),
    p_mid_in_d(vec2dr_in[12]),
    pflx_d(vec2dr_in[13]),
    q_mid_in_d(vec2dr_in[14]),
    q_star_d(vec2dr_in[15]),
    ql_d(vec2dr_in[16]),
    qtnd_d(vec2dr_in[17]),
    rprd_d(vec2dr_in[18]),
    t_mid_d(vec2dr_in[19]),
    t_star_d(vec2dr_in[20]),
    z_int_in_d(vec2dr_in[21]),
    z_mid_in_d(vec2dr_in[22]),
    zdu_d(vec2dr_in[23]);

  view1di_d
    gather_index_d(vec1di_in[0]),
    jcbot_d(vec1di_in[1]),
    jctop_d(vec1di_in[2]),
    jt_d(vec1di_in[3]),
    msemax_klev_g_d(vec1di_in[4]);

  const auto policy = ekat::TeamPolicyFactory<ExeSpace>::get_default_team_policy(d.pcols, d.pver);

  // unpack data scalars because we do not want the lambda to capture d
  const Real time_step = d.time_step;
  const Int pver = d.pver;
  const Int pverp = d.pverp;
  const bool is_first_step = d.is_first_step;
  view0di_d lengath_d("lengath_h");
  auto      lengath_h = Kokkos::create_mirror_view(lengath_d);

  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    // Get single-column subviews of all inputs, shouldn't need any i-indexing
    // after this.
    const auto t_mid_c = ekat::subview(t_mid_d, i);
    const auto q_mid_in_c = ekat::subview(q_mid_in_d, i);
    const auto omega_c = ekat::subview(omega_d, i);
    const auto p_mid_in_c = ekat::subview(p_mid_in_d, i);
    const auto p_int_in_c = ekat::subview(p_int_in_d, i);
    const auto p_del_in_c = ekat::subview(p_del_in_d, i);
    const auto z_mid_in_c = ekat::subview(z_mid_in_d, i);
    const auto z_int_in_c = ekat::subview(z_int_in_d, i);
    const auto t_star_c = ekat::subview(t_star_d, i);
    const auto q_star_c = ekat::subview(q_star_d, i);
    const auto heat_c = ekat::subview(heat_d, i);
    const auto qtnd_c = ekat::subview(qtnd_d, i);
    const auto mcon_c = ekat::subview(mcon_d, i);
    const auto pflx_c = ekat::subview(pflx_d, i);
    const auto zdu_c = ekat::subview(zdu_d, i);
    const auto mflx_up_c = ekat::subview(mflx_up_d, i);
    const auto entr_up_c = ekat::subview(entr_up_d, i);
    const auto detr_up_c = ekat::subview(detr_up_d, i);
    const auto mflx_dn_c = ekat::subview(mflx_dn_d, i);
    const auto entr_dn_c = ekat::subview(entr_dn_d, i);
    const auto p_del_c = ekat::subview(p_del_d, i);
    const auto ql_c = ekat::subview(ql_d, i);
    const auto rprd_c = ekat::subview(rprd_d, i);
    const auto dlf_c = ekat::subview(dlf_d, i);

    ZMF::zm_conv_main(
      team,
      pver,
      pverp,
      is_first_step,
      time_step,
      t_mid_c,
      q_mid_in_c,
      omega_c,
      p_mid_in_c,
      p_int_in_c,
      p_del_in_c,
      geos_d(i),
      z_mid_in_c,
      z_int_in_c,
      pbl_hgt_d(i),
      tpert_d(i),
      landfrac_d(i),
      t_star_c,
      q_star_c,
      lengath_d(),
      gather_index_d(i),
      msemax_klev_g_d(i),
      jctop_d(i),
      jcbot_d(i),
      jt_d(i),
      prec_d(i),
      heat_c,
      qtnd_c,
      cape_d(i),
      dcape_d(i),
      mcon_c,
      pflx_c,
      zdu_c,
      mflx_up_c,
      entr_up_c,
      detr_up_c,
      mflx_dn_c,
      entr_dn_c,
      p_del_c,
      dsubcld_d(i),
      ql_c,
      rliq_d(i),
      rprd_c,
      dlf_c);
  });

  // Get outputs back, start with scalars
  Kokkos::deep_copy(lengath_h, lengath_d);
  d.lengath = lengath_h();

  // Now get arrays
  std::vector<view1dr_d> vec1dr_out = {cape_d, dcape_d, dsubcld_d, prec_d, rliq_d};
  ekat::device_to_host({d.cape, d.dcape, d.dsubcld, d.prec, d.rliq}, d.pcols, vec1dr_out);

  std::vector<view2dr_d> vec2dr_out = {detr_up_d, dlf_d, entr_dn_d, entr_up_d, heat_d, mcon_d, mflx_dn_d, mflx_up_d, p_del_d, pflx_d, ql_d, qtnd_d, rprd_d, zdu_d};
  std::vector<int> vec2dr_out_0_sizes = {d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols};
  std::vector<int> vec2dr_out_1_sizes = {d.pver, d.pver, d.pver, d.pver, d.pver, d.pverp, d.pver, d.pver, d.pver, d.pverp, d.pver, d.pver, d.pver, d.pver};
  ekat::device_to_host({d.detr_up, d.dlf, d.entr_dn, d.entr_up, d.heat, d.mcon, d.mflx_dn, d.mflx_up, d.p_del, d.pflx, d.ql, d.qtnd, d.rprd, d.zdu}, vec2dr_out_0_sizes, vec2dr_out_1_sizes, vec2dr_out);

  std::vector<view1di_d> vec1di_out = {gather_index_d, jcbot_d, jctop_d, jt_d, msemax_klev_g_d};
  ekat::device_to_host({d.gather_index, d.jcbot, d.jctop, d.jt, d.msemax_klev_g}, d.pcols, vec1di_out);

  zm_finalize_cxx();
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
// end glue impls

} // namespace zm
} // namespace scream
