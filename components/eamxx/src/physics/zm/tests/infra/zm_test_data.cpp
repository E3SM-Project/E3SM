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
using view1di_d = ZMF::view_1d<Int>;
using view1db_d = ZMF::view_1d<bool>;
using view1dr_d = ZMF::view_1d<Real>;
using view2dr_d = ZMF::view_2d<Real>;
using view3dr_d = ZMF::view_3d<Real>;

using WSM = typename ZMF::WorkspaceManager;

extern "C" {

void zm_common_init_bridge_f();

void zm_common_finalize_bridge_f();

void zm_find_mse_max_f(Int pcols, Int ncol, Int pver, Int num_msg, Int *msemax_top_k, bool pergro_active, Real *temperature, Real *zmid, Real *sp_humidity, Int *msemax_klev, Real *mse_max_val);

void ientropy_bridge_f(Real s, Real p, Real qt, Real* t, Real* qst, Real tfg);

void entropy_bridge_f(Real tk, Real p, Real qtot, Real* entropy);

void zm_transport_tracer_bridge_f(Int pcols, Int pver, bool* doconvtran, Real* q, Int ncnst, Real* mu, Real* md, Real* du, Real* eu, Real* ed, Real* dp, Int* jt, Int* mx, Int* ideep, Int il1g, Int il2g, Real* fracis, Real* dqdt, Real* dpdry, Real dt);

void zm_transport_momentum_bridge_f(Int pcols, Int ncol, Int pver, Int pverp, Real* wind_in, Int nwind, Real* mu, Real* md, Real* du, Real* eu, Real* ed, Real* dp, Int* jt, Int* mx, Int* ideep, Int il1g, Int il2g, Real* wind_tend, Real* pguall, Real* pgdall, Real* icwu, Real* icwd, Real dt, Real* seten);

void compute_dilute_cape_bridge_f(Int pcols, Int ncol, Int pver, Int pverp, Int num_cin, Int num_msg, Real* sp_humidity_in, Real* temperature_in, Real* zmid, Real* pmid, Real* pint, Int* pblt, Real* tpert, Real* parcel_temp, Real* parcel_qsat, Int* msemax_klev, Real* lcl_temperature, Int* lcl_klev, Int* eql_klev, Real* cape, bool calc_msemax_klev, Int* prev_msemax_klev, bool use_input_tq_mx, Real* q_mx, Real* t_mx);

void find_mse_max_bridge_f(Int pcols, Int ncol, Int pver, Int num_msg, Int* msemax_top_k, bool pergro_active, Real* temperature, Real* zmid, Real* sp_humidity, Int* msemax_klev, Real* mse_max_val);

void compute_dilute_parcel_bridge_f(Int pcols, Int ncol, Int pver, Int num_msg, Int* klaunch, Real* pmid, Real* temperature, Real* sp_humidity, Real* tpert, Int* pblt, Real* parcel_temp, Real* parcel_vtemp, Real* parcel_qsat, Real* lcl_pmid, Real* lcl_temperature, Int* lcl_klev);

void compute_cape_from_parcel_bridge_f(Int pcols, Int ncol, Int pver, Int pverp, Int num_cin, Int num_msg, Real* temperature, Real* tv, Real* zmid, Real* sp_humidity, Real* pint, Int* msemax_klev, Real* lcl_pmid, Int* lcl_klev, Real* parcel_qsat, Real* parcel_temp, Real* parcel_vtemp, Int* eql_klev, Real* cape);
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

void zm_find_mse_max(zm_data_find_mse_max& d){
  d.transition<ekat::TransposeDirection::c2f>();
  zm_find_mse_max_f( d.pcols,
                     d.ncol,
                     d.pver,
                     d.num_msg,
                     d.msemax_top_k,
                     d.pergro_active,
                     d.temperature,
                     d.zmid,
                     d.sp_humidity,
                     d.msemax_klev,
                     d.mse_max_val );
  d.transition<ekat::TransposeDirection::f2c>();
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

  WSM wsm(d.pverp, 12, policy);

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
  compute_cape_from_parcel_bridge_f(d.pcols, d.ncol, d.pver, d.pverp, d.num_cin, d.num_msg, d.temperature, d.tv, d.zmid, d.sp_humidity, d.pint, d.msemax_klev, d.lcl_pmid, d.lcl_klev, d.parcel_qsat, d.parcel_temp, d.parcel_vtemp, d.eql_klev, d.cape);
  zm_common_finalize_f();
  d.transition<ekat::TransposeDirection::f2c>();
}

void compute_cape_from_parcel(ComputeCapeFromParcelData& d)
{
  zm_common_init();

  // create device views and copy
  std::vector<view1dr_d> vec1dr_in(2);
  ekat::host_to_device({d.cape, d.lcl_pmid}, d.pcols, vec1dr_in);

  std::vector<view2dr_d> vec2dr_in(8);
  std::vector<int> vec2dr_in_0_sizes = {d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols, d.pcols};
  std::vector<int> vec2dr_in_1_sizes = {d.pver, d.pver, d.pver, d.pverp, d.pver, d.pver, d.pver, d.pver};
  ekat::host_to_device({d.parcel_qsat, d.parcel_temp, d.parcel_vtemp, d.pint, d.sp_humidity, d.temperature, d.tv, d.zmid}, vec2dr_in_0_sizes, vec2dr_in_1_sizes, vec2dr_in);

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
    tv_d(vec2dr_in[6]),
    zmid_d(vec2dr_in[7]);

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
    const auto zmid_c = ekat::subview(zmid_d, i);
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
      zmid_c,
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
// end glue impls

} // namespace zm
} // namespace scream
