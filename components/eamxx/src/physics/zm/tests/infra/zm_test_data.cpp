#include "zm_test_data.hpp"
#include "zm_functions.hpp"

#include <ekat_pack_kokkos.hpp>

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

void zm_transport_tracer_bridge_f(Int pcols, Int ncol, Int pver, bool* doconvtran, Real* q, Int ncnst, Real* mu, Real* md, Real* du, Real* eu, Real* ed, Real* dp, Int* jt, Int* mx, Int* ideep, Int il1g, Int il2g, Real* fracis, Real* dqdt, Real* dpdry, Real dt);
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
  zm_common_init_f(); // Might need more specific init
  entropy_bridge_f(d.tk, d.p, d.qtot, &d.entropy);
  zm_common_finalize_f();
  d.transition<ekat::TransposeDirection::f2c>();
}

void entropy(EntropyData& d)
{
  zm_common_init(); // Might need more specific init

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
  zm_common_init_f(); // Might need more specific init
  zm_transport_tracer_f(d.pcols, d.ncol, d.pver, d.doconvtran, d.q, d.ncnst, d.mu, d.md, d.du, d.eu, d.ed, d.dp, d.jt, d.mx, d.ideep, d.il1g, d.il2g, d.fracis, d.dqdt, d.dpdry, d.dt);
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

  // unpack data scalars because we do not want the lambda to capture d
  const Real dt = d.dt;
  const Int il1g = d.il1g;
  const Int il2g = d.il2g;
  const Int ncnst = d.ncnst;
  const Int ncol = d.ncol;
  const Int pver = d.pver;

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
      pver,
      doconvtran_c,
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
      ideep_d(i),
      il1g,
      il2g,
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
// end glue impls

} // namespace zm
} // namespace scream
