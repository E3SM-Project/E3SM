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


// Wrapper around gw_init for cxx
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
    ZMF::ientropy(
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
// end glue impls

} // namespace zm
} // namespace scream
