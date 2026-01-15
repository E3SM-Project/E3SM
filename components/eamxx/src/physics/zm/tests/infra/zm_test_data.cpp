#include "zm_test_data.hpp"

#include <ekat_pack_kokkos.hpp>

#include <random>

using scream::Real;
using scream::Int;

//
// A C++ interface to ZM fortran calls and vice versa
//

namespace scream {
namespace zm {

extern "C" {

void zm_find_mse_max_f( Int  pcols,
                        Int  ncol,
                        Int  pver,
                        Int  num_msg,
                        Int  *msemax_top_k,
                        bool pergro_active,
                        Real *temperature,
                        Real *zmid,
                        Real *sp_humidity,
                        Int *msemax_klev,
                        Real *mse_max_val );

void ientropy_bridge_f(Int rcall, Real s, Real p, Real qt, Real* t, Real* qst, Real tfg);
} // extern "C" : end _f decls

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
  zm_common_init();
  ientropy_c(d.rcall, d.s, d.p, d.qt, &d.t, &d.qst, d.tfg);
  d.transition<ekat::TransposeDirection::f2c>();
}

void ientropy(IentropyData& d)
{
  zm_common_init(); // Might need more specific init

  // create device views and copy
  const auto policy = ekat::TeamPolicyFactory<ExeSpace>::get_default_team_policy(1, 1);

  // unpack data scalars because we do not want the lambda to capture d
  const Real p = d.p;
  const Real qt = d.qt;
  const Real s = d.s;
  const Real tfg = d.tfg;
  const Int rcall = d.rcall;
  view0dr_h qst_h("qst_h");
  view0dr_d qst_d = Kokkos::create_mirror_view_and_copy(DefaultDevice(), qst_h);
  view0dr_h t_h("t_h");
  view0dr_d t_d = Kokkos::create_mirror_view_and_copy(DefaultDevice(), t_h);

  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    // Get single-column subviews of all inputs, shouldn't need any i-indexing
    // after this.

    SHF::ientropy(
      team,
      rcall,
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

// end glue impls

} // namespace zm
} // namespace scream
