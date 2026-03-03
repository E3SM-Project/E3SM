#include "p3_functions.hpp" // for ETI only but harmless for GPU
#include "share/physics/physics_functions.hpp" // also for ETI not on GPUs
#include "share/physics/physics_saturation_impl.hpp"

#include <ekat_subview_utils.hpp>
#include <ekat_team_policy_utils.hpp>

namespace scream {
namespace p3 {

/*
 * Implementation of p3 main function. Clients should NOT #include
 * this file, #include p3_functions.hpp instead.
 */

template <>
void Functions<Real,DefaultDevice>
::p3_main_part1_disp(
  const Int& nj,
  const Int& nk,
  const bool& predictNc,
  const bool& prescribedCCN,
  const Scalar& dt,
  const uview_2d<const Pack>& pres,
  const uview_2d<const Pack>& dpres,
  const uview_2d<const Pack>& dz,
  const uview_2d<const Pack>& nc_nuceat_tend,
  const uview_2d<const Pack>& nccn_prescribed,
  const uview_2d<const Pack>& inv_exner,
  const uview_2d<const Pack>& exner,
  const uview_2d<const Pack>& inv_cld_frac_l,
  const uview_2d<const Pack>& inv_cld_frac_i,
  const uview_2d<const Pack>& inv_cld_frac_r,
  const uview_2d<Pack>& T_atm,
  const uview_2d<Pack>& rho,
  const uview_2d<Pack>& inv_rho,
  const uview_2d<Pack>& qv_sat_l,
  const uview_2d<Pack>& qv_sat_i,
  const uview_2d<Pack>& qv_supersat_i,
  const uview_2d<Pack>& rhofacr,
  const uview_2d<Pack>& rhofaci,
  const uview_2d<Pack>& acn,
  const uview_2d<Pack>& qv,
  const uview_2d<Pack>& th_atm,
  const uview_2d<Pack>& qc,
  const uview_2d<Pack>& nc,
  const uview_2d<Pack>& qr,
  const uview_2d<Pack>& nr,
  const uview_2d<Pack>& qi,
  const uview_2d<Pack>& ni,
  const uview_2d<Pack>& qm,
  const uview_2d<Pack>& bm,
  const uview_2d<Pack>& qc_incld,
  const uview_2d<Pack>& qr_incld,
  const uview_2d<Pack>& qi_incld,
  const uview_2d<Pack>& qm_incld,
  const uview_2d<Pack>& nc_incld,
  const uview_2d<Pack>& nr_incld,
  const uview_2d<Pack>& ni_incld,
  const uview_2d<Pack>& bm_incld,
  const uview_1d<bool>& nucleationPossible,
  const uview_1d<bool>& hydrometeorsPresent,
  const P3Runtime& runtime_options)
{
  using ExeSpace = typename KT::ExeSpace;
  using TPF      = ekat::TeamPolicyFactory<ExeSpace>;

  const Int nk_pack = ekat::npack<Pack>(nk);
  const auto policy = TPF::get_default_team_policy(nj, nk_pack);
  // p3_cloud_sedimentation loop
  Kokkos::parallel_for("p3_main_part1",
      policy, KOKKOS_LAMBDA(const MemberType& team) {

    const Int i = team.league_rank();

    p3_main_part1(
      team, nk, predictNc, prescribedCCN, dt,
      ekat::subview(pres, i), ekat::subview(dpres, i), ekat::subview(dz, i), ekat::subview(nc_nuceat_tend, i),
      ekat::subview(nccn_prescribed, i), ekat::subview(inv_exner, i), ekat::subview(exner, i), ekat::subview(inv_cld_frac_l, i),
      ekat::subview(inv_cld_frac_i, i), ekat::subview(inv_cld_frac_r, i),
      ekat::subview(T_atm, i), ekat::subview(rho, i), ekat::subview(inv_rho, i), ekat::subview(qv_sat_l, i),
      ekat::subview(qv_sat_i, i), ekat::subview(qv_supersat_i, i), ekat::subview(rhofacr, i), ekat::subview(rhofaci, i), ekat::subview(acn, i),
      ekat::subview(qv, i), ekat::subview(th_atm, i), ekat::subview(qc, i), ekat::subview(nc, i), ekat::subview(qr, i), ekat::subview(nr, i), ekat::subview(qi, i),
      ekat::subview(ni, i), ekat::subview(qm, i), ekat::subview(bm, i), ekat::subview(qc_incld, i), ekat::subview(qr_incld, i), ekat::subview(qi_incld, i),
      ekat::subview(qm_incld, i), ekat::subview(nc_incld, i), ekat::subview(nr_incld, i), ekat::subview(ni_incld, i), ekat::subview(bm_incld, i),
      nucleationPossible(i), hydrometeorsPresent(i), runtime_options);

  });
}

} // namespace p3
} // namespace scream

