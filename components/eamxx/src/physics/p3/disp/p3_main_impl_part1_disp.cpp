
#include "p3_functions.hpp" // for ETI only but harmless for GPU
#include "physics/share/physics_functions.hpp" // also for ETI not on GPUs
#include "physics/share/physics_saturation_impl.hpp"

#include "ekat/kokkos/ekat_subview_utils.hpp"

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
  const uview_2d<const Spack>& pres,
  const uview_2d<const Spack>& dpres,
  const uview_2d<const Spack>& dz,
  const uview_2d<const Spack>& nc_nuceat_tend,
  const uview_2d<const Spack>& nccn_prescribed,
  const uview_2d<const Spack>& inv_exner,
  const uview_2d<const Spack>& exner,
  const uview_2d<const Spack>& inv_cld_frac_l,
  const uview_2d<const Spack>& inv_cld_frac_i,
  const uview_2d<const Spack>& inv_cld_frac_r,
  const uview_2d<Spack>& T_atm,
  const uview_2d<Spack>& rho,
  const uview_2d<Spack>& inv_rho,
  const uview_2d<Spack>& qv_sat_l,
  const uview_2d<Spack>& qv_sat_i,
  const uview_2d<Spack>& qv_supersat_i,
  const uview_2d<Spack>& rhofacr,
  const uview_2d<Spack>& rhofaci,
  const uview_2d<Spack>& acn,
  const uview_2d<Spack>& qv,
  const uview_2d<Spack>& th_atm,
  const uview_2d<Spack>& qc,
  const uview_2d<Spack>& nc,
  const uview_2d<Spack>& qr,
  const uview_2d<Spack>& nr,
  const uview_2d<Spack>& qi,
  const uview_2d<Spack>& ni,
  const uview_2d<Spack>& qm,
  const uview_2d<Spack>& bm,
  const uview_2d<Spack>& qc_incld,
  const uview_2d<Spack>& qr_incld,
  const uview_2d<Spack>& qi_incld,
  const uview_2d<Spack>& qm_incld,
  const uview_2d<Spack>& nc_incld,
  const uview_2d<Spack>& nr_incld,
  const uview_2d<Spack>& ni_incld,
  const uview_2d<Spack>& bm_incld,
  const uview_1d<bool>& nucleationPossible,
  const uview_1d<bool>& hydrometeorsPresent,
  const P3Runtime& runtime_options)
{
  using ExeSpace = typename KT::ExeSpace;
  const Int nk_pack = ekat::npack<Spack>(nk);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(nj, nk_pack);
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

