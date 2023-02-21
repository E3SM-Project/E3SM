#ifndef P3_MAIN_IMPL_PART_1_HPP
#define P3_MAIN_IMPL_PART_1_HPP

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

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::p3_main_part1(
  const MemberType& team,
  const Int& nk,
  const bool& predictNc,
  const bool& do_prescribed_CCN,
  const Scalar& dt,
  const uview_1d<const Spack>& pres,
  const uview_1d<const Spack>& dpres,
  const uview_1d<const Spack>& dz,
  const uview_1d<const Spack>& nc_nuceat_tend,
  const uview_1d<const Spack>& nccn_prescribed,
  const uview_1d<const Spack>& inv_exner,
  const uview_1d<const Spack>& exner,
  const uview_1d<const Spack>& inv_cld_frac_l,
  const uview_1d<const Spack>& inv_cld_frac_i,
  const uview_1d<const Spack>& inv_cld_frac_r,
  const uview_1d<const Spack>& latent_heat_vapor,
  const uview_1d<const Spack>& latent_heat_sublim,
  const uview_1d<const Spack>& latent_heat_fusion,
  const uview_1d<Spack>& T_atm,
  const uview_1d<Spack>& rho,
  const uview_1d<Spack>& inv_rho,
  const uview_1d<Spack>& qv_sat_l,
  const uview_1d<Spack>& qv_sat_i,
  const uview_1d<Spack>& qv_supersat_i,
  const uview_1d<Spack>& rhofacr,
  const uview_1d<Spack>& rhofaci,
  const uview_1d<Spack>& acn,
  const uview_1d<Spack>& qv,
  const uview_1d<Spack>& th_atm,
  const uview_1d<Spack>& qc,
  const uview_1d<Spack>& nc,
  const uview_1d<Spack>& qr,
  const uview_1d<Spack>& nr,
  const uview_1d<Spack>& qi,
  const uview_1d<Spack>& ni,
  const uview_1d<Spack>& qm,
  const uview_1d<Spack>& bm,
  const uview_1d<Spack>& qc_incld,
  const uview_1d<Spack>& qr_incld,
  const uview_1d<Spack>& qi_incld,
  const uview_1d<Spack>& qm_incld,
  const uview_1d<Spack>& nc_incld,
  const uview_1d<Spack>& nr_incld,
  const uview_1d<Spack>& ni_incld,
  const uview_1d<Spack>& bm_incld,
  bool& nucleationPossible,
  bool& hydrometeorsPresent)
{
  // Get access to saturation functions
  using physics = scream::physics::Functions<Scalar, Device>;

  // load constants into local vars
  constexpr Scalar g            = C::gravit;
  constexpr Scalar rho_1000mb   = C::RHO_1000MB;
  constexpr Scalar rho_600mb    = C::RHO_600MB;
  constexpr Scalar rho_h2o      = C::RHO_H2O;
  constexpr Scalar nccnst       = C::NCCNST;
  constexpr Scalar T_zerodegc   = C::T_zerodegc;
  constexpr Scalar qsmall       = C::QSMALL;
  constexpr Scalar inv_cp       = C::INV_CP;

  nucleationPossible = false;
  hydrometeorsPresent = false;
  team.team_barrier();

  const Int nk_pack = ekat::npack<Spack>(nk);

  //
  // calculate some time-varying atmospheric variables
  // AaronDonahue - changed "rho" to be defined on nonhydrostatic
  // assumption, consistent with pressure based coordinate system
  //              - moved latent heat calculation to above.  Latent
  // heat is determined by calling a p3_util function so that it
  // can be made consistent with E3SM definition of latent heat
  //
  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, nk_pack), [&] (Int k) {

    const auto range_pack = ekat::range<IntSmallPack>(k*Spack::n);
    const auto range_mask = range_pack < nk;

    rho(k)          = dpres(k)/dz(k) / g;
    inv_rho(k)      = 1 / rho(k);
    qv_sat_l(k)     = physics::qv_sat(T_atm(k), pres(k), false, range_mask, physics::MurphyKoop, "p3::p3_main_part1 (liquid)");
    qv_sat_i(k)     = physics::qv_sat(T_atm(k), pres(k), true,  range_mask, physics::MurphyKoop, "p3::p3_main_part1 (ice)");

    qv_supersat_i(k) = qv(k) / qv_sat_i(k) - 1;

    rhofacr(k) = pow(rho_1000mb * inv_rho(k), sp(.54));
    rhofaci(k) = pow(rho_600mb * inv_rho(k), sp(.54));
    Spack dum  = sp(1.496e-6) * pow(T_atm(k), sp(1.5)) / (T_atm(k) + 120); // this is mu
    acn(k)     = g * rho_h2o / (18 * dum); // 'a' parameter for droplet fallspeed (Stokes' law)

    if ( (T_atm(k) < T_zerodegc && qv_supersat_i(k) >= -0.05).any() ) {
      nucleationPossible = true;
    }

    // apply mass clipping if dry and mass is sufficiently small
    // (implying all mass is expected to evaporate/sublimate in one time step)
    auto drymass = qc(k) < qsmall;
    auto not_drymass = !drymass && range_mask;
    qv(k).set(drymass, qv(k) + qc(k));
    th_atm(k).set(drymass, th_atm(k) - inv_exner(k) * qc(k) * latent_heat_vapor(k) * inv_cp);
    qc(k).set(drymass, 0);
    nc(k).set(drymass, 0);
    if ( not_drymass.any() ) {
      hydrometeorsPresent = true; // updated further down
      // Apply droplet activation here (before other microphysical processes) for consistency with qc increase by saturation
      // adjustment already applied in macrophysics. If prescribed drop number is used, this is also a good place to
      // prescribe that value

      if (do_prescribed_CCN) {
         nc(k).set(not_drymass, max(nc(k), nccn_prescribed(k)/inv_cld_frac_l(k)));
      } else if (predictNc) {
         nc(k).set(not_drymass, max(nc(k) + nc_nuceat_tend(k) * dt, 0.0));
      } else {
         // nccnst is in units of #/m3 so needs to be converted.
         nc(k).set(not_drymass, nccnst*inv_rho(k));
      }

    }

    drymass = qr(k) < qsmall;
    not_drymass = !drymass && range_mask;
    qv(k).set(drymass, qv(k) + qr(k));
    th_atm(k).set(drymass, th_atm(k) - inv_exner(k) * qr(k) * latent_heat_vapor(k) * inv_cp);
    qr(k).set(drymass, 0);
    nr(k).set(drymass, 0);
    if ( not_drymass.any() ) {
      hydrometeorsPresent = true; // updated further down
    }

    drymass = (qi(k) < qsmall || (qi(k) < 1.e-8 && qv_supersat_i(k) < -0.1));
    not_drymass = !drymass && range_mask;
    qv(k).set(drymass, qv(k) + qi(k));
    th_atm(k).set(drymass, th_atm(k) - inv_exner(k) * qi(k) * latent_heat_sublim(k) * inv_cp);
    qi(k).set(drymass, 0);
    ni(k).set(drymass, 0);
    qm(k).set(drymass, 0);
    bm(k).set(drymass, 0);
    if ( not_drymass.any() ) {
      hydrometeorsPresent = true; // final update
    }

    drymass = (qi(k) >= qsmall && qi(k) < 1.e-8 && T_atm(k) >= T_zerodegc);
    qr(k).set(drymass, qr(k) + qi(k));
    th_atm(k).set(drymass, th_atm(k) - inv_exner(k) * qi(k) * latent_heat_fusion(k) * inv_cp);
    qi(k).set(drymass, 0);
    ni(k).set(drymass, 0);
    qm(k).set(drymass, 0);
    bm(k).set(drymass, 0);

    T_atm(k) = th_atm(k) * exner(k);

    calculate_incloud_mixingratios(
      qc(k), qr(k), qi(k), qm(k), nc(k), nr(k), ni(k), bm(k),
      inv_cld_frac_l(k), inv_cld_frac_i(k), inv_cld_frac_r(k),
      qc_incld(k), qr_incld(k), qi_incld(k), qm_incld(k), nc_incld(k), nr_incld(k), ni_incld(k), bm_incld(k));
  });
  team.team_barrier();
}

} // namespace p3
} // namespace scream

#endif // P3_MAIN_IMPL_PART_1_HPP
