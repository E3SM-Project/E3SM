#ifndef ZM_ZM_CONV_MCSP_TEND_IMPL_HPP
#define ZM_ZM_CONV_MCSP_TEND_IMPL_HPP

#include "zm_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace zm {

/*
 * Implementation of zm zm_conv_mcsp_tend. Clients should NOT
 * #include this file, but include zm_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::zm_conv_mcsp_tend(
  // Inputs
  const MemberType& team,
  const Workspace& workspace,
  const ZmRuntimeOpt& runtime_opt,
  const Int& pver, // number of mid-point vertical levels
  const Int& pverp, // number of interface vertical levels
  const Real& ztodt, // 2x physics time step
  const Int& jctop, // cloud top level indices
  const uview_1d<const Real>& state_pmid, // physics state mid-point pressure
  const uview_1d<const Real>& state_pint, // physics state interface pressure
  const uview_1d<const Real>& state_pdel, // physics state pressure thickness
  const uview_1d<const Real>& state_s, // physics state dry energy
  const uview_1d<const Real>& state_q, // physics state specific humidity
  const uview_1d<const Real>& state_u, // physics state u momentum
  const uview_1d<const Real>& state_v, // physics state v momentum
  const uview_1d<const Real>& ptend_zm_s, // input ZM tendency for dry energy (DSE)
  const uview_1d<const Real>& ptend_zm_q, // input ZM tendency for specific humidity (qv)
  // Inputs/Outputs
  const uview_1d<Real>& ptend_s, // output tendency of DSE
  const uview_1d<Real>& ptend_q, // output tendency of qv
  const uview_1d<Real>& ptend_u, // output tendency of u-wind
  const uview_1d<Real>& ptend_v, // output tendency of v-wind
  // Outputs
  const uview_1d<Real>& mcsp_dt_out, // final MCSP tendency for DSE
  const uview_1d<Real>& mcsp_dq_out, // final MCSP tendency for qv
  const uview_1d<Real>& mcsp_du_out, // final MCSP tendency for u wind
  const uview_1d<Real>& mcsp_dv_out, // final MCSP tendency for v wind
  Real& mcsp_freq, // MSCP frequency for output
  Real& mcsp_shear, // shear used to check against threshold
  Real& zm_depth) // pressure depth of ZM heating
{
  //----------------------------------------------------------------------------
  // Purpose: perform MCSP tendency calculations
  //----------------------------------------------------------------------------

  if (!runtime_opt.mcsp_enabled) return;

  //----------------------------------------------------------------------------
  // initialize variables

  const bool do_mcsp_t = (runtime_opt.mcsp_t_coeff > 0);
  const bool do_mcsp_q = (runtime_opt.mcsp_q_coeff > 0);
  const bool do_mcsp_u = (runtime_opt.mcsp_u_coeff > 0);
  const bool do_mcsp_v = (runtime_opt.mcsp_v_coeff > 0);

  // Allocate temporary arrays
  uview_1d<Real> mcsp_tend_s, mcsp_tend_q, mcsp_tend_u, mcsp_tend_v;
  workspace.template take_many_contiguous_unsafe<4>(
    {"mcsp_tend_s", "mcsp_tend_q", "mcsp_tend_u", "mcsp_tend_v"},
    {&mcsp_tend_s, &mcsp_tend_q, &mcsp_tend_u, &mcsp_tend_v});

  Real zm_avg_tend_s = 0;     // mass weighted column average DSE tendency from ZM
  Real zm_avg_tend_q = 0;     // mass weighted column average qv tendency from ZM
  Real pdel_sum = 0;          // column integrated pressure thickness
  Real mcsp_avg_tend_s = 0;   // mass weighted column average MCSP tendency of DSE
  Real mcsp_avg_tend_q = 0;   // mass weighted column average MCSP tendency of qv
  Real mcsp_avg_tend_k = 0;   // mass weighted column average MCSP tendency of kinetic energy

  // Initialize arrays
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, pver), [&] (const Int& k) {
    mcsp_tend_s(k) = 0;
    mcsp_tend_q(k) = 0;
    mcsp_tend_u(k) = 0;
    mcsp_tend_v(k) = 0;
    mcsp_dt_out(k) = 0;
    mcsp_dq_out(k) = 0;
    mcsp_du_out(k) = 0;
    mcsp_dv_out(k) = 0;
  });
  team.team_barrier();

  //----------------------------------------------------------------------------
  // calculate shear

  zm_conv_mcsp_calculate_shear(team, pver, state_pmid, state_u, mcsp_shear);

  //----------------------------------------------------------------------------
  // calculate mass weighted column average tendencies from ZM

  zm_depth = 0;
  if (jctop != pver - 1) {
    // integrate pressure and ZM tendencies over column using parallel reduction
    Kokkos::parallel_reduce(Kokkos::TeamVectorRange(team, jctop, pver),
      [&] (const Int& k, Real& tend_s_sum, Real& tend_q_sum, Real& pdel_tot) {
        tend_s_sum += ptend_zm_s(k) * state_pdel(k);
        tend_q_sum += ptend_zm_q(k) * state_pdel(k);
        pdel_tot += state_pdel(k);
      },
      zm_avg_tend_s, zm_avg_tend_q, pdel_sum);
    team.team_barrier();

    // normalize integrated ZM tendencies by total mass
    zm_avg_tend_s /= pdel_sum;
    zm_avg_tend_q /= pdel_sum;
    // calculate diagnostic zm_depth
    zm_depth = state_pint(pver) - state_pmid(jctop);
  }
  team.team_barrier();

  //----------------------------------------------------------------------------
  // Note: To conserve total energy we need to account for the kinteic energy tendency
  // which we can obtain from the velocity tendencies based on the following:
  //   KE_new = (u_new^2 + v_new^2)/2
  //          = [ (u_old+du)^2 + (v_old+dv)^2 ]/2
  //          = [ ( u_old^2 + 2*u_old*du + du^2 ) + ( v_old^2 + 2*v_old*dv + dv^2 ) ]/2
  //          = ( u_old^2 + v_old^2 )/2 + ( 2*u_old*du + du^2 + 2*v_old*dv + dv^2 )/2
  //          = KE_old + [ 2*u_old*du + du^2 + 2*v_old*dv + dv^2 ] /2

  //----------------------------------------------------------------------------
  // calculate MCSP tendencies

  // check that ZM produced tendencies over a depth that exceeds the threshold
  if (zm_depth >= ZMC::MCSP_conv_depth_min) {
    // check that ZM provided a non-zero column total heating
    if (zm_avg_tend_s > 0) {
      // check that there is sufficient wind shear to justify coherent organization
      if (std::abs(mcsp_shear) >= ZMC::MCSP_shear_min &&
          std::abs(mcsp_shear) < ZMC::MCSP_shear_max) {
        // Calculate tendencies and integrate them using parallel reduce
        Kokkos::parallel_reduce(Kokkos::TeamVectorRange(team, jctop, pver),
          [&] (const Int& k, Real& avg_s, Real& avg_q, Real& avg_k) {

            // See eq 7-8 of Moncrieff et al. (2017) - also eq (5) of Moncrieff & Liu (2006)
            const Real pdepth_mid_k = state_pint(pver) - state_pmid(k);
            const Real pdepth_total = state_pint(pver) - state_pmid(jctop);

            // specify the assumed vertical structure
            if (do_mcsp_t) mcsp_tend_s(k) = -1 * runtime_opt.mcsp_t_coeff * std::sin(2 * PC::Pi * (pdepth_mid_k / pdepth_total));
            if (do_mcsp_q) mcsp_tend_q(k) = -1 * runtime_opt.mcsp_q_coeff * std::sin(2 * PC::Pi * (pdepth_mid_k / pdepth_total));
            if (do_mcsp_u) mcsp_tend_u(k) = runtime_opt.mcsp_u_coeff * (std::cos(PC::Pi * (pdepth_mid_k / pdepth_total)));
            if (do_mcsp_v) mcsp_tend_v(k) = runtime_opt.mcsp_v_coeff * (std::cos(PC::Pi * (pdepth_mid_k / pdepth_total)));

            // scale the vertical structure by the ZM heating/drying tendencies
            if (do_mcsp_t) mcsp_tend_s(k) = zm_avg_tend_s * mcsp_tend_s(k);
            if (do_mcsp_q) mcsp_tend_q(k) = zm_avg_tend_q * mcsp_tend_q(k);

            // integrate the DSE/qv tendencies for energy/mass fixer
            if (do_mcsp_t) avg_s += mcsp_tend_s(k) * state_pdel(k) / pdel_sum;
            if (do_mcsp_q) avg_q += mcsp_tend_q(k) * state_pdel(k) / pdel_sum;

            // integrate the change in kinetic energy (KE) for energy fixer
            if (do_mcsp_u || do_mcsp_v) {
              const Real tend_k = (2 * mcsp_tend_u(k) * ztodt * state_u(k) + mcsp_tend_u(k) * mcsp_tend_u(k) * ztodt * ztodt
                                  + 2 * mcsp_tend_v(k) * ztodt * state_v(k) + mcsp_tend_v(k) * mcsp_tend_v(k) * ztodt * ztodt) / 2 / ztodt;
              avg_k += tend_k * state_pdel(k) / pdel_sum;
            }
          },
          mcsp_avg_tend_s, mcsp_avg_tend_q, mcsp_avg_tend_k);
        team.team_barrier();
      } // shear threshold
    } // zm_avg_tend_s > 0
  } // zm_depth >= MCSP_conv_depth_min

  //----------------------------------------------------------------------------
  // calculate final output tendencies

  mcsp_freq = 0;
  bool any_tend = false;

  // Calculate final tendencies and check frequency in a single parallel_reduce
  Kokkos::parallel_reduce(Kokkos::TeamVectorRange(team, jctop, pver),
    [&] (const Int& k, bool& local_freq) {
      // subtract mass weighted average tendencies for energy/mass conservation
      mcsp_dt_out(k) = mcsp_tend_s(k) - mcsp_avg_tend_s;
      mcsp_dq_out(k) = mcsp_tend_q(k) - mcsp_avg_tend_q;
      mcsp_du_out(k) = mcsp_tend_u(k);
      mcsp_dv_out(k) = mcsp_tend_v(k);

      // make sure kinetic energy correction is added to DSE tendency
      // to conserve total energy whenever momentum tendencies are calculated
      if (do_mcsp_u || do_mcsp_v) {
        mcsp_dt_out(k) = mcsp_dt_out(k) - mcsp_avg_tend_k;
      }

      // update output tendencies
      if (do_mcsp_t) ptend_s(k) = ptend_s(k) + mcsp_dt_out(k);
      if (do_mcsp_q) ptend_q(k) = ptend_q(k) + mcsp_dq_out(k);
      if (do_mcsp_u) ptend_u(k) = ptend_u(k) + mcsp_du_out(k);
      if (do_mcsp_v) ptend_v(k) = ptend_v(k) + mcsp_dv_out(k);

      // adjust units for diagnostic outputs
      if (do_mcsp_t) mcsp_dt_out(k) = mcsp_dt_out(k) / PC::Cpair.value;

      // update frequency if MCSP contributes any tendency in the column
      if (std::abs(mcsp_tend_s(k)) > 0 || std::abs(mcsp_tend_q(k)) > 0 ||
          std::abs(mcsp_tend_u(k)) > 0 || std::abs(mcsp_tend_v(k)) > 0) {
        local_freq = true;
      }
    },
    Kokkos::LOr<bool>(any_tend));
  team.team_barrier();

  if (any_tend) mcsp_freq = 1;

  workspace.template release_many_contiguous<4>(
    {&mcsp_tend_s, &mcsp_tend_q, &mcsp_tend_u, &mcsp_tend_v});
}

} // namespace zm
} // namespace scream

#endif
