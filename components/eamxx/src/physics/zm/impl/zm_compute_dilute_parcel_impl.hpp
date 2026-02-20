#ifndef ZM_COMPUTE_DILUTE_PARCEL_IMPL_HPP
#define ZM_COMPUTE_DILUTE_PARCEL_IMPL_HPP

#include "zm_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace zm {

/*
 * Implementation of zm compute_dilute_parcel. Clients should NOT
 * #include this file, but include zm_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::compute_dilute_parcel(
  // Inputs
  const MemberType& team,
  const Workspace& workspace,
  const ZmRuntimeOpt& runtime_opt,
  const Int& pver, // number of mid-point vertical levels
  const Int& num_msg, // number of missing moisture levels at the top of model
  const Int& klaunch, // index of parcel launch level based on max MSE
  const uview_1d<const Real>& pmid, // ambient env pressure at cell center
  const uview_1d<const Real>& temperature, // ambient env temperature at cell center
  const uview_1d<const Real>& sp_humidity, // ambient env specific humidity at cell center
  const Real& tpert, // PBL temperature perturbation
  const Int& pblt, // index of pbl depth
  // Inputs/Outputs
  const uview_1d<Real>& parcel_temp, // Parcel temperature
  const uview_1d<Real>& parcel_vtemp, // Parcel virtual temperature
  const uview_1d<Real>& parcel_qsat, // Parcel water vapour (sat value above lcl)
  Real& lcl_pmid, // lifting condensation level (LCL) pressure
  Real& lcl_temperature, // lifting condensation level (LCL) temperature
  Int& lcl_klev) // lifting condensation level (LCL) vertical index
{
  // Allocate temporary arrays
  uview_1d<Real> tmix, qtmix, qsmix, smix, xsh2o, ds_xsh2o, ds_freeze;
  workspace.template take_many_contiguous_unsafe<7>(
    {"tmix", "qtmix", "qsmix", "smix", "xsh2o", "ds_xsh2o", "ds_freeze"},
    {&tmix, &qtmix, &qsmix, &smix, &xsh2o, &ds_xsh2o, &ds_freeze});

  //----------------------------------------------------------------------------
  // The original ZM scheme only treated PBL-rooted convection. A PBL temperature
  // perturbation (tpert) was then used to increase the parcel temperatue at launch
  // level, which is in PBL. The dcape_ull or ull trigger enables ZM scheme to treat
  // elevated convection with launch level above PBL. If parcel launch level is
  // above PBL top, tempeature perturbation in PBL should not be able to influence
  // it. In this situation, the temporary variable tpert_loc is reset to zero.
  Real tpert_loc = tpert;
  if (runtime_opt.tpert_fix && klaunch < pblt) {
    tpert_loc = 0.0;
  }

  // Initialize arrays
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, pver), [&] (const Int& k) {
    tmix(k) = 0.0;
    qtmix(k) = 0.0;
    qsmix(k) = 0.0;
    smix(k) = 0.0;
    xsh2o(k) = 0.0;
    ds_xsh2o(k) = 0.0;
    ds_freeze(k) = 0.0;
  });
  team.team_barrier();

  // Due to sequential dependencies in both major loops, none of this can be parallelized
  Kokkos::single(Kokkos::PerTeam(team), [&] {
    // Scalar variables for parcel properties
    Real mp = 0.0;
    Real qtp = 0.0;
    Real sp = 0.0;
    Real sp0 = 0.0;
    Real qtp0 = 0.0;
    Real mp0 = 0.0;

    // Entrainment loop (from launch level upward)
    for (Int k = pver - 1; k >= num_msg; --k) {
      if (k == klaunch) {
        // Initialize values at launch level
        mp0 = 1.0; // initial relative mass - value of 1.0 does not change for undilute (dmpdp=0)
        qtp0 = sp_humidity(k); // initial total water - assuming subsaturated
        sp0 = entropy(temperature(k), pmid(k), qtp0);
        smix(k) = sp0;
        qtmix(k) = qtp0;
        const Real tfguess = temperature(k);
        invert_entropy(team, smix(k), pmid(k), qtmix(k), tfguess, tmix(k), qsmix(k));

      } else if (k < klaunch) {
        // Set environmental values for this level
        const Real dp = pmid(k) - pmid(k + 1);
        const Real qtenv = ZMC::half * (sp_humidity(k) + sp_humidity(k + 1));
        const Real tenv = ZMC::half * (temperature(k) + temperature(k + 1));
        const Real penv = ZMC::half * (pmid(k) + pmid(k + 1));
        const Real senv = entropy(tenv, penv, qtenv);

        // Determine fractional entrainment rate 1/pa given value 1/m
        const Real dpdz = -(penv * PC::gravit.value) / (PC::Rair.value * tenv); // [mb/m]
        const Real dzdp = 1 / dpdz; // [m/mb]
        const Real dmpdp = runtime_opt.dmpdz * dzdp; // Fractional entrainment [1/mb]

        // sum entrainment to current level - entrain q,s out of intervening dp layers,
        // assuming linear variation (i.e. entrain the mean of the 2 stored values)
        sp  -= dmpdp * dp * senv;
        qtp -= dmpdp * dp * qtenv;
        mp  -= dmpdp * dp;

        // Entrain s and qt to next level
        smix(k) = (sp0 + sp) / (mp0 + mp);
        qtmix(k) = (qtp0 + qtp) / (mp0 + mp);

        // Invert entropy from s and q to determine T and saturation-capped q of mixture
        // temperature(i,k) used as a first guess so that it converges faster
        Real tfguess = tmix(k + 1);
        invert_entropy(team, smix(k), pmid(k), qtmix(k), tfguess, tmix(k), qsmix(k));

        // determine if we are at the LCL if this is first level where qsmix<=qtmix on ascending
        if (qsmix(k) <= qtmix(k) && qsmix(k + 1) > qtmix(k + 1)) {
          lcl_klev = k;
          const Real qxsk = qtmix(k) - qsmix(k);
          const Real qxskp1 = qtmix(k + 1) - qsmix(k + 1);
          const Real dqxsdp = (qxsk - qxskp1) / dp;
          lcl_pmid = pmid(k + 1) - qxskp1 / dqxsdp;
          const Real dsdp = (smix(k) - smix(k + 1)) / dp;
          const Real dqtdp = (qtmix(k) - qtmix(k + 1)) / dp;
          const Real slcl = smix(k + 1) + dsdp * (lcl_pmid - pmid(k + 1));
          const Real qtlcl = qtmix(k + 1) + dqtdp * (lcl_pmid - pmid(k + 1));
          Real qslcl;
          tfguess = tmix(k);
          invert_entropy(team, slcl, lcl_pmid, qtlcl, tfguess, lcl_temperature, qslcl);
        }
      }
    }

    //----------------------------------------------------------------------------
    // We now have a profile of entropy and total water of the entraining parcel
    // Varying with height from the launch level klaunch parcel=environment. To the
    // top allowed level for the existence of convection. If we stop now it will
    // provide some estimate of buoyancy without the effects of freezing/condensation.
    //
    // Instead, we will adjust these values such that the water held in vapor is
    // <=qsmix. We assume that the cloud holds a certain amount of condensate (lwmax)
    // and the rest is rained out (xsh2o). This provides latent heating to the
    // mixed parcel and so this has to be added back to it.

    // Precipitation/freezing loop - iterate twice for accuracy
    for (Int k = pver - 1; k >= num_msg; --k) {
      if (k == klaunch) {
        // Initialize values at launch level - assume no liquid water
        parcel_temp(k) = tmix(k);
        parcel_qsat(k) = sp_humidity(k);
        parcel_vtemp(k) = (parcel_temp(k) + runtime_opt.tpert_fac * tpert_loc) *
          (1 + PC::ZVIR * parcel_qsat(k)) / (1 + parcel_qsat(k));

      } else if (k < klaunch) {
        // Iterate nit_lheat times for s,qt changes
        for (Int ii = 0; ii < ZMC::nit_lheat; ++ii) {
          // Rain (xsh2o) is excess condensate, bar lwmax (accumulated loss from qtmix)
          xsh2o(k) = ekat::impl::max(0.0, qtmix(k) - qsmix(k) - ZMC::lwmax);

          // Contribution to ds from precip loss of condensate (accumulated change from smix)
          ds_xsh2o(k) = ds_xsh2o(k + 1) - PC::CpLiq.value * std::log(tmix(k) / PC::Tmelt.value) *
            ekat::impl::max(0.0, xsh2o(k) - xsh2o(k + 1));

          // Calculate entropy of freezing => ( latice x amount of water involved ) / T

          // One off freezing of condensate
          if (tmix(k) <= PC::Tmelt.value && ds_freeze(k + 1) == 0.0) {
            // entropy change from latent heat
            ds_freeze(k) = (PC::LatIce.value / tmix(k)) *
              ekat::impl::max(0.0, qtmix(k) - qsmix(k) - xsh2o(k));
          }

          if (tmix(k) <= PC::Tmelt.value && ds_freeze(k + 1) != 0.0) {
            // Continual freezing of additional condensate
            ds_freeze(k) = ds_freeze(k + 1) + (PC::LatIce.value / tmix(k)) *
              ekat::impl::max(0.0, qsmix(k + 1) - qsmix(k));
          }

          // Adjust entropy accordingly to sum of ds
          const Real new_s = smix(k) + ds_xsh2o(k) + ds_freeze(k);

          // Adjust liquid water accordingly to xsh2o
          const Real new_q = qtmix(k) - xsh2o(k);

          // Invert entropy to get updated Tmix and qsmix of parcel
          Real tfguess = tmix(k);
          invert_entropy(team, new_s, pmid(k), new_q, tfguess, tmix(k), qsmix(k));
        }

        // Parcel temp is temp of mixture
        // parcel virtual temp should be density temp with new_q total water
        parcel_temp(k) = tmix(k);

        // parcel_vtemp=tprho in the presence of condensate (i.e. when new_q > qsmix)
        const Real new_q = qtmix(k) - xsh2o(k);
        if (new_q > qsmix(k)) { // super-saturated so condensate present - reduces buoyancy
          parcel_qsat(k) = qsmix(k);
        } else { // just saturated/sub-saturated - no condensate virtual effects
          parcel_qsat(k) = new_q;
        }

        parcel_vtemp(k) = (parcel_temp(k) + runtime_opt.tpert_fac * tpert_loc) *
          (1 + PC::ZVIR * parcel_qsat(k)) / (1 + new_q);
      }
    }
  });

  workspace.template release_many_contiguous<7>(
    {&tmix, &qtmix, &qsmix, &smix, &xsh2o, &ds_xsh2o, &ds_freeze});
}

} // namespace zm
} // namespace scream

#endif
