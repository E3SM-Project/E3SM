#ifndef ZM_COMPUTE_DILUTE_CAPE_IMPL_HPP
#define ZM_COMPUTE_DILUTE_CAPE_IMPL_HPP

#include "zm_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace zm {

/*
 * Implementation of zm compute_dilute_cape. Clients should NOT
 * #include this file, but include zm_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::compute_dilute_cape(
  // Inputs
  const MemberType& team,
  const Workspace& workspace,
  const ZmRuntimeOpt& runtime_opt,
  const Int& pver, // number of mid-point vertical levels
  const Int& pverp, // number of interface vertical levels
  const Int& num_cin, // num of negative buoyancy regions that are allowed before the conv. top and CAPE calc are completed
  const Int& num_msg, // index of highest level convection is allowed
  const uview_1d<const Real>& sp_humidity_in, // specific humidity
  const uview_1d<const Real>& temperature_in, // temperature
  const uview_1d<const Real>& zmid, // altitude/height at mid-levels
  const uview_1d<const Real>& pmid, // pressure at mid-levels
  const uview_1d<const Real>& pint, // pressure at interfaces
  const Int& pblt, // index of pbl top used as upper limit index of max MSE search
  const Real& tpert, // perturbation temperature by pbl processes
  const bool& calc_msemax_klev, // true for normal procedure, otherwise use prev_msemax_klev from 1st call
  const Int& prev_msemax_klev, // values of msemax_klev from previous call for dcape closure
  const bool& use_input_tq_mx, // if .true., use input values of prev_msemax_klev, q_mx, t_mx in the CAPE calculation
  // Inputs/Outputs
  const uview_1d<Real>& parcel_qsat, // parcel saturation mixing ratio
  Int& msemax_klev, // index of max MSE at parcel launch level
  Int& lcl_klev, // index of lifting condensation level (i.e. cloud bottom)
  Int& eql_klev, // index of equilibrium level (i.e. cloud top)
  Real& cape, // convective available potential energy
  Real& q_mx, // specified sp humidity to apply at level of max MSE if use_input_tq_mx=.true.
  Real& t_mx, // specified temperature to apply at level of max MSE if use_input_tq_mx=.true.)
  // Outputs
  const uview_1d<Real>& parcel_temp, // parcel temperature
  Real& lcl_temperature) // lifting condensation level (LCL) temperature
{
  //----------------------------------------------------------------------------
  // Purpose: calculate convective available potential energy (CAPE), lifting
  //          condensation level (LCL), and convective top with dilute parcel ascent
  // Method: parcel temperature based on a plume model with constant entrainment
  // Original Author: Richard Neale - September 2004
  // References:
  //   Raymond, D. J., and A. M. Blyth, 1986: A Stochastic Mixing Model for
  //     Nonprecipitating Cumulus Clouds. J. Atmos. Sci., 43, 2708–2718
  //   Raymond, D. J., and A. M. Blyth, 1992: Extension of the Stochastic Mixing
  //     Model to Cumulonimbus Clouds. J. Atmos. Sci., 49, 1968–1983
  //----------------------------------------------------------------------------
#ifdef PERGRO
  const bool pergro_active = true;
#else
  const bool pergro_active = false;
#endif

  // Allocate temporary arrays
  uview_1d<Real> sp_humidity, temperature, tv, parcel_vtemp;
  workspace.template take_many_contiguous_unsafe<4>(
    {"sp_humidity", "temperature", "tv", "parcel_vtemp"},
    {&sp_humidity, &temperature, &tv, &parcel_vtemp});

  Real lcl_pmid = 0.0;
  Real mse_max_val = 0.0;
  Int pblt_ull = 0;
  Int msemax_top_k = 0;

  //----------------------------------------------------------------------------
  // Copy the incoming temperature and specific humidity values to local arrays
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, pver), [&] (const Int& k) {
    temperature(k) = temperature_in(k);
    sp_humidity(k) = sp_humidity_in(k);
  });
  team.team_barrier();

  //----------------------------------------------------------------------------
  // initialize msemax_klev and potentially modify T/q
  if (use_input_tq_mx) {
    // note - in this case we expect:
    // (1) the incoming array prev_msemax_klev contains prev identified launching level index, and
    // (2) the arrays q_mx and t_mx contain q and T values at the old launching level
    //     at the time when the old launching level was identified.
    // Copy the old values to work arrays for calculations in the rest of this subroutine
    msemax_klev = prev_msemax_klev;
    sp_humidity(msemax_klev) = q_mx;
    temperature(msemax_klev) = t_mx;
  } else {
    msemax_klev = pver - 1;
  }

  //----------------------------------------------------------------------------
  // calculate virtual temperature
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, pver), [&] (const Int& k) {
    tv(k) = temperature(k) * (1 + PC::ZVIR * sp_humidity(k)) / (1 + sp_humidity(k));
  });
  team.team_barrier();

  //----------------------------------------------------------------------------
  // Initialize parcel properties
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, pver), [&] (const Int& k) {
    parcel_temp(k) = temperature(k);
    parcel_qsat(k) = sp_humidity(k);
    parcel_vtemp(k) = tv(k);
  });
  team.team_barrier();

  //----------------------------------------------------------------------------
  // Find new upper bound for parcel starting level - unrestricted launch level (ULL)
  if (runtime_opt.trig_ull) {
    pblt_ull = 0;
    Kokkos::single(Kokkos::PerTeam(team), [&] () {
      for (Int k = pver - 2; k >= num_msg; --k) {
        if ((pmid(k)    <= ZMC::ull_upper_launch_pressure) &&
            (pmid(k + 1) > ZMC::ull_upper_launch_pressure)) {
          pblt_ull = k;
        }
      }
    });
  }

  //----------------------------------------------------------------------------
  // Set level of max moist static energy for parcel initialization
  if (runtime_opt.trig_dcape && (!calc_msemax_klev)) {
    // Use max moist static energy level that is passed in
    msemax_klev = prev_msemax_klev;
  } else if (!use_input_tq_mx) {
    if (runtime_opt.trig_ull) {
      msemax_top_k = pblt_ull;
    } else {
      msemax_top_k = pblt;
    }
    find_mse_max(team, runtime_opt, pver, num_msg, msemax_top_k, pergro_active,
                 temperature, zmid, sp_humidity,
                 msemax_klev, mse_max_val);
  }

  //----------------------------------------------------------------------------
  // Save launching level T, q for output
  if (!use_input_tq_mx) {
    q_mx = sp_humidity(msemax_klev);
    t_mx = temperature(msemax_klev);
  }
  // save LCL values for compute_dilute_parcel()
  lcl_klev = msemax_klev;
  lcl_pmid = pmid(msemax_klev);
  lcl_temperature = temperature(msemax_klev);

  //----------------------------------------------------------------------------
  // entraining parcel calculation
  compute_dilute_parcel(team, workspace, runtime_opt, pver, num_msg, msemax_klev,
                        pmid, temperature, sp_humidity, tpert, pblt,
                        parcel_temp, parcel_vtemp, parcel_qsat,
                        lcl_pmid, lcl_temperature, lcl_klev);

  //----------------------------------------------------------------------------
  // calculate CAPE
  compute_cape_from_parcel(team, workspace, runtime_opt, pver, pverp, num_cin, num_msg,
                          temperature, tv, zmid, sp_humidity, pint,
                          msemax_klev, lcl_pmid, lcl_klev,
                          parcel_qsat, parcel_temp, parcel_vtemp,
                          eql_klev, cape);

  workspace.template release_many_contiguous<4>(
    {&sp_humidity, &temperature, &tv, &parcel_vtemp});
}

} // namespace zm
} // namespace scream

#endif
