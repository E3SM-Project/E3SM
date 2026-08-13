#ifndef ZM_COMPUTE_CAPE_FROM_PARCEL_IMPL_HPP
#define ZM_COMPUTE_CAPE_FROM_PARCEL_IMPL_HPP

#include "zm_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace zm {

/*
 * Implementation of zm compute_cape_from_parcel. Clients should NOT
 * #include this file, but include zm_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::compute_cape_from_parcel(
  // Inputs
  const MemberType& team,
  const Workspace& workspace,
  const ZmRuntimeOpt& runtime_opt,
  const Int& pver, // number of mid-point vertical levels
  const Int& pverp, // number of interface vertical levels
  const Int& num_cin, // num of negative buoyancy regions that are allowed before the conv. top and CAPE calc are completed
  const Int& num_msg, // number of missing moisture levels at the top of model
  const uview_1d<const Real>& temperature, // temperature
  const uview_1d<const Real>& tv, // virtual temperature
  const uview_1d<const Real>& sp_humidity, // specific humidity
  const uview_1d<const Real>& pint, // pressure at interfaces
  const Int& msemax_klev, // index of max MSE at parcel launch level
  const Real& lcl_pmid, // lifting condensation level (LCL) pressure
  const Int& lcl_klev, // lifting condensation level (LCL) index
  // Inputs/Outputs
  const uview_1d<Real>& parcel_qsat, // parcel saturation mixing ratio
  const uview_1d<Real>& parcel_temp, // parcel temperature
  const uview_1d<Real>& parcel_vtemp, // parcel virtual temperature
  Int& eql_klev, // index of equilibrium level (i.e. cloud top)
  Real& cape) // convective available potential energy
{
  //----------------------------------------------------------------------------
  // Purpose: calculate convective available potential energy (CAPE)
  //          from parcel thermodynamic properties from compute_dilute_parcel()
  //----------------------------------------------------------------------------

  // Allocate temporary arrays from workspace:
  //   buoyancy:        parcel buoyancy at each level (size pver)
  //   cape_tmp:        provisional CAPE value for each negative buoyancy region (size num_cin, held in pver-sized slot)
  //   eql_klev_tmp_r:  provisional equilibrium level index per region, stored as Real (size num_cin, held in pver-sized slot)
  uview_1d<Real> buoyancy, cape_tmp, eql_klev_tmp_r;
  workspace.template take_many_contiguous_unsafe<3>(
    {"buoyancy", "cape_tmp", "eql_klev_tmp"},
    {&buoyancy, &cape_tmp, &eql_klev_tmp_r});

  //----------------------------------------------------------------------------
  // Initialize output scalars and temporary arrays
  Kokkos::single(Kokkos::PerTeam(team), [&]() {
    eql_klev = pver - 1;
    cape     = 0;
  });
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, num_cin), [&](Int n) {
    eql_klev_tmp_r(n) = static_cast<Real>(pver - 1);
    cape_tmp(n)       = 0;
  });
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, pver), [&](Int k) {
    buoyancy(k) = 0;
  });
  team.team_barrier();

  //----------------------------------------------------------------------------
  // Calculate buoyancy at each level from launch level to model top.
  // Levels where the parcel is not active (above launch level or LCL pressure
  // too low) are reset to environmental values.
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, num_msg, pver), [&](Int k) {
    if (k <= msemax_klev && lcl_pmid >= ZMC::lcl_pressure_threshold) {
      buoyancy(k) = parcel_vtemp(k) - tv(k) + runtime_opt.tiedke_add;
    } else {
      parcel_qsat(k)  = sp_humidity(k);
      parcel_temp(k)  = temperature(k);
      parcel_vtemp(k) = tv(k);
    }
  });
  team.team_barrier();

  //----------------------------------------------------------------------------
  // Find the convective equilibrium level for each allowed negative-buoyancy
  // region.  This loop has a sequential dependency on neg_buoyancy_cnt, so it
  // must run on a single thread.
  //
  // Index mapping (Fortran 1-based → C++ 0-based):
  //   Fortran k = num_msg+2 .. pver  →  C++ k = num_msg+1 .. pver-1
  //   buoyancy(k+1) in Fortran       →  buoyancy(k+1) in C++ (same offset)
  //   The guard k < lcl_klev ensures k+1 <= lcl_klev <= pver-1, staying in bounds.
  Kokkos::single(Kokkos::PerTeam(team), [&]() {
    Int neg_buoyancy_cnt = 0;
    for (Int k = num_msg + 1; k <= pver - 1; ++k) {
      if (k < lcl_klev && lcl_pmid >= ZMC::lcl_pressure_threshold) {
        if (buoyancy(k + 1) > 0 && buoyancy(k) <= 0) {
          neg_buoyancy_cnt = Kokkos::min(num_cin, neg_buoyancy_cnt + 1);
          // Store 0-based C++ level index into the Real array slot
          eql_klev_tmp_r(neg_buoyancy_cnt - 1) = static_cast<Real>(k);
        }
      }
    }
  });
  team.team_barrier();

  //----------------------------------------------------------------------------
  // Integrate buoyancy to obtain provisional CAPE for each negative-buoyancy
  // region.  The outer loop over n (num_cin regions) is serial; the inner
  // reduction over vertical levels is parallel.
  for (Int n = 0; n < num_cin; ++n) {
    const Int eql_klev_n = static_cast<Int>(eql_klev_tmp_r(n));
    Real local_cape = 0;
    Kokkos::parallel_reduce(Kokkos::TeamVectorRange(team, num_msg, pver),
      [&](Int k, Real& sum) {
        if (lcl_pmid >= ZMC::lcl_pressure_threshold &&
            k <= msemax_klev && k > eql_klev_n) {
          sum += PC::Rair.value * buoyancy(k) * std::log(pint(k + 1) / pint(k));
        }
      }, local_cape);
    cape_tmp(n) = local_cape;
    team.team_barrier();
  }

  //----------------------------------------------------------------------------
  // Select the maximum CAPE from all provisional values and set the
  // corresponding equilibrium level.  Apply a non-negative limiter.
  Kokkos::single(Kokkos::PerTeam(team), [&]() {
    for (Int n = 0; n < num_cin; ++n) {
      if (cape_tmp(n) > cape) {
        cape     = cape_tmp(n);
        eql_klev = static_cast<Int>(eql_klev_tmp_r(n));
      }
    }
    cape = Kokkos::max(cape, Real(0));
  });
  team.team_barrier();

  workspace.template release_many_contiguous<3>(
    {&buoyancy, &cape_tmp, &eql_klev_tmp_r});
}

} // namespace zm
} // namespace scream

#endif
