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
  const uview_1d<const Real>& zmid, // height/altitude at mid-levels
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
  // Allocate temporary arrays
  uview_1d<Real> buoyancy, cape_tmp, eql_klev_tmpf;
  workspace.template take_many_contiguous_unsafe<3>(
    {"buoyancy", "cape_tmp", "eql_klev_tmp"},
    {&buoyancy, &cape_tmp, &eql_klev_tmpf});

  constexpr Real Rair = PC::Rair.value;

  uview_1d<Int> eql_klev_tmp(reinterpret_cast<Int*>(eql_klev_tmpf.data()), num_cin);

  // Initialize variables
  eql_klev = pver - 1;
  cape = 0.0;

  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, num_cin), [&] (const Int& n) {
    eql_klev_tmp(n) = pver - 1;
    cape_tmp(n) = 0.0;
  });
  team.team_barrier();

  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, pver), [&] (const Int& k) {
    buoyancy(k) = 0.0;
  });
  team.team_barrier();

  // Calculate buoyancy
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, num_msg, pver), [&] (const Int& k) {
    // Define buoyancy from launch level to equilibrium level
    if (k <= msemax_klev && lcl_pmid >= ZMC::lcl_pressure_threshold) {
      buoyancy(k) = parcel_vtemp(k) - tv(k) + runtime_opt.tiedke_add;
    } else {
      parcel_qsat(k) = sp_humidity(k);
      parcel_temp(k) = temperature(k);
      parcel_vtemp(k) = tv(k);
    }
  });
  team.team_barrier();

  // Find convective equilibrium level accounting for negative buoyancy levels
  Kokkos::single(Kokkos::PerTeam(team), [&] {
    Int neg_buoyancy_cnt = 0;
    for (Int k = num_msg + 1; k < pver; ++k) {
      if (k < lcl_klev && lcl_pmid >= ZMC::lcl_pressure_threshold) {
        if (buoyancy(k + 1) > 0.0 && buoyancy(k) <= 0.0) {
          neg_buoyancy_cnt = ekat::impl::min(num_cin, neg_buoyancy_cnt + 1);
          eql_klev_tmp(neg_buoyancy_cnt - 1) = k;
        }
      }
    }
  });

  // Integrate buoyancy to obtain possible CAPE values
  for (Int n = 0; n < num_cin; ++n) {
    Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, num_msg, pver),
      [&] (const Int& k, Real& cape_n) {
        if (lcl_pmid >= ZMC::lcl_pressure_threshold &&
            k <= msemax_klev && k > eql_klev_tmp(n)) {
          cape_n += Rair * buoyancy(k) * std::log(pint(k + 1) / pint(k));
        }
      }, cape_tmp(n));
  }
  team.team_barrier();

  Kokkos::single(Kokkos::PerTeam(team), [&] {
    // Find maximum cape from all possible tentative CAPE values
    for (Int n = 0; n < num_cin; ++n) {
      if (cape_tmp(n) > cape) {
        cape = cape_tmp(n);
        eql_klev = eql_klev_tmp(n);
      }
    }

    // Apply limiter to ensure CAPE is positive
    cape = ekat::impl::max(cape, 0.0);
  });

  workspace.template release_many_contiguous<3>(
    {&buoyancy, &cape_tmp, &eql_klev_tmpf});
}

} // namespace zm
} // namespace scream

#endif
