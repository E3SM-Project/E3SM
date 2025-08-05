#include "tms_test_data.hpp"

#include "share/util/eamxx_deep_copy.hpp"

#include <ekat_team_policy_utils.hpp>

#include <random>

using scream::Real;

namespace scream {
namespace tms {

//
// _f function definitions. These expect data in C layout
//
void compute_tms_f(int ncols, int nlevs,
                   Real *u_wind, Real *v_wind, Real *t_mid, Real *p_mid, Real *exner, Real *z_mid,
                   Real *sgh, Real *landfrac, Real *ksrf, Real *taux, Real *tauy)
{
  using TMSFunc  = Functions<Real, DefaultDevice>;

  using Scalar     = typename TMSFunc::Scalar;
  using Spack      = ekat::Pack<Scalar,SCREAM_PACK_SIZE>;
  using view_1d    = typename TMSFunc::view_1d<Scalar>;
  using view_2d    = typename TMSFunc::view_2d<Spack>;
  using view_2d_s  = typename TMSFunc::view_2d<Scalar>;
  using view_3d    = typename TMSFunc::view_3d<Spack>;
  using ExeSpace   = typename TMSFunc::KT::ExeSpace;
  using MemberType = typename TMSFunc::KT::MemberType;
  using TPF        = ekat::TeamPolicyFactory<ExeSpace>;

  // Initialize Kokkos views, sync to device
  std::vector<view_1d> temp_d_1d(2);
  std::vector<view_2d> temp_d_2d(6);
  ScreamDeepCopy::copy_to_device({sgh, landfrac}, ncols, temp_d_1d);
  ekat::host_to_device({u_wind, v_wind, t_mid, p_mid, exner, z_mid},
                        ncols, nlevs, temp_d_2d, true);

  view_1d
    sgh_d     (temp_d_1d[0]),
    landfrac_d(temp_d_1d[1]),
    ksrf_d    ("ksrf_d", ncols),
    taux_d    ("taux_d", ncols),
    tauy_d    ("tauy_d", ncols);

  view_2d
    u_wind_d(temp_d_2d[0]),
    v_wind_d(temp_d_2d[1]),
    t_mid_d (temp_d_2d[2]),
    p_mid_d (temp_d_2d[3]),
    exner_d (temp_d_2d[4]),
    z_mid_d (temp_d_2d[5]);

  // calculate_tms treats u/v_wind and taux/y as multiple component arrays.
  const auto nlev_packs = ekat::npack<Spack>(nlevs);
  view_3d   horiz_wind_d("horiz_wind_d", ncols, 2, nlev_packs);
  view_2d_s tau_d("tau_d", ncols, 2);
  const auto policy = TPF::get_default_team_policy(ncols, nlev_packs);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const int i = team.league_rank();
    tau_d(i, 0) = taux_d(i);
    tau_d(i, 1) = tauy_d(i);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev_packs), [&] (const Int& k) {
      horiz_wind_d(i,0,k) = u_wind_d(i,k);
      horiz_wind_d(i,1,k) = v_wind_d(i,k);
    });
  });

  // C++ compute_tms function implementation
  TMSFunc::compute_tms(ncols, nlevs,
                       ekat::scalarize(horiz_wind_d),
                       ekat::scalarize(t_mid_d),
                       ekat::scalarize(p_mid_d),
                       ekat::scalarize(exner_d),
                       ekat::scalarize(z_mid_d),
                       sgh_d, landfrac_d, ksrf_d, tau_d);

  // Transfer data back to individual arrays (only for output variables)
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const int i = team.league_rank();
    taux_d(i) = tau_d(i, 0);
    tauy_d(i) = tau_d(i, 1);
  });

  // Sync back to host
  std::vector<view_1d> output_data = {ksrf_d, taux_d, tauy_d};
  ScreamDeepCopy::copy_to_host({ksrf, taux, tauy}, ncols, output_data);
}

void compute_tms(ComputeTMSData& d)
{
  compute_tms_f(d.ncols, d.nlevs, d.u_wind, d.v_wind, d.t_mid, d.p_mid, d.exner,
                d.z_mid, d.sgh, d.landfrac, d.ksrf, d.taux, d.tauy);
}

} // namespace tms
} // namespace scream
