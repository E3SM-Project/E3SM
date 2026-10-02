#include "rayleigh_friction_test_data.hpp"

#include <ekat_pack_kokkos.hpp>
#include <ekat_team_policy_utils.hpp>

#include <vector>

using scream::Real;

namespace scream {
namespace rayleigh_friction {

using RFFunc     = Functions<Real, DefaultDevice>;
using view_1d    = typename RFFunc::view_1d<Real>;
using view_2d    = typename RFFunc::view_2d<Real>;
using view_3d    = typename RFFunc::view_3d<Real>;
using ExeSpace   = typename RFFunc::KT::ExeSpace;
using MemberType = typename RFFunc::KT::MemberType;
using TPF        = ekat::TeamPolicyFactory<ExeSpace>;

void rayleigh_friction_init(RayleighFrictionInitData& d)
{
  view_1d otau_d("otau_d", d.nlevs);

  RFFunc::rayleigh_friction_init(d.nlevs, d.rayk0, d.raykrange, d.raytau0, otau_d);

  // Sync back to host
  std::vector<view_1d> output_data = {otau_d};
  ekat::device_to_host({d.otau}, d.nlevs, output_data);
}

void rayleigh_friction_tend(RayleighFrictionTendData& d)
{
  const int ncols = d.ncols;
  const int nlevs = d.nlevs;

  // Initialize Kokkos views, sync to device
  std::vector<view_1d> temp_d_1d(1);
  std::vector<view_2d> temp_d_2d(3);
  ekat::host_to_device({d.otau}, nlevs, temp_d_1d);
  ekat::host_to_device({d.u_wind, d.v_wind, d.t_mid}, ncols, nlevs, temp_d_2d);

  view_1d otau_d(temp_d_1d[0]);
  view_2d
    u_wind_d(temp_d_2d[0]),
    v_wind_d(temp_d_2d[1]),
    t_mid_d (temp_d_2d[2]);

  // rayleigh_friction_tend treats u/v_wind as a multiple component array
  view_3d horiz_winds_d("horiz_winds_d", ncols, 2, nlevs);
  const auto policy = TPF::get_default_team_policy(ncols, nlevs);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const int i = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlevs), [&] (const int& k) {
      horiz_winds_d(i,0,k) = u_wind_d(i,k);
      horiz_winds_d(i,1,k) = v_wind_d(i,k);
    });
  });

  RFFunc::rayleigh_friction_tend(ncols, nlevs, d.dt, otau_d, horiz_winds_d, t_mid_d);

  // Transfer data back to individual arrays
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const int i = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlevs), [&] (const int& k) {
      u_wind_d(i,k) = horiz_winds_d(i,0,k);
      v_wind_d(i,k) = horiz_winds_d(i,1,k);
    });
  });

  // Sync back to host
  std::vector<view_2d> inout_data = {u_wind_d, v_wind_d, t_mid_d};
  ekat::device_to_host({d.u_wind, d.v_wind, d.t_mid}, ncols, nlevs, inout_data);
}

} // namespace rayleigh_friction
} // namespace scream
