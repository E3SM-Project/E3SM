#include <catch2/catch.hpp>

#include <Kokkos_Core.hpp>
#include <ekat_team_policy_utils.hpp>

#include <cmath>
#include <type_traits>

#include "physics/mam/mam_coupling.hpp"

namespace {

using Real = scream::Real;
using HostSpace = Kokkos::HostSpace;
using ExecSpace = scream::DefaultDevice::execution_space;
using view_2d = scream::mam_coupling::view_2d;
using MemberType = scream::mam_coupling::Team;

bool nearly_equal(const Real a, const Real b,
                  const Real rtol, const Real atol = 1e-14) {
  return std::abs(a - b) <= atol + rtol * std::abs(b);
}

}  // namespace

TEST_CASE("mam_layer_geometry_uses_wet_atmosphere", "[mam4][geometry][kokkos]") {
  namespace mc = scream::mam_coupling;
  using C = scream::physics::Constants<Real>;
  const Real qv = GENERATE(Real(0), Real(0.005), Real(0.01), Real(0.02));
  const bool optics_dry_mass = GENERATE(false, true);
  constexpr int nlev = mam4::nlev;
  constexpr Real dp = 1000, p = 80000, temperature = 280, omega = -0.1;
  const Real dry_dp = optics_dry_mass ? dp * (1 - qv) : dp;
  INFO("qv=" << qv << ", optics_dry_mass=" << optics_dry_mass);

  view_2d wet_qv("wet_qv", 1, nlev), wet_dp("wet_dp", 1, nlev);
  view_2d dry_qv("dry_qv", 1, nlev), aerosol_dp("aerosol_dp", 1, nlev);
  view_2d pmid("pmid", 1, nlev), temp("temp", 1, nlev);
  view_2d pressure_velocity("omega", 1, nlev), dz_before("dz_before", 1, nlev);
  Kokkos::deep_copy(wet_qv, qv);
  Kokkos::deep_copy(wet_dp, dp);
  Kokkos::deep_copy(dry_qv, qv / (1 - qv));
  Kokkos::deep_copy(aerosol_dp, dry_dp);
  Kokkos::deep_copy(pmid, p);
  Kokkos::deep_copy(temp, temperature);
  Kokkos::deep_copy(pressure_velocity, omega);

  mc::WetAtmosphere wet_atm;
  wet_atm.qv = wet_qv;
  wet_atm.p_del = wet_dp;
  wet_atm.dz = view_2d("dz", 1, nlev);
  wet_atm.z_iface = view_2d("z_iface", 1, nlev + 1);
  wet_atm.z_mid = view_2d("z_mid", 1, nlev);
  mc::DryAtmosphere dry_atm;
  dry_atm.z_surf = 0;
  dry_atm.qv = dry_qv;
  dry_atm.p_del = aerosol_dp;
  dry_atm.p_mid = pmid;
  dry_atm.T_mid = temp;
  dry_atm.omega = pressure_velocity;
  dry_atm.w_updraft = view_2d("w_updraft", 1, nlev);

  using TPF = ekat::TeamPolicyFactory<ExecSpace>;
  const auto policy = TPF::get_thread_range_parallel_scan_team_policy(1, nlev);
  Kokkos::parallel_for("mam_layer_geometry", policy,
    KOKKOS_LAMBDA(const MemberType& team) {
      mc::compute_vertical_layer_heights(team, wet_atm, dry_atm, 0);
      team.team_barrier();
      Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev),
        [&](const int k) { dz_before(0, k) = wet_atm.dz(0, k); });
      team.team_barrier();
      mc::compute_updraft_velocities(team, wet_atm, dry_atm, 0);
    });
  const auto dz = Kokkos::create_mirror_view_and_copy(HostSpace(), wet_atm.dz);
  const auto before = Kokkos::create_mirror_view_and_copy(HostSpace(), dz_before);
  const auto zi = Kokkos::create_mirror_view_and_copy(HostSpace(), wet_atm.z_iface);
  const auto zm = Kokkos::create_mirror_view_and_copy(HostSpace(), wet_atm.z_mid);
  const auto w = Kokkos::create_mirror_view_and_copy(HostSpace(), dry_atm.w_updraft);
  const auto mass = Kokkos::create_mirror_view_and_copy(HostSpace(), aerosol_dp);
  const Real expected_dz =
      temperature * (C::RD.value * (1 - qv) + C::RV.value * qv) * dp /
      (C::gravit.value * p);
  constexpr Real rtol = std::is_same<Real, float>::value ? 1e-5 : 1e-8;
  REQUIRE(zi(0, nlev) == 0);
  for (int k = 0; k < nlev; ++k) {
    REQUIRE(nearly_equal(dz(0, k), expected_dz, rtol));
    REQUIRE(dz(0, k) == before(0, k));
    REQUIRE(nearly_equal(zi(0, k), (nlev - k) * expected_dz, rtol));
    REQUIRE(nearly_equal(zm(0, k), (nlev - k - Real(0.5)) * expected_dz, rtol));
    REQUIRE(nearly_equal(w(0, k), -omega * expected_dz / dp, rtol));
    REQUIRE(mass(0, k) == dry_dp);
  }
}
