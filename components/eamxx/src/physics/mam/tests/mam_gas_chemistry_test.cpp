#include <catch2/catch.hpp>

#include <mam4xx/gas_phase_chemistry.hpp>

#include <cmath>
#include <limits>
#include <string>

TEST_CASE("mam_gas_chemistry_strided_view", "[mam4][gas_chem][kokkos]")
{
  using mam4::Real;
  using namespace mam4::gas_chemistry;

  // Check the selected production macro against the independent CMake option.
#if defined(MAM4XX_GAS_CHEMISTRY_SOLVER_LEGACY)
  REQUIRE(std::string(EAMXX_TEST_MAM_GAS_CHEMISTRY_SOLVER) == "legacy");
#else
  REQUIRE(std::string(EAMXX_TEST_MAM_GAS_CHEMISTRY_SOLVER) == "direct");
#endif

  Real removal_rate = 0;
  SECTION("zero forcing and rates preserve the state") {}
  SECTION("first order removal follows backward Euler")
  {
    removal_rate = Real(1e-5); // s^-1
  }

  constexpr int ncol  = 2;
  constexpr int nlev  = 3;
  constexpr int steps = 3;
  constexpr Real dt   = 60; // s
  using Device        = typename mam4::DeviceType::Device;
  // Species entries of each column/level subview have stride ncol*nlev.
  // This exercises the view interface used by EAMxx rather than a local array.
  using View = Kokkos::View<Real ***, Kokkos::LayoutLeft, Device>;
  View state("gas_chemistry_state", ncol, nlev, gas_pcnst);
  View history("gas_chemistry_history", (steps + 1) * ncol, nlev, gas_pcnst);
  using ResultView = Kokkos::View<ImpSolResult ***, Device>;
  ResultView results("gas_chemistry_results", steps, ncol, nlev);

  Kokkos::parallel_for(
      "gas_chemistry_strided_view", Kokkos::RangePolicy<mam4::ExecutionSpace>(0, ncol * nlev),
      KOKKOS_LAMBDA(const int i) {
        const int icol                      = i / nlev;
        const int lev                       = i % nlev;
        const auto qq                       = Kokkos::subview(state, icol, lev, Kokkos::ALL());
        Real heterogeneous_rates[gas_pcnst] = {};
        for (int j = 0; j < gas_pcnst; ++j) {
          // Synthetic values retain each species' existing work-array basis.
          qq[j]                  = Real((icol + 1) * (lev + 1) * (j + 1)) * Real(1e-9);
          heterogeneous_rates[j] = j == o3_idx ? Real(0) : removal_rate;
          history(icol, lev, j)  = qq[j];
        }
        Real photo_rates[mam4::mo_photo::phtcnt] = {};
        Real forcing[extcnt]                     = {};
        // Prescribed number densities [molecules cm^-3]. Zero oxidants turn
        // off chemical reactions, leaving only the specified removal.
        Real invariants[nfs] = {Real(2.5e19), Real(1.95e19), Real(5.25e18), Real(1e17),
                                Real(0),      Real(0),       Real(0),       Real(0)};

        for (int step = 1; step <= steps; ++step) {
          results(step - 1, icol, lev) = mam4::microphysics::gas_phase_chemistry(
              Real(300), dt, photo_rates, forcing, invariants, heterogeneous_rates, qq);
          for (int j = 0; j < gas_pcnst; ++j) {
            history(step * ncol + icol, lev, j) = qq[j];
          }
        }
      });
  Kokkos::fence();

  const auto host         = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), history);
  const auto host_results = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), results);
  const Real tolerance    = Real(64) * std::numeric_limits<Real>::epsilon();
  for (int icol = 0; icol < ncol; ++icol) {
    for (int lev = 0; lev < nlev; ++lev) {
      for (int step = 0; step < steps; ++step) {
        CAPTURE(icol, lev, step, removal_rate);
        const auto result = host_results(step, icol, lev);
        REQUIRE(result.success());
        REQUIRE(result.outcome == ImpSolOutcome::Converged);
        REQUIRE(result.failed_attempts == 0);
        REQUIRE(result.cut_count == 0);
        REQUIRE(result.accepted_steps == 1);
        REQUIRE(result.requested_interval == dt);
        REQUIRE(result.accepted_interval == dt);
      }
      for (int j = 0; j < gas_pcnst; ++j) {
        Real expected = host(icol, lev, j);
        for (int step = 1; step <= steps; ++step) {
          CAPTURE(icol, lev, j, step, removal_rate);
          const Real actual = host(step * ncol + icol, lev, j);
          REQUIRE(std::isfinite(actual));
          REQUIRE(actual >= 0);
          if (j == o3_idx) {
            REQUIRE(actual == host(icol, lev, j));
          } else if (removal_rate == 0) {
            // Legacy Newton arithmetic can introduce roundoff even for dq/dt=0.
            REQUIRE(std::abs(actual - expected) <= tolerance * std::abs(expected));
          } else {
            // Independent scalar reference for dq/dt = -removal_rate*q.
            expected /= Real(1) + dt * removal_rate;
            REQUIRE(std::abs(actual - expected) <= tolerance * std::abs(expected));
            REQUIRE(actual < host((step - 1) * ncol + icol, lev, j));
          }
        }
      }
    }
  }
}
