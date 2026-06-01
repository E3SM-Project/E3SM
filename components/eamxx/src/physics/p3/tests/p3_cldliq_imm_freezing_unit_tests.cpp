#include "catch2/catch.hpp"

#include "p3_functions.hpp"
#include "p3_test_data.hpp"
#include "p3_unit_tests_common.hpp"

#include "share/core/eamxx_types.hpp"

#include <cmath>
#include <array>
#include <algorithm>
#include <limits>

namespace scream {
namespace p3 {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestCldliqImmersionFreezing : public UnitWrap::UnitTest<D>::Base {

  // These identities involve exp, tgamma, and ratio formation, so use a
  // modestly loose epsilon-scaled tolerance that remains robust in SP builds.
  static constexpr Scalar identity_tol =
    1000 * std::numeric_limits<Scalar>::epsilon();
  static constexpr Scalar zero_tol = Scalar(1e-30);

  static_assert(max_pack_size >= 8,
                "This test assumes at least 8 scenario slots.");
  static_assert(max_pack_size % Pack::n == 0,
                "max_pack_size must be divisible by Pack::n.");

  struct ImmFreezeResult {
    Scalar mass;
    Scalar number;
  };

  struct ImmFreezeLane {
    Scalar T_atm, lamc, mu_c, cdist1, qc_incld, inv_qc_relvar;
    bool context;
    Scalar mass_in, number_in;
    Scalar mass_out, number_out;
  };

  static bool rel_close(const Scalar got,
                        const Scalar expected,
                        const Scalar rtol = identity_tol,
                        const Scalar atol = zero_tol)
  {
    const auto scale = std::max(std::abs(got), std::abs(expected));
    return std::abs(got - expected) <= atol + rtol * scale;
  }

  static Scalar cube_host(const Scalar x) { return x * x * x; }

  static Scalar expected_mass_number_ratio(const Scalar mu_c,
                                           const Scalar lamc)
  {
    return (C::CONS6 / C::CONS5)
         * (mu_c + 4)
         * (mu_c + 5)
         * (mu_c + 6)
         / cube_host(lamc);
  }

  ImmFreezeResult run_case(
    const Scalar T_atm,
    const Scalar lamc,
    const Scalar mu_c,
    const Scalar cdist1,
    const Scalar qc_incld,
    const Scalar inv_qc_relvar,
    const Scalar exponent,
    const bool context_value,
    const Scalar mass_seed = 0,
    const Scalar number_seed = 0) const
  {
    auto lanes = make_lanes();

    lanes[0].T_atm = T_atm;
    lanes[0].lamc = lamc;
    lanes[0].mu_c = mu_c;
    lanes[0].cdist1 = cdist1;
    lanes[0].qc_incld = qc_incld;
    lanes[0].inv_qc_relvar = inv_qc_relvar;
    lanes[0].context = context_value;
    lanes[0].mass_in = mass_seed;
    lanes[0].number_in = number_seed;

    for (Int s = 1; s < max_pack_size; ++s) {
      lanes[s].context = false;
      lanes[s].mass_in = Scalar(0.0);
      lanes[s].number_in = Scalar(0.0);
      lanes[s].T_atm = C::T_rainfrz.value - Scalar(1.0);
      lanes[s].lamc = Scalar(5.0);
      lanes[s].mu_c = Scalar(2.0);
      lanes[s].cdist1 = Scalar(1.0);
      lanes[s].qc_incld = Scalar(2.0) * C::QSMALL;
      lanes[s].inv_qc_relvar = Scalar(2.0);
    }

    run_kernel(lanes, exponent);

    return {lanes[0].mass_out, lanes[0].number_out};
  }

  std::array<ImmFreezeLane, max_pack_size> make_lanes() const
  {
    std::array<ImmFreezeLane, max_pack_size> lanes;
    for (Int s = 0; s < max_pack_size; ++s) {
      lanes[s].T_atm = C::T_rainfrz.value - (1.0 + 0.5 * s);
      lanes[s].lamc = 2.5 + 0.15 * s;
      lanes[s].mu_c = 1.0 + 0.2 * s;
      lanes[s].cdist1 = 0.25 + 0.1 * s;
      lanes[s].qc_incld = 5.0 * C::QSMALL;
      lanes[s].inv_qc_relvar = 2.0 + 0.25 * s;
      lanes[s].context = true;
      lanes[s].mass_in = 0.0;
      lanes[s].number_in = 0.0;
      lanes[s].mass_out = 0.0;
      lanes[s].number_out = 0.0;
    }
    return lanes;
  }

  void run_kernel(std::array<ImmFreezeLane, max_pack_size>& lanes,
                  const Scalar exponent) const
  {
    using KTH = KokkosTypes<HostDevice>;

    KTH::view_1d<ImmFreezeLane> lanes_host("lanes_host", max_pack_size);
    view_1d<ImmFreezeLane> lanes_device("lanes_device", max_pack_size);
    std::copy(lanes.begin(), lanes.end(), lanes_host.data());
    Kokkos::deep_copy(lanes_device, lanes_host);

    Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& i) {
      const Int offset = i * Pack::n;
      Pack T_atm, lamc, mu_c, cdist1, qc_incld, inv_qc_relvar;
      Pack qc2qi_hetero_freeze_tend, nc2ni_immers_freeze_tend;
      Mask context;

      for (Int s = 0; s < Pack::n; ++s) {
        const Int vs = offset + s;
        T_atm[s] = lanes_device(vs).T_atm;
        lamc[s] = lanes_device(vs).lamc;
        mu_c[s] = lanes_device(vs).mu_c;
        cdist1[s] = lanes_device(vs).cdist1;
        qc_incld[s] = lanes_device(vs).qc_incld;
        inv_qc_relvar[s] = lanes_device(vs).inv_qc_relvar;
        qc2qi_hetero_freeze_tend[s] = lanes_device(vs).mass_in;
        nc2ni_immers_freeze_tend[s] = lanes_device(vs).number_in;
        context.set(s, lanes_device(vs).context);
      }

      typename Functions::P3Runtime runtime_options;
      runtime_options.immersion_freezing_exponent = exponent;

      Functions::cldliq_immersion_freezing(
        T_atm, lamc, mu_c, cdist1, qc_incld, inv_qc_relvar,
        qc2qi_hetero_freeze_tend, nc2ni_immers_freeze_tend,
        runtime_options, context);

      for (Int s = 0; s < Pack::n; ++s) {
        const Int vs = offset + s;
        lanes_device(vs).mass_out = qc2qi_hetero_freeze_tend[s];
        lanes_device(vs).number_out = nc2ni_immers_freeze_tend[s];
      }
    });

    Kokkos::deep_copy(lanes_host, lanes_device);
    std::copy(lanes_host.data(), lanes_host.data() + max_pack_size, lanes.begin());
  }

  // For active lanes, cloud-liquid immersion freezing computes two moments of
  // the same cloud droplet size distribution:
  //
  //   M = CONS6 * cdist1 * Gamma(mu_c + 7)
  //     * exp(a * (T_zerodegc - T_atm)) * lamc^-6
  //   N = CONS5 * cdist1 * Gamma(mu_c + 4)
  //     * exp(a * (T_zerodegc - T_atm)) * lamc^-3
  //
  // The property tests below avoid duplicating the full production formula.
  // They instead check scaling identities, activation behavior, and the
  // moment-ratio identity:
  //
  //   M / N = (CONS6 / CONS5)
  //         * (mu_c + 4)(mu_c + 5)(mu_c + 6) / lamc^3.
  void run_phys()
  {
    const auto require_rel_close = [&](const Scalar got, const Scalar expected,
                                       const Scalar rtol = identity_tol,
                                       const Scalar atol = zero_tol) {
      INFO("got=" << got << " expected=" << expected
           << " rtol=" << rtol << " atol=" << atol);
      REQUIRE(rel_close(got, expected, rtol, atol));
    };

    const auto require_preserved = [&](const Scalar got, const Scalar expected) {
      INFO("got=" << got << " expected=" << expected);
      REQUIRE(got == expected);
    };

    const auto require_positive_finite = [&](const ImmFreezeResult& result) {
      INFO("mass=" << result.mass << " number=" << result.number);
      REQUIRE(std::isfinite(result.mass));
      REQUIRE(std::isfinite(result.number));
      REQUIRE(result.mass > 0);
      REQUIRE(result.number > 0);
    };

    SECTION("activation_and_thresholds") {
      const Scalar lamc = Scalar(5.0);
      const Scalar mu_c = Scalar(2.0);
      const Scalar cdist1 = Scalar(0.75);
      const Scalar inv_qc_relvar = Scalar(2.0);
      const Scalar exponent = Scalar(0.65);
      const Scalar seed_mass = Scalar(-123.0);
      const Scalar seed_number = Scalar(-456.0);
      const Scalar t_cold = C::T_rainfrz.value - Scalar(5.0);
      const Scalar t_threshold = C::T_rainfrz.value;
      const Scalar t_warm = C::T_rainfrz.value + Scalar(1.0);
      const Scalar qc_small = Scalar(0.9) * C::QSMALL;
      const Scalar qc_edge = C::QSMALL;
      const Scalar qc_active = Scalar(2.0) * C::QSMALL;

      const auto below_qsmall = run_case(t_cold, lamc, mu_c, cdist1, qc_small,
                                         inv_qc_relvar, exponent, true,
                                         seed_mass, seed_number);
      require_preserved(below_qsmall.mass, seed_mass);
      require_preserved(below_qsmall.number, seed_number);

      const auto at_qsmall = run_case(t_cold, lamc, mu_c, cdist1, qc_edge,
                                      inv_qc_relvar, exponent, true,
                                      seed_mass, seed_number);
      require_positive_finite(at_qsmall);

      const auto active = run_case(t_cold, lamc, mu_c, cdist1, qc_active,
                                   inv_qc_relvar, exponent, true,
                                   seed_mass, seed_number);
      require_positive_finite(active);

      const auto at_freezing = run_case(t_threshold, lamc, mu_c, cdist1, qc_active,
                                        inv_qc_relvar, exponent, true,
                                        seed_mass, seed_number);
      require_positive_finite(at_freezing);

      const auto warm = run_case(t_warm, lamc, mu_c, cdist1, qc_active,
                                 inv_qc_relvar, exponent, true,
                                 seed_mass, seed_number);
      require_preserved(warm.mass, seed_mass);
      require_preserved(warm.number, seed_number);

      const auto masked = run_case(t_cold, lamc, mu_c, cdist1, qc_active,
                                   inv_qc_relvar, exponent, false,
                                   seed_mass, seed_number);
      require_preserved(masked.mass, seed_mass);
      require_preserved(masked.number, seed_number);
    }

    SECTION("temperature_supercooling_scaling") {
      const Scalar lamc = Scalar(5.0);
      const Scalar mu_c = Scalar(2.0);
      const Scalar cdist1 = Scalar(0.75);
      const Scalar qc_incld = Scalar(2.0) * C::QSMALL;
      const Scalar inv_qc_relvar = Scalar(2.0);
      const Scalar exponent = Scalar(0.65);
      const Scalar t_warm = C::T_rainfrz.value;
      const Scalar t_cold = C::T_rainfrz.value - Scalar(10.0);

      const auto warm = run_case(t_warm, lamc, mu_c, cdist1, qc_incld,
                                 inv_qc_relvar, exponent, true);
      const auto cold = run_case(t_cold, lamc, mu_c, cdist1, qc_incld,
                                 inv_qc_relvar, exponent, true);
      const auto expected_ratio = std::exp(exponent * (t_warm - t_cold));

      REQUIRE(cold.mass > warm.mass);
      REQUIRE(cold.number > warm.number);
      require_rel_close(cold.mass / warm.mass, expected_ratio);
      require_rel_close(cold.number / warm.number, expected_ratio);
    }

    SECTION("runtime_exponent_sensitivity") {
      const Scalar T_atm = C::T_rainfrz.value - Scalar(5.0);
      const Scalar lamc = Scalar(5.0);
      const Scalar mu_c = Scalar(2.0);
      const Scalar cdist1 = Scalar(0.75);
      const Scalar qc_incld = Scalar(2.0) * C::QSMALL;
      const Scalar inv_qc_relvar = Scalar(2.0);
      const Scalar a1 = Scalar(0.2);
      const Scalar a2 = Scalar(1.1);

      const auto r1 = run_case(T_atm, lamc, mu_c, cdist1, qc_incld,
                               inv_qc_relvar, a1, true);
      const auto r2 = run_case(T_atm, lamc, mu_c, cdist1, qc_incld,
                               inv_qc_relvar, a2, true);
      const auto theta = C::T_zerodegc.value - T_atm;
      const auto expected_ratio = std::exp((a2 - a1) * theta);

      REQUIRE(r2.mass > r1.mass);
      REQUIRE(r2.number > r1.number);
      require_rel_close(r2.mass / r1.mass, expected_ratio);
      require_rel_close(r2.number / r1.number, expected_ratio);

      const auto z1 = run_case(C::T_rainfrz.value, lamc, mu_c, cdist1, qc_incld,
               inv_qc_relvar, Scalar(0.0), true);
      const auto z2 = run_case(C::T_rainfrz.value - Scalar(10.0), lamc, mu_c,
               cdist1, qc_incld, inv_qc_relvar, Scalar(0.0),
               true);
      REQUIRE(z1.mass == z2.mass);
      REQUIRE(z1.number == z2.number);
    }

    SECTION("lambda_power_law_scaling") {
      const Scalar T_atm = C::T_rainfrz.value - Scalar(5.0);
      const Scalar mu_c = Scalar(2.0);
      const Scalar cdist1 = Scalar(0.75);
      const Scalar qc_incld = Scalar(2.0) * C::QSMALL;
      const Scalar inv_qc_relvar = Scalar(2.0);
      const Scalar exponent = Scalar(0.65);
      const Scalar lam1 = Scalar(2.5);
      const Scalar lam2 = Scalar(5.0);

      const auto r1 = run_case(T_atm, lam1, mu_c, cdist1, qc_incld,
                               inv_qc_relvar, exponent, true);
      const auto r2 = run_case(T_atm, lam2, mu_c, cdist1, qc_incld,
                               inv_qc_relvar, exponent, true);

      require_rel_close(r2.mass / r1.mass, std::pow(lam1 / lam2, 6));
      require_rel_close(r2.number / r1.number, std::pow(lam1 / lam2, 3));
    }

    SECTION("mass_number_moment_identity") {
      constexpr std::array<Scalar, 3> mus = {
        Scalar(0.0), Scalar(2.0), Scalar(5.0)
      };
      constexpr std::array<Scalar, 3> lams = {
        Scalar(2.0), Scalar(5.0), Scalar(10.0)
      };
      const Scalar T_atm = C::T_rainfrz.value - Scalar(5.0);
      const Scalar cdist1 = Scalar(0.75);
      const Scalar qc_incld = Scalar(2.0) * C::QSMALL;
      const Scalar inv_qc_relvar = Scalar(2.0);
      const Scalar exponent = Scalar(0.65);

      for (const auto mu : mus) {
        for (const auto lam : lams) {
          const auto result = run_case(T_atm, lam, mu, cdist1, qc_incld,
                                       inv_qc_relvar, exponent, true);
          require_positive_finite(result);
          require_rel_close(result.mass / result.number,
                            expected_mass_number_ratio(mu, lam));
        }
      }
    }

    SECTION("volume_nucleation_prefers_larger_than_mean_droplets") {
      constexpr std::array<Scalar, 3> mus = {
        Scalar(0.0), Scalar(2.0), Scalar(5.0)
      };
      constexpr std::array<Scalar, 3> lams = {
        Scalar(2.0), Scalar(5.0), Scalar(10.0)
      };
      const Scalar T_atm = C::T_rainfrz.value - Scalar(5.0);
      const Scalar cdist1 = Scalar(0.75);
      const Scalar qc_incld = Scalar(2.0) * C::QSMALL;
      const Scalar inv_qc_relvar = Scalar(2.0);
      const Scalar exponent = Scalar(0.65);

      for (const auto mu : mus) {
        for (const auto lam : lams) {
          const auto result = run_case(T_atm, lam, mu, cdist1, qc_incld,
                                       inv_qc_relvar, exponent, true);
          require_positive_finite(result);
          const Scalar frozen_mean_mass = result.mass / result.number;
          const Scalar cloud_mean_mass =
            (C::CONS6 / C::CONS5)
            * (mu + 1) * (mu + 2) * (mu + 3)
            / cube_host(lam);
          const Scalar expected_ratio =
            ((mu + 4) * (mu + 5) * (mu + 6))
            / ((mu + 1) * (mu + 2) * (mu + 3));

          REQUIRE(frozen_mean_mass > cloud_mean_mass);
          require_rel_close(frozen_mean_mass / cloud_mean_mass,
                            expected_ratio);
        }
      }
    }

    SECTION("distribution_prefactor_scaling") {
      const Scalar T_atm = C::T_rainfrz.value - Scalar(5.0);
      const Scalar lamc = Scalar(5.0);
      const Scalar mu_c = Scalar(2.0);
      const Scalar qc_incld = Scalar(2.0) * C::QSMALL;
      const Scalar inv_qc_relvar = Scalar(2.0);
      const Scalar exponent = Scalar(0.65);
      const Scalar c1 = Scalar(0.25);
      const Scalar c2 = Scalar(1.25);

      const auto r1 = run_case(T_atm, lamc, mu_c, c1, qc_incld,
                               inv_qc_relvar, exponent, true);
      const auto r2 = run_case(T_atm, lamc, mu_c, c2, qc_incld,
                               inv_qc_relvar, exponent, true);

      require_rel_close(r2.mass / r1.mass, c2 / c1);
      require_rel_close(r2.number / r1.number, c2 / c1);
    }

    SECTION("zero_distribution_prefactor_gives_zero") {
      const Scalar T_atm = C::T_rainfrz.value - Scalar(5.0);
      const Scalar lamc = Scalar(5.0);
      const Scalar mu_c = Scalar(2.0);
      const Scalar cdist1 = Scalar(0.0);
      const Scalar qc_incld = Scalar(2.0) * C::QSMALL;
      const Scalar inv_qc_relvar = Scalar(2.0);
      const Scalar exponent = Scalar(0.65);

      const auto result = run_case(T_atm, lamc, mu_c, cdist1, qc_incld,
                                   inv_qc_relvar, exponent, true,
                                   Scalar(-123.0), Scalar(-456.0));

      REQUIRE(result.mass == 0);
      REQUIRE(result.number == 0);
    }

    SECTION("scalar_multiplier_cancellation_in_mass_number_ratio") {
      const Scalar T_base = C::T_rainfrz.value - Scalar(5.0);
      const Scalar T_cold = C::T_rainfrz.value - Scalar(10.0);
      const Scalar lamc = Scalar(5.0);
      const Scalar mu_c = Scalar(2.0);
      const Scalar qc_incld = Scalar(2.0) * C::QSMALL;
      const Scalar inv_qc_relvar = Scalar(2.0);
      const Scalar c_base = Scalar(0.5);
      const Scalar c_bigger = Scalar(1.4);
      const Scalar a_base = Scalar(0.65);
      const Scalar a_bigger = Scalar(1.1);

      const auto base = run_case(T_base, lamc, mu_c, c_base, qc_incld,
                                 inv_qc_relvar, a_base, true);
      const auto colder = run_case(T_cold, lamc, mu_c, c_base, qc_incld,
                                   inv_qc_relvar, a_base, true);
      const auto bigger_a = run_case(T_base, lamc, mu_c, c_base, qc_incld,
                                     inv_qc_relvar, a_bigger, true);
      const auto bigger_c = run_case(T_base, lamc, mu_c, c_bigger, qc_incld,
                                     inv_qc_relvar, a_base, true);
      const auto base_ratio = base.mass / base.number;

      require_rel_close(colder.mass / colder.number, base_ratio);
      require_rel_close(bigger_a.mass / bigger_a.number, base_ratio);
      require_rel_close(bigger_c.mass / bigger_c.number, base_ratio);
    }

    SECTION("qc_gate_only_behavior") {
      const Scalar T_atm = C::T_rainfrz.value - Scalar(5.0);
      const Scalar lamc = Scalar(5.0);
      const Scalar mu_c = Scalar(2.0);
      const Scalar cdist1 = Scalar(0.75);
      const Scalar inv_qc_relvar = Scalar(2.0);
      const Scalar exponent = Scalar(0.65);
      const Scalar qc1 = C::QSMALL;
      const Scalar qc2 = Scalar(10.0) * C::QSMALL;

      const auto r1 = run_case(T_atm, lamc, mu_c, cdist1, qc1,
                               inv_qc_relvar, exponent, true);
      const auto r2 = run_case(T_atm, lamc, mu_c, cdist1, qc2,
                               inv_qc_relvar, exponent, true);

      REQUIRE(r1.mass == r2.mass);
      REQUIRE(r1.number == r2.number);
    }

    SECTION("inv_qc_relvar_currently_inactive") {
      const Scalar T_atm = C::T_rainfrz.value - Scalar(5.0);
      const Scalar lamc = Scalar(5.0);
      const Scalar mu_c = Scalar(2.0);
      const Scalar cdist1 = Scalar(0.75);
      const Scalar qc_incld = Scalar(2.0) * C::QSMALL;
      const Scalar exponent = Scalar(0.65);

      const auto low_r = run_case(T_atm, lamc, mu_c, cdist1, qc_incld,
                                  Scalar(0.5), exponent, true);
      const auto high_r = run_case(T_atm, lamc, mu_c, cdist1, qc_incld,
                                   Scalar(10.0), exponent, true);

      REQUIRE(low_r.mass == high_r.mass);
      REQUIRE(low_r.number == high_r.number);
    }

    SECTION("context_mask_preserves_inactive_lanes") {
      auto lanes = make_lanes();
      const Scalar exponent = Scalar(0.65);
      for (Int s = 0; s < max_pack_size; ++s) {
        lanes[s].context = false;
        lanes[s].mass_in = Scalar(-1000.0) - s;
        lanes[s].number_in = Scalar(-2000.0) - s;
      }

      // These are scenario slots in a max_pack_size host array, not assumptions
      // about Pack::n. run_kernel processes them in chunks of Pack::n, so this
      // remains valid when SCREAM_PACK_SIZE == 1.
      lanes[0].context = true; // active
      lanes[1].context = true; // active

      lanes[2].context = true; // warm but context true: threshold-inactive
      lanes[2].T_atm = C::T_rainfrz.value + Scalar(1.0);

      lanes[3].context = true; // below QSMALL but context true: threshold-inactive
      lanes[3].qc_incld = Scalar(0.9) * C::QSMALL;

      lanes[4].T_atm = C::T_rainfrz.value - Scalar(5.0); // context false but otherwise active
      lanes[4].qc_incld = Scalar(5.0) * C::QSMALL;

      lanes[5].context = true; // active
      lanes[7].context = true; // active at inclusive temperature threshold
      lanes[7].T_atm = C::T_rainfrz.value;

      run_kernel(lanes, exponent);

      for (Int s = 0; s < max_pack_size; ++s) {
        const bool active = lanes[s].context
          && lanes[s].qc_incld >= C::QSMALL
          && lanes[s].T_atm <= C::T_rainfrz.value;
        INFO("lane=" << s << " active=" << active
             << " mass_out=" << lanes[s].mass_out
             << " number_out=" << lanes[s].number_out);
        if (active) {
          REQUIRE(std::isfinite(lanes[s].mass_out));
          REQUIRE(std::isfinite(lanes[s].number_out));
          REQUIRE(lanes[s].mass_out > 0);
          REQUIRE(lanes[s].number_out > 0);
        } else {
          REQUIRE(lanes[s].mass_out == lanes[s].mass_in);
          REQUIRE(lanes[s].number_out == lanes[s].number_in);
        }
      }
    }

    SECTION("finite_nonnegative_outputs") {
      constexpr std::array<Scalar, 3> mus = {
        Scalar(0.0), Scalar(2.0), Scalar(5.0)
      };
      constexpr std::array<Scalar, 3> lams = {
        Scalar(2.0), Scalar(5.0), Scalar(10.0)
      };
      const Scalar exponent = Scalar(0.65);
      const Scalar cdist1 = Scalar(0.75);
      const Scalar qc_incld = Scalar(2.0) * C::QSMALL;
      const Scalar inv_qc_relvar = Scalar(2.0);
      const Scalar T_atm = C::T_rainfrz.value - Scalar(5.0);

      for (Int i = 0; i < static_cast<Int>(mus.size()); ++i) {
        const auto result = run_case(T_atm, lams[i], mus[i], cdist1, qc_incld,
                                     inv_qc_relvar, exponent, true);
        INFO("mass=" << result.mass << " number=" << result.number);
        REQUIRE(std::isfinite(result.mass));
        REQUIRE(std::isfinite(result.number));
        REQUIRE(result.mass >= 0);
        REQUIRE(result.number >= 0);
        REQUIRE(result.mass > 0);
        REQUIRE(result.number > 0);
      }
    }
}

void run_bfb()
{
  // This BFB path uses the default all-true context, so cloud-liquid
  // immersion freezing is active only when cloud liquid is present above
  // QSMALL and the temperature is at or below the rain-freezing threshold.
  constexpr Scalar qsmall = C::QSMALL;

  constexpr Scalar t_freezing = 0.9 * C::T_rainfrz.value,
                   t_not_freezing = 2.0 * C::T_rainfrz.value;
  constexpr Scalar qc_incld_small = 0.9 * qsmall;
  constexpr Scalar qc_incld_not_small = 2.0 * qsmall;
  constexpr Scalar lamc1 = 0.1, lamc2 = 0.2, lamc3 = 0.3, lamc4 = 0.4;
  constexpr Scalar mu_c1 = 0.2, mu_c2 = 0.4, mu_c3 = 0.6, mu_c4 = 0.8;
  constexpr Scalar cdist11 = 0.25, cdist12 = 0.5, cdist13 = 0.75, cdist14 = 1.0;
  constexpr Scalar inv_qc_relvar_val = 1;

  CldliqImmersionFreezingData cldliq_imm_freezing_data[max_pack_size] = {
    // T_atm, lamc, mu_c, cdist1, qc_incld, inv_qc_relvar
    {t_not_freezing, lamc1, mu_c1, cdist11, qc_incld_small, inv_qc_relvar_val},
    {t_not_freezing, lamc2, mu_c2, cdist12, qc_incld_small, inv_qc_relvar_val},
    {t_not_freezing, lamc3, mu_c3, cdist13, qc_incld_small, inv_qc_relvar_val},
    {t_not_freezing, lamc4, mu_c4, cdist14, qc_incld_small, inv_qc_relvar_val},

    {t_not_freezing, lamc1, mu_c1, cdist11, qc_incld_not_small, inv_qc_relvar_val},
    {t_not_freezing, lamc2, mu_c2, cdist12, qc_incld_not_small, inv_qc_relvar_val},
    {t_not_freezing, lamc3, mu_c3, cdist13, qc_incld_not_small, inv_qc_relvar_val},
    {t_not_freezing, lamc4, mu_c4, cdist14, qc_incld_not_small, inv_qc_relvar_val},

    {t_freezing, lamc1, mu_c1, cdist11, qc_incld_small, inv_qc_relvar_val},
    {t_freezing, lamc2, mu_c2, cdist12, qc_incld_small, inv_qc_relvar_val},
    {t_freezing, lamc3, mu_c3, cdist13, qc_incld_small, inv_qc_relvar_val},
    {t_freezing, lamc4, mu_c4, cdist14, qc_incld_small, inv_qc_relvar_val},

    {t_freezing, lamc1, mu_c1, cdist11, qc_incld_not_small, inv_qc_relvar_val},
    {t_freezing, lamc2, mu_c2, cdist12, qc_incld_not_small, inv_qc_relvar_val},
    {t_freezing, lamc3, mu_c3, cdist13, qc_incld_not_small, inv_qc_relvar_val},
    {t_freezing, lamc4, mu_c4, cdist14, qc_incld_not_small, inv_qc_relvar_val}
  };

  // Sync to device
  view_1d<CldliqImmersionFreezingData> device_data("cldliq_imm_freezing", max_pack_size);
  const auto host_data = Kokkos::create_mirror_view(device_data);
  std::copy(&cldliq_imm_freezing_data[0], &cldliq_imm_freezing_data[0] + max_pack_size,
            host_data.data());
  Kokkos::deep_copy(device_data, host_data);

  // Read baseline data
  if (this->m_baseline_action == COMPARE) {
    for (Int i = 0; i < max_pack_size; ++i) {
      cldliq_imm_freezing_data[i].read(Base::m_ifile);
    }
  }

  // Run the lookup from a kernel and copy results back to host
  Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& i) {
    const Int offset = i * Pack::n;

    // Init pack inputs
    Pack T_atm, lamc, mu_c, cdist1, qc_incld, inv_qc_relvar;
    for (Int s = 0, vs = offset; s < Pack::n; ++s, ++vs) {
      T_atm[s]        = device_data(vs).T_atm;
      lamc[s]         = device_data(vs).lamc;
      mu_c[s]         = device_data(vs).mu_c;
      cdist1[s]       = device_data(vs).cdist1;
      qc_incld[s]     = device_data(vs).qc_incld;
      inv_qc_relvar[s] = device_data(vs).inv_qc_relvar;
    }

    Pack qc2qi_hetero_freeze_tend{0.0};
    Pack nc2ni_immers_freeze_tend{0.0};

    typename Functions::P3Runtime runtime_options;
    Functions::cldliq_immersion_freezing(
      T_atm, lamc, mu_c, cdist1, qc_incld, inv_qc_relvar,
      qc2qi_hetero_freeze_tend, nc2ni_immers_freeze_tend,
      runtime_options);

    // Copy results back into views
    for (Int s = 0, vs = offset; s < Pack::n; ++s, ++vs) {
      device_data(vs).qc2qi_hetero_freeze_tend  = qc2qi_hetero_freeze_tend[s];
      device_data(vs).nc2ni_immers_freeze_tend  = nc2ni_immers_freeze_tend[s];
    }
  });

  // Sync back to host.
  Kokkos::deep_copy(host_data, device_data);

  // Validate results.
  if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
    for (Int s = 0; s < max_pack_size; ++s) {
      REQUIRE(cldliq_imm_freezing_data[s].qc2qi_hetero_freeze_tend == host_data[s].qc2qi_hetero_freeze_tend);
      REQUIRE(cldliq_imm_freezing_data[s].nc2ni_immers_freeze_tend == host_data[s].nc2ni_immers_freeze_tend);
    }
  }
  else if (this->m_baseline_action == GENERATE) {
    for (Int s = 0; s < max_pack_size; ++s) {
      host_data(s).write(Base::m_ofile);
    }
  }
}

};

}
}
}

namespace {

TEST_CASE("p3_cldliq_immersion_freezing", "[p3_functions]")
{
  using T = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestCldliqImmersionFreezing;

  T t;
  SECTION("properties") {
    t.run_phys();
  }

  SECTION("bfb") {
    t.run_bfb();
  }
}

} // namespace
