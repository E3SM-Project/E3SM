#include "catch2/catch.hpp"

#include "p3_functions.hpp"
#include "p3_unit_tests_common.hpp"

#include "share/core/eamxx_types.hpp"

namespace scream {
namespace p3 {
namespace unit_test {

/*
 * Property tests of the merge of emulated warm-rain rates into P3's
 * (warm_rain_emulator_merge).
 */

template <typename D>
struct UnitWrap::UnitTest<D>::TestWarmRainEmulatorMerge : public UnitWrap::UnitTest<D>::Base {

  void run_property()
  {
    int nerr = 0;
    Kokkos::parallel_reduce("warm_rain_emulator_merge", 1, KOKKOS_LAMBDA(const Int&, int& err) {
      const Scalar dt = 100, inv_dt = 1/dt, scale = 0.36;
      // P3's own rates (au, nc2nr, ncautr, nc_sc, ac, nc_ac, nr_sc), and the emulated ones
      const Scalar p3[7]  = {1e-9, 2e2, 3e1, 0, 4e-9, 5e2, 6};
      const Scalar emu[7] = {7e-9, 8e2, 9e1, -1e3, 1e-8, 2e3, 3};
      const Scalar nc0 = 5e7, nr0 = 1e4;

      // The masks, as lane-dependent as packs allow: every combination shows up
      for (int combo = 0; combo < 4; ++combo) {
        for (int self_collection = 0; self_collection < 2; ++self_collection) {
          for (int in_context = 0; in_context < 2; ++in_context) {
            Pack use_cloud, use_rain;
            for (int s = 0; s < Pack::n; ++s) {
              const int c = (combo + s) % 4;
              use_cloud[s] = c % 2;
              use_rain[s]  = c / 2;
            }
            Pack r[7];
            for (int n = 0; n < 7; ++n) r[n] = p3[n];
            Functions::warm_rain_emulator_merge(
              Pack(nc0), Pack(nr0), inv_dt, scale, self_collection == 1,
              Pack(emu[0]), Pack(emu[1]), Pack(emu[2]), Pack(emu[3]), Pack(emu[4]), Pack(emu[5]),
              Pack(emu[6]), use_cloud, use_rain,
              r[0], r[1], r[2], r[3], r[4], r[5], r[6], Mask(in_context == 1));

            for (int s = 0; s < Pack::n; ++s) {
              const bool uc = use_cloud[s] > 0.5, ur = use_rain[s] > 0.5;
              Scalar expected[7];
              for (int n = 0; n < 7; ++n) expected[n] = p3[n];
              if (in_context == 1) {
                if (uc) {
                  // The emulator's cloud rates
                  for (int n = 0; n < 6; ++n) expected[n] = emu[n];
                  if (self_collection == 0) expected[3] = 0;
                } else {
                  // P3, with scaled autoconversion
                  for (int n = 0; n < 3; ++n) expected[n] = p3[n] * scale;
                }
                if (ur) expected[6] = emu[6];
              }
              for (int n = 0; n < 7; ++n) {
                if (r[n][s] != expected[n]) ++err;
              }
            }
          }
        }
      }

      // Limiters: emulated rates cannot break what the conservation checks assume
      {
        Pack r[7];
        for (int n = 0; n < 7; ++n) r[n] = p3[n];
        // Sources that are negative, and self-collection that removes more than there is
        Functions::warm_rain_emulator_merge(
          Pack(nc0), Pack(nr0), inv_dt, scale, true,
          Pack(-1), Pack(-1), Pack(-1), Pack(-1e30), Pack(-1), Pack(-1), Pack(1e30),
          Pack(1), Pack(1), r[0], r[1], r[2], r[3], r[4], r[5], r[6]);
        for (int s = 0; s < Pack::n; ++s) {
          if (r[0][s] != 0 || r[1][s] != 0 || r[2][s] != 0 || r[4][s] != 0 || r[5][s] != 0) ++err;
          if (r[3][s] != -nc0 * inv_dt) ++err;
          // (ncautr is zero here)
          if (r[6][s] != nr0 * inv_dt) ++err;
        }

        // Self-collection of the wrong sign, and negative number concentrations
        for (int n = 0; n < 7; ++n) r[n] = p3[n];
        Functions::warm_rain_emulator_merge(
          Pack(-nc0), Pack(-nr0), inv_dt, scale, true,
          Pack(emu[0]), Pack(emu[1]), Pack(emu[2]), Pack(1e3), Pack(emu[4]), Pack(emu[5]), Pack(-5),
          Pack(1), Pack(1), r[0], r[1], r[2], r[3], r[4], r[5], r[6]);
        for (int s = 0; s < Pack::n; ++s) {
          if (r[3][s] != 0 || r[6][s] != 0) ++err;
        }

        // Rain self-collection can also remove what autoconversion adds
        for (int n = 0; n < 7; ++n) r[n] = p3[n];
        Functions::warm_rain_emulator_merge(
          Pack(nc0), Pack(nr0), inv_dt, scale, true,
          Pack(emu[0]), Pack(emu[1]), Pack(emu[2]), Pack(emu[3]), Pack(emu[4]), Pack(emu[5]), Pack(1e30),
          Pack(1), Pack(1), r[0], r[1], r[2], r[3], r[4], r[5], r[6]);
        for (int s = 0; s < Pack::n; ++s) {
          // (to round-off: the compiler may contract the sum)
          const Scalar expected = nr0 * inv_dt + emu[2];
          if (Kokkos::abs(r[6][s] - expected) > 1e-12 * expected) ++err;
        }
      }
    }, nerr);
    Kokkos::fence();
    REQUIRE(nerr == 0);
  }

}; // TestWarmRainEmulatorMerge

} // namespace unit_test
} // namespace p3
} // namespace scream

namespace {

TEST_CASE("p3_warm_rain_emulator_merge", "[p3_warm_rain_emulator]"){
  using T = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestWarmRainEmulatorMerge;

  T t;
  t.run_property();
}

} // namespace
