// Catch2 v2 single header
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "fpe_guard.hpp"

#include <cfenv>

namespace emulator {
namespace inference {
namespace test {

#ifdef __GLIBC__

TEST_CASE("FpeGuard suspends traps and restores the prior set", "[fpe]") {
  // The set an E3SM debug build enables
  const int debug_traps = FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW;

  fedisableexcept(FE_ALL_EXCEPT);
  feenableexcept(debug_traps);

  {
    FpeGuard guard;
    REQUIRE(fegetexcept() == 0);

    // Would trap without the guard
    volatile double zero = 0.0;
    volatile double inf = 1.0 / zero;
    (void)inf;
  }

  REQUIRE(fegetexcept() == debug_traps);
  // A flag left set would trap as soon as the traps came back
  REQUIRE(fetestexcept(FE_DIVBYZERO) == 0);

  fedisableexcept(FE_ALL_EXCEPT);
}

TEST_CASE("FpeGuards nest", "[fpe]") {
  fedisableexcept(FE_ALL_EXCEPT);
  feenableexcept(FE_INVALID);
  {
    FpeGuard outer;
    {
      FpeGuard inner;
      REQUIRE(fegetexcept() == 0);
    }
    REQUIRE(fegetexcept() == 0);
  }
  REQUIRE(fegetexcept() == FE_INVALID);

  fedisableexcept(FE_ALL_EXCEPT);
}

#endif

TEST_CASE("FpeGuard with no traps enabled changes nothing", "[fpe]") {
  FpeGuard guard;
  SUCCEED();
}

} // namespace test
} // namespace inference
} // namespace emulator
