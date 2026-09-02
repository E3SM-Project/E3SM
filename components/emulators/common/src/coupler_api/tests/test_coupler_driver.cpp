#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "fake_components.hpp"
#include <coupler.hpp>
#include <coupler_driver.hpp>

#include <filesystem>

using namespace e3sm::coupler;
using namespace e3sm::coupler::test;

static_assert(E3SMComponent<FakeAtmosphere>);
static_assert(E3SMComponent<FakeLand>);

TEST_CASE("CouplerDriver registers and couples fake components") {
  constexpr std::size_t ncols = 4;

  FakeAtmosphere atm{ncols};
  FakeLand lnd{ncols};

  CouplerDriver driver;

  driver.add_component(atm);
  driver.add_component(lnd);

  const std::string config = TEST_DATA_DIR "/fake_coupling.yaml";

  driver.initialize(config);

  REQUIRE(lnd.temperature().size() == ncols);
  REQUIRE(lnd.surface_flux().size() == ncols);
  REQUIRE(lnd.precipitation().size() == ncols);

  REQUIRE(atm.temperature().size() == ncols);
  REQUIRE(atm.surface_flux().size() == ncols);
  REQUIRE(atm.precipitation().size() == ncols);

  SECTION("atmosphere fields are exchanged to land") {
    driver.run_component(atm.name());
    for (std::size_t i = 0; i < ncols; ++i) {
      REQUIRE(atm.temperature()[i] == 280.0 + static_cast<double>(i));
      REQUIRE(atm.precipitation()[i] == 0.1 * static_cast<double>(i));
      REQUIRE(atm.surface_flux()[i] == 0.0);
    }

    driver.run_component(lnd.name());
    for (std::size_t i = 0; i < ncols; ++i) {
      REQUIRE(lnd.temperature()[i] == atm.temperature()[i]);
      REQUIRE(lnd.precipitation()[i] == atm.precipitation()[i]);
      REQUIRE(lnd.surface_flux()[i] == 2.0 * lnd.precipitation()[i]);
    }
    driver.run_component(atm.name());
    for (std::size_t i = 0; i < ncols; ++i) {
      REQUIRE(atm.surface_flux()[i] == lnd.surface_flux()[i]);
    }
  }
}
