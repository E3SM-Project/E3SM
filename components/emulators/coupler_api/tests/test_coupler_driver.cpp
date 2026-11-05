#include <memory>
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "fake_components.hpp"
#include <coupler.hpp>
#include <coupler_driver.hpp>

using namespace e3sm::coupler;
using namespace e3sm::coupler::test;

static_assert(CoupledProcess<FakeAtmosphere>);
static_assert(CoupledProcess<FakeLand>);

namespace {
// Registers nothing, so it can be added to a driver without being in the plan
class UnplannedComponent {
public:
  static constexpr std::string_view name() { return "unplanned"; }
  void populate_registry(FieldRegistry&) {}
  void configure_coupling(const ActiveCouplingFields& plan) {
    coupling_.configure(plan);
  }
  void run(const ProcessRunOpts) {}
  void export_fields(std::span<ExportBuffer> buffers) const {
    coupling_.export_fields(buffers);
  }
  void import_fields(std::span<const ImportBuffer> buffers) {
    coupling_.import_fields(buffers);
  }
  std::span<const int> get_local_gids() const {
    return coupling_.get_local_gids();
  }

private:
  CouplingAdapter coupling_;
};
static_assert(CoupledProcess<UnplannedComponent>);

// Has `state` as both an input and an output, and advances it by one each run
class StatefulComponent {
public:
  explicit StatefulComponent(std::size_t ncols)
      : input_(ncols, 0.0), output_(ncols, 0.0) {}

  static constexpr std::string_view name() { return "stateful"; }
  void populate_registry(FieldRegistry& registry) {
    coupling_.register_coupled_field(registry,
                                     {.role = FieldRole::Import,
                                      .component = "stateful",
                                      .name = "state",
                                      .size = input_.size()},
                                     input_);
    coupling_.register_coupled_field(registry,
                                     {.role = FieldRole::Export,
                                      .component = "stateful",
                                      .name = "state",
                                      .size = output_.size()},
                                     output_);
  }
  void configure_coupling(const ActiveCouplingFields& plan) {
    coupling_.configure(plan);
  }
  void run(const ProcessRunOpts) {
    for (std::size_t i = 0; i < input_.size(); ++i) {
      output_[i] = input_[i] + 1.0;
    }
  }
  void export_fields(std::span<ExportBuffer> buffers) const {
    coupling_.export_fields(buffers);
  }
  void import_fields(std::span<const ImportBuffer> buffers) {
    coupling_.import_fields(buffers);
  }
  std::span<const int> get_local_gids() const {
    return coupling_.get_local_gids();
  }

private:
  std::vector<double> input_;
  std::vector<double> output_;
  CouplingAdapter coupling_;
};
static_assert(CoupledProcess<StatefulComponent>);

// Binds one caller-owned span as `state` in both directions
class ExternalComponent {
public:
  explicit ExternalComponent(std::span<double> state) : state_(state) {}

  static constexpr std::string_view name() { return "external"; }
  void populate_registry(FieldRegistry& registry) {
    coupling_.register_coupled_field(registry,
                                     {.role = FieldRole::Import,
                                      .component = "external",
                                      .name = "state",
                                      .size = state_.size()},
                                     state_);
    coupling_.register_coupled_field(registry,
                                     {.role = FieldRole::Export,
                                      .component = "external",
                                      .name = "state",
                                      .size = state_.size()},
                                     state_);
  }
  void configure_coupling(const ActiveCouplingFields& plan) {
    coupling_.configure(plan);
  }
  void run(const ProcessRunOpts) {}
  void export_fields(std::span<ExportBuffer> buffers) const {
    coupling_.export_fields(buffers);
  }
  void import_fields(std::span<const ImportBuffer> buffers) {
    coupling_.import_fields(buffers);
  }
  std::span<const int> get_local_gids() const {
    return coupling_.get_local_gids();
  }

private:
  std::span<double> state_;
  CouplingAdapter coupling_;
};
static_assert(CoupledProcess<ExternalComponent>);
} // namespace

TEST_CASE("CouplerDriver rejects an unknown merge type") {
  constexpr std::size_t ncols = 4;

  CouplerDriver driver{};

  driver.add_component(std::make_unique<FakeAtmosphere>(ncols));
  driver.add_component(std::make_unique<FakeLand>(ncols));

  const std::string config = TEST_DATA_DIR "/invalid_merge_type.yaml";

  REQUIRE_THROWS_WITH(driver.initialize(config),
                      Catch::Matchers::Contains("Unknown merge type"));
}

TEST_CASE("A prognostic field advances through external storage") {
  std::vector<double> state{10.0, 20.0};

  CouplerDriver driver{};
  driver.add_component(std::make_unique<StatefulComponent>(state.size()));
  driver.add_component(std::make_unique<ExternalComponent>(state));
  driver.initialize(TEST_DATA_DIR "/prognostic_plan.yaml");

  const ProcessRunOpts opts{3600.0, 1};
  for (int step = 1; step <= 3; ++step) {
    CAPTURE(step);
    driver.export_component("external");
    driver.run_component("stateful", opts);
    driver.import_component("external");

    REQUIRE(state[0] == 10.0 + step);
    REQUIRE(state[1] == 20.0 + step);
  }
}

TEST_CASE("Two routes cannot deliver to the same input") {
  using Catch::Matchers::Contains;
  std::vector<double> state(2, 0.0);

  CouplerDriver driver{};
  driver.add_component(std::make_unique<StatefulComponent>(state.size()));
  driver.add_component(std::make_unique<ExternalComponent>(state));

  REQUIRE_THROWS_WITH(driver.initialize(TEST_DATA_DIR "/double_delivery_plan.yaml"),
                      Contains("state_again") && Contains("state_in") &&
                          Contains("'stateful'"));
}

TEST_CASE("FieldRegistry identifies a field by component, name, and role") {
  using Catch::Matchers::Contains;
  FieldRegistry registry;

  const auto import_id = registry.register_field(
      {.role = FieldRole::Import, .component = "c", .name = "f", .size = 1});
  const auto export_id = registry.register_field(
      {.role = FieldRole::Export, .component = "c", .name = "f", .size = 2});

  REQUIRE(import_id != export_id);
  REQUIRE(registry.get_id("c", "f", FieldRole::Import) == import_id);
  REQUIRE(registry.get_id("c", "f", FieldRole::Export) == export_id);
  REQUIRE(registry.get("c", "f", FieldRole::Export).size == 2);

  REQUIRE_THROWS_WITH(
      registry.register_field(
          {.role = FieldRole::Export, .component = "c", .name = "f", .size = 2}),
      Contains("'f'") && Contains("Export") && Contains("twice"));
}

TEST_CASE("CouplerDriver errors name the problem") {
  using Catch::Matchers::Contains;
  constexpr std::size_t ncols = 4;

  CouplerDriver driver{};
  driver.add_component(std::make_unique<FakeAtmosphere>(ncols));
  driver.add_component(std::make_unique<FakeLand>(ncols));

  SECTION("component in the driver but not the plan") {
    driver.add_component(std::make_unique<UnplannedComponent>());
    REQUIRE_THROWS_WITH(driver.initialize(TEST_DATA_DIR "/fake_coupling.yaml"),
                        Contains("'unplanned'") &&
                            Contains("not in the coupling plan") &&
                            Contains("atm"));
  }

  SECTION("source that imports the field") {
    REQUIRE_THROWS_WITH(
        driver.initialize(TEST_DATA_DIR "/source_imports_plan.yaml"),
        Contains("'surface_flux'") && Contains("'atm'") &&
            Contains("Sources must export"));
  }

  SECTION("destination that exports the field") {
    REQUIRE_THROWS_WITH(
        driver.initialize(TEST_DATA_DIR "/destination_exports_plan.yaml"),
        Contains("'temperature'") && Contains("'atm'") &&
            Contains("Destinations must import"));
  }

  SECTION("component listed twice for a field") {
    REQUIRE_THROWS_WITH(
        driver.initialize(TEST_DATA_DIR "/duplicate_destination_plan.yaml"),
        Contains("'temperature'") && Contains("'lnd'") &&
            Contains("more than once"));
  }

  SECTION("field with no destinations") {
    REQUIRE_THROWS_WITH(
        driver.initialize(TEST_DATA_DIR "/no_destinations_plan.yaml"),
        Contains("'temperature'") && Contains("no destinations"));
  }

  SECTION("route without a field name") {
    REQUIRE_THROWS_WITH(
        driver.initialize(TEST_DATA_DIR "/missing_name_plan.yaml"),
        Contains("missing_name_plan.yaml") && Contains("name"));
  }

  SECTION("field the component did not register") {
    REQUIRE_THROWS_WITH(
        driver.initialize(TEST_DATA_DIR "/unknown_field_plan.yaml"),
        Contains("'not_a_field'") && Contains("'atm'") &&
            Contains("temperature"));
  }

  SECTION("component that was never added") {
    REQUIRE_THROWS_WITH(
        driver.initialize(TEST_DATA_DIR "/unknown_component_plan.yaml"),
        Contains("'ocn'") && Contains("added to the driver"));
  }

  SECTION("missing plan file") {
    REQUIRE_THROWS_WITH(driver.initialize("no_such_plan.yaml"),
                        Contains("Coupling plan file not found") &&
                            Contains("no_such_plan.yaml"));
  }

  SECTION("exchange before initialize") {
    REQUIRE_THROWS_WITH(driver.export_component("atm"),
                        Contains("before initialize()"));
    REQUIRE_THROWS_WITH(driver.import_component("atm"),
                        Contains("before initialize()"));
  }

  SECTION("unknown component name") {
    driver.initialize(TEST_DATA_DIR "/fake_coupling.yaml");
    REQUIRE_THROWS_WITH(driver.get_component("ocn"),
                        Contains("'ocn'") && Contains("atm") &&
                            Contains("lnd"));
  }
}

TEST_CASE("CouplerDriver rejects destinations of a different size") {
  using Catch::Matchers::Contains;

  CouplerDriver driver{};
  driver.add_component(std::make_unique<FakeAtmosphere>(4));
  driver.add_component(std::make_unique<FakeLand>(5));

  REQUIRE_THROWS_WITH(driver.initialize(TEST_DATA_DIR "/fake_coupling.yaml"),
                      Contains("value(s)") && Contains("'atm'") &&
                          Contains("'lnd'"));
}

TEST_CASE("CouplerDriver registers and couples fake components") {
  constexpr std::size_t ncols = 4;

  CouplerDriver driver{};

  driver.add_component(std::make_unique<FakeAtmosphere>(ncols));
  driver.add_component(std::make_unique<FakeLand>(ncols));

  const std::string config = TEST_DATA_DIR "/fake_coupling.yaml";

  driver.initialize(config);

  FakeLand& lnd = driver.get_component_as<FakeLand>("lnd");
  FakeAtmosphere& atm = driver.get_component_as<FakeAtmosphere>("atm");

  REQUIRE(lnd.temperature().size() == ncols);
  REQUIRE(lnd.surface_flux().size() == ncols);
  REQUIRE(lnd.precipitation().size() == ncols);

  REQUIRE(atm.temperature().size() == ncols);
  REQUIRE(atm.surface_flux().size() == ncols);
  REQUIRE(atm.precipitation().size() == ncols);

  // Create Run Opt
  const ProcessRunOpts opts{3600.0,1};

  SECTION("atmosphere fields are exchanged to land") {
    driver.run_component(atm.name(), opts);
    for (std::size_t i = 0; i < ncols; ++i) {
      REQUIRE(atm.temperature()[i] == 280.0 + static_cast<double>(i));
      REQUIRE(atm.precipitation()[i] == 0.1 * static_cast<double>(i));
      REQUIRE(atm.surface_flux()[i] == 0.0);
    }

    driver.run_component(lnd.name(), opts);
    for (std::size_t i = 0; i < ncols; ++i) {
      REQUIRE(lnd.temperature()[i] == atm.temperature()[i]);
      REQUIRE(lnd.precipitation()[i] == atm.precipitation()[i]);
      REQUIRE(lnd.surface_flux()[i] == 2.0 * lnd.precipitation()[i]);
    }
    driver.run_component(atm.name(), opts);
    for (std::size_t i = 0; i < ncols; ++i) {
      REQUIRE(atm.surface_flux()[i] == lnd.surface_flux()[i]);
    }
  }
}
