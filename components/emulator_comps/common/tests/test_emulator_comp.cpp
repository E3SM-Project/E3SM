// Catch2 v2 single header
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "emulator_comp.hpp"
#include "emulator_context.hpp"

namespace emulator {
namespace test {

// Concrete implementation for testing
class TestComp : public EmulatorComp {
public:
  TestComp() : EmulatorComp(CompType::ATM) {}

  // Track calls for verification
  bool init_called = false;
  bool run_called = false;
  bool final_called = false;
  int last_dt = 0;

protected:
  void init_impl() override { init_called = true; }
  void run_impl(int dt) override {
    run_called = true;
    last_dt = dt;
  }
  void final_impl() override { final_called = true; }
};

// Test components for different CompTypes
class TestOcnComp : public EmulatorComp {
public:
  TestOcnComp() : EmulatorComp(CompType::OCN) {}

protected:
  void init_impl() override {}
  void run_impl(int) override {}
  void final_impl() override {}
};

class TestIceComp : public EmulatorComp {
public:
  TestIceComp() : EmulatorComp(CompType::ICE) {}

protected:
  void init_impl() override {}
  void run_impl(int) override {}
  void final_impl() override {}
};

class TestLndComp : public EmulatorComp {
public:
  TestLndComp() : EmulatorComp(CompType::LND) {}

protected:
  void init_impl() override {}
  void run_impl(int) override {}
  void final_impl() override {}
};

TEST_CASE("EmulatorComp construction", "[emulator_comp]") {
  TestComp comp;
  REQUIRE(comp.type() == CompType::ATM);
  REQUIRE(comp.comp_id() == -1);
  REQUIRE(comp.name().empty()); // name is empty before create_instance
  REQUIRE_FALSE(comp.is_initialized());
  REQUIRE(comp.step_count() == 0);
}

TEST_CASE("EmulatorComp different types", "[emulator_comp]") {
  TestComp atm;
  TestOcnComp ocn;
  TestIceComp ice;
  TestLndComp lnd;

  REQUIRE(atm.type() == CompType::ATM);
  REQUIRE(ocn.type() == CompType::OCN);
  REQUIRE(ice.type() == CompType::ICE);
  REQUIRE(lnd.type() == CompType::LND);
}

TEST_CASE("EmulatorComp create_instance", "[emulator_comp]") {
  TestComp comp;
  comp.create_instance(42, "test_atm");

  REQUIRE(comp.comp_id() == 42);
  REQUIRE(comp.name() == "test_atm");
}

TEST_CASE("EmulatorComp lifecycle", "[emulator_comp]") {
  TestComp comp;
  comp.create_instance(1, "test");

  SECTION("initialize calls init_impl") {
    REQUIRE_FALSE(comp.init_called);
    comp.initialize();
    REQUIRE(comp.init_called);
    REQUIRE(comp.is_initialized());
  }

  SECTION("run calls run_impl and increments step count") {
    comp.initialize();
    REQUIRE(comp.step_count() == 0);

    comp.run(3600);
    REQUIRE(comp.run_called);
    REQUIRE(comp.last_dt == 3600);
    REQUIRE(comp.step_count() == 1);

    comp.run(1800);
    REQUIRE(comp.step_count() == 2);
  }

  SECTION("finalize calls final_impl") {
    comp.initialize();
    REQUIRE_FALSE(comp.final_called);
    comp.finalize();
    REQUIRE(comp.final_called);
    REQUIRE_FALSE(comp.is_initialized());
  }
}

TEST_CASE("EmulatorComp error handling", "[emulator_comp]") {
  TestComp comp;
  comp.create_instance(1, "test");

  SECTION("run before initialize throws") {
    REQUIRE_THROWS_AS(comp.run(100), std::runtime_error);
  }

  SECTION("double initialize throws") {
    comp.initialize();
    REQUIRE_THROWS_AS(comp.initialize(), std::runtime_error);
  }

  SECTION("finalize without initialize is safe") {
    REQUIRE_NOTHROW(comp.finalize());
  }

  SECTION("re-initialization after finalize works") {
    comp.initialize();
    REQUIRE(comp.is_initialized());
    comp.run(100);
    REQUIRE(comp.step_count() == 1);

    comp.finalize();
    REQUIRE_FALSE(comp.is_initialized());

    // Re-initialize and verify it works
    comp.initialize();
    REQUIRE(comp.is_initialized());
    comp.run(200);
    REQUIRE(comp.step_count() == 2); // step count persists
  }
}

TEST_CASE("EmulatorComp with EmulatorContext", "[emulator_comp][integration]") {
  auto &ctx = EmulatorContext::singleton();
  ctx.clean_up();

  // Create component via context
  auto &comp = ctx.create<TestComp>();

  REQUIRE(ctx.has<TestComp>());
  REQUIRE(comp.type() == CompType::ATM);

  // Use component from context
  comp.create_instance(99, "context_test");
  comp.initialize();
  comp.run(100);

  // Verify through context access
  const auto &ref = ctx.get<TestComp>();
  REQUIRE(ref.comp_id() == 99);
  REQUIRE(ref.step_count() == 1);

  ctx.clean_up();
}

} // namespace test
} // namespace emulator
