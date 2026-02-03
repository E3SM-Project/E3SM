// Catch2 v2 single header
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "emulator.hpp"
#include "emulator_registry.hpp"

namespace emulator {
namespace test {

// Concrete implementation for testing
class TestEmulator : public Emulator {
public:
  TestEmulator(int id = -1, const std::string &name = "")
      : Emulator(EmulatorType::ATM_COMP, id, name) {}

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

// Test emulators for different EmulatorTypes
class TestOcnEmulator : public Emulator {
public:
  TestOcnEmulator(int id = -1, const std::string &name = "")
      : Emulator(EmulatorType::OCN_COMP, id, name) {}

protected:
  void init_impl() override {}
  void run_impl(int) override {}
  void final_impl() override {}
};

class TestIceEmulator : public Emulator {
public:
  TestIceEmulator(int id = -1, const std::string &name = "")
      : Emulator(EmulatorType::ICE_COMP, id, name) {}

protected:
  void init_impl() override {}
  void run_impl(int) override {}
  void final_impl() override {}
};

class TestLndEmulator : public Emulator {
public:
  TestLndEmulator(int id = -1, const std::string &name = "")
      : Emulator(EmulatorType::LND_COMP, id, name) {}

protected:
  void init_impl() override {}
  void run_impl(int) override {}
  void final_impl() override {}
};

TEST_CASE("Emulator construction", "[emulator]") {
  TestEmulator emu;
  REQUIRE(emu.type() == EmulatorType::ATM_COMP);
  REQUIRE(emu.id() == -1);
  REQUIRE(emu.name().empty());
  REQUIRE_FALSE(emu.is_initialized());
  REQUIRE(emu.step_count() == 0);
}

TEST_CASE("Emulator construction with args", "[emulator]") {
  TestEmulator emu(42, "test_atm");
  REQUIRE(emu.id() == 42);
  REQUIRE(emu.name() == "test_atm");
}

TEST_CASE("Emulator different types", "[emulator]") {
  TestEmulator atm;
  TestOcnEmulator ocn;
  TestIceEmulator ice;
  TestLndEmulator lnd;

  REQUIRE(atm.type() == EmulatorType::ATM_COMP);
  REQUIRE(ocn.type() == EmulatorType::OCN_COMP);
  REQUIRE(ice.type() == EmulatorType::ICE_COMP);
  REQUIRE(lnd.type() == EmulatorType::LND_COMP);
}

TEST_CASE("Emulator lifecycle", "[emulator]") {
  TestEmulator emu(1, "test");

  SECTION("initialize calls init_impl") {
    REQUIRE_FALSE(emu.init_called);
    emu.initialize();
    REQUIRE(emu.init_called);
    REQUIRE(emu.is_initialized());
  }

  SECTION("run calls run_impl and increments step count") {
    emu.initialize();
    REQUIRE(emu.step_count() == 0);

    emu.run(3600);
    REQUIRE(emu.run_called);
    REQUIRE(emu.last_dt == 3600);
    REQUIRE(emu.step_count() == 1);

    emu.run(1800);
    REQUIRE(emu.step_count() == 2);
  }

  SECTION("finalize calls final_impl") {
    emu.initialize();
    REQUIRE_FALSE(emu.final_called);
    emu.finalize();
    REQUIRE(emu.final_called);
    REQUIRE_FALSE(emu.is_initialized());
  }
}

TEST_CASE("Emulator error handling", "[emulator]") {
  TestEmulator emu(1, "test");

  SECTION("run before initialize throws") {
    REQUIRE_THROWS_AS(emu.run(100), std::runtime_error);
  }

  SECTION("double initialize throws") {
    emu.initialize();
    REQUIRE_THROWS_AS(emu.initialize(), std::runtime_error);
  }

  SECTION("finalize without initialize is safe") {
    REQUIRE_NOTHROW(emu.finalize());
  }

  SECTION("re-initialization after finalize works") {
    emu.initialize();
    REQUIRE(emu.is_initialized());
    emu.run(100);
    REQUIRE(emu.step_count() == 1);

    emu.finalize();
    REQUIRE_FALSE(emu.is_initialized());

    // Re-initialize and verify it works
    emu.initialize();
    REQUIRE(emu.is_initialized());
    emu.run(200);
    REQUIRE(emu.step_count() == 2); // step count persists
  }
}

TEST_CASE("Emulator with EmulatorRegistry", "[emulator][integration]") {
  auto &reg = EmulatorRegistry::instance();
  reg.clean_up();

  // Create emulator via registry with name
  auto &emu = reg.create<TestEmulator>("test_emu", 99, "registry_test");

  REQUIRE(reg.has("test_emu"));
  REQUIRE(emu.type() == EmulatorType::ATM_COMP);
  REQUIRE(emu.id() == 99);

  // Use emulator from registry
  emu.initialize();
  emu.run(100);

  // Verify through registry access
  const auto &ref = reg.get<TestEmulator>("test_emu");
  REQUIRE(ref.id() == 99);
  REQUIRE(ref.step_count() == 1);

  // Test get_mut
  auto &mut_ref = reg.get_mut<TestEmulator>("test_emu");
  mut_ref.run(200);
  REQUIRE(ref.step_count() == 2);

  reg.clean_up();
}

TEST_CASE("EmulatorRegistry get_mut throws for unknown name",
          "[emulator_registry]") {
  auto &reg = EmulatorRegistry::instance();
  reg.clean_up();

  REQUIRE_THROWS_AS(reg.get_mut<TestEmulator>("nonexistent"),
                    std::runtime_error);

  reg.clean_up();
}

} // namespace test
} // namespace emulator
