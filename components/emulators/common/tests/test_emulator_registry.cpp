// Catch2 v2 single header - define CATCH_CONFIG_MAIN in one translation unit
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "emulator_registry.hpp"

#include <string>

namespace emulator {
namespace test {

// Simple test type for the registry
struct TestComponent {
  int value;
  std::string name;

  TestComponent() : value(0), name("default") {}
  TestComponent(int v, const std::string &n) : value(v), name(n) {}
};

// Another test type to verify multiple types work
struct AnotherComponent {
  double data;
  explicit AnotherComponent(double d = 0.0) : data(d) {}
};

TEST_CASE("EmulatorRegistry instance access", "[emulator_registry]") {
  auto &reg1 = EmulatorRegistry::instance();
  auto &reg2 = EmulatorRegistry::instance();
  REQUIRE(&reg1 == &reg2);
}

TEST_CASE("EmulatorRegistry create and retrieve", "[emulator_registry]") {
  auto &reg = EmulatorRegistry::instance();
  reg.clean_up(); // Start fresh

  SECTION("Create default-constructed object") {
    auto &comp = reg.create<TestComponent>("default_comp");
    REQUIRE(comp.value == 0);
    REQUIRE(comp.name == "default");
  }

  SECTION("Create with arguments") {
    auto &comp = reg.create<TestComponent>("arg_comp", 42, "test");
    REQUIRE(comp.value == 42);
    REQUIRE(comp.name == "test");
  }

  reg.clean_up();
}

TEST_CASE("EmulatorRegistry get methods", "[emulator_registry]") {
  auto &reg = EmulatorRegistry::instance();
  reg.clean_up();

  reg.create<TestComponent>("getter_comp", 100, "getter_test");

  SECTION("get() returns const reference") {
    const auto &comp = reg.get<TestComponent>("getter_comp");
    REQUIRE(comp.value == 100);
    REQUIRE(comp.name == "getter_test");
  }

  SECTION("get_mut() returns mutable reference") {
    auto &comp = reg.get_mut<TestComponent>("getter_comp");
    comp.value = 200;
    comp.name = "modified";

    const auto &check = reg.get<TestComponent>("getter_comp");
    REQUIRE(check.value == 200);
    REQUIRE(check.name == "modified");
  }

  reg.clean_up();
}

TEST_CASE("cleanup_emulator_registry free function", "[emulator_registry]") {
  auto &reg = EmulatorRegistry::instance();
  reg.clean_up();

  reg.create<TestComponent>("cleanup_comp", 999, "cleanup_test");
  REQUIRE(reg.has("cleanup_comp"));

  cleanup_emulator_registry(); // Test the convenience function
  REQUIRE_FALSE(reg.has("cleanup_comp"));
}

TEST_CASE("EmulatorRegistry has() method", "[emulator_registry]") {
  auto &reg = EmulatorRegistry::instance();
  reg.clean_up();

  REQUIRE_FALSE(reg.has("my_comp"));

  reg.create<TestComponent>("my_comp");
  REQUIRE(reg.has("my_comp"));

  reg.clean_up();
  REQUIRE_FALSE(reg.has("my_comp"));
}

TEST_CASE("EmulatorRegistry multiple instances of same type",
          "[emulator_registry]") {
  auto &reg = EmulatorRegistry::instance();
  reg.clean_up();

  // Create multiple instances of the same type with different names
  auto &comp1 = reg.create<TestComponent>("comp1", 1, "first");
  auto &comp2 = reg.create<TestComponent>("comp2", 2, "second");

  REQUIRE(reg.has("comp1"));
  REQUIRE(reg.has("comp2"));

  REQUIRE(comp1.value == 1);
  REQUIRE(comp2.value == 2);

  // Retrieve by name
  const auto &ref1 = reg.get<TestComponent>("comp1");
  const auto &ref2 = reg.get<TestComponent>("comp2");

  REQUIRE(ref1.name == "first");
  REQUIRE(ref2.name == "second");

  reg.clean_up();
}

TEST_CASE("EmulatorRegistry multiple types", "[emulator_registry]") {
  auto &reg = EmulatorRegistry::instance();
  reg.clean_up();

  reg.create<TestComponent>("test_comp", 1, "first");
  reg.create<AnotherComponent>("another_comp", 3.14);

  REQUIRE(reg.has("test_comp"));
  REQUIRE(reg.has("another_comp"));

  const auto &tc = reg.get<TestComponent>("test_comp");
  const auto &ac = reg.get<AnotherComponent>("another_comp");

  REQUIRE(tc.value == 1);
  REQUIRE(ac.data == 3.14);

  reg.clean_up();
}

TEST_CASE("EmulatorRegistry error handling", "[emulator_registry]") {
  auto &reg = EmulatorRegistry::instance();
  reg.clean_up();

  SECTION("get() throws when object not found") {
    REQUIRE_THROWS_AS(reg.get<TestComponent>("nonexistent"),
                      std::runtime_error);
  }

  SECTION("get_mut() throws when object not found") {
    REQUIRE_THROWS_AS(reg.get_mut<TestComponent>("nonexistent"),
                      std::runtime_error);
  }

  SECTION("create() throws on duplicate name") {
    reg.create<TestComponent>("dup_comp");
    REQUIRE_THROWS_AS(reg.create<TestComponent>("dup_comp"),
                      std::runtime_error);
  }

  reg.clean_up();
}

} // namespace test
} // namespace emulator
