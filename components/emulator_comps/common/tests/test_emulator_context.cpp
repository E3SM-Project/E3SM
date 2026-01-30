// Catch2 v2 single header - define CATCH_CONFIG_MAIN in one translation unit
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "emulator_context.hpp"

#include <string>

namespace emulator {
namespace test {

// Simple test type for the context
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

TEST_CASE("EmulatorContext singleton access", "[emulator_context]") {
  auto &ctx1 = EmulatorContext::singleton();
  auto &ctx2 = EmulatorContext::singleton();
  REQUIRE(&ctx1 == &ctx2);
}

TEST_CASE("EmulatorContext create and retrieve", "[emulator_context]") {
  auto &ctx = EmulatorContext::singleton();
  ctx.clean_up(); // Start fresh

  SECTION("Create default-constructed object") {
    auto &comp = ctx.create<TestComponent>();
    REQUIRE(comp.value == 0);
    REQUIRE(comp.name == "default");
  }

  SECTION("Create with arguments") {
    auto &comp = ctx.create<TestComponent>(42, "test");
    REQUIRE(comp.value == 42);
    REQUIRE(comp.name == "test");
  }

  ctx.clean_up();
}

TEST_CASE("EmulatorContext get methods", "[emulator_context]") {
  auto &ctx = EmulatorContext::singleton();
  ctx.clean_up();

  ctx.create<TestComponent>(100, "getter_test");

  SECTION("get() returns const reference") {
    const auto &comp = ctx.get<TestComponent>();
    REQUIRE(comp.value == 100);
    REQUIRE(comp.name == "getter_test");
  }

  SECTION("get_non_const() returns mutable reference") {
    auto &comp = ctx.get_non_const<TestComponent>();
    comp.value = 200;
    comp.name = "modified";

    const auto &check = ctx.get<TestComponent>();
    REQUIRE(check.value == 200);
    REQUIRE(check.name == "modified");
  }

  ctx.clean_up();
}

TEST_CASE("cleanup_emulator_context free function", "[emulator_context]") {
  auto &ctx = EmulatorContext::singleton();
  ctx.clean_up();

  ctx.create<TestComponent>(999, "cleanup_test");
  REQUIRE(ctx.has<TestComponent>());

  cleanup_emulator_context(); // Test the convenience function
  REQUIRE_FALSE(ctx.has<TestComponent>());
}

TEST_CASE("EmulatorContext has() method", "[emulator_context]") {
  auto &ctx = EmulatorContext::singleton();
  ctx.clean_up();

  REQUIRE_FALSE(ctx.has<TestComponent>());

  ctx.create<TestComponent>();
  REQUIRE(ctx.has<TestComponent>());

  ctx.clean_up();
  REQUIRE_FALSE(ctx.has<TestComponent>());
}

TEST_CASE("EmulatorContext multiple types", "[emulator_context]") {
  auto &ctx = EmulatorContext::singleton();
  ctx.clean_up();

  ctx.create<TestComponent>(1, "first");
  ctx.create<AnotherComponent>(3.14);

  REQUIRE(ctx.has<TestComponent>());
  REQUIRE(ctx.has<AnotherComponent>());

  const auto &tc = ctx.get<TestComponent>();
  const auto &ac = ctx.get<AnotherComponent>();

  REQUIRE(tc.value == 1);
  REQUIRE(ac.data == 3.14);

  ctx.clean_up();
}

TEST_CASE("EmulatorContext error handling", "[emulator_context]") {
  auto &ctx = EmulatorContext::singleton();
  ctx.clean_up();

  SECTION("get() throws when object not found") {
    REQUIRE_THROWS_AS(ctx.get<TestComponent>(), std::runtime_error);
  }

  SECTION("get_non_const() throws when object not found") {
    REQUIRE_THROWS_AS(ctx.get_non_const<TestComponent>(), std::runtime_error);
  }

  SECTION("create() throws on duplicate") {
    ctx.create<TestComponent>();
    REQUIRE_THROWS_AS(ctx.create<TestComponent>(), std::runtime_error);
  }

  ctx.clean_up();
}

} // namespace test
} // namespace emulator
