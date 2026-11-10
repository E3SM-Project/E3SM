#include <catch2/catch_test_macros.hpp>
#include <dexpr/ast.hpp>
#include <dexpr/lexer.hpp>
#include <dexpr/parser.hpp>
#include <dexpr/supported_functions.hpp>

#include <string>

namespace dexpr {

namespace { // anonymous

// Stands in for a component's vocabulary, covering each shape validate_calls
// has to reason about.
FunctionRegistry test_registry() {
  FunctionRegistry reg;
  reg.add({.name = "mean",
           .desc = "average over a dimension",
           .min_positional = 1,
           .max_positional = 1,
           .keywords = {{"weights", false}},
           .example = "T_mid.mean('lev')"});
  reg.add({.name = "prev",
           .desc = "value at the previous step",
           .min_positional = 0,
           .max_positional = 0,
           .keywords = {},
           .example = "T_mid.prev()"});
  reg.add({.name = "clamp",
           .desc = "bound a field",
           .min_positional = 1,
           .max_positional = 2,
           .keywords = {{"fill", true}},
           .example = "clamp(T_mid, fill=0)"});
  return reg;
}

// Returns "" if `input` validates, else the ValidationError message.
std::string validation_error(const std::string& input,
                             const FunctionRegistry& reg) {
  parser::Parser parser{Lexer{input}};
  const auto expr = parser.parse();
  REQUIRE(expr != nullptr);
  try {
    validate_calls(*expr, reg);
  } catch (const ValidationError& e) {
    return e.what();
  }
  return "";
}

bool contains(const std::string& haystack, const std::string& needle) {
  return haystack.find(needle) != std::string::npos;
}

// Returns "" if `reg` describes itself consistently, else the error message.
std::string registry_error(const FunctionRegistry& reg) {
  try {
    validate_registry(reg);
  } catch (const ValidationError& e) {
    return e.what();
  }
  return "";
}

// Every field spelled out: -Wmissing-field-initializers requires it.
FunctionSpec spec(std::string name, int min_pos, int max_pos,
                  std::vector<ParamSpec> kws, std::string example) {
  return FunctionSpec{std::move(name), "a description", min_pos, max_pos,
                      std::move(kws),  std::move(example)};
}

} // anonymous namespace

TEST_CASE("registry_add_and_find") {
  FunctionRegistry reg;
  CHECK(reg.find("mean") == nullptr);
  CHECK_FALSE(reg.contains("mean"));

  reg.add(spec("mean", 0, 0, {}, "x.mean()"));
  REQUIRE(reg.find("mean") != nullptr);
  CHECK(reg.find("mean")->desc == "a description");
  CHECK(reg.contains("mean"));

  // Lookups are case sensitive, like identifiers everywhere else.
  CHECK(reg.find("Mean") == nullptr);
}

TEST_CASE("registry_rejects_duplicates") {
  FunctionRegistry reg;
  reg.add(spec("mean", 0, 0, {}, "x.mean()"));
  CHECK_THROWS_AS(reg.add(spec("mean", 0, 0, {}, "x.mean()")),
                  std::invalid_argument);
  // The first registration is the one that stands.
  REQUIRE(reg.find("mean") != nullptr);
  CHECK(reg.find("mean")->example == "x.mean()");

  CHECK_THROWS_AS(reg.add(spec("", 0, 0, {}, "x()")),
                  std::invalid_argument);
}

TEST_CASE("registry_names_are_ordered") {
  const auto reg = test_registry();
  const auto names = reg.names();
  REQUIRE(names.size() == 3);
  CHECK(names[0] == "clamp");
  CHECK(names[1] == "mean");
  CHECK(names[2] == "prev");
}

TEST_CASE("builtin_functions_are_seeded") {
  const auto& reg = builtin_functions();
  CHECK(reg.contains("where"));
  CHECK(reg.contains("sum"));
  CHECK(reg.contains("derivative"));
  CHECK(reg.contains("tend"));
  // Nothing consults the builtins implicitly; an empty registry knows nothing.
  CHECK_FALSE(FunctionRegistry{}.contains("where"));
}

TEST_CASE("validate_accepts_well_formed_calls") {
  const auto reg = test_registry();

  CHECK(validation_error("T_mid.mean('lev')", reg) == "");
  CHECK(validation_error("T_mid.mean('lev', weights='dp')", reg) == "");
  CHECK(validation_error("T_mid.prev()", reg) == "");
  CHECK(validation_error("clamp(T_mid, fill=0)", reg) == "");
  CHECK(validation_error("clamp(T_mid, 1, fill=0)", reg) == "");
  // `f(x)` and `x.f()` both validate: the receiver is not positional. A
  // component that implements only one spelling rejects the other itself.
  CHECK(validation_error("T_mid.clamp(1, fill=0)", reg) == "");

  // Expressions with no calls at all are trivially valid.
  CHECK(validation_error("(qc+qv)*p_mid", reg) == "");
  CHECK(validation_error("['col', 'lev']", reg) == "");
}

TEST_CASE("validate_rejects_unknown_function") {
  const auto reg = test_registry();

  const auto msg = validation_error("T_mid.meen('lev')", reg);
  CHECK(contains(msg, "unknown function 'meen'"));
  // The message lists what is available.
  CHECK(contains(msg, "clamp, mean, prev"));

  CHECK(contains(validation_error("nope(x)", FunctionRegistry{}),
                 "no functions are registered"));
}

TEST_CASE("validate_rejects_bad_arity") {
  const auto reg = test_registry();

  CHECK(contains(validation_error("T_mid.mean()", reg),
                 "'mean' takes 1 positional argument(s), got 0"));
  CHECK(contains(validation_error("T_mid.mean('lev', 'col')", reg),
                 "'mean' takes 1 positional argument(s), got 2"));
  CHECK(contains(validation_error("T_mid.prev(1)", reg),
                 "'prev' takes 0 positional argument(s), got 1"));
  CHECK(contains(validation_error("clamp(fill=0)", reg),
                 "'clamp' takes 1 to 2 positional argument(s), got 0"));

  // Keyword arguments do not count towards positional arity.
  CHECK(validation_error("T_mid.mean('lev', weights='dp')", reg) == "");
}

TEST_CASE("validate_rejects_bad_keywords") {
  const auto reg = test_registry();

  const auto msg = validation_error("T_mid.mean('lev', weight='dp')", reg);
  CHECK(contains(msg, "'mean' has no argument 'weight'"));
  CHECK(contains(msg, "accepts: weights"));

  CHECK(contains(validation_error("T_mid.prev(x=1)", reg),
                 "'prev' has no argument 'x'"));
  CHECK(contains(validation_error("T_mid.mean('lev', weights='dp', weights='dz')", reg),
                 "got argument 'weights' more than once"));
  CHECK(contains(validation_error("clamp(T_mid)", reg),
                 "'clamp' requires argument 'fill'"));
}

TEST_CASE("validate_recurses_into_operands_and_arguments") {
  const auto reg = test_registry();

  // Receiver of a method call.
  CHECK(contains(validation_error("T_mid.meen('lev').mean('col')", reg),
                 "unknown function 'meen'"));
  // Inside a positional argument.
  CHECK(contains(validation_error("clamp(a.meen('lev'), fill=0)", reg),
                 "unknown function 'meen'"));
  // Inside a keyword's value; the keyword name itself is not looked up.
  CHECK(contains(validation_error("T_mid.mean('lev', weights=a.meen('c'))", reg),
                 "unknown function 'meen'"));
  // Across operators, unary operators, and array elements.
  CHECK(contains(validation_error("qc + a.meen('lev')", reg),
                 "unknown function 'meen'"));
  CHECK(contains(validation_error("not a.meen('lev')", reg),
                 "unknown function 'meen'"));
  CHECK(contains(validation_error("[a.meen('lev')]", reg),
                 "unknown function 'meen'"));
}

TEST_CASE("validate_collects_every_problem") {
  const auto reg = test_registry();

  // All three in one pass.
  const auto msg = validation_error("T_mid.meen('lev') + qc.prev(1) + a.mean()", reg);
  CHECK(contains(msg, "unknown function 'meen'"));
  CHECK(contains(msg, "'prev' takes 0 positional argument(s), got 1"));
  CHECK(contains(msg, "'mean' takes 1 positional argument(s), got 0"));
}

TEST_CASE("validate_rejects_non_name_callee") {
  const auto reg = test_registry();
  CHECK(contains(validation_error("(a+b)(c)", reg),
                 "call target is not a function name"));
}

TEST_CASE("registry_rejects_malformed_specs") {
  FunctionRegistry reg;
  // Arity that can never be satisfied
  CHECK_THROWS_AS(reg.add(spec("a", 2, 1, {}, "x.a(1)")), std::invalid_argument);
  CHECK_THROWS_AS(reg.add(spec("b", -1, -1, {}, "x.b()")), std::invalid_argument);
  // Keyword arguments that cannot be addressed
  CHECK_THROWS_AS(reg.add(spec("c", 0, 0, {{"", true}}, "x.c()")),
                  std::invalid_argument);
  CHECK_THROWS_AS(reg.add(spec("e", 0, 0, {{"dims", true}, {"dims", false}}, "x.e()")),
                  std::invalid_argument);
  // None of them landed
  CHECK(reg.names().empty());
}

TEST_CASE("validate_registry_accepts_a_self_consistent_vocabulary") {
  CHECK(registry_error(test_registry()) == "");
  // The builtins must describe themselves correctly too
  CHECK(registry_error(builtin_functions()) == "");
}

TEST_CASE("validate_registry_catches_a_spec_its_example_contradicts") {
  // Declared arity contradicts the example.
  FunctionRegistry wrong_arity;
  wrong_arity.add(spec("mean", 0, 0, {}, "T_mid.mean('lev')"));
  CHECK(contains(registry_error(wrong_arity), "does not match the spec"));

  // A keyword the spec never declared
  FunctionRegistry wrong_kw;
  wrong_kw.add(spec("mean", 1, 1, {}, "T_mid.mean('lev', weights='dp')"));
  CHECK(contains(registry_error(wrong_kw), "has no argument 'weights'"));

  // A required keyword the example forgets
  FunctionRegistry missing_kw;
  missing_kw.add(spec("isel", 0, 0, {{"lev", true}}, "T_mid.isel()"));
  CHECK(contains(registry_error(missing_kw), "requires argument 'lev'"));

  // An example that is not an expression at all
  FunctionRegistry unparseable;
  unparseable.add(spec("mean", 1, 1, {}, "T_mid.mean('lev'"));
  CHECK(contains(registry_error(unparseable), "does not parse"));

  // An example demonstrating somebody else's function
  FunctionRegistry wrong_function;
  wrong_function.add(spec("mean", 1, 1, {}, "T_mid.mean('lev')"));
  wrong_function.add(spec("prev", 0, 0, {}, "T_mid.mean('lev')"));
  CHECK(contains(registry_error(wrong_function), "never calls 'prev'"));

  // No example at all: nothing to check the spec against
  FunctionRegistry no_example;
  no_example.add(spec("mean", 1, 1, {}, ""));
  CHECK(contains(registry_error(no_example), "no example"));
}

TEST_CASE("validate_registry_reports_every_bad_spec_at_once") {
  FunctionRegistry reg;
  reg.add(spec("a", 0, 0, {}, "x.a("));
  reg.add(spec("b", 0, 0, {}, ""));
  const auto msg = registry_error(reg);
  CHECK(contains(msg, "'a'"));
  CHECK(contains(msg, "'b'"));
}

TEST_CASE("function_spec_to_string") {
  FunctionSpec s{.name = "mean",
                 .desc = "average over a dimension",
                 .min_positional = 1,
                 .max_positional = 2,
                 .keywords = {{"weights", false}, {"dim", true}},
                 .example = "T_mid.mean('lev', dim='x')"};
  const auto str = s.to_string();
  CHECK(contains(str, "mean("));
  CHECK(contains(str, "<arg>"));
  CHECK(contains(str, "[<arg>]"));
  CHECK(contains(str, "[weights=..]"));
  CHECK(contains(str, "dim=.."));
  CHECK(contains(str, "average over a dimension"));

  // `dexpr functions` shows how to write the call, not just its shape.
  CHECK(contains(str, "e.g. T_mid.mean('lev', dim='x')"));
}

} // namespace dexpr
