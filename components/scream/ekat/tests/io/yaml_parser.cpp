#include <catch2/catch.hpp>
#include "ekat/scream_parse_yaml_file.hpp"

namespace {

TEST_CASE ("yaml_parser","") {
  using namespace scream;

  std::string fname = "input.yaml";
  ParameterList params("parameters");
  REQUIRE_NOTHROW ( parse_yaml_file(fname,params) );

  // Check some of the loaded parameters.
  // NOTE: if you change input.yaml, you may have to
  //       change some of these checks as well.

  REQUIRE (params.isSublist("Constants"));
  REQUIRE (params.isSublist("Options"));

  auto& constants = params.sublist("Constants");
  auto& options   = params.sublist("Options");

  REQUIRE (constants.isParameter("Two Logicals"));
  REQUIRE (constants.isParameter("One String"));
  REQUIRE (options.isParameter("My Int"));
  REQUIRE (options.isParameter("My Bool"));

  std::vector<int> logicals;
  REQUIRE_NOTHROW(logicals = constants.get<std::vector<int>>("Two Logicals"));

  REQUIRE (logicals.size()==2);
  REQUIRE (logicals[0] == 1);
  REQUIRE (logicals[1] == 0);

  std::string str = "multiple words string";
  REQUIRE (constants.get<std::string>("One String") == str);

  REQUIRE (options.get<int>("My Int") == -2);
  REQUIRE (options.get<bool>("My Bool") == false);

  REQUIRE (constants.isSublist("Nested Sublist"));
  auto& sl = constants.sublist("Nested Sublist");
  REQUIRE (sl.get<int>("The Answer") == 42);
}

} // anonymous namespace
