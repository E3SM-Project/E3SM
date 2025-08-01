#include <catch2/catch.hpp>

#include "share/io/eamxx_io_utils.hpp"

#include "ekat/ekat_parameter_list.hpp"
#include "ekat/ekat_assert.hpp"

#include <vector>
#include <string>
#include <map>

namespace scream {

TEST_CASE("io_field_aliasing") {
  SECTION("parse_field_alias") {
    // Test normal case with alias
    auto [alias1, field1] = parse_field_alias("LWP:=:LiqWaterPath");
    REQUIRE(alias1 == "LWP");
    REQUIRE(field1 == "LiqWaterPath");
    
    // Test with spaces
    auto [alias2, field2] = parse_field_alias(" SWP :=: SolidWaterPath ");
    REQUIRE(alias2 == "SWP");
    REQUIRE(field2 == "SolidWaterPath");
    
    // Test no alias
    auto [alias3, field3] = parse_field_alias("Temperature");
    REQUIRE(alias3 == "Temperature");
    REQUIRE(field3 == "Temperature");
    
    // Test complex field names
    auto [alias4, field4] = parse_field_alias("RH:=:RelativeHumidity");
    REQUIRE(alias4 == "RH");
    REQUIRE(field4 == "RelativeHumidity");
    
    // Test error cases
    REQUIRE_THROWS(parse_field_alias("LWP:=:"));
    REQUIRE_THROWS(parse_field_alias(":=:LiqWaterPath"));
    REQUIRE_THROWS(parse_field_alias("  :=:  "));
  }
  
  SECTION("process_field_aliases") {
    std::vector<std::string> field_specs = {
      "T",                          // No alias
      "LWP:=:LiqWaterPath",        // With alias
      "SWP:=:SolidWaterPath",      // With alias
      "qv",                        // No alias
      "RH:=:RelativeHumidity"      // With alias
    };
    
    auto [alias_map, alias_names] = process_field_aliases(field_specs);
    
    REQUIRE(alias_map.size() == 5);
    REQUIRE(alias_names.size() == 5);
    
    // Check mappings
    REQUIRE(alias_map["T"] == "T");
    REQUIRE(alias_map["LWP"] == "LiqWaterPath");
    REQUIRE(alias_map["SWP"] == "SolidWaterPath");
    REQUIRE(alias_map["qv"] == "qv");
    REQUIRE(alias_map["RH"] == "RelativeHumidity");
    
    // Check order preservation
    REQUIRE(alias_names[0] == "T");
    REQUIRE(alias_names[1] == "LWP");
    REQUIRE(alias_names[2] == "SWP");
    REQUIRE(alias_names[3] == "qv");
    REQUIRE(alias_names[4] == "RH");
  }
  
  SECTION("duplicate_alias_detection") {
    std::vector<std::string> duplicate_specs = {
      "LWP:=:LiqWaterPath",
      "LWP:=:SolidWaterPath"  // Duplicate alias
    };
    REQUIRE_THROWS(process_field_aliases(duplicate_specs));
  }
  
  SECTION("mixed_aliases_and_regular_fields") {
    std::vector<std::string> mixed_specs = {
      "temperature",
      "LWP:=:LiqWaterPath",
      "pressure",
      "RH:=:RelativeHumidity"
    };
    
    auto [alias_map, alias_names] = process_field_aliases(mixed_specs);
    
    REQUIRE(alias_map.size() == 4);
    REQUIRE(alias_map["temperature"] == "temperature");
    REQUIRE(alias_map["LWP"] == "LiqWaterPath");
    REQUIRE(alias_map["pressure"] == "pressure");
    REQUIRE(alias_map["RH"] == "RelativeHumidity");
  }
  
  SECTION("empty_field_list") {
    std::vector<std::string> empty_specs;
    auto [alias_map, alias_names] = process_field_aliases(empty_specs);
    
    REQUIRE(alias_map.empty());
    REQUIRE(alias_names.empty());
  }
}

} // namespace scream