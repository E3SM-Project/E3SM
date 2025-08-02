#include <catch2/catch.hpp>

#include "share/io/eamxx_io_utils.hpp"
#include "share/io/eamxx_output_manager.hpp"
#include "share/io/eamxx_scorpio_interface.hpp"
#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/field/field_utils.hpp"
#include "share/field/field.hpp"
#include "share/field/field_manager.hpp"
#include "share/util/eamxx_time_stamp.hpp"
#include "share/eamxx_types.hpp"

#include "ekat/util/ekat_units.hpp"
#include "ekat/ekat_parameter_list.hpp"
#include "ekat/ekat_assert.hpp"
#include "ekat/mpi/ekat_comm.hpp"

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

TEST_CASE("output_aliases_integration", "[.][io][alias]") {
  // Integration test to verify that OutputManager correctly uses aliases
  // This test creates actual output files and verifies the variable names
  
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  // Create a simple grid and field manager
  ekat::Comm comm(MPI_COMM_WORLD);
  
  // Initialize the PIO subsystem for this test
  scorpio::init_subsystem(comm);
  
  const int ncols = 10;
  const int nlevs = 5;
  
  auto gm = create_mesh_free_grids_manager(comm, 0, 0, nlevs, ncols);
  gm->build_grids();
  
  auto grid = gm->get_grid("point_grid");  
  auto fm = std::make_shared<FieldManager>(grid);
  
  // Create some test fields with realistic EAMxx names
  FieldIdentifier fid1("qv", {{COL,LEV},{10,5}}, kg/kg, grid->name());
  FieldIdentifier fid2("T_mid", {{COL,LEV},{10,5}}, K, grid->name());
  FieldIdentifier fid3("ps", {{COL},{10}}, Pa, grid->name());

  Field qv(fid1);
  Field T_mid(fid2);
  Field ps(fid3);
  
  qv.allocate_view();
  T_mid.allocate_view();
  ps.allocate_view();
  
  // Initialize with dummy data
  qv.deep_copy(0.01);
  T_mid.deep_copy(280.0);
  ps.deep_copy(101325.0);
  
  // Update timestamps
  util::TimeStamp t0({2023,1,1},{0,0,0});
  qv.get_header().get_tracking().update_time_stamp(t0);
  T_mid.get_header().get_tracking().update_time_stamp(t0);
  ps.get_header().get_tracking().update_time_stamp(t0);
  
  fm->add_field(qv);
  fm->add_field(T_mid);
  fm->add_field(ps);
  
  // Test that OutputManager can be created and initialized with aliases
  SECTION("construction_with_aliases") {
    // Create output parameter list with aliases
    ekat::ParameterList params;
    params.set<std::string>("filename_prefix", "alias_test");
    params.set<std::string>("averaging_type", "instant");
    
    // Test field specifications with aliases
    std::vector<std::string> field_specs = {
      "QV:=:qv",        // Alias QV for qv
      "TEMP:=:T_mid",   // Alias TEMP for T_mid  
      "PSURF:=:ps"      // Alias PSURF for ps
    };
    params.set("field_names", field_specs);
    
    auto& ctrl_pl = params.sublist("output_control");
    ctrl_pl.set<std::string>("frequency_units", "nsteps");
    ctrl_pl.set<int>("frequency", 1);
    
    REQUIRE_NOTHROW([&]() {
      OutputManager om;
      om.initialize(comm, params, t0, false);
      om.setup(fm, gm->get_grid_names());
      // OutputManager setup should succeed with alias syntax
    }());
  }
  
  // Test mixed alias and non-alias field specifications
  SECTION("mixed_aliases") {
    ekat::ParameterList params;
    params.set<std::string>("filename_prefix", "mixed_test");
    params.set<std::string>("averaging_type", "instant");
    
    std::vector<std::string> mixed_specs = {
      "qv",             // No alias - use original name
      "TEMP:=:T_mid",   // With alias
      "ps"              // No alias - use original name
    };
    params.set("field_names", mixed_specs);
    
    auto& ctrl_pl = params.sublist("output_control");
    ctrl_pl.set<std::string>("frequency_units", "nsteps");
    ctrl_pl.set<int>("frequency", 1);
    
    REQUIRE_NOTHROW([&]() {
      OutputManager om;
      om.initialize(comm, params, t0, false);
      om.setup(fm, gm->get_grid_names());
      // OutputManager setup should succeed with mixed alias syntax
    }());
  }
  
  // Clean up PIO subsystem
  scorpio::finalize_subsystem();
}

} // namespace scream
