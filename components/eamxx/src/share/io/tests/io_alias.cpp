#include <catch2/catch.hpp>

// Test file for EAMxx field aliasing functionality
// This file tests the field aliasing feature that allows users to specify
// alternative names for fields in output files using the syntax "ALIAS:=FIELDNAME"
//
// Tests include:
// 1. Parsing of alias specifications with various formats
// 2. Processing of mixed alias and non-alias field lists
// 3. Integration with OutputManager to create files with aliased names
// 4. Verification that aliased names appear correctly in NetCDF output files

#include "share/eamxx_types.hpp"
#include "share/field/field.hpp"
#include "share/field/field_manager.hpp"
#include "share/field/field_utils.hpp"
#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/io/eamxx_io_utils.hpp"
#include "share/io/eamxx_output_manager.hpp"
#include "share/io/eamxx_scorpio_interface.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/util/eamxx_time_stamp.hpp"

#include <ekat_assert.hpp>
#include <ekat_comm.hpp>
#include <ekat_parameter_list.hpp>
#include <ekat_units.hpp>

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace scream {

// Helper function to create test fields and field manager
std::shared_ptr<FieldManager>
create_test_field_manager(const std::shared_ptr<const AbstractGrid> &grid,
                          const util::TimeStamp &t0, int ncols, int nlevs) {

  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  auto fm = std::make_shared<FieldManager>(grid);

  // Create some test fields with realistic EAMxx names
  FieldIdentifier fid1("qv", {{COL, LEV}, {ncols, nlevs}}, kg / kg, grid->name());
  FieldIdentifier fid2("T_mid", {{COL, LEV}, {ncols, nlevs}}, K, grid->name());
  FieldIdentifier fid3("ps", {{COL}, {ncols}}, Pa, grid->name());

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
  qv.get_header().get_tracking().update_time_stamp(t0);
  T_mid.get_header().get_tracking().update_time_stamp(t0);
  ps.get_header().get_tracking().update_time_stamp(t0);

  fm->add_field(qv);
  fm->add_field(T_mid);
  fm->add_field(ps);

  return fm;
}

TEST_CASE("io_field_aliasing") {
  SECTION("parse_field_alias") {
    // Test normal case with alias
    auto [alias1, field1] = parse_field_alias("LWP:=LiqWaterPath");
    REQUIRE(alias1 == "LWP");
    REQUIRE(field1 == "LiqWaterPath");

    // Test with spaces
    auto [alias2, field2] = parse_field_alias(" SWP := SolidWaterPath ");
    REQUIRE(alias2 == "SWP");
    REQUIRE(field2 == "SolidWaterPath");

    // Test no alias
    auto [alias3, field3] = parse_field_alias("Temperature");
    REQUIRE(alias3 == "Temperature");
    REQUIRE(field3 == "Temperature");

    // Test complex field names
    auto [alias4, field4] = parse_field_alias("RH:=RelativeHumidity");
    REQUIRE(alias4 == "RH");
    REQUIRE(field4 == "RelativeHumidity");

    // Test error cases
    REQUIRE_THROWS(parse_field_alias("LWP:="));
    REQUIRE_THROWS(parse_field_alias(":=LiqWaterPath"));
    REQUIRE_THROWS(parse_field_alias("  :=  "));
  }

  SECTION("process_field_aliases") {
    std::vector<std::string> field_specs = {
        "T",                   // No alias
        "LWP:=LiqWaterPath",   // With alias
        "SWP:=SolidWaterPath", // With alias
        "qv",                  // No alias
        "RH:=RelativeHumidity" // With alias
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
        "LWP:=LiqWaterPath",
        "LWP:=SolidWaterPath" // Duplicate alias
    };
    REQUIRE_THROWS(process_field_aliases(duplicate_specs));
  }

  SECTION("mixed_aliases_and_regular_fields") {
    std::vector<std::string> mixed_specs = {"temperature", "LWP:=LiqWaterPath", "pressure",
                                            "RH:=RelativeHumidity"};

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

// Separate test case with PIO initialization to avoid re-initialization issues
TEST_CASE("output_aliases_file_verification", "[io][alias]") {
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
  util::TimeStamp t0({2023, 1, 1}, {0, 0, 0});
  auto fm = create_test_field_manager(grid, t0, ncols, nlevs);

  // Test mixed alias and non-alias field specifications with actual file writing
  SECTION("mixed_aliases_file_output") {
    ekat::ParameterList params;
    params.set<std::string>("filename_prefix", "io_alias_mixed");
    params.set<std::string>("averaging_type", "INSTANT");

    std::vector<std::string> mixed_specs = {
        "qv",          // No alias - use original name
        "TEMP:=T_mid", // With alias
        "PSURF:=ps"    // With alias
    };
    params.set("field_names", mixed_specs);

    auto &ctrl_pl = params.sublist("output_control");
    ctrl_pl.set<std::string>("frequency_units", "nsteps");
    ctrl_pl.set<int>("frequency", 1);
    ctrl_pl.set<bool>("save_grid_data", false);

    // Set max snapshots to ensure file gets written and closed
    params.set<int>("max_snapshots_per_file", 2);
    params.set<std::string>("floating_point_precision", "single");

    OutputManager om;
    om.initialize(comm, params, t0, false);
    om.setup(fm, gm->get_grid_names());

    // Run enough timesteps to fill the file and force it to be written
    auto dt = 1; // 1 second timestep for nsteps frequency
    auto t  = t0;
    const int nsteps =
        3; // Run 3 steps with frequency=1 to get 3 outputs, exceeding max_snapshots=2

    for (int n = 0; n < nsteps; ++n) {
      om.init_timestep(t, dt);
      t += dt;

      // Update field values to make them different each timestep
      // Get the fields from the field manager to ensure we're updating the right ones
      auto qv_field    = fm->get_field("qv");
      auto T_mid_field = fm->get_field("T_mid");
      auto ps_field    = fm->get_field("ps");

      qv_field.deep_copy(0.01 + n * 0.001);
      T_mid_field.deep_copy(280.0 + n * 1.0);
      ps_field.deep_copy(101325.0 + n * 100.0);

      // Update field timestamps
      qv_field.get_header().get_tracking().update_time_stamp(t);
      T_mid_field.get_header().get_tracking().update_time_stamp(t);
      ps_field.get_header().get_tracking().update_time_stamp(t);

      om.run(t);
    }

    // Check that the file was closed since we exceeded max_snapshots
    const auto &file_specs = om.output_file_specs();
    REQUIRE(not file_specs.is_open);

    om.finalize();

    // Verify the NetCDF file was created with the expected filename
    std::string expected_filename = "io_alias_mixed.INSTANT.nsteps_x1.np" +
                                    std::to_string(comm.size()) + "." + t0.to_string() + ".nc";

    // First check if the file actually exists
    std::ifstream file_check(expected_filename);
    if (file_check.good()) {
      file_check.close();
    } else {
      REQUIRE(false); // Fail the test if file doesn't exist
    }

    // Now verify the file was created and contains the correct aliased field names
    // Create a new field manager with fields using aliased names for reading
    auto read_fm = std::make_shared<FieldManager>(grid);

    // Create fields with aliased names for reading back
    FieldIdentifier qv_read_id("qv", {{COL, LEV}, {ncols, nlevs}}, kg / kg, grid->name());
    FieldIdentifier temp_read_id("TEMP", {{COL, LEV}, {ncols, nlevs}}, K, grid->name());
    FieldIdentifier psurf_read_id("PSURF", {{COL}, {ncols}}, Pa, grid->name());

    Field qv_read(qv_read_id);
    Field temp_read(temp_read_id);
    Field psurf_read(psurf_read_id);

    qv_read.allocate_view();
    temp_read.allocate_view();
    psurf_read.allocate_view();

    qv_read.get_header().get_tracking().update_time_stamp(t0);
    temp_read.get_header().get_tracking().update_time_stamp(t0);
    psurf_read.get_header().get_tracking().update_time_stamp(t0);

    read_fm->add_field(qv_read);
    read_fm->add_field(temp_read);
    read_fm->add_field(psurf_read);

    // Set up reader parameter list
    ekat::ParameterList reader_pl;
    reader_pl.set("filename", expected_filename);

    // Use the aliased names in the field list for reading
    std::vector<std::string> aliased_names = {"qv", "TEMP", "PSURF"};
    reader_pl.set("field_names", aliased_names);

    // Try to read the file - this will fail if the aliases weren't written correctly
    REQUIRE_NOTHROW([&]() {
      AtmosphereInput reader(reader_pl, read_fm);
      reader.read_variables(0); // Read first timestep

      // Verify that we can access the fields we just read
      auto qv_field_read    = read_fm->get_field("qv");
      auto temp_field_read  = read_fm->get_field("TEMP");
      auto psurf_field_read = read_fm->get_field("PSURF");

      // Check that fields have reasonable data
      REQUIRE(qv_field_read.get_header().get_alloc_properties().get_num_scalars() > 0);
      REQUIRE(temp_field_read.get_header().get_alloc_properties().get_num_scalars() > 0);
      REQUIRE(psurf_field_read.get_header().get_alloc_properties().get_num_scalars() > 0);

      // The fact that we reached here means:
      // 1. The NetCDF file was created successfully
      // 2. The file contains variables with the aliased names (qv, TEMP, PSURF)
      // 3. The aliases were properly applied during output
      // 4. We can successfully read back the data using the aliased names
    }());

    // Clean up PIO subsystem
    scorpio::finalize_subsystem();
  }
}

} // namespace scream
