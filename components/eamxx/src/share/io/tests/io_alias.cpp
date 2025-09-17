#include <catch2/catch.hpp>

#include "share/core/eamxx_types.hpp"
#include "share/field/field.hpp"
#include "share/manager/field_manager.hpp"
#include "share/field/field_utils.hpp"
#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/io/eamxx_io_utils.hpp"
#include "share/io/eamxx_output_manager.hpp"
#include "share/io/eamxx_scorpio_interface.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/util/eamxx_time_stamp.hpp"
#include "diagnostics/register_diagnostics.hpp"

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

Field col_iota (const std::shared_ptr<const AbstractGrid> &grid)
{
  const auto nondim = ekat::units::Units::nondimensional();
  const auto layout = grid->get_vertical_layout(true);
  const int nlevs = grid->get_num_vertical_levels();

  // A helper field to set values depending on lev
  FieldIdentifier fid1("col", layout, nondim, grid->name());
  Field col(fid1);
  col.allocate_view();
  auto col_h = col.get_view<Real*,Host>();
  for (int i=0; i<nlevs; ++i) {
    col_h[i] = i+1;
  }
  col.sync_to_dev();
  return col;
}

void calc_fields (Field& qv, Field& T_mid, Field& ps, const Field& col, int n)
{
  using namespace ShortFieldTagsNames;
  int ncols = ps.get_header().get_identifier().get_layout().dim(0);
  qv.deep_copy(1);
  T_mid.deep_copy(1);
  for (int icol=0; icol<ncols; ++icol) {
    qv.subfield(COL,icol).update(col,icol,n);
    T_mid.subfield(COL,icol).update(col,-icol,n);
    ps.subfield(COL,icol).deep_copy(icol+n);
  }
}

// Helper function to create test fields and field manager
std::shared_ptr<FieldManager>
create_test_field_manager(const std::shared_ptr<const AbstractGrid> &grid,
                          const util::TimeStamp &t0)
{

  using namespace ekat::units;

  const auto layout2d = grid->get_2d_scalar_layout();
  const auto layout3d = grid->get_3d_scalar_layout(true);
  auto fm = std::make_shared<FieldManager>(grid);

  // Create some test fields with realistic EAMxx names
  FieldIdentifier fid1("qv", layout3d, kg / kg, grid->name());
  FieldIdentifier fid2("T_mid", layout3d, K, grid->name());
  FieldIdentifier fid3("ps", layout2d, Pa, grid->name());

  Field qv(fid1);
  Field T_mid(fid2);
  Field ps(fid3);

  qv.allocate_view();
  T_mid.allocate_view();
  ps.allocate_view();

  // Update timestamps
  qv.get_header().get_tracking().update_time_stamp(t0);
  T_mid.get_header().get_tracking().update_time_stamp(t0);
  ps.get_header().get_tracking().update_time_stamp(t0);

  fm->add_field(qv);
  fm->add_field(T_mid);
  fm->add_field(ps);

  // Initialize fields
  calc_fields(qv,T_mid,ps,col_iota(grid),0);

  return fm;
}

// Separate test case with PIO initialization to avoid re-initialization issues
TEST_CASE("io_with_aliases") {
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;
  using strvec_t = std::vector<std::string>;
  
  // This test uses diagnostics, so make sure they are registered
  register_diagnostics();

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
  auto fm = create_test_field_manager(grid, t0);

  ekat::ParameterList params;
  params.set<std::string>("filename_prefix", "io_alias_mixed");
  params.set<std::string>("averaging_type", "INSTANT");

  std::vector<std::string> mixed_specs = {
      "qv",          // No alias - use original name
      "TEMP:=T_mid", // With alias
      "PSURF:=ps",   // With alias
      "surf_temp:=TEMP_at_model_bot"
  };
  params.set("field_names", mixed_specs);

  auto &ctrl_pl = params.sublist("output_control");
  ctrl_pl.set<std::string>("frequency_units", "nsteps");
  ctrl_pl.set<int>("frequency", 1);
  ctrl_pl.set<bool>("save_grid_data", false);

  // Set max snapshots to ensure file gets written and closed
  params.set<std::string>("floating_point_precision", "single");

  OutputManager om;
  om.initialize(comm, params, t0, false);
  om.setup(fm, gm->get_grid_names());

  // Run enough timesteps to fill the file and force it to be written
  auto dt = 1; // 1 second timestep for nsteps frequency
  auto t  = t0;
  const int nsteps =
      3; // Run 3 steps with frequency=1 to get 3 outputs, exceeding max_snapshots=2

  auto col = col_iota(fm->get_grid());
  for (int n = 1; n <= nsteps; ++n) {
    om.init_timestep(t, dt);
    t += dt;

    // Update field values to make them different each timestep
    // Get the fields from the field manager to ensure we're updating the right ones
    auto qv_field    = fm->get_field("qv");
    auto T_mid_field = fm->get_field("T_mid");
    auto ps_field    = fm->get_field("ps");

    calc_fields (qv_field, T_mid_field, ps_field, col, n);

    // Update field timestamps
    qv_field.get_header().get_tracking().update_time_stamp(t);
    T_mid_field.get_header().get_tracking().update_time_stamp(t);
    ps_field.get_header().get_tracking().update_time_stamp(t);

    om.run(t);
  }

  om.finalize();

  // Verify the NetCDF file was created with the expected filename
  std::string expected_filename = "io_alias_mixed.INSTANT.nsteps_x1.np" +
                                  std::to_string(comm.size()) + "." + t0.to_string() + ".nc";

  // First check if the file actually exists
  std::ifstream file_check(expected_filename);
  REQUIRE (file_check.good());
  file_check.close();

  // Now verify the file was created and contains the correct aliased field names
  // Create a new field manager with fields using aliased names for reading
  auto read_fm = std::make_shared<FieldManager>(grid);

  // Create fields with aliased names for reading back
  FieldIdentifier qv_read_id("qv", {{COL, LEV}, {ncols, nlevs}}, kg / kg, grid->name());
  FieldIdentifier temp_read_id("TEMP", {{COL, LEV}, {ncols, nlevs}}, K, grid->name());
  FieldIdentifier psurf_read_id("PSURF", {{COL}, {ncols}}, Pa, grid->name());
  FieldIdentifier surf_temp_id("surf_temp", {{COL}, {ncols}}, K, grid->name());

  Field qv_read(qv_read_id);
  Field temp_read(temp_read_id);
  Field psurf_read(psurf_read_id);
  Field surf_temp_read(surf_temp_id);

  qv_read.allocate_view();
  temp_read.allocate_view();
  psurf_read.allocate_view();
  surf_temp_read.allocate_view();

  qv_read.get_header().get_tracking().update_time_stamp(t0);
  temp_read.get_header().get_tracking().update_time_stamp(t0);
  psurf_read.get_header().get_tracking().update_time_stamp(t0);
  surf_temp_read.get_header().get_tracking().update_time_stamp(t0);

  read_fm->add_field(qv_read);
  read_fm->add_field(temp_read);
  read_fm->add_field(psurf_read);
  read_fm->add_field(surf_temp_read);

  // Set up reader parameter list
  ekat::ParameterList reader_pl;
  reader_pl.set("filename", expected_filename);

  // Use the aliased names in the field list for reading
  reader_pl.set<strvec_t>("field_names", {"qv", "TEMP", "PSURF", "surf_temp"});

  // Try to read the file - this will fail if the aliases weren't written correctly
  AtmosphereInput reader(reader_pl, read_fm);

  reader.set_logger(console_logger(ekat::logger::LogLevel::trace));
  auto qv_field_read        = read_fm->get_field("qv");
  auto temp_field_read      = read_fm->get_field("TEMP");
  auto psurf_field_read     = read_fm->get_field("PSURF");
  auto surf_temp_field_read = read_fm->get_field("surf_temp");

  auto qv_tgt    = qv_field_read.clone();
  auto T_mid_tgt = temp_field_read.clone();
  auto ps_tgt    = psurf_field_read.clone();
  auto T_mid_bot_tgt = T_mid_tgt.subfield(LEV,nlevs-1);
  for (int n=0; n<=nsteps; ++n) {
    reader.read_variables(n); // Read n-th timestamp

    calc_fields (qv_tgt, T_mid_tgt, ps_tgt, col, n);

    if (not views_are_equal(qv_tgt,qv_field_read)) {
      print_field_hyperslab(qv_tgt.alias("qv_tgt"));
      print_field_hyperslab(qv_field_read.alias("qv_field_read"));
    }
    REQUIRE (views_are_equal(qv_tgt,qv_field_read));
    REQUIRE (views_are_equal(T_mid_tgt,temp_field_read));
    REQUIRE (views_are_equal(ps_tgt,psurf_field_read));
    if (not views_are_equal(T_mid_bot_tgt,surf_temp_field_read)) {
      print_field_hyperslab(T_mid_bot_tgt.alias("T_mid_bot_tgt"));
      print_field_hyperslab(surf_temp_read.alias("surf_temp_read"));
    }
    REQUIRE (views_are_equal(T_mid_bot_tgt,surf_temp_field_read));
  }
  reader.finalize();

  // Clean up PIO subsystem
  scorpio::finalize_subsystem();
}

} // namespace scream
