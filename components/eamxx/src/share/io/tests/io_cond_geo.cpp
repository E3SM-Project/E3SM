#include <catch2/catch.hpp>

#include "share/io/eamxx_output_manager.hpp"

#include "share/data_managers/mesh_free_grids_manager.hpp"

#include "share/field/field_utils.hpp"
#include "share/field/field.hpp"
#include "share/data_managers/field_manager.hpp"

#include "share/util/eamxx_time_stamp.hpp"
#include "share/core/eamxx_types.hpp"

#include <ekat_units.hpp>
#include <ekat_parameter_list.hpp>
#include <ekat_assert.hpp>
#include <ekat_comm.hpp>

namespace scream {

namespace {

// Helper: build a simple grids manager with one "point_grid"
std::shared_ptr<const GridsManager>
get_gm (const ekat::Comm& comm)
{
  const int ngcols = std::max(static_cast<int>(comm.size())-1, 1);
  const int nlevs  = 4;
  auto gm = create_mesh_free_grids_manager(comm, 0, 0, nlevs, ngcols);
  gm->build_grids();
  return gm;
}

// Create a field manager with:
//   - f_nodim : layout (COL)         -- no swband dim
//   - f_swband: layout (COL, swband) -- has the swband dim
void
add_fields (const std::shared_ptr<const AbstractGrid>& grid,
            const std::shared_ptr<FieldManager>& fm,
            const util::TimeStamp& t0,
            const int nswbands)
{
  using FL  = FieldLayout;
  using FID = FieldIdentifier;
  using namespace ShortFieldTagsNames;

  const int nlcols = grid->get_num_local_dofs();

  // A plain field with no sw-band dimension
  {
    FL  fl({COL}, {nlcols});
    FID fid("f_nodim", fl, ekat::units::none, grid->name());
    Field f(fid);
    f.allocate_view();
    f.deep_copy(1.0);
    f.get_header().get_tracking().update_time_stamp(t0);
    fm->add_field(f);
  }

  // A field with the swband dimension (uses the grid's vector layout helper
  // so the dim name is properly set to "swband")
  {
    FL  fl = grid->get_2d_vector_layout(nswbands, "swband");
    FID fid("f_swband", fl, ekat::units::none, grid->name());
    Field f(fid);
    f.allocate_view();
    f.deep_copy(2.0);
    f.get_header().get_tracking().update_time_stamp(t0);
    fm->add_field(f);
  }
}

// Add conditional geo data (swband and swband_bounds) to the grid.
// These carry the "io_output_if_dim_exists" flag so they should only appear
// in output files where the "swband" dimension is already registered.
void
add_conditional_geo (const std::shared_ptr<const AbstractGrid>& grid,
                     const int nswbands)
{
  using namespace ShortFieldTagsNames;

  // swband: 1-D coord variable (dim name = "swband")
  FieldLayout bnd_layout({CMP}, {nswbands}, {"swband"});
  Field swband(FieldIdentifier("swband", bnd_layout, ekat::units::none, grid->name()));
  swband.allocate_view();
  swband.deep_copy(1.0);
  swband.get_header().set_extra_data("io_output_if_dim_exists", std::string("swband"));
  grid->set_geometry_data(swband);

  // swband_bounds: 2-D bounds variable (shape: nswbands x 2)
  FieldLayout bnds_layout({CMP}, {nswbands}, {"swband"});
  bnds_layout.append_dim(CMP, 2);
  Field swband_bounds(FieldIdentifier("swband_bounds", bnds_layout, ekat::units::none, grid->name()));
  swband_bounds.allocate_view();
  swband_bounds.deep_copy(0.5);
  swband_bounds.get_header().set_extra_data("io_output_if_dim_exists", std::string("swband"));
  grid->set_geometry_data(swband_bounds);
}

// Build an OutputManager parameter list for a single-field instant output
ekat::ParameterList
make_om_params (const std::string& prefix,
                const std::vector<std::string>& field_names,
                const bool save_grid_data)
{
  ekat::ParameterList pl;
  pl.set("filename_prefix", prefix);
  pl.set("field_names", field_names);
  pl.set("averaging_type", std::string("INSTANT"));
  auto& ctrl = pl.sublist("output_control");
  ctrl.set("frequency_units", std::string("nsteps"));
  ctrl.set("frequency", 1);
  ctrl.set("save_grid_data", save_grid_data);
  return pl;
}

} // anonymous namespace

TEST_CASE ("io_conditional_geo_data")
{
  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);

  const int nswbands = 3;
  auto t0 = util::TimeStamp({2023, 2, 17}, {0, 0, 0});
  const int dt = 1;

  auto gm   = get_gm(comm);
  auto grid = gm->get_grid("point_grid");

  // Add conditional geo data to the grid
  add_conditional_geo(grid, nswbands);

  // Build a field manager with both kinds of fields
  auto fm = std::make_shared<FieldManager>(grid);
  add_fields(grid, fm, t0, nswbands);

  // -----------------------------------------------------------------------
  // Case 1: output only the field WITHOUT the swband dim.
  //         The conditional geo data should NOT appear in the output file.
  // -----------------------------------------------------------------------
  {
    auto pl = make_om_params("io_cond_geo_no_swband",
                             {"f_nodim"}, /*save_grid_data=*/true);

    OutputManager om;
    om.initialize(comm, pl, t0, false);
    om.setup(fm, gm->get_grid_names());

    auto ts = t0;
    ts += dt;
    om.init_timestep(t0, dt);
    om.run(ts);
    om.finalize();

    // Verify: "swband" dim should NOT be in the file
    std::string fname = "io_cond_geo_no_swband.INSTANT.nsteps_x1.np"
                        + std::to_string(comm.size())
                        + "." + t0.to_string() + ".nc";
    scorpio::register_file(fname, scorpio::FileMode::Read);
    REQUIRE_FALSE (scorpio::has_dim(fname, "swband"));
    REQUIRE_FALSE (scorpio::has_var(fname, "swband"));
    REQUIRE_FALSE (scorpio::has_var(fname, "swband_bounds"));
    scorpio::release_file(fname);
  }

  // -----------------------------------------------------------------------
  // Case 2: output the field WITH the swband dim.
  //         The conditional geo data SHOULD appear in the output file.
  // -----------------------------------------------------------------------
  {
    auto pl = make_om_params("io_cond_geo_with_swband",
                             {"f_swband"}, /*save_grid_data=*/true);

    OutputManager om;
    om.initialize(comm, pl, t0, false);
    om.setup(fm, gm->get_grid_names());

    auto ts = t0;
    ts += dt;
    om.init_timestep(t0, dt);
    om.run(ts);
    om.finalize();

    // Verify: "swband" dim should be in the file and both geo vars present
    std::string fname = "io_cond_geo_with_swband.INSTANT.nsteps_x1.np"
                        + std::to_string(comm.size())
                        + "." + t0.to_string() + ".nc";
    scorpio::register_file(fname, scorpio::FileMode::Read);
    REQUIRE (scorpio::has_dim(fname, "swband"));
    REQUIRE (scorpio::has_var(fname, "swband"));
    REQUIRE (scorpio::has_var(fname, "swband_bounds"));
    scorpio::release_file(fname);
  }

  scorpio::finalize_subsystem();
}

} // namespace scream
