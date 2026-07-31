#include <catch2/catch.hpp>

#include "share/io/eamxx_output_manager.hpp"

#include "share/data_managers/mesh_free_grids_manager.hpp"

#include "share/field/field_reader.hpp"
#include "share/field/field_utils.hpp"
#include "share/field/field.hpp"
#include "share/data_managers/field_manager.hpp"

#include "share/util/eamxx_universal_constants.hpp"
#include "share/core/eamxx_setup_random_test.hpp"
#include "share/util/eamxx_time_stamp.hpp"
#include "share/core/eamxx_types.hpp"

#include <ekat_units.hpp>
#include <ekat_parameter_list.hpp>
#include <ekat_assert.hpp>
#include <ekat_comm.hpp>

#include <iomanip>
#include <memory>

namespace scream {

void add (const Field& f, const double v) {
  auto data = f.get_internal_view_data<Real,Host>();
  auto nscalars = f.get_header().get_alloc_properties().get_num_scalars();
  for (int i=0; i<nscalars; ++i) {
    data[i] += v;
  }
  f.sync_to_dev();
}

util::TimeStamp get_t0 () {
  return util::TimeStamp({2000,1,15},{0,0,0});
}

std::shared_ptr<const GridsManager>
get_gm (const ekat::Comm& comm)
{
  // For 2+ ranks tests, this will check IO works correctly
  // even if one rank owns 0 dofs
  const int ngcols = std::max(comm.size()-1,1);
  const int nlevs = 4;
  auto gm = create_mesh_free_grids_manager(comm,0,0,nlevs,ngcols);
  gm->build_grids();
  return gm;
}

std::shared_ptr<FieldManager>
get_fm (const std::shared_ptr<const AbstractGrid>& grid,
        const util::TimeStamp& t0, int seed)
{
  using FL  = FieldLayout;
  using FID = FieldIdentifier;
  using namespace ShortFieldTagsNames;

  // Note: we use a discrete set of random values, so we can
  // check answers without risk of non-bfb diffs due to ops order
  std::vector<Real> values;
  for (int i=0; i<=100; ++i)
    values.push_back(static_cast<Real>(i));

  const int nlcols = grid->get_num_local_dofs();
  const int nlevs  = grid->get_num_vertical_levels();

  std::vector<FL> layouts =
  {
    FL({COL         }, {nlcols        }),
    FL({COL,     LEV}, {nlcols,  nlevs}),
    FL({COL,CMP,ILEV}, {nlcols,2,nlevs+1})
  };

  auto fm = std::make_shared<FieldManager>(grid);

  int count=0;
  for (const auto& fl : layouts) {
    FID fid("f_"+std::to_string(count),fl,ekat::units::none,grid->name());
    Field f(fid);
    f.allocate_view();
    randomize_discrete (f,seed++,values);
    f.get_header().get_tracking().update_time_stamp(t0);
    fm->add_field(f);
    ++count;
  }

  return fm;
}

// Returns fields after initialization
void write (const int seed, const ekat::Comm& comm)
{
  // Create grid
  auto gm = get_gm(comm);
  auto grid = gm->get_grid("point_grid");

  // Time advance parameters
  auto t0 = get_t0();
  const int dt = 86400*30; // 30 days

  // Create some fields
  auto fm = get_fm(grid,t0,seed);
  std::vector<std::string> fnames;
  for (auto it : fm->get_repo()) {
    fnames.push_back(it.second->name());
  }

  // Create output params
  ekat::ParameterList om_pl;
  om_pl.set("filename_prefix",std::string("io_monthly"));
  om_pl.set("field_names",fnames);
  om_pl.set("averaging_type", std::string("instant"));
  om_pl.set("file_max_storage_type",std::string("one_month"));
  om_pl.set("floating_point_precision",std::string("single"));
  auto& ctrl_pl = om_pl.sublist("output_control");
  ctrl_pl.set("frequency_units",std::string("nsteps"));
  ctrl_pl.set("frequency",1);
  ctrl_pl.set("save_grid_data",false);

  // Create Output manager
  OutputManager om;
  om.initialize(comm,om_pl,t0,false);
  om.setup(fm,gm->get_grid_names());

  // Time loop: do 11 steps, since we already did Jan output at t0
  const int nsteps = 11;
  auto t = t0;
  for (int n=0; n<nsteps; ++n) {
    // Update time
    t += dt;

    om.init_timestep(t,dt);

    // Add 1 to all fields entries
    for (const auto& name : fnames) {
      auto f = fm->get_field(name);
      add(f,1);
    }

    // Run output manager
    om.run (t);
  }

  // Close file and cleanup
  om.finalize();
}

void read (const int seed, const ekat::Comm& comm)
{
  // Time quantities
  auto t0 = get_t0();
  int dt = 86400*30;

  // Get gm
  auto gm = get_gm (comm);
  auto grid = gm->get_grid("point_grid");
  auto gids = grid->get_partitioned_dim_gids();

  // Get initial fields. Use wrong seed for fm, so fields are not
  // inited with right data (avoid getting right answer without reading).
  auto fm0 = get_fm(grid,t0,seed);
  std::vector<Field> fields;
  for (auto it : fm0->get_repo()) {
    fields.push_back(it.second->clone());
  }

  // Get filename from timestamp
  std::string casename = "io_monthly";
  auto get_filename = [&](const util::TimeStamp& t) {
    auto t_str = t.to_string().substr(0,7);
    std::string fname = casename
                      + ".INSTANT.nsteps_x1"
                      + ".np" + std::to_string(comm.size())
                      + "." + t_str
                      + ".nc";
    return fname;
  };

  for (int n=0; n<12; ++n) {
    auto t = t0 + n*dt;
    auto filename = get_filename(t);

    // There should be just one time snapshot per file
    REQUIRE(scorpio::get_dimlen(filename,"time")==1);

    read_fields (filename,fields,gids,comm);

    for (const auto& f : fields) {
      auto f0 = fm0->get_field(f.name()).clone(CloneFlags::CopyData);
      add(f0,n);
      REQUIRE (views_are_equal(f,f0));
    }
  }
}

// Write 12 months of data using the given averaging type with monthly storage.
// This tests that the timestamp index in the file is taken from last_write_ts
// (start of the averaging window) rather than next_write_ts (end of window)
// for non-Instant averaging types.
void write_avg_type (const std::string& avg_type, const int seed, const ekat::Comm& comm)
{
  auto gm = get_gm(comm);
  auto grid = gm->get_grid("point_grid");

  auto t0 = get_t0();
  const int dt = 86400*30; // 30 days

  auto fm = get_fm(grid,t0,seed);
  std::vector<std::string> fnames;
  for (auto it : fm->get_repo()) {
    fnames.push_back(it.second->name());
  }

  ekat::ParameterList om_pl;
  om_pl.set("filename_prefix",std::string("io_monthly_avg"));
  om_pl.set("field_names",fnames);
  om_pl.set("averaging_type", avg_type);
  om_pl.set("file_max_storage_type",std::string("one_month"));
  om_pl.set("floating_point_precision",std::string("single"));
  auto& ctrl_pl = om_pl.sublist("output_control");
  ctrl_pl.set("frequency_units",std::string("nsteps"));
  ctrl_pl.set("frequency",1);
  ctrl_pl.set("save_grid_data",false);

  OutputManager om;
  om.initialize(comm,om_pl,t0,false);
  om.setup(fm,gm->get_grid_names());

  // Run 12 steps: one per month
  const int nsteps = 12;
  auto t = t0;
  for (int n=0; n<nsteps; ++n) {
    om.init_timestep(t,dt);
    t += dt;
    for (const auto& name : fnames) {
      auto f = fm->get_field(name);
      add(f,1);
    }
    om.run(t);
  }
  om.finalize();
}

// Verify that for AVERAGE/MIN/MAX averaging with monthly storage:
//  - 12 separate files were created (one per month)
//  - each file contains exactly 1 time snapshot
//  - each file contains the correct averaged/min/max data
void read_avg_type (const std::string& avg_type, const int seed, const ekat::Comm& comm)
{
  auto t0 = get_t0();
  int dt = 86400*30;

  auto gm = get_gm(comm);
  auto grid = gm->get_grid("point_grid");
  auto gids = grid->get_partitioned_dim_gids();

  auto fm0 = get_fm(grid,t0,seed);
  std::vector<Field> fields;
  for (auto it : fm0->get_repo()) {
    fields.push_back(it.second->clone());
  }

  // For non-Instant averaging, the filename is derived from last_write_ts
  // (the start of the averaging window). With freq=1 nstep, the n-th window
  // starts at t0 + n*dt and ends at t0 + (n+1)*dt.
  // So the n-th output file is named using the timestamp t0 + n*dt.
  std::string casename = "io_monthly_avg";
  auto get_filename = [&](const util::TimeStamp& t) {
    auto t_str = t.to_string().substr(0,7); // YYYY-MM
    std::string fname = casename
                      + "." + avg_type + ".nsteps_x1"
                      + ".np" + std::to_string(comm.size())
                      + "." + t_str
                      + ".nc";
    return fname;
  };

  for (int n=0; n<12; ++n) {
    auto window_start = t0 + n*dt;
    auto filename = get_filename(window_start);

    // Each file must hold exactly one snapshot (one per month).
    // Before the bug fix, setup_file used next_write_ts (end of window) instead
    // of last_write_ts (start of window), so the file's time_idx was set to the
    // wrong month. This caused snapshot_fits to return true for an extra month,
    // leaving the file open and allowing multiple snapshots to accumulate.
    REQUIRE(scorpio::get_dimlen(filename,"time")==1);

    read_fields(filename,fields,gids,comm);

    // With 1 sample per averaging window, AVERAGE=MIN=MAX equals the single sample.
    // At window n (0-indexed), the field was incremented n+1 times before the write.
    for (const auto& f : fields) {
      auto f0 = fm0->get_field(f.name()).clone(CloneFlags::CopyData);
      add(f0,n+1);
      REQUIRE(views_are_equal(f,f0));
    }
  }
}

TEST_CASE ("io_monthly") {
  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);

  auto seed = get_random_test_seed(&comm);

  if (comm.am_i_root()) {
    std::cout << "   -> Testing monthly output with INSTANT averaging ...\n";
  }
  write(seed,comm);
  read (seed,comm);
  if (comm.am_i_root()) {
    std::cout << "   -> Testing monthly output with INSTANT averaging ... PASS\n";
  }

  // Test that for AVERAGE/MIN/MAX averaging with monthly storage, each month's
  // output goes into a separate correctly-named file (one snapshot per file).
  // This tests the fix for the timestamp handling bug in setup_file, where
  // set_time_idx now uses last_write_ts instead of next_write_ts for non-Instant
  // averaging types.
  for (const auto& avg_type : {"AVERAGE","MIN","MAX"}) {
    if (comm.am_i_root()) {
      std::cout << "   -> Testing monthly output with " << avg_type << " averaging ...\n";
    }
    write_avg_type(avg_type,seed,comm);
    read_avg_type (avg_type,seed,comm);
    if (comm.am_i_root()) {
      std::cout << "   -> Testing monthly output with " << avg_type << " averaging ... PASS\n";
    }
  }

  scorpio::finalize_subsystem();
}

} // anonymous namespace
