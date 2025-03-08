#include <catch2/catch.hpp>

#include "share/io/eamxx_output_manager.hpp"
#include "share/io/scorpio_input.hpp"

#include "share/grid/mesh_free_grids_manager.hpp"

#include "share/field/field_utils.hpp"
#include "share/field/field.hpp"
#include "share/field/field_manager.hpp"

#include "share/util/eamxx_universal_constants.hpp"
#include "share/util/eamxx_setup_random_test.hpp"
#include "share/util/eamxx_time_stamp.hpp"
#include "share/eamxx_types.hpp"

#include "ekat/util/ekat_units.hpp"
#include "ekat/ekat_parameter_list.hpp"
#include "ekat/ekat_assert.hpp"
#include "ekat/mpi/ekat_comm.hpp"
#include "ekat/util/ekat_test_utils.hpp"

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
        const util::TimeStamp& t0, const int seed)
{
  using FL  = FieldLayout;
  using FID = FieldIdentifier;
  using namespace ShortFieldTagsNames;

  // Random number generation stuff
  // NOTES
  //  - Use integers, so we can check answers without risk of
  //    non bfb diffs due to different order of sums.
  //  - Uniform_int_distribution returns an int, and the randomize
  //    util checks that return type matches the Field data type.
  //    So wrap the int pdf in a lambda, that does the cast.
  std::mt19937_64 engine(seed);
  auto my_pdf = [&](std::mt19937_64& engine) -> Real {
    std::uniform_int_distribution<int> pdf (0,100);
    Real v = pdf(engine);
    return v;
  };

  const int nlcols = grid->get_num_local_dofs();
  const int nlevs  = grid->get_num_vertical_levels();

  std::vector<FL> layouts =
  {
    FL({COL         }, {nlcols        }),
    FL({COL,     LEV}, {nlcols,  nlevs}),
    FL({COL,CMP,ILEV}, {nlcols,2,nlevs+1})
  };

  auto fm = std::make_shared<FieldManager>(grid);

  const auto units = ekat::units::Units::nondimensional();
  int count=0;
  for (const auto& fl : layouts) {
    FID fid("f_"+std::to_string(count),fl,units,grid->name());
    Field f(fid);
    f.allocate_view();
    randomize (f,engine,my_pdf);
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
  auto grid = gm->get_grid("Point Grid");

  // Time advance parameters
  auto t0 = get_t0();
  const int dt = 86400*30; // 30 days

  // Create some fields
  auto fm = get_fm(grid,t0,seed);
  std::vector<std::string> fnames;
  for (auto it : *fm) {
    fnames.push_back(it.second->name());
  }

  // Create output params
  ekat::ParameterList om_pl;
  om_pl.set("filename_prefix",std::string("io_monthly"));
  om_pl.set("Field Names",fnames);
  om_pl.set("Averaging Type", std::string("Instant"));
  om_pl.set("file_max_storage_type",std::string("one_month"));
  om_pl.set("Floating Point Precision",std::string("single"));
  auto& ctrl_pl = om_pl.sublist("output_control");
  ctrl_pl.set("frequency_units",std::string("nsteps"));
  ctrl_pl.set("Frequency",1);
  ctrl_pl.set("save_grid_data",false);

  // Create Output manager
  OutputManager om;
  om.initialize(comm,om_pl,t0,false);
  om.setup(fm,gm);

  // Time loop: do 11 steps, since we already did Jan output at t0
  const int nsteps = 11;
  auto t = t0;
  for (int n=0; n<nsteps; ++n) {
    om.init_timestep(t,dt);
    // Update time
    t += dt;

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
  auto grid = gm->get_grid("Point Grid");

  // Get initial fields. Use wrong seed for fm, so fields are not
  // inited with right data (avoid getting right answer without reading).
  auto fm0 = get_fm(grid,t0,seed);
  auto fm  = get_fm(grid,t0,-seed-1);
  std::vector<std::string> fnames;
  for (auto it : *fm) {
    fnames.push_back(it.second->name());
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

  // Create reader pl
  ekat::ParameterList reader_pl;
  reader_pl.set("Field Names",fnames);

  for (int n=0; n<12; ++n) {
    auto t = t0 + n*dt;
    auto filename = get_filename(t);

    // There should be just one time snapshot per file
    REQUIRE(scorpio::get_dimlen(filename,"time")==1);

    reader_pl.set("Filename",filename);
    AtmosphereInput reader(reader_pl,fm);
    reader.read_variables();

    for (const auto& fn : fnames) {
      auto f0 = fm0->get_field(fn).clone();
      auto f  = fm->get_field(fn);
      add(f0,n);
      REQUIRE (views_are_equal(f,f0));
    }
  }
}

TEST_CASE ("io_monthly") {
  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);

  auto seed = get_random_test_seed(&comm);

  if (comm.am_i_root()) {
    std::cout << "   -> Testing output with one file per month ...\n";
  }
  write(seed,comm);
  read (seed,comm);
  if (comm.am_i_root()) {
    std::cout << "   -> Testing output with one file per month ... PASS\n";
  }
  scorpio::finalize_subsystem();
}

} // anonymous namespace
