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

constexpr int num_output_steps = 5;

void add (const Field& f, const double v) {
  auto data = f.get_internal_view_data<Real,Host>();
  auto nscalars = f.get_header().get_alloc_properties().get_num_scalars();
  for (int i=0; i<nscalars; ++i) {
    data[i] += v;
  }
  f.sync_to_dev();
}

int get_dt (const std::string& freq_units) {
  int dt;
  if (freq_units=="nsteps") {
    dt = 1;
  } else if (freq_units=="nsecs") {
    dt = 1;
  } else if (freq_units=="nmins") {
    dt = 60;
  } else if (freq_units=="nhours") {
    dt = 60*60;
  } else if (freq_units=="ndays") {
    dt = 60*60*24;
  } else {
    EKAT_ERROR_MSG ("Error! Unsupported freq_units\n"
        " - freq_units: " + freq_units + "\n"
        " - valid units: nsteps, nsecs, nmins, nhours, ndays\n");
  }
  return dt;
}

util::TimeStamp get_t0 () {
  return util::TimeStamp({2023,2,17},{0,0,0});
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
  using stratts_t = std::map<std::string,std::string>;
  for (const auto& fl : layouts) {
    FID fid("f_"+std::to_string(count),fl,units,grid->name());
    Field f(fid);
    f.allocate_view();
    auto& str_atts = f.get_header().get_extra_data<stratts_t>("io: string attributes");
    str_atts["test"] = f.name();
    randomize (f,engine,my_pdf);
    f.get_header().get_tracking().update_time_stamp(t0);
    fm->add_field(f);
    ++count;
  }

  return fm;
}

// Returns fields after initialization
void write (const std::string& avg_type, const std::string& freq_units,
            const int freq, const int seed, const ekat::Comm& comm)
{
  // Create grid
  auto gm = get_gm(comm);
  auto grid = gm->get_grid("Point Grid");

  // Time advance parameters
  auto t0 = get_t0();
  const int dt = get_dt(freq_units);

  // Create some fields
  auto fm = get_fm(grid,t0,seed);
  std::vector<std::string> fnames;
  for (auto it : *fm) {
    fnames.push_back(it.second->name());
  }

  // Create output params
  ekat::ParameterList om_pl;
  om_pl.set("filename_prefix",std::string("io_basic"));
  om_pl.set("Field Names",fnames);
  om_pl.set("Averaging Type", avg_type);
  auto& ctrl_pl = om_pl.sublist("output_control");
  ctrl_pl.set("frequency_units",freq_units);
  ctrl_pl.set("Frequency",freq);
  ctrl_pl.set("save_grid_data",false);

  // While setting this is in practice irrelevant (we would close
  // the file anyways at the end of the run), we can test that the OM closes
  // the file AS SOON as it's full (before calling finalize)
  int max_snaps = num_output_steps;
  if (avg_type=="INSTANT") {
    ++max_snaps;
  }
  om_pl.set("Max Snapshots Per File", max_snaps);

  // Create Output manager
  OutputManager om;

  // Attempt to use invalid fp precision string
  om_pl.set("Floating Point Precision",std::string("triple"));
  om.initialize(comm,om_pl,t0,false);
  REQUIRE_THROWS (om.setup(fm,gm));
  om_pl.set("Floating Point Precision",std::string("single"));
  om.initialize(comm,om_pl,t0,false);
  om.setup(fm,gm);

  // Time loop: ensure we always hit 3 output steps
  const int nsteps = num_output_steps*freq;
  auto t = t0;
  for (int n=0; n<nsteps; ++n) {
    om.init_timestep(t,dt);
    // Update time
    t += dt;

    // Add 1 to all fields entries
    for (const auto& name : fnames) {
      auto f = fm->get_field(name);
      add(f,1.0);
    }

    // Run output manager
    om.run (t);
  }

  // Check that the file was closed, since we reached full capacity
  const auto& file_specs = om.output_file_specs();
  REQUIRE (not file_specs.is_open);

  // Close file and cleanup
  om.finalize();
}

void read (const std::string& avg_type, const std::string& freq_units,
           const int freq, const int seed, const ekat::Comm& comm)
{
  // Only INSTANT writes at t=0
  bool instant = avg_type=="INSTANT";

  // Time quantities
  auto t0 = get_t0();
  int num_writes = num_output_steps + (instant ? 1 : 0);

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

  // Create reader pl
  ekat::ParameterList reader_pl;
  std::string casename = "io_basic";
  auto filename = casename
    + "." + avg_type
    + "." + freq_units
    + "_x" + std::to_string(freq)
    + ".np" + std::to_string(comm.size())
    + "." + t0.to_string()
    + ".nc";
  reader_pl.set("Filename",filename);
  reader_pl.set("Field Names",fnames);
  AtmosphereInput reader(reader_pl,fm);

  // We added 1.0 to the input fields for each timestep
  // Hence, at output step N, we should get
  //  avg=INSTANT: output = f(N) = f(0) + N*freq
  //  avg=MAX:     output = f(N) = f(0) + N*freq
  //  avg=MIN:     output = f(N*freq+dt)
  //  avg=AVERAGE: output = f(0) + N*freq + (freq+1)/2
  // The last one comes from
  //   (a+1 + a+2 +..+a+freq)/freq =
  //   a + sum(i)/freq = a + (freq(freq+1)/2)/freq
  //   = a + (freq+1)/2
  double delta = (freq+1)/2.0;

  for (int n=0; n<num_writes; ++n) {
    reader.read_variables(n);
    for (const auto& fn : fnames) {
      auto f0 = fm0->get_field(fn).clone();
      auto f  = fm->get_field(fn);
      if (avg_type=="MIN") {
        // The 1st snap in the avg window (the smallest)
        // is one past window_start=n*freq
        add(f0,n*freq+1);
        REQUIRE (views_are_equal(f,f0));
      } else if (avg_type=="MAX") {
        add(f0,(n+1)*freq);
        REQUIRE (views_are_equal(f,f0));
      } else if (avg_type=="INSTANT") {
        add(f0,n*freq);
        REQUIRE (views_are_equal(f,f0));
      } else {
        add(f0,n*freq+delta);
        REQUIRE (views_are_equal(f,f0));
      }
    }
  }

  // Check that the expected metadata was appropriately set for each variable
  for (const auto& fn: fnames) {
    auto att_fill = scorpio::get_attribute<float>(filename,fn,"_FillValue");
    REQUIRE(att_fill==constants::DefaultFillValue<float>().value);

    auto att_str = scorpio::get_attribute<std::string>(filename,fn,"test");
    REQUIRE (att_str==fn);
  }
}

TEST_CASE ("io_basic") {
  std::vector<std::string> freq_units = {
    "nsteps",
    "nsecs",
    "nmins",
    "nhours",
    "ndays"
  };
  std::vector<std::string> avg_type = {
    "INSTANT",
    "MAX",
    "MIN",
    "AVERAGE"
  };

  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);

  auto seed = get_random_test_seed(&comm);

  const int freq = 5;
  auto print = [&] (const std::string& s, int line_len = -1) {
    if (comm.am_i_root()) {
      if (line_len<0) {
        std::cout << s;
      } else {
        std::cout << std::left << std::setw(line_len) << std::setfill('.') << s;
      }
    }
  };

  for (const auto& units : freq_units) {
    print ("-> Output frequency: " + units + "\n");
    for (const auto& avg : avg_type) {
      print("   -> Averaging type: " + avg + " ", 40);
      write(avg,units,freq,seed,comm);
      read (avg,units,freq,seed,comm);
      print(" PASS\n");
    }
  }
  scorpio::finalize_subsystem();
}

} // anonymous namespace
