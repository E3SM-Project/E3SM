#include <catch2/catch.hpp>

#include "share/io/eamxx_output_manager.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/io/eamxx_io_utils.hpp"

#include "share/data_managers/mesh_free_grids_manager.hpp"

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

constexpr int num_output_steps = 5;
constexpr Real fill_value = constants::fill_value<Real>;
constexpr Real fill_threshold = 0.5;

int get_dt (const std::string& freq_units) {
  int dt = 0;
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
  const int nlcols = 3;
  const int nlevs = 4;
  const int ngcols = nlcols*comm.size();
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
  for (const auto& fl : layouts) {
    FID fid("f_"+std::to_string(fl.size()),fl,units,grid->name());
    Field f(fid,true);
    f.deep_copy(0.0);
    f.get_header().get_tracking().update_time_stamp(t0);
    f.create_mask();

    fm->add_field(f);
  }

  return fm;
}

// Returns fields after initialization
void write (const std::string& avg_type, const std::string& freq_units,
            const int freq, const int seed, const ekat::Comm& comm)
{
  // Create grid
  auto gm = get_gm(comm);
  auto grid = gm->get_grid("point_grid");

  // Time advance parameters
  auto t0 = get_t0();
  const int dt = get_dt(freq_units);

  // Create some fields
  auto fm = get_fm(grid,t0,seed);
  std::vector<std::string> fnames;
  for (auto it : fm->get_repo()) {
    fnames.push_back(it.second->name());
  }

  // Create output params
  ekat::ParameterList om_pl;
  om_pl.set("filename_prefix",std::string("io_filled"));
  om_pl.set("field_names",fnames);
  om_pl.set("averaging_type", avg_type);
  om_pl.set<Real>("fill_threshold",fill_threshold);
  om_pl.set("track_avg_cnt",true);
  auto& ctrl_pl = om_pl.sublist("output_control");
  ctrl_pl.set("frequency_units",freq_units);
  ctrl_pl.set("frequency",freq);
  ctrl_pl.set("save_grid_data",false);

  // Create Output manager
  OutputManager om;
  om.initialize(comm,om_pl,t0,false);
  om.setup(fm,gm->get_grid_names());

  // Time loop: ensure we always hit 3 output steps
  const int nsteps = num_output_steps*freq;
  auto t = t0;
  for (int n=1; n<=nsteps; ++n) {
    om.init_timestep(t,dt);
    // Update time
    t += dt;

    // Set fields to n+1 or the fill_value, depending on step:
    //  - n if n is odd
    //  - fill_value if n is even
    bool valid_step = n % 2 == 1; 
    Real setval = valid_step ? Real(n) : fill_value;
    for (const auto& fn : fnames) {
      auto f = fm->get_field(fn);
      f.deep_copy(setval);
      auto& m = f.get_mask();
      m.deep_copy(valid_step ? 1 : 0);
    }

    // Run output manager
    om.run (t);
  }

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
  int first_snap = instant ? 0 : 1;
  int last_snap = num_output_steps;

  // Get gm
  auto gm = get_gm (comm);
  auto grid = gm->get_grid("point_grid");

  // Get initial fields. Use wrong seed for fm, so fields are not
  // inited with right data (avoid getting right answer without reading).
  auto fm0 = get_fm(grid,t0,seed);
  auto fm  = get_fm(grid,t0,-seed-1);
  std::vector<std::string> fnames;
  for (auto it : fm->get_repo()) {
    fnames.push_back(it.second->name());
  }

  // Create reader pl
  ekat::ParameterList reader_pl;
  std::string casename = "io_filled";
  auto filename = casename
    + "." + avg_type
    + "." + freq_units
    + "_x" + std::to_string(freq)
    + ".np" + std::to_string(comm.size())
    + "." + t0.to_string()
    + ".nc";
  reader_pl.set("filename",filename);
  reader_pl.set("field_names",fnames);
  AtmosphereInput reader(reader_pl,fm);

  // We set the value n (the step number) to each input field for each odd valued timestep
  // and fill_value for each even valued timestep
  // Hence, at output step N = snap*freq (snap=1,...,num_writes), we should get
  //  avg=INSTANT: output = N if (N%2=1), else Fillvalue (that's b/c we have t0 output, which is FV)
  //  avg=MAX:     output = N-1 if (N%2=0), else N
  //  avg=MIN:     output = N + 1, where n is the first timesnap of the Nth output step.
  //                        we add + 1 more in cases where (N%2=0) because that means the first snap was filled.
  //  avg=AVERAGE: output = sum([i if i%2==1]) divided by the list length if > 0.5*freq, else FillValue
  for (int n=first_snap; n<=last_snap; ++n) {
    reader.read_variables(n-first_snap);
    for (const auto& fn : fnames) {
      auto f0 = fm0->get_field(fn).clone();
      auto f  = fm->get_field(fn);
      if (avg_type=="MIN") {
        int first_avg_step = (n-1)*freq + 1;
        Real test_val = first_avg_step%2 == 0 ? first_avg_step+1 : first_avg_step;
        f0.deep_copy(test_val);
        REQUIRE (views_are_equal(f,f0));
      } else if (avg_type=="MAX") {
        int last_avg_step = n*freq;
        Real test_val = last_avg_step%2 == 0 ? last_avg_step-1 : last_avg_step;
        f0.deep_copy(test_val);
        REQUIRE (views_are_equal(f,f0));
      } else if (avg_type=="INSTANT") {
        Real test_val = (n*freq%2==0) ? fill_value : n*freq;
        f0.deep_copy(test_val);
        REQUIRE (views_are_equal(f,f0));
      } else { // Is avg_type = AVERAGE
        // Note, for AVERAGE type output with filling we need to check that the
        // number of contributing fill steps surpasses the fill_threshold, if not
        // then we know that the snap will reflect the fill value.
        int first_avg_step = (n-1)*freq + 1;
        int last_avg_step = n*freq;
        Real test_val = 0;
        int nvalid = 0;
        for (int i=first_avg_step; i<=last_avg_step; ++i) {
          if (i%2==1) {
            test_val += i;
            ++nvalid;
          }
        }
        test_val = nvalid > fill_threshold*freq ? test_val/nvalid : fill_value;

        f0.deep_copy(test_val);
        REQUIRE (views_are_equal(f,f0));
      }
    }
  }

  // Check that the fill value gets appropriately set for each variable
  for (const auto& fn: fnames) {
    // NOTE: use float, since default fp_precision for I/O is 'single'
    auto att_fill = scorpio::get_attribute<float>(filename,fn,"_FillValue");
    REQUIRE(att_fill==constants::fill_value<Real>);
  }
}

TEST_CASE ("io_filled") {
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
    print ("-> Output frequency: freq=" + std::to_string(freq) + ", units=" + units + "\n");
    for (const auto& avg : avg_type) {
      print("   -> Averaging type: " + avg + " ", 40);
      write(avg,units,freq,seed,comm);
      read(avg,units,freq,seed,comm);
      print(" PASS\n");
    }
  }
  scorpio::finalize_subsystem();
}

} // anonymous namespace
