#include <catch2/catch.hpp>

#include "share/io/scream_output_manager.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/io/scream_io_utils.hpp"

#include "share/grid/mesh_free_grids_manager.hpp"

#include "share/field/field_utils.hpp"
#include "share/field/field.hpp"
#include "share/field/field_manager.hpp"

#include "share/util/scream_universal_constants.hpp"
#include "share/util/scream_setup_random_test.hpp"
#include "share/util/scream_time_stamp.hpp"
#include "share/scream_types.hpp"

#include "ekat/util/ekat_units.hpp"
#include "ekat/ekat_parameter_list.hpp"
#include "ekat/ekat_assert.hpp"
#include "ekat/mpi/ekat_comm.hpp"
#include "ekat/util/ekat_test_utils.hpp"

#include <iomanip>
#include <memory>

namespace scream {

constexpr int num_output_steps = 5;
constexpr Real FillValue = constants::DefaultFillValue<float>().value;
constexpr Real fill_threshold = 0.5;

void set (const Field& f, const double v) {
  auto data = f.get_internal_view_data<Real,Host>();
  auto nscalars = f.get_header().get_alloc_properties().get_num_scalars();
  for (int i=0; i<nscalars; ++i) {
    data[i] = v;
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
    Field f(fid);
    f.allocate_view();
    f.deep_copy(0.0); // For the "filled" field we start with a filled value.
    f.get_header().get_tracking().update_time_stamp(t0);
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
  om_pl.set("filename_prefix",std::string("io_filled"));
  om_pl.set("Field Names",fnames);
  om_pl.set("Averaging Type", avg_type);
  om_pl.set<double>("fill_value",FillValue);
  om_pl.set<Real>("fill_threshold",fill_threshold);
  om_pl.set("track_avg_cnt",true);
  auto& ctrl_pl = om_pl.sublist("output_control");
  ctrl_pl.set("frequency_units",freq_units);
  ctrl_pl.set("Frequency",freq);
  ctrl_pl.set("save_grid_data",false);

  // Create Output manager
  OutputManager om;
  om.initialize(comm,om_pl,t0,false);
  om.setup(fm,gm);

  // Time loop: ensure we always hit 3 output steps
  const int nsteps = num_output_steps*freq;
  auto t = t0;
  for (int n=0; n<nsteps; ++n) {
    om.init_timestep(t,dt);
    // Update time
    t += dt;

    // Set fields to n or the FillValue, depending on timesnap
    Real setval = ((n+1) % 2 == 0) ? 1.0*(n+1) : FillValue;
    for (const auto& n : fnames) {
      auto f = fm->get_field(n);
      set(f,setval);
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
  std::string casename = "io_filled";
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

  // We set the value n to each input field for each odd valued timestep and FillValue for each even valued timestep
  // Hence, at output step N = snap*freq, we should get
  //  avg=INSTANT: output = N if (N%2=0), else Fillvalue
  //  avg=MAX:     output = N if (N%2=0), else N-1
  //  avg=MIN:     output = N + 1, where n is the first timesnap of the Nth output step.
  //                        we add + 1 more in cases where (N%2=0) because that means the first snap was filled.
  //  avg=AVERAGE: output = a + M+1 = a + M*(M+1)/M
  // The last one comes from
  //   a + 2*(1 + 2 +..+M)/M =
  //   a + 2*sum(i)/M = a + 2*(M(M+1)/2)/M,
  //         where M = freq/2 + ( N%2=0 ? 0 : 1 ),
  //               a = floor(N/freq)*freq + ( N%2=0 ? 0 : -1)
  for (int n=0; n<num_writes; ++n) {
    reader.read_variables(n);
    for (const auto& fn : fnames) {
      auto f0 = fm0->get_field(fn).clone();
      auto f  = fm->get_field(fn);
      if (avg_type=="MIN") {
        Real test_val = ((n+1)*freq%2==0) ? n*freq+1 : n*freq+2;
        set(f0,test_val);
        REQUIRE (views_are_equal(f,f0));
      } else if (avg_type=="MAX") {
        Real test_val = ((n+1)*freq%2==0) ? (n+1)*freq : (n+1)*freq-1;
        set(f0,test_val);
        REQUIRE (views_are_equal(f,f0));
      } else if (avg_type=="INSTANT") {
        Real test_val = (n*freq%2==0) ? n*freq : FillValue;
        set(f0,test_val);
        REQUIRE (views_are_equal(f,f0));
      } else { // Is avg_type = AVERAGE
        // Note, for AVERAGE type output with filling we need to check that the
        // number of contributing fill steps surpasses the fill_threshold, if not
        // then we know that the snap will reflect the fill value.
        Real test_val;
        Real M = freq/2 + (n%2==0 ? 0.0 :  1.0);
        Real a = n*freq + (n%2==0 ? 0.0 : -1.0);
        test_val = (M/freq > fill_threshold) ? a + (M+1.0) : FillValue;
        set(f0,test_val);
        REQUIRE (views_are_equal(f,f0));
      }
    }
  }

  // Check that the fill value gets appropriately set for each variable
  for (const auto& fn: fnames) {
    // NOTE: use float, since default fp_precision for I/O is 'single'
    auto att_fill = scorpio::get_attribute<float>(filename,fn,"_FillValue");
    REQUIRE(att_fill==constants::DefaultFillValue<float>().value);
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
    print ("-> Output frequency: " + units + "\n");
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
