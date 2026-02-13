#include <catch2/catch.hpp>

#include "share/io/eamxx_output_manager.hpp"
#include "share/io/scorpio_input.hpp"

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

constexpr int num_output_steps = 3;

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
  const int ngcols = std::max(comm.size(),1);
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

  const auto units = ekat::units::Units::nondimensional();
  int count=0;
  using stratts_t = std::map<std::string,std::string>;
  for (const auto& fl : layouts) {
    FID fid("f_"+std::to_string(count),fl,units,grid->name());
    Field f(fid);
    f.allocate_view();
    auto& str_atts = f.get_header().get_extra_data<stratts_t>("io: string attributes");
    str_atts["test"] = f.name();
    randomize_discrete (f,seed++,values);
    f.get_header().get_tracking().update_time_stamp(t0);
    fm->add_field(f);
    ++count;
  }

  return fm;
}

// Write output with or without transpose
void write (const std::string& avg_type, const std::string& freq_units,
            const int freq, const int seed, const ekat::Comm& comm,
            const bool transpose)
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
  std::string prefix = transpose ? "io_transpose_T" : "io_transpose_N";
  om_pl.set("filename_prefix", prefix);
  om_pl.set("field_names",fnames);
  om_pl.set("averaging_type", avg_type);
  if (transpose) {
    om_pl.set("transpose", true);  // Enable transposed output
  }
  auto& ctrl_pl = om_pl.sublist("output_control");
  ctrl_pl.set("frequency_units",freq_units);
  ctrl_pl.set("frequency",freq);
  ctrl_pl.set("save_grid_data",false);

  int max_snaps = num_output_steps;
  if (avg_type=="INSTANT") {
    ++max_snaps;
  }
  om_pl.set("max_snapshots_per_file", max_snaps);

  // Create Output manager
  OutputManager om;
  om_pl.set("floating_point_precision",std::string("single"));
  om.initialize(comm,om_pl,t0,false);
  om.setup(fm,gm->get_grid_names());

  // Time loop
  const int nsteps = num_output_steps*freq;
  auto t = t0;
  for (int n=0; n<nsteps; ++n) {
    om.init_timestep(t,dt);
    t += dt;

    // Add 1 to all fields entries
    for (const auto& name : fnames) {
      auto f = fm->get_field(name);
      add(f,1.0);
    }

    om.run (t);
  }

  // Check that the file was closed
  const auto& file_specs = om.output_file_specs();
  REQUIRE (not file_specs.is_open);

  om.finalize();
}

// Verify transposed data matches non-transposed data
void verify_transpose(const std::string& avg_type, const ekat::Comm& comm) {
  using namespace scorpio;
  
  const int freq = 3;
  std::string prefix_N = "io_transpose_N";
  std::string prefix_T = "io_transpose_T";
  std::string suffix = avg_type + ".nsteps_x3.np" + std::to_string(comm.size()) + ".2023-02-17-00000.nc";
  std::string file_N = prefix_N + "." + suffix;
  std::string file_T = prefix_T + "." + suffix;
  
  // Read dimensions from normal file
  register_file(file_N, Read);
  int ncol_N = get_dimlen(file_N, "ncol");
  int lev_N = get_dimlen(file_N, "lev");
  int time_N = get_dimlen(file_N, "time");
  
  // Read dimensions from transposed file
  register_file(file_T, Read);
  int ncol_T = get_dimlen(file_T, "ncol");
  int lev_T = get_dimlen(file_T, "lev");
  int time_T = get_dimlen(file_T, "time");
  
  // Verify dimensions are swapped correctly
  REQUIRE(ncol_N == ncol_T);
  REQUIRE(lev_N == lev_T);
  REQUIRE(time_N == time_T);
  
  // For field f_1 which has layout (ncol, lev), verify data is correctly transposed
  // Read a subset of data from both files and verify correspondence
  std::string varname = "f_1";
  
  // Get the variable dimensions
  auto dims_N = get_var_dims(file_N, varname);
  auto dims_T = get_var_dims(file_T, varname);
  
  // Normal file should be (time, ncol, lev)
  // Transposed file should be (time, lev, ncol)
  REQUIRE(dims_N.size() == 3);
  REQUIRE(dims_T.size() == 3);
  REQUIRE(dims_N[0] == "time");
  REQUIRE(dims_N[1] == "ncol");
  REQUIRE(dims_N[2] == "lev");
  REQUIRE(dims_T[0] == "time");
  REQUIRE(dims_T[1] == "lev");
  REQUIRE(dims_T[2] == "ncol");
  
  // Read data from first timestep to verify
  const int t = 0;
  std::vector<Real> data_N(ncol_N * lev_N);
  std::vector<Real> data_T(lev_T * ncol_T);
  
  // Set decomposition for normal file (ncol decomposed)
  set_dim_decomp(file_N, "ncol", ncol_N, 0, comm.rank() * ncol_N, comm.size() == 1 ? ncol_N : 1);
  
  // Set decomposition for transposed file (ncol decomposed, but it's last dimension)
  set_dim_decomp(file_T, "ncol", ncol_T, 0, comm.rank() * ncol_T, comm.size() == 1 ? ncol_T : 1);
  
  // Read the data
  read_var(file_N, varname, data_N.data());
  read_var(file_T, varname, data_T.data());
  
  // Verify that data_N[col, lev] == data_T[lev, col]
  // For C-order storage: data_N[col * lev_N + lev] == data_T[lev * ncol_T + col]
  for (int col = 0; col < (comm.size() == 1 ? ncol_N : 1); ++col) {
    for (int lev = 0; lev < lev_N; ++lev) {
      Real val_N = data_N[col * lev_N + lev];
      Real val_T = data_T[lev * ncol_T + col];
      REQUIRE(val_N == val_T);
    }
  }
  
  release_file(file_N);
  release_file(file_T);
}

TEST_CASE ("io_transpose") {
  std::vector<std::string> avg_type = {
    "INSTANT",
    "MAX",
    "MIN",
    "AVERAGE"
  };

  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);

  auto seed = get_random_test_seed(&comm);

  const int freq = 3;
  const std::string freq_units = "nsteps";
  
  auto print = [&] (const std::string& s, int line_len = -1) {
    if (comm.am_i_root()) {
      if (line_len<0) {
        std::cout << s;
      } else {
        std::cout << std::left << std::setw(line_len) << std::setfill('.') << s;
      }
    }
  };

  print ("Testing transposed output\n");
  for (const auto& avg : avg_type) {
    print("   -> Averaging type: " + avg + " ", 40);
    
    // Write both transposed and non-transposed versions
    write(avg,freq_units,freq,seed,comm,false);  // Non-transposed
    write(avg,freq_units,freq,seed,comm,true);   // Transposed
    
    // Verify the transposed data is correct
    verify_transpose(avg, comm);
    
    print(" PASS\n");
  }
  
  scorpio::finalize_subsystem();
}

} // namespace scream
