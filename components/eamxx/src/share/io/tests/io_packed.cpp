#include <catch2/catch.hpp>

#include "share/io/scream_output_manager.hpp"
#include "share/io/scorpio_input.hpp"

#include "share/grid/mesh_free_grids_manager.hpp"

#include "share/field/field_utils.hpp"
#include "share/field/field.hpp"
#include "share/field/field_manager.hpp"

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

void add (const Field& f, const double v) {
  auto data = f.get_internal_view_data<Real,Host>();
  auto nscalars = f.get_header().get_alloc_properties().get_num_scalars();
  for (int i=0; i<nscalars; ++i) {
    data[i] += v;
  }
}

util::TimeStamp get_t0 () {
  return util::TimeStamp({2023,2,17},{0,0,0});
}

std::shared_ptr<const GridsManager>
get_gm (const ekat::Comm& comm)
{
  const int nlcols = 3;
  const int nlevs = 16;
  const int ngcols = nlcols*comm.size();
  ekat::ParameterList gm_params;
  gm_params.set("number_of_global_columns",ngcols);
  gm_params.set("number_of_vertical_levels",nlevs);
  auto gm = create_mesh_free_grids_manager(comm,gm_params);
  gm->build_grids();
  return gm;
}

std::shared_ptr<FieldManager>
get_fm (const std::shared_ptr<const AbstractGrid>& grid,
        const util::TimeStamp& t0, const int seed, const int ps)
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
    FL({COL,     LEV}, {nlcols,  nlevs}),
    FL({COL,CMP,ILEV}, {nlcols,2,nlevs+1})
  };

  auto fm = std::make_shared<FieldManager>(grid);
  fm->registration_begins();
  fm->registration_ends();
  
  const auto units = ekat::units::Units::nondimensional();
  for (const auto& fl : layouts) {
    FID fid("f_"+std::to_string(fl.size()),fl,units,grid->name());
    Field f(fid);
    f.get_header().get_alloc_properties().request_allocation(ps);
    f.allocate_view();
    randomize (f,engine,my_pdf);
    f.get_header().get_tracking().update_time_stamp(t0);
    fm->add_field(f);
  }

  return fm;
}

// Returns fields after initialization
void write (const int freq, const int seed, const int ps, const ekat::Comm& comm)
{
  // Create grid
  auto gm = get_gm(comm);
  auto grid = gm->get_grid("Point Grid");

  // Time advance parameters
  auto t0 = get_t0();

  // Create some fields
  auto fm = get_fm(grid,t0,seed,ps);
  std::vector<std::string> fnames;
  for (auto it : *fm) {
    fnames.push_back(it.second->name());
  }

  // Create output params
  ekat::ParameterList om_pl;
  om_pl.set("MPI Ranks in Filename",true);
  om_pl.set("filename_prefix","io_packed_ps"+std::to_string(ps));
  om_pl.set("Field Names",fnames);
  om_pl.set("Averaging Type", std::string("INSTANT"));
  auto& ctrl_pl = om_pl.sublist("output_control");
  ctrl_pl.set("frequency_units",std::string("nsteps"));
  ctrl_pl.set("Frequency",freq);
  ctrl_pl.set("MPI Ranks in Filename",true);
  ctrl_pl.set("save_grid_data",false);

  // Create Output manager
  OutputManager om;
  om.setup(comm,om_pl,fm,gm,t0,t0,false);

  // Run output manager
  om.run (t0);

  // Close file and cleanup
  om.finalize();
}

void read (const int freq, const int seed, const int ps_write, const int ps_read, const ekat::Comm& comm)
{
  // Time quantities
  auto t0 = get_t0();

  // Get gm
  auto gm = get_gm (comm);
  auto grid = gm->get_grid("Point Grid");

  // Get initial fields. Use wrong seed for fm, so fields are not
  // inited with right data (avoid getting right answer without reading).
  auto fm0 = get_fm(grid,t0,seed,ps_read);
  auto fm  = get_fm(grid,t0,-seed-1,ps_read);
  std::vector<std::string> fnames;
  for (auto it : *fm) {
    fnames.push_back(it.second->name());
  }

  // Create reader pl
  ekat::ParameterList reader_pl;
  std::string casename = "io_packed_ps"+std::to_string(ps_write);
  auto filename = casename
    + ".INSTANT.nsteps"
    + "_x" + std::to_string(freq)
    + ".np" + std::to_string(comm.size())
    + "." + t0.to_string()
    + ".nc";
  reader_pl.set("Filename",filename);
  reader_pl.set("Field Names",fnames);
  AtmosphereInput reader(reader_pl,fm);

  reader.read_variables();
  for (const auto& fn : fnames) {
    auto f0 = fm0->get_field(fn);
    auto f  = fm->get_field(fn);
    if (not views_are_equal(f,f0)) {
      print_field_hyperslab(f,{},{},std::cout);
      print_field_hyperslab(f0,{},{},std::cout);
    }
    REQUIRE (views_are_equal(f,f0));
  }
}

TEST_CASE ("io_packs") {
  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::eam_init_pio_subsystem(comm);

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

  for (const auto ps_write : {1,2,4,8}) {
    print ("-> Pack size write: " + std::to_string(ps_write) + "\n");
    write(freq,seed,ps_write,comm);
    for (const auto ps_read : {1,2,4,8,16}) {
      print ("  -> Pack size read: " + std::to_string(ps_read) + " ",40);
      read(freq,seed,ps_write,ps_read,comm);
      print(" PASS\n");
    }
  }
  scorpio::eam_pio_finalize();
}

} // anonymous namespace
