#include <catch2/catch.hpp>

#include "share/atm_process/atmosphere_diagnostic.hpp"

#include "share/io/eamxx_output_manager.hpp"
#include "share/io/scorpio_input.hpp"

#include "share/grid/mesh_free_grids_manager.hpp"

#include "share/field/field_utils.hpp"
#include "share/field/field.hpp"
#include "share/field/field_manager.hpp"

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

class MyDiag : public AtmosphereDiagnostic
{
public:
  MyDiag (const ekat::Comm& comm, const ekat::ParameterList&params)
    : AtmosphereDiagnostic(comm,params)
  {
    //Do nothing
  }

  std::string name() const override { return "MyDiag"; }

  void set_grids (const std::shared_ptr<const GridsManager> gm) override {
    using namespace ekat::units;
    using namespace ShortFieldTagsNames;
    using FL = FieldLayout;

    const auto grid = gm->get_grid("Point Grid");
    const auto& grid_name = grid->name();
    m_num_cols  = grid->get_num_local_dofs(); // Number of columns on this rank
    m_num_levs  = grid->get_num_vertical_levels();  // Number of levels per column

    std::vector<FieldTag> tag_2d = {COL,LEV};
    std::vector<Int>     dims_2d = {m_num_cols,m_num_levs};
    FL lt( tag_2d, dims_2d );

    m_f_in = "f_"+std::to_string(lt.size());

    auto units = Units::nondimensional();
    add_field<Required>(m_f_in,lt,units,grid_name);

    // We have to initialize the m_diagnostic_output:
    FieldIdentifier fid (name(), lt, units, grid_name);
    m_diagnostic_output = Field(fid);
    m_diagnostic_output.allocate_view();
    m_one = m_diagnostic_output.clone("one");
    m_one.deep_copy(1.0);
  }

  void init_timestep (const util::TimeStamp& start_of_step) override {
    m_t_beg = start_of_step;
  }

protected:

  void compute_diagnostic_impl () override {
    const auto& f_in  = get_field_in(m_f_in);

    const auto& t = f_in.get_header().get_tracking().get_time_stamp();
    const double dt = t - m_t_beg;

    m_diagnostic_output.deep_copy(f_in);
    m_diagnostic_output.update(m_one,dt,2.0);
  }

  void initialize_impl (const RunType /* run_type */ ) override {
    m_diagnostic_output.get_header().get_tracking().update_time_stamp(start_of_step_ts());
  }

  // Clean up
  void finalize_impl ( /* inputs */ ) override {}

  // Internal variables
  int m_num_cols, m_num_levs;

  std::string m_f_in;

  util::TimeStamp m_t_beg;
  Field m_one;
};

util::TimeStamp get_t0 () {
  return util::TimeStamp({2023,2,17},{0,0,0});
}

constexpr double get_dt () {
  return 10;
};

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
        const util::TimeStamp& t0, const int seed,
        const bool add_diag_field = false)
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

  auto fm = std::make_shared<FieldManager>(grid);

  const auto units = ekat::units::Units::nondimensional();
  FL fl ({COL,LEV}, {nlcols,nlevs});

  FID fid("f_"+std::to_string(fl.size()),fl,units,grid->name());
  Field f(fid);
  f.allocate_view();
  randomize (f,engine,my_pdf);
  f.get_header().get_tracking().update_time_stamp(t0);
  fm->add_field(f);

  if (add_diag_field) {
    FID diag_id("MyDiag",fl,units,grid->name());
    Field diag(diag_id);
    diag.allocate_view();
    randomize (diag,engine,my_pdf);
    fm->add_field(diag);
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
  auto dt = get_dt();

  // Create some fields
  auto fm = get_fm(grid,t0,seed);
  std::vector<std::string> fnames;
  for (auto it : fm->get_repo()) {
    const auto& fn = it.second->name();
    fnames.push_back(fn);
  }
  fnames.push_back("MyDiag");

  // Create output params
  ekat::ParameterList om_pl;
  om_pl.set("filename_prefix",std::string("io_diags"));
  om_pl.set("Field Names",fnames);
  om_pl.set("Averaging Type", std::string("INSTANT"));
  auto& ctrl_pl = om_pl.sublist("output_control");
  ctrl_pl.set("frequency_units",std::string("nsteps"));
  ctrl_pl.set("Frequency",1);
  ctrl_pl.set("save_grid_data",false);

  // Create Output manager
  OutputManager om;
  om.initialize(comm, om_pl, t0, false);
  om.setup(fm,gm->get_grid_names());

  // Run output manager
  for (auto it : fm->get_repo()) {
    auto& f = *it.second;
    Field one = f.clone("one");
    one.deep_copy(1.0);
    f.get_header().get_tracking().update_time_stamp(t0+dt);
    f.update(one,1.0,1.0);
  }
  om.init_timestep(t0,dt);
  om.run (t0+dt);

  // Close file and cleanup
  om.finalize();
}

void read (const int seed, const ekat::Comm& comm)
{
  // Time quantities
  auto t0 = get_t0();

  // Get gm
  auto gm = get_gm (comm);
  auto grid = gm->get_grid("Point Grid");

  // Get initial fields
  auto fm0 = get_fm(grid,t0,seed);
  auto fm  = get_fm(grid,t0,seed,true);

  std::vector<std::string> fnames;
  std::string f_name;
  for (auto it : fm->get_repo()) {
    const auto& fn = it.second->name();
    fnames.push_back(fn);
    if (fn!="MyDiag") {
      f_name = fn;
    }
  }

  auto f0 = fm0->get_field(f_name).clone();
  auto f  = fm->get_field(f_name);

  // Sanity check
  REQUIRE (f_name!="");

  // Create reader pl
  ekat::ParameterList reader_pl;
  std::string casename = "io_diags";
  auto filename = casename
    + ".INSTANT.nsteps_x1"
    + ".np" + std::to_string(comm.size())
    + "." + t0.to_string()
    + ".nc";
  reader_pl.set("Filename",filename);
  reader_pl.set("Field Names",fnames);
  AtmosphereInput reader(reader_pl,fm);

  Field one = f0.clone("one");
  one.deep_copy(1.0);
  for (int i=0; i<2; ++i) {
    reader.read_variables(i);

    // Check regular field is correct
    REQUIRE (views_are_equal(f,f0));

    // Check diag field is correct
    const auto t = t0+i*get_dt();
    const double dt = t-t0;
    auto d = fm->get_field("MyDiag");
    auto d0 = f0.clone();
    d0.update(one,dt,2.0);
    REQUIRE (views_are_equal(d,d0));

    // Update f0
    f0.update(one,1.0,1.0);
  }
}

TEST_CASE ("io_diags") {
  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);

  // Make MyDiag available via diag factory
  auto& diag_factory = AtmosphereDiagnosticFactory::instance();
  diag_factory.register_product("MyDiag",&create_atmosphere_diagnostic<MyDiag>);

  auto seed = get_random_test_seed(&comm);

  auto print = [&] (const std::string& s, int line_len = -1) {
    if (comm.am_i_root()) {
      if (line_len<0) {
        std::cout << s;
      } else {
        std::cout << std::left << std::setw(line_len) << std::setfill('.') << s;
      }
    }
  };

  print ("-> Write diagnostic output ", 40);
  write(seed,comm);
  read(seed,comm);
  print(" PASS\n");
  scorpio::finalize_subsystem();
}

} // anonymous namespace
