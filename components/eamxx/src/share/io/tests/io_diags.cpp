#include <catch2/catch.hpp>

#include "share/atm_process/atmosphere_diagnostic.hpp"

#include "share/io/eamxx_output_manager.hpp"
#include "share/io/scorpio_input.hpp"

#include "share/grid/mesh_free_grids_manager.hpp"

#include "share/field/field_utils.hpp"
#include "share/field/field.hpp"
#include "share/manager/field_manager.hpp"

#include "eamxx_setup_random_test.hpp"
#include "share/util/eamxx_time_stamp.hpp"
#include "share/core/eamxx_types.hpp"

#include <ekat_units.hpp>
#include <ekat_parameter_list.hpp>
#include <ekat_assert.hpp>
#include <ekat_comm.hpp>

#include <iomanip>
#include <memory>

namespace scream {

// Small class to expose protected members of Atm output
class OutputTester : public AtmosphereOutput
{
public:
  OutputTester(const ekat::Comm& comm, const ekat::ParameterList& params,
               const std::shared_ptr<const fm_type>& field_mgr, const std::string& grid_name)
    : AtmosphereOutput (comm,params,field_mgr,grid_name)
  { /* Nothing to do here */ }

  std::list<diag_ptr_type> get_diags () const { return m_diagnostics; }
};

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
    const auto grid = gm->get_grid("point_grid");
    const auto& grid_name = grid->name();
    auto units = ekat::units::Units::nondimensional();
    auto layout = grid->get_3d_scalar_layout(true);
    add_field<Required>("my_f",layout,units,grid_name);

    // We have to initialize the m_diagnostic_output:
    FieldIdentifier fid ("MyDiag", layout, units, grid_name);
    m_diagnostic_output = Field(fid);
    m_diagnostic_output.allocate_view();
    m_one = m_diagnostic_output.clone("one");
    m_one.deep_copy(1.0);
  }

  void init_timestep (const util::TimeStamp& start_of_step) override {
    m_t_beg = start_of_step;
  }

  int get_num_evaluations () const { return m_num_evaluations; }
protected:

  void compute_diagnostic_impl () override {
    const auto& f_in  = get_field_in("my_f");

    const auto& t = f_in.get_header().get_tracking().get_time_stamp();
    const double dt = t - m_t_beg;

    m_diagnostic_output.deep_copy(f_in);
    m_diagnostic_output.update(m_one,dt,2.0);

    ++m_num_evaluations;
  }

  void initialize_impl (const RunType /* run_type */ ) override {
    m_diagnostic_output.get_header().get_tracking().update_time_stamp(start_of_step_ts());
  }

  // Clean up
  void finalize_impl ( /* inputs */ ) override {}

  util::TimeStamp m_t_beg;
  Field m_one;

  int m_num_evaluations = 0;
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

  FID fid("my_f",fl,units,grid->name());
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
  auto grid = gm->get_grid("point_grid");

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
  om_pl.set("field_names",fnames);
  om_pl.set("averaging_type", std::string("instant"));
  auto& ctrl_pl = om_pl.sublist("output_control");
  ctrl_pl.set("frequency_units",std::string("nsteps"));
  ctrl_pl.set("frequency",1);
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
  auto grid = gm->get_grid("point_grid");

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
  reader_pl.set("filename",filename);
  reader_pl.set("field_names",fnames);
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

TEST_CASE ("diags_sharing_check") {
  // NOTE: this test does NOT do any scorpio call. We create AtmosphereOutput objects
  // rather than OutputManager objects (which WOULD do scorpio calls), and set up the
  // output specs so that when we call the run method ONCE, we don't trigger a write.
  // All we want to do is verify that the two writers share the same diag object,
  // and that the diag object only runs once, skipping the compute call for the second stream
  ekat::Comm comm(MPI_COMM_WORLD);

  // Make MyDiag available via diag factory
  auto& diag_factory = AtmosphereDiagnosticFactory::instance();
  diag_factory.register_product("MyDiag",&create_atmosphere_diagnostic<MyDiag>);

  // Time stamp
  auto time = get_t0();

  // Create grids and field managers
  auto gm = get_gm (comm);
  auto grid = gm->get_grid("point_grid");
  auto fm = get_fm(grid,time,123);

  // Create output params
  ekat::ParameterList params;
  params.set("averaging_type", std::string("average"));
  params.sublist("fields").sublist(grid->name()).set<std::vector<std::string>>("field_names",{"MyDiag"});

  // Create Output testers
  OutputTester out1(comm,params,fm,grid->name());
  OutputTester out2(comm,params,fm,grid->name());

  // Check the diags in the two streams are the SAME instance
  auto diags1 = out1.get_diags();
  auto diags2 = out2.get_diags();

  REQUIRE (diags1.size()==1);
  REQUIRE (diags2.size()==1);
  REQUIRE (diags1.front()==diags2.front());
  auto d = std::dynamic_pointer_cast<MyDiag>(diags1.front());

  // Run outputs and verify the diag was evaluated only once
  out1.init_timestep(time);
  out2.init_timestep(time);

  out1.run("UNUSED", false, false, 0, false);
  out2.run("UNUSED", false, false, 0, false);

  REQUIRE (d->get_num_evaluations()==1);

  // Update diag input, then run again, and verify the diag was evaluated one more time
  time += 1.0;
  fm->get_field("my_f").get_header().get_tracking().update_time_stamp(time);
  out1.run("UNUSED", false, false, 0, false);
  out2.run("UNUSED", false, false, 0, false);

  REQUIRE (d->get_num_evaluations()==2);
}

} // anonymous namespace
