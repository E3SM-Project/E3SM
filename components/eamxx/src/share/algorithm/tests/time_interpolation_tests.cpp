#include <catch2/catch.hpp>

#include "share/data_managers/mesh_free_grids_manager.hpp"

#include "share/field/field_utils.hpp"
#include "share/field/field.hpp"
#include "share/data_managers/field_manager.hpp"

#include "share/algorithm/eamxx_time_interpolation.hpp"
#include "share/core/eamxx_setup_random_test.hpp"
#include "share/util/eamxx_time_stamp.hpp"

#include "share/io/eamxx_output_manager.hpp"
#include "share/core/eamxx_config.hpp"
#include "share/io/eamxx_io_utils.hpp"
#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>

#include <ekat_parameter_list.hpp>
/*-----------------------------------------------------------------------------------------------
 * Test TimeInterpolation class
 *-----------------------------------------------------------------------------------------------*/

namespace scream {

// Test Constants
constexpr Real tol            = std::numeric_limits<Real>::epsilon()*1e5;
constexpr int  dt             = 100;
constexpr int  total_snaps    = 10;
constexpr int  snap_freq      = 4;
constexpr int  slope_freq     = snap_freq; // We will change the slope every slope_freq steps to ensure that the data is not all on the same line.
constexpr int  snaps_per_file = 3;

// Functions needed to set the test up
std::shared_ptr<const GridsManager> get_gm (const ekat::Comm& comm, const int ncols, const int nlevs);
std::shared_ptr<FieldManager> get_fm (const std::shared_ptr<const AbstractGrid>& grid, const util::TimeStamp& t0, const int seed);
util::TimeStamp init_timestamp();
std::vector<std::string> create_test_data_files(const ekat::Comm& comm, const std::shared_ptr<const GridsManager>& gm,const util::TimeStamp& t0, const int seed);

// Functions needed to run the test
void update_field_data(const Real slope, const Real dt, Field& field);
bool views_are_approx_equal(const Field& f0, const Field& f1, const Real tol);

// Helper randomizer
Real my_pdf(std::mt19937_64& engine) {
  std::uniform_int_distribution<int> pdf (1,100);
  Real v = pdf(engine);
  return v;
};

// Wrapper for output manager that also extracts the list of files
class OutputManager4Test : public scream::OutputManager
{
public:
  OutputManager4Test()
    : OutputManager()
  {
    // Do Nothing
  }

  void runme(const util::TimeStamp& ts) {
    run(ts);
    update_file_list();
  }

  std::vector<std::string> get_list_of_files() { return m_list_of_files; }
private:
  void update_file_list() {
    if (m_output_file_specs.is_open and std::find(m_list_of_files.begin(),m_list_of_files.end(), m_output_file_specs.filename) == m_list_of_files.end()) {
      m_list_of_files.push_back(m_output_file_specs.filename);
    }
  }
  std::vector<std::string> m_list_of_files;
};
/*-----------------------------------------------------------------------------------------------*/
TEST_CASE ("eamxx_time_interpolation_simple") {
  printf("TimeInterpolation - Simple Case...\n\n\n");
  // Setup basic test params
  ekat::Comm comm(MPI_COMM_WORLD);
  auto seed = get_random_test_seed(&comm);
  std::mt19937_64 engine(seed);
  const auto t0 = init_timestamp();

  const int nlevs  = SCREAM_PACK_SIZE*2+1;
  const int ncols  = comm.size()*2 + 1;

  // Get a grids manager for the test
  auto grids_man = get_gm(comm, ncols, nlevs);
  const auto& grid = grids_man->get_grid("point_grid");
  // Now create a fields manager to store initial data for testing.
  auto fields_man_t0 = get_fm(grid, t0, seed);
  // Construct a time interpolation object and add all of the fields to it.
  printf(  "Constructing a time interpolation object ...\n");
  util::TimeInterpolation time_interpolator(grid);
  for (auto ff_pair : fields_man_t0->get_repo()) {
    const auto ff   = ff_pair.second;
    time_interpolator.add_field(*ff);
    time_interpolator.initialize_data_from_field(*ff);
  }
  time_interpolator.initialize_timestamps(t0);
  printf(  "Constructing a time interpolation object ... DONE\n");

  // Set the field data for a step 10 dt in the future to be used for interpolation.  We
  // create a new field manager so we can continue to have the intial condition for reference.
  printf(  "Setting data for other end of time interpolation, at t = 10 dt ...\n");
  auto slope = my_pdf(engine);
  auto fields_man_tf = get_fm(grid, t0, seed);
  auto t1 = t0 + 10*dt;
  for (auto ff_pair : fields_man_tf->get_repo()) {
    auto  ff = ff_pair.second;
    update_field_data(slope, t1.seconds_from(t0), *ff);
    time_interpolator.update_data_from_field(*ff);
  }
  time_interpolator.update_timestamp(t1);
  printf(  "Setting data for other end of time interpolation, at t = 10 dt ...DONE\n");

  // Now check that the interpolator is working as expected.  Should be able to
  // match the interpolated fields against the results of update_field_data at any
  // time between 0 and 10 dt.
  printf(  "Testing all timesteps ... slope = %f\n",slope);
  auto fields_man_test = get_fm(grid, t0, seed);
  t1 = t0;
  for (int tt = 1; tt<=10; tt++) {
    t1 += dt;
    printf("                        ... t = %s\n",t1.to_string().c_str());
    time_interpolator.perform_time_interpolation(t1);
    for (const auto& [name,ff] : fields_man_test->get_repo()) {
      update_field_data(slope, dt, *ff);
      REQUIRE(views_are_approx_equal(*ff,time_interpolator.get_field(name),tol));
    }
  }
  printf("                        ... DONE\n");
  printf("TimeInterpolation - Simple Case...DONE\n\n\n");
} // TEST_CASE eamxx_time_interpolation_simple
/*-----------------------------------------------------------------------------------------------*/
TEST_CASE ("eamxx_time_interpolation_data_from_file") {
  printf("TimeInterpolation - From File Case...\n\n\n");
  // Setup basic test
  printf("   - Test Basics...\n");
  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);
  auto seed = get_random_test_seed(&comm);
  std::mt19937_64 engine(seed);
  const auto t0 = init_timestamp();

  const int nlevs  = SCREAM_PACK_SIZE*2+1;
  const int ncols  = comm.size()*2 + 1;
  printf("   - Test Basics...DONE\n");

  // Get a grids manager for the test
  printf("   - Grids Manager...\n");
  auto grids_man = get_gm(comm, ncols, nlevs);
  const auto& grid = grids_man->get_grid("point_grid");
  printf("   - Grids Manager...DONE\n");
  // Now create a fields manager to store initial data for testing.
  printf("   - Fields Manager...\n");
  auto fields_man_t0 = get_fm(grid, t0, seed);
  auto fields_man_deep = get_fm(grid, t0, seed);  // A field manager for checking deep copies.
  std::vector<std::string> fnames;
  for (auto it : fields_man_t0->get_repo()) {
    fnames.push_back(it.second->name());
  }
  printf("   - Fields Manager...DONE\n");
  // Construct the files of interpolation data
  printf("   - create test data files...\n");
  auto list_of_files = create_test_data_files(comm, grids_man, t0, seed);
  printf("   - create test data files...DONE\n");

  // Construct a time interpolation object using the list of files with the data
  printf(  "Constructing a time interpolation object ...\n");
  util::TimeInterpolation time_interpolator(grid,list_of_files);
  util::TimeInterpolation time_interpolator_deep(grid,list_of_files);
  for (auto name : fnames) {
    auto ff      = fields_man_t0->get_field(name);
    auto ff_deep = fields_man_deep->get_field(name);
    time_interpolator.add_field(ff);
    time_interpolator_deep.add_field(ff_deep,true);
  }
  time_interpolator.initialize_data_from_files();
  time_interpolator_deep.initialize_data_from_files();
  printf(  "Constructing a time interpolation object ... DONE\n");

  // Now check that the interpolator is working as expected.  Should be able to
  // match the interpolated fields against the results of update_field_data at any
  // time between 0 and 10 dt.
  printf(  "Testing all timesteps ...\n");
  const int max_steps      = snap_freq*total_snaps;
  // The strategy is to update the model state following the same linear updates
  // that we used to generate the time interpolation data.  The fields produced
  // by the time interpolator should match.
  auto ts = t0;
  for (int nn=0; nn<=max_steps; nn++) {
    // Update Time and slope
    if (nn > 0) {
      ts += dt;
    }
    printf("                        ... t = %s\n",ts.to_string().c_str());
    // Perform time interpolation
    if (nn > 0) {
      const Real slope = ((nn-1) / slope_freq) + 1;
      for (auto name : fnames) {
        auto field = fields_man_t0->get_field(name);
        update_field_data(slope,dt,field);
	// We set the deep copy fields to wrong values to stress test that everything still works.
	auto field_deep = fields_man_deep->get_field(name);
	field_deep.deep_copy(-9999.0);
      }
    }
    time_interpolator.perform_time_interpolation(ts);
    time_interpolator_deep.perform_time_interpolation(ts);
    // Now compare the interp_fields to the fields in the field manager which should be updated.
    for (auto name : fnames) {
      auto field      = fields_man_t0->get_field(name);
      auto field_deep = fields_man_deep->get_field(name);
      // Check that the shallow copies match the expected values
      REQUIRE(views_are_approx_equal(field,time_interpolator.get_field(name),tol));
      // Check that the deep fields which were not updated directly are the same as the ones stored in the time interpolator.
      REQUIRE(views_are_equal(field_deep,time_interpolator_deep.get_field(name)));
      // Check that the deep and shallow fields match showing that both approaches got the correct answer.
      REQUIRE(views_are_equal(field,field_deep));
    }

  }


  time_interpolator.finalize();
  time_interpolator_deep.finalize();
  printf("                        ... DONE\n");

  // All done with IO
  scorpio::finalize_subsystem();

  printf("TimeInterpolation - From File Case...DONE\n\n\n");
} // TEST_CASE eamxx_time_interpolation_data_from_file
/*-----------------------------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------------------------*/
bool views_are_approx_equal(const Field& f0, const Field& f1, const Real tol)
{
  const auto& l0    = f0.get_header().get_identifier().get_layout();
  const auto& l1    = f1.get_header().get_identifier().get_layout();
  EKAT_REQUIRE_MSG(l0==l1,"Error! views_are_approx_equal - the two fields don't have matching layouts.");
  // Take advantage of field utils update, min and max to assess the max difference between the two fields
  // simply.
  auto ft = f0.clone();
  ft.update(f1,1.0,-1.0);
  auto d_min = field_min(ft).as<Real>();
  auto d_max = field_max(ft).as<Real>();
  if (std::abs(d_min) > tol or std::abs(d_max) > tol) {
    printf("The two copies of (%16s) are NOT approx equal within a tolerance of %e.\n     The min and max errors are %e and %e respectively.\n",f0.name().c_str(),tol,d_min,d_max);
    return false;
  } else {
    return true;
  }

}
/*-----------------------------------------------------------------------------------------------*/
void update_field_data(const Real slope, const Real dt, Field& field)
{
  const auto& l1    = field.get_header().get_identifier().get_layout();
  const auto& dims  = l1.dims();
  switch (l1.rank()) {
    case 1:
      {
        auto view = field.template get_view<Real*,Host>();
        for (int ii=0; ii<dims[0]; ++ii) {
          view(ii) += slope * dt;
        }
      }
      break;
    case 2:
      {
        auto view = field.template get_view<Real**,Host>();
        for (int ii=0; ii<dims[0]; ++ii) {
          for (int jj=0; jj<dims[1]; ++jj) {
            view(ii,jj) += slope * dt;
          }
        }
      }
      break;
    case 3:
      {
        auto view = field.template get_view<Real***,Host>();
        for (int ii=0; ii<dims[0]; ++ii) {
          for (int jj=0; jj<dims[1]; ++jj) {
            for (int kk=0; kk<dims[2]; ++kk) {
              view(ii,jj,kk) += slope * dt;
            }
          }
        }
      }
      break;
    default:
      EKAT_ERROR_MSG("Error! Unsupported field rank.\n");
  }
}
/*-----------------------------------------------------------------------------------------------*/
util::TimeStamp init_timestamp()
{
  return util::TimeStamp({2023,5,16},{0,0,0});
}
/*-----------------------------------------------------------------------------------------------*/
/* Create a grids manager for the test */
std::shared_ptr<const GridsManager> get_gm (const ekat::Comm& comm, const int ncols, const int nlevs)
{
  using vos_t = std::vector<std::string>;
  ekat::ParameterList gm_params;
  gm_params.set("grids_names",vos_t{"point_grid"});
  auto& pl = gm_params.sublist("point_grid");
  pl.set<std::string>("type","point_grid");
  pl.set("aliases",vos_t{"physics"});
  pl.set<int>("number_of_global_columns", ncols);
  pl.set<int>("number_of_vertical_levels", nlevs);
  auto gm = create_mesh_free_grids_manager(comm,gm_params);
  gm->build_grids();
  return gm;
}
/*-----------------------------------------------------------------------------------------------*/
/* Create a fields manager for the test */
std::shared_ptr<FieldManager>
get_fm (const std::shared_ptr<const AbstractGrid>& grid, const util::TimeStamp& t0, int seed)
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

  auto fm = std::make_shared<FieldManager>(grid,RepoState::Closed);

  const auto units = ekat::units::Units::nondimensional();
  std::vector<Real> values;
  for (int i=1; i<=100; ++i)
    values.push_back(static_cast<Real>(i));
  for (const auto& fl : layouts) {
    int gl_size = fl.size();
    grid->get_comm().all_reduce(&gl_size,1,MPI_SUM);
    FID fid("f_"+std::to_string(gl_size),fl,units,grid->name());
    Field f(fid);
    f.allocate_view();
    randomize_discrete (f,seed++,values);
    f.get_header().get_tracking().update_time_stamp(t0);
    fm->add_field(f);
  }

  return fm;
}
/*-----------------------------------------------------------------------------------------------*/
/* Construct data for multiple time snaps and write the data to file, to be used for testing
 * the capability of TimeInterpolation to handle data read from multiple files.
 */
std::vector<std::string> create_test_data_files(
		const ekat::Comm& comm,
		const std::shared_ptr<const GridsManager>& gm,
		const util::TimeStamp& t0,
	  const int seed)
{
  // We initialize a local field manager to use for output
  auto fm = get_fm(gm->get_grid("point_grid"), t0, seed);
  // We will write data for 10 snaps for this test.  We set the max snaps per file to 3 to
  // ensure that a) there is more than 1 file and b) at least one file has fewer snap then
  // the others.
  const int max_steps      = snap_freq*total_snaps;
  // Gather the set of fields from the field manager
  std::vector<std::string> fnames;
  for (auto it : fm->get_repo()) {
    fnames.push_back(it.second->name());
  }
  // Create the output parameters
  ekat::ParameterList om_pl;
  om_pl.set("filename_prefix",std::string("source_data_for_time_interpolation"));
  om_pl.set("field_names",fnames);
  om_pl.set("averaging_type", std::string("instant"));
  om_pl.set("max_snapshots_per_file",snaps_per_file);
  auto& ctrl_pl = om_pl.sublist("output_control");
  ctrl_pl.set("frequency_units",std::string("nsteps"));
  ctrl_pl.set("frequency",snap_freq);
  ctrl_pl.set("save_grid_data",false);
  // Create an output manager, note we use a subclass defined in this test so we can extract
  // the list of files created by the output manager.
  OutputManager4Test om;
  om.initialize(comm,om_pl,t0,false);
  om.setup(fm,gm->get_grid_names());

  // Time loop to create and write data
  auto tw = t0;
  for (int nn=0; nn<=max_steps; nn++) {
    // Update Time and slope
    tw += dt;
    const Real slope = (nn / slope_freq) + 1;
    // Update Fields
    for (auto ff : fnames) {
      auto field = fm->get_field(ff);
      update_field_data(slope,dt,field);
    }
    // Run output manager
    om.runme(tw);
  }

  // Now that all the data is written we finalize everything and return the list of files.
  om.finalize();

  // Jumble the order of files to also test the sorting algorithm
  auto list_of_files_tmp = om.get_list_of_files();
  std::vector<std::string> list_of_files;
  for (size_t ii = 1; ii<list_of_files_tmp.size(); ii++) {
    list_of_files.push_back(list_of_files_tmp[ii]);
  }
  list_of_files.push_back(list_of_files_tmp[0]);
  EKAT_REQUIRE(list_of_files.size()==list_of_files_tmp.size());
  return list_of_files;

}
/*-----------------------------------------------------------------------------------------------*/

class NoleapInputSession {
public:
  explicit NoleapInputSession(const ekat::Comm& comm, const bool reverse = false)
    : m_comm(comm), m_leap(use_leap_year()), m_policy(use_noleap_input()), m_reverse(use_gregorian_input())
  {
    scorpio::init_subsystem(comm);
    set_use_noleap_input(false);
    set_use_gregorian_input(false);
    set_use_leap_year(not reverse);
    set_use_noleap_input(not reverse);
    set_use_gregorian_input(reverse);
  }

  ~NoleapInputSession() {
    for (const auto& filename : m_files) {
      while (scorpio::is_file_open(filename)) {
        scorpio::release_file(filename);
      }
    }
    scorpio::finalize_subsystem();
    m_comm.barrier();
    if (m_comm.am_i_root()) {
      for (const auto& filename : m_files) {
        std::remove(filename.c_str());
      }
    }
    m_comm.barrier();
    set_use_noleap_input(false);
    set_use_gregorian_input(false);
    set_use_leap_year(m_leap);
    set_use_noleap_input(m_policy);
    set_use_gregorian_input(m_reverse);
  }

  std::vector<std::string> files(const int epoch_year, const int ncols,
      const std::string& tag, const std::vector<double>& offsets,
      const std::string& calendar = "noleap", const std::string& method = "time: point",
      const bool bounds = false, const bool declare_calendar = true,
      const std::vector<double>& record_values = {}, const std::string& epoch_date = "")
  {
    EKAT_REQUIRE_MSG(record_values.empty() or record_values.size()==offsets.size(),"Invalid fixture values");
    std::vector<std::string> result;
    for (int record = 0; record < static_cast<int>(offsets.size()); ++record) {
      const auto filename = "calendar_noleap_" + std::to_string(m_comm.size()) + "_" +
          std::to_string(epoch_year) + "_" + tag + "_" + std::to_string(record) + ".nc";
      m_files.push_back(filename);
      result.push_back(filename);
      scorpio::register_file(filename,scorpio::FileMode::Write);
      scorpio::define_dim(filename,"ncol",ncols);
      const auto date = epoch_date.empty() ? std::to_string(epoch_year)+"-01-01" : epoch_date;
      scorpio::define_time(filename,"days since " + date + " 00:00:00");
      scorpio::define_var(filename,"forcing",{"ncol"},"real",true);
      if (declare_calendar) {
        scorpio::set_attribute(filename,"time","calendar",calendar);
      }
      if (not method.empty()) {
        scorpio::set_attribute(filename,"forcing","cell_methods",method);
      }
      if (bounds) {
        scorpio::define_dim(filename,"nbnd",2);
        scorpio::define_var(filename,"time_bnds",{"nbnd"},"double",true);
        scorpio::set_attribute(filename,"time","bounds","time_bnds");
      }
      scorpio::set_attribute(filename,"GLOBAL","case_t0",date+"-00000");
      scorpio::enddef(filename);
      scorpio::set_dim_decomp(filename,"ncol");
      const int local_columns = scorpio::get_dimlen_local(filename,"ncol");
      const std::vector<Real> values(local_columns,record_values.empty() ? 10+4*record : record_values[record]);
      scorpio::update_time(filename,offsets[record]);
      scorpio::write_var(filename,"forcing",values.data());
      if (bounds) {
        const double endpoints[2] = {offsets[record]-0.125,offsets[record]+0.125};
        scorpio::write_var(filename,"time_bnds",endpoints);
      }
      scorpio::release_file(filename);
    }
    std::reverse(result.begin(),result.end());
    return result;
  }

private:
  ekat::Comm m_comm;
  bool m_leap;
  bool m_policy;
  bool m_reverse;
  std::vector<std::string> m_files;
};

Field noleap_input_field(const std::shared_ptr<const AbstractGrid>& grid) {
  const FieldIdentifier identifier("forcing",
      FieldLayout({ShortFieldTagsNames::COL},{grid->get_num_local_dofs()}),
      ekat::units::Units::nondimensional(),grid->name());
  Field field(identifier);
  field.allocate_view();
  return field;
}

TEST_CASE("calendar_noleap_file_values", "[calendar_noleap]") {
  const ekat::Comm comm(MPI_COMM_WORLD);
  NoleapInputSession session(comm);
  const int ncols = 2*comm.size()+1;
  const auto gm = get_gm(comm,ncols,2);
  const auto grid = gm->get_grid("point_grid");
  auto field = noleap_input_field(grid);
  for (const int year : {1900,1999,2000,2004,2100,2400}) {
    const bool leap = year%4==0 and (year%100!=0 or year%400==0);
    for (const std::string frequency : {"daily","subdaily","monthly","multiyear"}) {
      INFO("year=" << year << " frequency=" << frequency);
      int epoch = year;
      std::vector<double> offsets = {58,59,60};
      auto query = leap ? util::TimeStamp(year,2,29,0,0,0) : util::TimeStamp(year,2,28,12,0,0);
      Real expected = 12;
      if (frequency=="subdaily") {
        offsets[0] = 58.75;
        if (not leap) query = util::TimeStamp(year,2,28,21,0,0);
        expected = leap ? 10.8 : 12;
      } else if (frequency=="monthly") {
        offsets = {31,59,90};
        query = util::TimeStamp(year,2,leap ? 29 : 28,0,0,0);
        expected = 10 + 4*Real(leap ? 28 : 27)/Real(leap ? 29 : 28);
      } else if (frequency=="multiyear") {
        --epoch;
        for (auto& offset : offsets) offset += 365;
      }
      const auto files = session.files(epoch,ncols,frequency,offsets);
      util::TimeInterpolation reader(grid,files);
      reader.add_field(field);
      reader.initialize_data_from_files();
      reader.perform_time_interpolation(query);
      auto expected_field = field.clone();
      expected_field.deep_copy(expected);
      REQUIRE(views_are_approx_equal(expected_field,reader.get_field("forcing"),tol));
      reader.perform_time_interpolation(util::TimeStamp(year,3,1,0,0,0));
      expected_field.deep_copy(14);
      REQUIRE(views_are_approx_equal(expected_field,reader.get_field("forcing"),tol));
      reader.finalize();
    }
  }
}

TEST_CASE("calendar_noleap_gregorian_controls", "[calendar_noleap]") {
  const ekat::Comm comm(MPI_COMM_WORLD);
  NoleapInputSession session(comm);
  const int ncols = 2*comm.size()+1;
  const auto gm = get_gm(comm,ncols,2);
  const auto grid = gm->get_grid("point_grid");
  auto field = noleap_input_field(grid);
  for (const bool enabled : {false,true}) {
    set_use_noleap_input(enabled);
    const auto files = session.files(2000,ncols,enabled ? "gregorian_on" : "gregorian_off",
        {58,60,61},"gregorian");
    util::TimeInterpolation reader(grid,files);
    reader.add_field(field);
    reader.initialize_data_from_files();
    reader.perform_time_interpolation(util::TimeStamp(2000,2,29,0,0,0));
    auto expected = field.clone();
    expected.deep_copy(12);
    REQUIRE(views_are_approx_equal(expected,reader.get_field("forcing"),tol));
    reader.finalize();
  }
}

TEST_CASE("calendar_historical_source_ranges", "[calendar_noleap]") {
  const ekat::Comm comm(MPI_COMM_WORLD);
  NoleapInputSession session(comm);
  const int ncols = 2*comm.size()+1;
  const auto gm = get_gm(comm,ncols,2);
  const auto grid = gm->get_grid("point_grid");
  auto field = noleap_input_field(grid);
  for (const bool enabled : {false,true}) {
    set_use_noleap_input(enabled);
    for (const std::string flavor : {"standard","gregorian"}) {
      const auto files = session.files(1,ncols,flavor+std::to_string(enabled),
          {730181,730182},flavor,"time: point",false,true,{},"0001-01-01");
      REQUIRE_THROWS_WITH(util::TimeInterpolation(grid,files),
          Catch::Matchers::Contains("unsupported historical calendar="));
    }
    const auto bounds = session.files(1582,ncols,"reform_bounds"+std::to_string(enabled),
        {0,1},"standard","time: point",true,true,{},"1582-10-15");
    REQUIRE_THROWS_WITH(util::TimeInterpolation(grid,bounds),
        Catch::Matchers::Contains("unsupported historical calendar="));
  }
  set_use_noleap_input(false);
  set_use_leap_year(false);
  set_use_gregorian_input(true);
  const auto proleptic = session.files(1,ncols,"historical_proleptic",{730181,730182},
      "proleptic_gregorian","time: point",false,true,{},"0001-01-01");
  util::TimeInterpolation reader(grid,proleptic);
  reader.add_field(field);
  reader.initialize_data_from_files();
  reader.perform_time_interpolation(util::TimeStamp(2000,3,3,0,0,0));
  auto expected = field.clone();
  expected.deep_copy(10);
  REQUIRE(views_are_approx_equal(expected,reader.get_field("forcing"),tol));
  reader.finalize();
}

TEST_CASE("calendar_noleap_metadata_rejection", "[calendar_noleap]") {
  const ekat::Comm comm(MPI_COMM_WORLD);
  NoleapInputSession session(comm);
  const int ncols = 2*comm.size()+1;
  const auto gm = get_gm(comm,ncols,2);
  const auto grid = gm->get_grid("point_grid");
  auto field = noleap_input_field(grid);
    const auto missing_calendar = session.files(2000,ncols,"missing_calendar",
      {58,59,60},"noleap","time: point",false,false);
    REQUIRE_THROWS_WITH(util::TimeInterpolation(grid,missing_calendar),
      Catch::Matchers::Contains("explicit time calendar metadata is required"));
  for (const std::string calendar : {"", "360_day", "julian"}) {
    INFO("calendar=" << calendar);
    const auto files = session.files(2000,ncols,"bad_calendar_"+calendar,{58,59,60},calendar);
    REQUIRE_THROWS_WITH(util::TimeInterpolation(grid,files),
      Catch::Matchers::Contains(calendar.empty() ? "NOLEAP_INPUT: empty calendar metadata" :
                            "NOLEAP_INPUT: unsupported calendar="));
  }
  for (const std::string method : {"time: mean", "time: sum", "time:", "time: point time: mean"}) {
    INFO("method=" << method);
    const auto files = session.files(2000,ncols,"bad_method_"+std::to_string(method.size()),
        {58,59,60},"365_day",method);
    auto initialize = [&]() {
      util::TimeInterpolation reader(grid,files);
      reader.add_field(field);
      reader.initialize_data_from_files();
    };
    REQUIRE_THROWS_WITH(initialize(),Catch::Matchers::Contains("expected a single time: point"));
  }
  const auto files = session.files(2000,ncols,"ambiguous_bounds",{58,59,60},"noleap","",true);
  auto initialize = [&]() {
    util::TimeInterpolation reader(grid,files);
    reader.add_field(field);
    reader.initialize_data_from_files();
  };
  REQUIRE_THROWS_WITH(initialize(),Catch::Matchers::Contains("time bounds without time: point"));
}

TEST_CASE("calendar_noleap_aliases_and_off_control", "[calendar_noleap]") {
  const ekat::Comm comm(MPI_COMM_WORLD);
  NoleapInputSession session(comm);
  const int ncols = 2*comm.size()+1;
  const auto gm = get_gm(comm,ncols,2);
  const auto grid = gm->get_grid("point_grid");
  auto field = noleap_input_field(grid);
  for (const std::string calendar : {"noleap", "no_leap", "365_day", "365day", "NO_LEAP"}) {
    const auto files = session.files(2000,ncols,"alias_"+calendar,{58,59,60},calendar);
    util::TimeInterpolation reader(grid,files);
    reader.add_field(field);
    reader.initialize_data_from_files();
    reader.perform_time_interpolation(util::TimeStamp(2000,2,29,12,0,0));
    auto expected = field.clone();
    expected.deep_copy(13);
    REQUIRE(views_are_approx_equal(expected,reader.get_field("forcing"),tol));
    reader.finalize();
  }
  set_use_noleap_input(false);
  const auto files = session.files(2000,ncols,"legacy_off",{58,59,60});
  util::TimeInterpolation reader(grid,files);
  reader.add_field(field);
  reader.initialize_data_from_files();
  reader.perform_time_interpolation(util::TimeStamp(2000,2,29,0,0,0));
  auto expected = field.clone();
  expected.deep_copy(14);
  REQUIRE(views_are_approx_equal(expected,reader.get_field("forcing"),tol));
  reader.finalize();
}

TEST_CASE("calendar_noleap_offset_and_sequence_rejection", "[calendar_noleap]") {
  const ekat::Comm comm(MPI_COMM_WORLD);
  NoleapInputSession session(comm);
  const int ncols = 2*comm.size()+1;
  const auto gm = get_gm(comm,ncols,2);
  const auto grid = gm->get_grid("point_grid");
  int example = 0;
  for (const double offset : {-1.0, std::numeric_limits<double>::quiet_NaN(),
                             std::numeric_limits<double>::infinity()}) {
    const auto files = session.files(2000,ncols,"invalid_offset_"+std::to_string(example++),
        {offset,59,60});
    REQUIRE_THROWS_WITH(util::TimeInterpolation(grid,files),
      Catch::Matchers::Contains(std::isfinite(offset) ? "expected finite nonnegative source time offset" :
                              "non-finite source time coordinate"));
  }
  const auto duplicate = session.files(2000,ncols,"duplicate",{58,59,59});
  REQUIRE_THROWS_WITH(util::TimeInterpolation(grid,duplicate),
      Catch::Matchers::Contains("same time step has been encountered more than once"));
  auto mixed = session.files(2000,ncols,"mixed_noleap",{58,59,60});
  const auto gregorian = session.files(2000,ncols,"mixed_gregorian",{61,62,63},"gregorian");
  mixed.insert(mixed.end(),gregorian.begin(),gregorian.end());
  REQUIRE_THROWS_WITH(util::TimeInterpolation(grid,mixed),
      Catch::Matchers::Contains("mixed-calendar input sequence"));
}

TEST_CASE("calendar_noleap_coverage_rejection", "[calendar_noleap]") {
  const ekat::Comm comm(MPI_COMM_WORLD);
  NoleapInputSession session(comm);
  const int ncols = 2*comm.size()+1;
  const auto gm = get_gm(comm,ncols,2);
  const auto grid = gm->get_grid("point_grid");
  const auto files = session.files(2000,ncols,"coverage",{58,59,60});
  auto field = noleap_input_field(grid);
  for (const bool lower : {true,false}) {
    util::TimeInterpolation reader(grid,files);
    reader.add_field(field);
    reader.initialize_data_from_files();
    const auto query = lower ? util::TimeStamp(2000,2,27,0,0,0) : util::TimeStamp(2000,3,3,0,0,0);
    REQUIRE_THROWS_WITH(reader.perform_time_interpolation(query),
        Catch::Matchers::Contains(lower ? "lower than the TimeInterpolation bounds" :
                                         "outside the bounds of the set of data files"));
    reader.finalize();
  }
}

TEST_CASE("calendar_gregorian_to_noleap_file_values", "[calendar_gregorian_to_noleap]") {
  const ekat::Comm comm(MPI_COMM_WORLD);
  NoleapInputSession session(comm,true);
  const int ncols = 2*comm.size()+1;
  const auto gm = get_gm(comm,ncols,2);
  const auto grid = gm->get_grid("point_grid");
  auto field = noleap_input_field(grid);
  for (const int year : {1900,1999,2000,2004,2100,2400}) {
    const bool leap = year%4==0 and (year%100!=0 or year%400==0);
    for (const std::string frequency : {"daily","sparse","subdaily","sparse_subdaily","monthly"}) {
      INFO("year=" << year << " frequency=" << frequency);
      std::vector<double> offsets = {58,59,60};
      std::vector<double> values = {10,100,14};
      auto query = util::TimeStamp(year,2,28,12,0,0);
      Real expected = 55;
      if (frequency=="sparse") {
        offsets = {58,leap ? 60.0 : 59.0,61};
        values = {10,14,18};
        expected = leap ? 11 : 12;
      } else if (frequency=="subdaily" or frequency=="sparse_subdaily") {
        offsets = {58.75,frequency=="sparse_subdaily" and leap ? 60.0 : 59.0,61};
        values = {10,14,18};
        query = util::TimeStamp(year,2,28,21,0,0);
        expected = frequency=="sparse_subdaily" and leap ? 10.4 : 12;
      } else if (frequency=="monthly") {
        offsets = {31,leap ? 60.0 : 59.0,90};
        values = {10,14,18};
        query = util::TimeStamp(year,2,28,0,0,0);
        expected = 10+4*Real(27)/Real(leap ? 29 : 28);
      }
      const auto files = session.files(year,ncols,frequency,offsets,"gregorian","time: point",false,true,values);
      util::TimeInterpolation reader(grid,files);
      reader.add_field(field);
      reader.initialize_data_from_files();
      reader.perform_time_interpolation(query);
      auto reference = field.clone();
      reference.deep_copy(expected);
      REQUIRE(views_are_approx_equal(reference,reader.get_field("forcing"),tol));
      if (frequency=="daily") {
        reader.perform_time_interpolation(util::TimeStamp(year,3,1,0,0,0));
        reference.deep_copy(leap ? 14 : 100);
        REQUIRE(views_are_approx_equal(reference,reader.get_field("forcing"),tol));
      }
      reader.finalize();
    }
  }
  REQUIRE_THROWS(util::TimeStamp(2000,2,29,0,0,0));
}

TEST_CASE("calendar_gregorian_to_noleap_epoch_and_controls", "[calendar_gregorian_to_noleap]") {
  const ekat::Comm comm(MPI_COMM_WORLD);
  NoleapInputSession session(comm,true);
  const int ncols = 2*comm.size()+1;
  const auto gm = get_gm(comm,ncols,2);
  const auto grid = gm->get_grid("point_grid");
  auto field = noleap_input_field(grid);
  for (const bool leap_epoch : {false,true}) {
    const auto files = leap_epoch ?
        session.files(2000,ncols,"leap_epoch",{0,1,2},"gregorian","time: point",false,true,{},"2000-02-29") :
        session.files(2000,ncols,"matching_noleap",{58,59,60});
    util::TimeInterpolation reader(grid,files);
    reader.add_field(field);
    reader.initialize_data_from_files();
    reader.perform_time_interpolation(util::TimeStamp(2000,3,1,12,0,0));
    auto reference = field.clone();
    reference.deep_copy(16);
    REQUIRE(views_are_approx_equal(reference,reader.get_field("forcing"),tol));
    reader.finalize();
  }
  const auto invalid = session.files(2001,ncols,"invalid_epoch",{0,1},"gregorian","time: point",false,true,{},"2001-02-29");
  REQUIRE_THROWS_WITH(util::TimeInterpolation(grid,invalid),Catch::Matchers::Contains("invalid source date"));
  REQUIRE_THROWS_WITH(set_use_noleap_input(true),Catch::Matchers::Contains("mutually exclusive"));
}

TEST_CASE("calendar_gregorian_to_noleap_coverage", "[calendar_gregorian_to_noleap]") {
  const ekat::Comm comm(MPI_COMM_WORLD);
  NoleapInputSession session(comm,true);
  const int ncols = 2*comm.size()+1;
  const auto gm = get_gm(comm,ncols,2);
  const auto grid = gm->get_grid("point_grid");
  auto field = noleap_input_field(grid);
  const auto files = session.files(2000,ncols,"reverse_coverage",{58,59,60},"gregorian");
  util::TimeInterpolation reader(grid,files);
  reader.add_field(field);
  reader.initialize_data_from_files();
  REQUIRE_THROWS_WITH(reader.perform_time_interpolation(util::TimeStamp(2000,2,27,0,0,0)),
      Catch::Matchers::Contains("lower than the TimeInterpolation bounds"));
  REQUIRE_THROWS_WITH(reader.perform_time_interpolation(util::TimeStamp(2000,3,2,0,0,0)),
      Catch::Matchers::Contains("outside the bounds of the set of data files"));
  reader.finalize();
}

TEST_CASE("calendar_gregorian_to_noleap_seconds", "[calendar_gregorian_to_noleap]") {
  REQUIRE(gregorian_input_seconds(20000101,0)==0);
  REQUIRE(gregorian_input_seconds(20000229,0)==59.0*86400);
  REQUIRE(gregorian_input_seconds(20000301,0)==60.0*86400);
  REQUIRE(gregorian_input_seconds(19990101,0)==-365.0*86400);
  REQUIRE(gregorian_input_seconds(19000101,0)==-36524.0*86400);
  REQUIRE(gregorian_input_seconds(21000101,0)==36525.0*86400);
  REQUIRE(gregorian_input_seconds(24000101,0)==146097.0*86400);
  REQUIRE(gregorian_input_seconds(20000101,0.5)==0.5);
  for (const int invalid : {19000229,21000229,20010229,20001301,20000100}) {
    REQUIRE_THROWS(gregorian_input_seconds(invalid,0));
  }
}

} // namespace scream
