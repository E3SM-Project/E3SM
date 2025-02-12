#include <catch2/catch.hpp>

#include "share/grid/mesh_free_grids_manager.hpp"

#include "share/field/field_utils.hpp"
#include "share/field/field.hpp"
#include "share/field/field_manager.hpp"

#include "share/util/eamxx_time_interpolation.hpp"
#include "share/util/eamxx_setup_random_test.hpp"
#include "share/util/eamxx_time_stamp.hpp"

#include "share/io/eamxx_output_manager.hpp"

#include "ekat/ekat_parameter_list.hpp"
/*-----------------------------------------------------------------------------------------------
 * Test TimeInterpolation class
 *-----------------------------------------------------------------------------------------------*/

namespace scream {

// Test Constants
constexpr Real tol            = std::numeric_limits<Real>::epsilon()*1e5;
constexpr int  slp_switch     = 4; // The frequency that we change the slope of the data written to file.
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
    if (std::find(m_list_of_files.begin(),m_list_of_files.end(), m_output_file_specs.filename) == m_list_of_files.end()) {
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
  const auto& grid = grids_man->get_grid("Point Grid");
  // Now create a fields manager to store initial data for testing.
  auto fields_man_t0 = get_fm(grid, t0, seed);
  // Construct a time interpolation object and add all of the fields to it.
  printf(  "Constructing a time interpolation object ...\n");
  util::TimeInterpolation time_interpolator(grid);
  for (auto ff_pair = fields_man_t0->begin(); ff_pair != fields_man_t0->end(); ff_pair++) {
    const auto ff   = ff_pair->second;
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
  for (auto ff_pair = fields_man_tf->begin(); ff_pair != fields_man_tf->end(); ff_pair++)
  {
    auto  ff = ff_pair->second;
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
    for (auto ff_pair = fields_man_test->begin(); ff_pair != fields_man_test->end(); ff_pair++)
    {
      const auto name = ff_pair->first;
      auto ff = ff_pair->second;
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
  const auto& grid = grids_man->get_grid("Point Grid");
  printf("   - Grids Manager...DONE\n");
  // Now create a fields manager to store initial data for testing.
  printf("   - Fields Manager...\n");
  auto fields_man_t0 = get_fm(grid, t0, seed);
  auto fields_man_deep = get_fm(grid, t0, seed);  // A field manager for checking deep copies.
  std::vector<std::string> fnames;
  for (auto it : *fields_man_t0) {
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
  auto d_min = field_min<Real>(ft);
  auto d_max = field_max<Real>(ft);
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
  gm_params.set("grids_names",vos_t{"Point Grid"});
  auto& pl = gm_params.sublist("Point Grid");
  pl.set<std::string>("type","point_grid");
  pl.set("aliases",vos_t{"Physics"});
  pl.set<int>("number_of_global_columns", ncols);
  pl.set<int>("number_of_vertical_levels", nlevs);
  auto gm = create_mesh_free_grids_manager(comm,gm_params);
  gm->build_grids();
  return gm;
}
/*-----------------------------------------------------------------------------------------------*/
/* Create a fields manager for the test */
std::shared_ptr<FieldManager> get_fm (const std::shared_ptr<const AbstractGrid>& grid, const util::TimeStamp& t0, const int seed)
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

  const int nlcols = grid->get_num_local_dofs();
  const int nlevs  = grid->get_num_vertical_levels();

  std::vector<FL> layouts =
  {
    FL({COL         }, {nlcols        }),
    FL({COL,     LEV}, {nlcols,  nlevs}),
    FL({COL,CMP,ILEV}, {nlcols,2,nlevs+1})
  };

  auto fm = std::make_shared<FieldManager>(grid);
  fm->registration_begins();
  fm->registration_ends();

  const auto units = ekat::units::Units::nondimensional();
  for (const auto& fl : layouts) {
    int gl_size = fl.size();
    grid->get_comm().all_reduce(&gl_size,1,MPI_SUM);
    FID fid("f_"+std::to_string(gl_size),fl,units,grid->name());
    Field f(fid);
    f.allocate_view();
    randomize (f,engine,my_pdf);
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
  auto fm = get_fm(gm->get_grid("Point Grid"), t0, seed);
  // We will write data for 10 snaps for this test.  We set the max snaps per file to 3 to
  // ensure that a) there is more than 1 file and b) at least one file has fewer snap then
  // the others.
  const int max_steps      = snap_freq*total_snaps;
  // Gather the set of fields from the field manager
  std::vector<std::string> fnames;
  for (auto it : *fm) {
    fnames.push_back(it.second->name());
  }
  // Create the output parameters
  ekat::ParameterList om_pl;
  om_pl.set("filename_prefix",std::string("source_data_for_time_interpolation"));
  om_pl.set("Field Names",fnames);
  om_pl.set("Averaging Type", std::string("INSTANT"));
  om_pl.set("Max Snapshots Per File",snaps_per_file);
  auto& ctrl_pl = om_pl.sublist("output_control");
  ctrl_pl.set("frequency_units",std::string("nsteps"));
  ctrl_pl.set("Frequency",snap_freq);
  ctrl_pl.set("save_grid_data",false);
  // Create an output manager, note we use a subclass defined in this test so we can extract
  // the list of files created by the output manager.
  OutputManager4Test om;
  om.initialize(comm,om_pl,t0,false);
  om.setup(fm,gm);

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

} // namespace scream
