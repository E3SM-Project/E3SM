#include <catch2/catch.hpp>

#include "share/grid/mesh_free_grids_manager.hpp"

#include "share/field/field_utils.hpp"
#include "share/field/field.hpp"
#include "share/field/field_manager.hpp"

#include "share/util/eamxx_time_interpolation.hpp"
#include "share/util/scream_setup_random_test.hpp"
#include "share/util/scream_time_stamp.hpp"

#include "ekat/ekat_parameter_list.hpp"
/*-----------------------------------------------------------------------------------------------
 * Test TimeInterpolation class
 *-----------------------------------------------------------------------------------------------*/

namespace scream {

// Test Constants
constexpr Real tol = std::numeric_limits<Real>::epsilon()*1e5;

// Functions needed to set the test up
std::shared_ptr<const GridsManager> get_gm (const ekat::Comm& comm, const int ncols, const int nlevs);
std::shared_ptr<FieldManager> get_fm (const std::shared_ptr<const AbstractGrid>& grid, const util::TimeStamp& t0, const int seed);
util::TimeStamp init_timestamp();

// Functions needed to run the test
void update_field_data(const Real slope, const Real dt, Field& field);
bool views_are_approx_equal(const Field& f0, const Field& f1, const Real tol);

// Helper randomizer
auto my_pdf = [&](std::mt19937_64& engine) -> Real {
  std::uniform_int_distribution<int> pdf (1,100);
  Real v = pdf(engine);
  return v;
};
/*-----------------------------------------------------------------------------------------------*/
TEST_CASE ("eamxx_time_interpolation_simple") {
  // Setup basic test params
  ekat::Comm comm(MPI_COMM_WORLD);
  auto seed = get_random_test_seed(&comm);
  std::mt19937_64 engine(seed); 
  const auto t0 = init_timestamp();

  const int nlevs  = SCREAM_PACK_SIZE*2+1;
  const int ncols  = comm.size()*2 + 1;
  const int dt     = 100;

  // Get a grids manager for the test
  auto grids_man = get_gm(comm, ncols, nlevs);
  const auto& grid = grids_man->get_grid("Point Grid");
  // Now create a fields manager to store initial data for testing.
  auto fields_man_t0 = get_fm(grid, t0, seed);
  // Construct a time interpolation object and add all of the fields to it.
  printf("Constructing a time interpolation object ...\n");
  util::TimeInterpolation time_interpolator(grid);
  for (auto ff_pair = fields_man_t0->begin(); ff_pair != fields_man_t0->end(); ff_pair++) {
    const auto ff   = ff_pair->second;
    time_interpolator.add_field(*ff);
    time_interpolator.initialize_data_from_field(*ff);
  }
  time_interpolator.initialize_timestamps(t0);
  printf("Constructing a time interpolation object ... DONE\n");

  // Set the field data for a step 10 dt in the future to be used for interpolation.  We
  // create a new field manager so we can continue to have the intial condition for reference.
  printf("Setting data for other end of time interpolation, at t = 10 dt ...\n");
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
  printf("Setting data for other end of time interpolation, at t = 10 dt ...DONE\n");

  // Now check that the interpolator is working as expected.  Should be able to
  // match the interpolated fields against the results of update_field_data at any
  // time between 0 and 10 dt.
  printf("Testing all timesteps ... slope = %f\n",slope);
  auto fields_man_test = get_fm(grid, t0, seed);
  t1 = t0;
  for (int tt = 1; tt<=10; tt++) {
    t1 += dt;
    printf("                      ... t = %s\n",t1.to_string().c_str());
    auto interp_fields = time_interpolator.perform_time_interpolation(t1);
    for (auto ff_pair = fields_man_test->begin(); ff_pair != fields_man_test->end(); ff_pair++)
    {
      const auto name = ff_pair->first;
      auto fcmp = interp_fields.at(name);
      auto ff = ff_pair->second;
      update_field_data(slope, dt, *ff);
      REQUIRE(views_are_approx_equal(*ff,interp_fields.at(name),tol));
    }
  }
  printf("                      ... DONE\n");


} // TEST_CASE eamxx_time_interpolation_simple
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
  ekat::ParameterList gm_params;
  gm_params.set("number_of_global_columns",ncols);
  gm_params.set("number_of_vertical_levels",nlevs);
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
    FID fid("f_"+std::to_string(fl.size()),fl,units,grid->name());
    Field f(fid);
    f.allocate_view();
    randomize (f,engine,my_pdf);
    f.get_header().get_tracking().update_time_stamp(t0);
    fm->add_field(f);
  }

  return fm;
}
/*-----------------------------------------------------------------------------------------------*/

} // namespace scream
