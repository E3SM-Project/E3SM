#include "catch2/catch.hpp"

#include "share/grid/point_grid.hpp"
#include "share/diagnostics/register_diagnostics.hpp"

#include "share/physics/physics_constants.hpp"

#include "share/core/eamxx_setup_random_test.hpp"
#include "share/field/field_utils.hpp"

#include <ekat_view_utils.hpp>

namespace scream {


//-----------------------------------------------------------------------------------------------//
template<typename DeviceT>
void run(std::mt19937_64& engine)
{
  using PC = scream::physics::Constants<Real>;
  using KT = ekat::KokkosTypes<DeviceT>;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // For input randomization
  int seed = get_random_test_seed(&comm);

  // Create a grids manager
  const int ncols = 13;
  auto grid = create_point_grid("physics",ncols,1,comm);

  // Create timestep
  const double dt=1800;

  // Initial time stamp
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Construct the Diagnostics
  register_diagnostics();
  ekat::ParameterList params;
  auto& diag_factory = DiagnosticFactory::instance();
  params.set<std::string>("precip_type","total");
  auto diag_total = diag_factory.create("precip_surf_mass_flux", comm, params, grid);
  params.set<std::string>("precip_type","ice");
  auto diag_ice   = diag_factory.create("precip_surf_mass_flux"  , comm, params, grid);
  params.set<std::string>("precip_type","liq");
  auto diag_liq   = diag_factory.create("precip_surf_mass_flux"  , comm, params, grid);
  // Set the required fields for the diagnostic.
  std::map<std::string,Field> input_fields;
  {
    using namespace ekat::units;
    auto scalar2d = grid->get_2d_scalar_layout();
    for (const auto& fname : diag_total->get_input_fields_names()) {
      FieldIdentifier fid(fname, scalar2d, kg/(m*m*s), grid->name());
      Field f(fid);
      f.allocate_view();
      f.get_header().get_tracking().update_time_stamp(t0);
      f.get_header().get_tracking().set_accum_start_time(t0);
      diag_total->set_input_field(f);
      if (fname=="precip_ice_surf_mass") {
        diag_ice->set_input_field(f);
      }
      if (fname=="precip_liq_surf_mass") {
        diag_liq->set_input_field(f);
      }
      input_fields.emplace(fname, f);
    }
  }

  // Initialize the diagnostic
  diag_total->initialize();
  diag_ice->initialize();
  diag_liq->initialize();

  // Run tests
  util::TimeStamp t = t0 + dt;
  
  // Construct random data to use for test
  // Get views of input data and set to random values
  auto precip_ice_surf_mass_f = input_fields["precip_ice_surf_mass"];
  auto precip_liq_surf_mass_f = input_fields["precip_liq_surf_mass"];

  randomize_uniform(precip_ice_surf_mass_f,seed++,0,1);
  randomize_uniform(precip_liq_surf_mass_f,seed++,0,1);

  precip_ice_surf_mass_f.get_header().get_tracking().update_time_stamp(t);
  precip_liq_surf_mass_f.get_header().get_tracking().update_time_stamp(t);
  
  // Run diagnostics and compare with manual calculation
  diag_total->compute(t0);
  diag_liq->compute(t0);
  diag_ice->compute(t0);

  Field precip_total_f = diag_total->get().clone();
  Field precip_liq_f   = diag_liq->get().clone();
  Field precip_ice_f   = diag_ice->get().clone();
  precip_total_f.deep_copy(0);
  precip_liq_f.deep_copy(0);
  precip_ice_f.deep_copy(0);
  auto precip_total_v = precip_total_f.get_view<Real*>();
  auto precip_liq_v   = precip_liq_f.get_view<Real*>();
  auto precip_ice_v   = precip_ice_f.get_view<Real*>();
  auto precip_ice_surf_mass_v = precip_ice_surf_mass_f.get_view<const Real*>();
  auto precip_liq_surf_mass_v = precip_liq_surf_mass_f.get_view<const Real*>();
  const Real rhodt = PC::RHO_H2O.value*dt;
  Kokkos::parallel_for("precip_total_surf_mass_flux_test",
                       typename KT::RangePolicy(0,ncols),
                       KOKKOS_LAMBDA(const int& icol) {
    precip_liq_v(icol)   = precip_liq_surf_mass_v(icol)/rhodt;
    precip_ice_v(icol)   = precip_ice_surf_mass_v(icol)/rhodt;
    precip_total_v(icol) = precip_ice_v(icol);
    precip_total_v(icol) += precip_liq_surf_mass_v(icol)/rhodt;
  });
  Kokkos::fence();

  REQUIRE(views_are_equal(diag_total->get(),precip_total_f));
  REQUIRE(views_are_equal(diag_liq->get(),precip_liq_f));
  REQUIRE(views_are_equal(diag_ice->get(),precip_ice_f));

} // run()

TEST_CASE("precip_total_surf_mass_flux_test", "precip_total_surf_mass_flux_test]"){
  using scream::Real;
  using Device = scream::DefaultDevice;

  constexpr int num_runs = 5;

  auto engine = scream::setup_random_test();

  printf(" -> Number of randomized runs: %d\n\n", num_runs);

  for (int irun=0; irun<num_runs; ++irun) {
    run<Device>(engine);
  }
  printf("ok!\n");

  printf("\n");

} // TEST_CASE

} // namespace
