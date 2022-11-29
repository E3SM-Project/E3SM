#include "catch2/catch.hpp"

#include "share/grid/mesh_free_grids_manager.hpp"
#include "diagnostics/precip_liq_surf_mass_flux.hpp"
#include "diagnostics/register_diagnostics.hpp"

#include "physics/share/physics_constants.hpp"

#include "share/util/scream_setup_random_test.hpp"
#include "share/field/field_utils.hpp"

#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/util/ekat_test_utils.hpp"

#include <iomanip>

namespace scream {

std::shared_ptr<GridsManager>
create_gm (const ekat::Comm& comm, const int ncols) {

  const int num_global_cols = ncols*comm.size();

  ekat::ParameterList gm_params;
  gm_params.set<int>("number_of_global_columns", num_global_cols);
  gm_params.set<int>("number_of_vertical_levels", 1);

  auto gm = create_mesh_free_grids_manager(comm,gm_params);
  gm->build_grids();

  return gm;
}

//-----------------------------------------------------------------------------------------------//
template<typename DeviceT>
void run(std::mt19937_64& engine)
{
  using PC = scream::physics::Constants<Real>;
  using KT = ekat::KokkosTypes<DeviceT>;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // Create a grids manager
  const int ncols = 13;
  auto gm = create_gm(comm,ncols);

  // Create timestep
  const int dt=1800;

  // Construct random input data
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF pdf_mass(0.0,1.0);

  // A time stamp
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Construct the Diagnostic
  ekat::ParameterList params;
  register_diagnostics();
  auto& diag_factory = AtmosphereDiagnosticFactory::instance();
  auto diag = diag_factory.create("PrecipLiqSurfMassFlux",comm,params);
  diag->set_grids(gm);

  // Set the required fields for the diagnostic.
  std::map<std::string,Field> input_fields;
  for (const auto& req : diag->get_required_field_requests()) {
    Field f(req.fid);
    auto & f_ap = f.get_header().get_alloc_properties();
    f_ap.request_allocation();
    f.allocate_view();
    const auto name = f.name();
    f.get_header().get_tracking().update_time_stamp(t0);
    diag->set_required_field(f.get_const());
    REQUIRE_THROWS(diag->set_computed_field(f));
    input_fields.emplace(name,f);
  }

  // Initialize the diagnostic
  diag->initialize(t0,RunType::Initial);

  // Run tests
  {
    // Construct random data to use for test
    // Get views of input data and set to random values
    const auto& precip_liq_surf_mass_f = input_fields["precip_liq_surf_mass"];
    const auto& precip_liq_surf_mass_v = precip_liq_surf_mass_f.get_view<Real*>();
    for (int icol=0;icol<ncols;icol++) {
      ekat::genRandArray(precip_liq_surf_mass_v, engine, pdf_mass);
    }

    // Run diagnostic and compare with manual calculation
    diag->compute_diagnostic(dt);
    const auto& diag_out = diag->get_diagnostic();
    Field theta_f = diag_out.clone();
    theta_f.deep_copy<double,Host>(0.0);
    theta_f.sync_to_dev();
    const auto& theta_v = theta_f.get_view<Real*>();
    const auto rho_h2o = PC::RHO_H2O;
    Kokkos::parallel_for("precip_liq_surf_mass_flux_test",
                         typename KT::RangePolicy(0,ncols),
                         KOKKOS_LAMBDA(const int& icol) {
      theta_v(icol) = precip_liq_surf_mass_v(icol)/rho_h2o/dt;
    });
    Kokkos::fence();
    REQUIRE(views_are_equal(diag_out,theta_f));
  }
 
  // Finalize the diagnostic
  diag->finalize(); 

} // run()

TEST_CASE("precip_liq_surf_mass_flux_test", "precip_liq_surf_mass_flux_test]"){
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
