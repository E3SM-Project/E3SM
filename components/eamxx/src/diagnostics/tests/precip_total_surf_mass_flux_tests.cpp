#include "catch2/catch.hpp"

#include "share/grid/mesh_free_grids_manager.hpp"
#include "diagnostics/precip_total_surf_mass_flux.hpp"
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

  // Initial time stamp
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Construct the Diagnostics
  ekat::ParameterList params;
  register_diagnostics();
  auto& diag_factory = AtmosphereDiagnosticFactory::instance();
  auto diag_total = diag_factory.create("PrecipTotalSurfMassFlux", comm, params);
  auto diag_ice   = diag_factory.create("PrecipIceSurfMassFlux",   comm, params);
  auto diag_liq   = diag_factory.create("PrecipLiqSurfMassFlux",   comm, params);
  diag_total->set_grids(gm);
  diag_ice->set_grids(gm);
  diag_liq->set_grids(gm);

  // Set the required fields for the diagnostic.
  std::map<std::string,Field> input_fields;
  for (const auto& req : diag_total->get_required_field_requests()) {
    Field f(req.fid);
    f.get_header().get_alloc_properties().request_allocation();
    f.allocate_view();
    const auto name = f.name();
    f.get_header().get_tracking().update_time_stamp(t0);
    f.get_header().get_tracking().set_accum_start_time(t0);
    diag_total->set_required_field(f.get_const());
    REQUIRE_THROWS(diag_total->set_computed_field(f));
    if (name=="precip_ice_surf_mass") {
      diag_ice->set_required_field(f.get_const());
      REQUIRE_THROWS(diag_ice->set_computed_field(f));
    }
    if (name=="precip_liq_surf_mass") {
      diag_liq->set_required_field(f.get_const());
      REQUIRE_THROWS(diag_liq->set_computed_field(f));
    }
    input_fields.emplace(name,f);
  }

  // Initialize the diagnostic
  diag_total->initialize(t0,RunType::Initial);
  diag_ice->initialize(t0,RunType::Initial);
  diag_liq->initialize(t0,RunType::Initial);

  // Run tests
  {
    util::TimeStamp t = t0 + dt;
    
    // Construct random data to use for test
    // Get views of input data and set to random values
    auto precip_ice_surf_mass_f = input_fields["precip_ice_surf_mass"];
    const auto& precip_ice_surf_mass_v = precip_ice_surf_mass_f.get_view<Real*>();
    auto precip_liq_surf_mass_f = input_fields["precip_liq_surf_mass"];
    const auto& precip_liq_surf_mass_v = precip_liq_surf_mass_f.get_view<Real*>();
    for (int icol=0;icol<ncols;icol++) {
      ekat::genRandArray(precip_ice_surf_mass_v, engine, pdf_mass);
      ekat::genRandArray(precip_liq_surf_mass_v, engine, pdf_mass);
    }

    precip_ice_surf_mass_f.get_header().get_tracking().update_time_stamp(t);
    precip_liq_surf_mass_f.get_header().get_tracking().update_time_stamp(t);
    
    // Run diagnostic and compare with manual calculation
    diag_total->compute_diagnostic();
    const auto& diag_total_out = diag_total->get_diagnostic();
    Field theta_f = diag_total_out.clone();
    theta_f.deep_copy<double,Host>(0.0);
    theta_f.sync_to_dev();
    const auto& theta_v = theta_f.get_view<Real*>();
    const auto rhodt = PC::RHO_H2O*dt;
    Kokkos::parallel_for("precip_total_surf_mass_flux_test",
                         typename KT::RangePolicy(0,ncols),
                         KOKKOS_LAMBDA(const int& icol) {
      theta_v(icol) = precip_ice_surf_mass_v(icol)/rhodt +
                      precip_liq_surf_mass_v(icol)/rhodt;
    });
    Kokkos::fence();
    REQUIRE(views_are_equal(diag_total_out,theta_f));

    // Test against sum of precip_ice/liq diagnostics
    diag_ice->compute_diagnostic();
    diag_liq->compute_diagnostic();

    const auto& diag_ice_out = diag_ice->get_diagnostic();
    const auto& diag_liq_out = diag_liq->get_diagnostic();

    Field sum_of_diags_f(diag_total_out.get_header().get_identifier());
    sum_of_diags_f.get_header().get_alloc_properties().request_allocation();
    sum_of_diags_f.allocate_view();
    const auto& sum_of_diags_v = sum_of_diags_f.get_view<Real*>();

    const auto& diag_ice_v     = diag_ice_out.get_view<const Real*>();
    const auto& diag_liq_v     = diag_liq_out.get_view<const Real*>();
    Kokkos::parallel_for("calculate_sum_of_diags",
                        typename KT::RangePolicy(0,ncols),
                        KOKKOS_LAMBDA(const int& icol) {
      sum_of_diags_v(icol) = diag_ice_v(icol) + diag_liq_v(icol);
    });
    Kokkos::fence();
    REQUIRE(views_are_equal(diag_total_out,sum_of_diags_f));
}

  // Finalize the diagnostic
  diag_total->finalize();
  diag_ice->finalize();
  diag_liq->finalize();
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
