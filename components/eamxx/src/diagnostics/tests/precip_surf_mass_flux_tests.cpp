#include "catch2/catch.hpp"

#include "share/grid/mesh_free_grids_manager.hpp"
#include "diagnostics/register_diagnostics.hpp"

#include "physics/share/physics_constants.hpp"

#include "share/util/eamxx_setup_random_test.hpp"
#include "share/field/field_utils.hpp"

#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/util/ekat_test_utils.hpp"

namespace scream {

std::shared_ptr<GridsManager>
create_gm (const ekat::Comm& comm, const int ncols) {

  const int num_global_cols = ncols*comm.size();

  using vos_t = std::vector<std::string>;
  ekat::ParameterList gm_params;
  gm_params.set("grids_names",vos_t{"Point Grid"});
  auto& pl = gm_params.sublist("Point Grid");
  pl.set<std::string>("type","point_grid");
  pl.set("aliases",vos_t{"Physics"});
  pl.set<int>("number_of_global_columns", num_global_cols);
  pl.set<int>("number_of_vertical_levels", 1);

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
  register_diagnostics();
  ekat::ParameterList params;
  auto& diag_factory = AtmosphereDiagnosticFactory::instance();
  params.set<std::string>("precip_type","total");
  auto diag_total = diag_factory.create("precip_surf_mass_flux", comm, params);
  params.set<std::string>("precip_type","ice");
  auto diag_ice   = diag_factory.create("precip_surf_mass_flux"  , comm, params);
  params.set<std::string>("precip_type","liq");
  auto diag_liq   = diag_factory.create("precip_surf_mass_flux"  , comm, params);
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
  util::TimeStamp t = t0 + dt;
  
  // Construct random data to use for test
  // Get views of input data and set to random values
  auto precip_ice_surf_mass_f = input_fields["precip_ice_surf_mass"];
  auto precip_ice_surf_mass_v = precip_ice_surf_mass_f.get_view<Real*>();
  auto precip_liq_surf_mass_f = input_fields["precip_liq_surf_mass"];
  auto precip_liq_surf_mass_v = precip_liq_surf_mass_f.get_view<Real*>();
  for (int icol=0;icol<ncols;icol++) {
    ekat::genRandArray(precip_ice_surf_mass_v, engine, pdf_mass);
    ekat::genRandArray(precip_liq_surf_mass_v, engine, pdf_mass);
  }

  precip_ice_surf_mass_f.get_header().get_tracking().update_time_stamp(t);
  precip_liq_surf_mass_f.get_header().get_tracking().update_time_stamp(t);
  
  // Run diagnostics and compare with manual calculation
  diag_total->compute_diagnostic();
  diag_liq->compute_diagnostic();
  diag_ice->compute_diagnostic();

  Field preicp_total_f = diag_total->get_diagnostic().clone();
  Field preicp_liq_f   = diag_liq->get_diagnostic().clone();
  Field preicp_ice_f   = diag_ice->get_diagnostic().clone();
  preicp_total_f.deep_copy<double>(0.0);
  preicp_liq_f.deep_copy<double>(0.0);
  preicp_ice_f.deep_copy<double>(0.0);
  auto precip_total_v = preicp_total_f.get_view<Real*>();
  auto precip_liq_v   = preicp_liq_f.get_view<Real*>();
  auto precip_ice_v   = preicp_ice_f.get_view<Real*>();
  const auto rhodt = PC::RHO_H2O*dt;
  Kokkos::parallel_for("precip_total_surf_mass_flux_test",
                       typename KT::RangePolicy(0,ncols),
                       KOKKOS_LAMBDA(const int& icol) {
    precip_liq_v(icol)   = precip_liq_surf_mass_v(icol)/rhodt;
    precip_ice_v(icol)   = precip_ice_surf_mass_v(icol)/rhodt;
    precip_total_v(icol) = precip_liq_v(icol) + precip_ice_v(icol);
  });
  Kokkos::fence();

  REQUIRE(views_are_equal(diag_total->get_diagnostic(),preicp_total_f));
  REQUIRE(views_are_equal(diag_liq->get_diagnostic(),preicp_liq_f));
  REQUIRE(views_are_equal(diag_ice->get_diagnostic(),preicp_ice_f));

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
