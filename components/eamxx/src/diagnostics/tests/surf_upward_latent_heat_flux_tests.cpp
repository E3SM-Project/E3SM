#include "catch2/catch.hpp"

#include "share/grid/mesh_free_grids_manager.hpp"
#include "diagnostics/surf_upward_latent_heat_flux.hpp"
#include "diagnostics/register_diagnostics.hpp"

#include "physics/share/physics_constants.hpp"

#include "share/util/eamxx_setup_random_test.hpp"
#include "share/field/field_utils.hpp"

#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/util/ekat_test_utils.hpp"
#include "ekat/logging/ekat_logger.hpp"

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
template<typename DeviceT, typename LoggerType>
void run(std::mt19937_64& engine, const ekat::Comm& comm, LoggerType& logger)
{
  using PC         = scream::physics::Constants<Real>;
  using KT         = ekat::KokkosTypes<DeviceT>;

  // Create a grids manager
  const int ncols = 13;
  auto gm = create_gm(comm,ncols);

  // Construct random input data
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF pdf_mass(0.0,1.0);

  // Initial time stamp
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Construct the Diagnostics
  register_diagnostics();
  ekat::ParameterList params;

  auto& diag_factory = AtmosphereDiagnosticFactory::instance();
  auto diag_latent_heat = diag_factory.create("surface_upward_latent_heat_flux", comm, params);
  diag_latent_heat->set_grids(gm);

  // Set the required fields for the diagnostic.
  std::map<std::string,Field> input_fields;
  for (const auto& req: diag_latent_heat->get_required_field_requests()) {
    Field f(req.fid);
    f.get_header().get_alloc_properties().request_allocation();
    f.allocate_view();
    const auto name = f.name();
    f.get_header().get_tracking().update_time_stamp(t0);
    diag_latent_heat->set_required_field(f.get_const());
    input_fields.emplace(name, f);
  }

  // Initialize the diagnostic
  diag_latent_heat->initialize(t0,RunType::Initial);
  const auto& surf_evap_f = input_fields["surf_evap"];
  const auto& surf_evap_v = surf_evap_f.get_view<Real*>();
  for (int icol=0; icol<ncols; ++icol) {
   ekat::genRandArray(surf_evap_v, engine, pdf_mass);
  }

  // Run diagnostic and compare with manual calculation
  diag_latent_heat->compute_diagnostic();
  const auto diag_latent_heat_out = diag_latent_heat->get_diagnostic();
  Field surf_lhf = diag_latent_heat_out.clone();
  surf_lhf.deep_copy(0);
  const auto& surf_lhf_v = surf_lhf.get_view<Real*>();
  constexpr auto latent_heat_evap = PC::LatVap; // [J/kg]
  Kokkos::parallel_for("surf_upward_latent_heat_flux_test",
    typename KT::RangePolicy(0, ncols),
    KOKKOS_LAMBDA (const int& icol) {
      surf_lhf_v(icol) = surf_evap_v(icol) * latent_heat_evap;
    });
  Kokkos::fence();

  if (!views_are_equal(diag_latent_heat_out, surf_lhf)) {
    // In case of failure, log additional info before aborting with
    // Catch2's REQUIRE macro
    logger.error("error: surf_lhf_v and diag_latent_heat_out are not passing the views_are_equal test.");
    auto surf_lhf_h = Kokkos::create_mirror_view(surf_lhf_v);
    diag_latent_heat_out.sync_to_host();
    auto diag_latent_heat_out_h = diag_latent_heat->get_diagnostic().get_view<Real*,Host>();
    Kokkos::deep_copy(surf_lhf_h, surf_lhf_v);
    for (int i=0; i<ncols; ++i) {
      logger.debug("\tat col {}: diag_latent_heat_out = {} surf_lhf = {}", i,
        diag_latent_heat_out_h(i), surf_lhf_h(i));
    }
  }
  REQUIRE(views_are_equal(diag_latent_heat_out, surf_lhf));

  // Finalize the diagnostic
  diag_latent_heat->finalize();
} // run

TEST_CASE("surf_upward_latent_heat_flux_test", "[surf_upward_latent_heat_flux_test]") {
  using scream::Real;
  using Device = scream::DefaultDevice;

  ekat::Comm comm(MPI_COMM_WORLD);
  ekat::logger::Logger<> logger("surf_upward_latent_heat_flux_test",
                                ekat::logger::LogLevel::debug, comm);
  constexpr int num_runs = 5;
  auto engine = scream::setup_random_test();
  logger.info(" -> Number of randomized runs: {}", num_runs);
  for (int irun=0; irun<num_runs; ++irun) {
    run<Device>(engine, comm, logger);
  }
  logger.info("Tests complete.");
  logger.info("ok!");

}

} // namespace scream
