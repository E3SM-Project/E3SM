#include "catch2/catch.hpp"

#include "share/grid/point_grid.hpp"
#include "share/diagnostics/register_diagnostics.hpp"
#include "share/physics/physics_constants.hpp"
#include "share/core/eamxx_setup_random_test.hpp"
#include "share/field/field_utils.hpp"

namespace scream {

template<typename DeviceT>
void run(std::mt19937_64& engine)
{
  using PC = physics::Constants<Real>;
  using KT = KokkosTypes<DeviceT>;
  using RP = typename KT::RangePolicy;

  ekat::Comm comm(MPI_COMM_WORLD);

  // Create a grids manager
  const int ncols = 13;
  auto grid = create_point_grid("physics",ncols*comm.size(),1,comm);

  // For randomized inputs
  int seed = get_random_test_seed(&comm);

  // Initial time stamp
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Construct the Diagnostics
  register_diagnostics();
  ekat::ParameterList params;

  auto& diag_factory = DiagnosticFactory::instance();
  auto diag_latent_heat = diag_factory.create("surface_upward_latent_heat_flux", comm, params, grid);
  // Set the required fields for the diagnostic.
  std::map<std::string,Field> input_fields;
  {
    using namespace ekat::units;
    auto scalar2d = grid->get_2d_scalar_layout();
    FieldIdentifier fid("surf_evap", scalar2d, kg/(m*m*s), grid->name());
    Field f(fid);
    f.allocate_view();
    f.get_header().get_tracking().update_time_stamp(t0);
    diag_latent_heat->set_input_field(f);
    input_fields.emplace("surf_evap", f);
  }

  // Initialize the diagnostic
  diag_latent_heat->initialize();
  const auto& surf_evap_f = input_fields["surf_evap"];
  randomize_uniform(surf_evap_f,seed++,0,1);

  // Run diagnostic and compare with manual calculation
  diag_latent_heat->compute(t0);
  const auto diag_latent_heat_out = diag_latent_heat->get();
  Field surf_lhf = diag_latent_heat_out.clone();
  surf_lhf.deep_copy(0);
  const auto& surf_lhf_v = surf_lhf.get_view<Real*>();
  constexpr Real latent_heat_evap = PC::LatVap.value; // [J/kg]
  const auto& surf_evap_v = surf_evap_f.get_view<Real*>();
  auto manual = KOKKOS_LAMBDA (const int& icol) {
    surf_lhf_v(icol) = surf_evap_v(icol) * latent_heat_evap;
  };
  Kokkos::parallel_for(RP(0, ncols),manual);
  Kokkos::fence();

  REQUIRE(views_are_equal(diag_latent_heat_out, surf_lhf));
}

TEST_CASE("surf_upward_latent_heat_flux_test", "[surf_upward_latent_heat_flux_test]") {
  using Device = scream::DefaultDevice;

  constexpr int num_runs = 5;
  auto engine = scream::setup_random_test();
  printf(" -> Number of randomized runs: %d\n", num_runs);
  for (int irun=0; irun<num_runs; ++irun) {
    run<Device>(engine);
  }
}

} // namespace scream
