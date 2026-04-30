#include "catch2/catch.hpp"

#include "share/grid/point_grid.hpp"
#include "share/diagnostics/register_diagnostics.hpp"

#include "share/physics/physics_constants.hpp"

#include "share/core/eamxx_setup_random_test.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"
#include "share/field/field_utils.hpp"

namespace scream {

template<typename DeviceT>
void run(std::mt19937_64& engine)
{
  using KT = ekat::KokkosTypes<DeviceT>;
  using RP = typename KT::RangePolicy;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // For input randomization
  int seed = get_random_test_seed(&comm);

  // Create a grid
  const int ncols = 1; //TODO should be set to the size of the communication group.
  const int nlevs = 33;
  auto grid = create_point_grid("physics",ncols,nlevs,comm);

  // A time stamp
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Construct the Diagnostic
  ekat::ParameterList params;
  register_diagnostics();
  auto& diag_factory = DiagnosticFactory::instance();
  auto diag = diag_factory.create("ShortwaveCloudForcing",comm,params, grid);
  // Set the required fields for the diagnostic.
  std::map<std::string,Field> input_fields;
  {
    using namespace ShortFieldTagsNames;
    using namespace ekat::units;
    auto m2 = (m*m).rename("m2");
    auto ilev_layout = grid->get_3d_scalar_layout(ILEV);
    for (const auto& fname : diag->get_input_fields_names()) {
      FieldIdentifier fid(fname, ilev_layout, W/m2, grid->name());
      Field f(fid,true);
      f.get_header().get_tracking().update_time_stamp(t0);
      diag->set_input_field(f);
      input_fields.emplace(fname, f);
    }
  }

  // Initialize the diagnostic
  diag->initialize();

  // Run tests
  {
    // Randomize inputs
    const auto& SW_flux_dn_f        = input_fields["SW_flux_dn"];
    const auto& SW_flux_up_f        = input_fields["SW_flux_up"];
    const auto& SW_clrsky_flux_dn_f = input_fields["SW_clrsky_flux_dn"];
    const auto& SW_clrsky_flux_up_f = input_fields["SW_clrsky_flux_up"];
    randomize_uniform(SW_flux_dn_f,seed++,0,400);
    randomize_uniform(SW_clrsky_flux_dn_f,seed++,0,400);

    // Run diagnostic and compare with manual calculation
    diag->compute(t0);
    const auto& diag_out = diag->get();
    auto SWCF_f = diag_out.clone();
    SWCF_f.deep_copy(0);

    auto SWCF_v = SWCF_f.get_view<Real*>();
    auto SW_flux_dn_v = SW_flux_dn_f.get_view<const Real**>();
    auto SW_flux_up_v = SW_flux_up_f.get_view<const Real**>();
    auto SW_clrsky_flux_dn_v = SW_clrsky_flux_dn_f.get_view<const Real**>();
    auto SW_clrsky_flux_up_v = SW_clrsky_flux_up_f.get_view<const Real**>();
    auto manual = KOKKOS_LAMBDA(const int icol) {
      SWCF_v(icol) = (SW_flux_dn_v(icol,0) - SW_flux_up_v(icol,0)) - (SW_clrsky_flux_dn_v(icol,0) - SW_clrsky_flux_up_v(icol,0));
    };
    RP policy(0,ncols);
    Kokkos::parallel_for("", policy, manual);
    Kokkos::fence();
    REQUIRE(views_are_equal(diag_out,SWCF_f));
  }
}

TEST_CASE("shortwave_cloud_forcing_test", "shortwave_cloud_forcing_test]"){
  // Run tests for both Real and Pack, and for (potentially) different pack sizes
  using scream::Real;
  using Device = scream::DefaultDevice;

  constexpr int num_runs = 5;

  auto engine = scream::setup_random_test();

  printf(" -> Number of randomized runs: %d\n\n", num_runs);

  printf(" -> Testing Pack<Real,%d> scalar type...",SCREAM_PACK_SIZE);
  for (int irun=0; irun<num_runs; ++irun) {
    run<Device>(engine);
  }
  printf("ok!\n");

  printf("\n");

} // TEST_CASE

} // namespace scream
