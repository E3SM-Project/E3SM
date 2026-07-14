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
  using PF = scream::PhysicsFunctions<DeviceT>;
  using PC = scream::physics::Constants<Real>;
  using KT = ekat::KokkosTypes<DeviceT>;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // For input randomization
  int seed = get_random_test_seed(&comm);

  // Create a grids manager - single column for these tests
  const int ncols = 1;
  const int nlevs = 33;

  auto grid = create_point_grid("physics",ncols,nlevs,comm);

  // A time stamp
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Construct the Diagnostic
  ekat::ParameterList params;
  register_diagnostics();
  auto& diag_factory = DiagnosticFactory::instance();
  auto diag = diag_factory.create("SeaLevelPressure",comm,params, grid);
  // Set the required fields for the diagnostic.
  std::map<std::string,Field> input_fields;
  {
    using namespace ShortFieldTagsNames;
    using namespace ekat::units;
    auto scalar3d = grid->get_3d_scalar_layout(LEV);
    auto scalar2d = grid->get_2d_scalar_layout();
    auto create_and_add = [&](const std::string& fname, const FieldLayout& layout, const ekat::units::Units& u) {
      FieldIdentifier fid(fname, layout, u, grid->name());
      Field f(fid);
      f.allocate_view();
      f.get_header().get_tracking().update_time_stamp(t0);
      diag->set_input_field(f);
      input_fields.emplace(fname, f);
    };
    create_and_add("T_mid", scalar3d, K);
    create_and_add("p_mid", scalar3d, Pa);
    create_and_add("phis",  scalar2d, pow(m,2)/pow(s,2));
  }

  // Initialize the diagnostic
  diag->initialize();

  // Run tests
  {
    // Construct random data to use for test
    // Get views of input data and set to random values
    const auto& T_mid_f       = input_fields["T_mid"];
    const auto& p_mid_f       = input_fields["p_mid"];
    const auto& phis_f        = input_fields["phis"];
    randomize_uniform(T_mid_f,seed++,200,400);
    randomize_uniform(p_mid_f,seed++,0,PC::P0.value);
    randomize_uniform(phis_f,seed++,100,400);
    const auto& T_mid_v       = T_mid_f.get_view<Real**>();
    const auto& p_mid_v       = p_mid_f.get_view<Real**>();
    const auto& phis_v        = phis_f.get_view<Real*>();

    // Run diagnostic and compare with manual calculation
    diag->compute(t0);
    const auto& diag_out = diag->get();
    Field p_sealevel_f = diag_out.clone();
    p_sealevel_f.deep_copy(0);

    const int surf_lev = nlevs - 1;
    const auto& p_sealevel_v = p_sealevel_f.get_view<Real*>();
    auto policy = typename KT::RangePolicy(0,ncols);
    Kokkos::parallel_for("", policy, KOKKOS_LAMBDA(const int icol) {
      p_sealevel_v(icol) = PF::calculate_psl(T_mid_v(icol,surf_lev),p_mid_v(icol,surf_lev),phis_v(icol));
    });
    Kokkos::fence();
    REQUIRE(views_are_equal(diag_out,p_sealevel_f));
  }
}

TEST_CASE("sea_level_pressure_test", "sea_level_pressure_test]"){
  // Run tests for both Real and Pack, and for (potentially) different pack sizes
  using scream::Real;
  using Device = scream::DefaultDevice;

  constexpr int num_runs = 5;

  auto engine = scream::setup_random_test();

  printf(" -> Number of randomized runs: %d\n\n", num_runs);

  for (int irun=0; irun<num_runs; ++irun) {
    printf(" -> run %d\n",irun);
    run<Device>(engine);
  }
} // TEST_CASE

} // namespace
