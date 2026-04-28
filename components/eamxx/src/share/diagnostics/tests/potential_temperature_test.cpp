#include "catch2/catch.hpp"

#include "share/grid/point_grid.hpp"
#include "share/diagnostics/register_diagnostics.hpp"

#include "share/physics/physics_constants.hpp"

#include "share/core/eamxx_setup_random_test.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"
#include "share/field/field_utils.hpp"

namespace scream {

template<typename DeviceT>
void run(std::mt19937_64& engine, int int_ptype)
{
  using PF      = PhysicsFunctions<DeviceT>;
  using PC      = physics::Constants<Real>;
  using KT      = ekat::KokkosTypes<DeviceT>;
  using MDRange = Kokkos::MDRangePolicy<typename KT::ExeSpace,Kokkos::Rank<2>>;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // For randomized inputs
  int seed = get_random_test_seed(&comm);

  // Create a grid
  const int ncols = 1;
  const int nlevs = 33;
  auto grid = create_point_grid("physics",ncols,nlevs,comm);

  // A time stamp
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Construct the Diagnostic
  ekat::ParameterList params;
  register_diagnostics();
  auto& diag_factory = DiagnosticFactory::instance();
  std::string ptype = int_ptype == 0 ? "Tot" : "Liq";
  params.set("temperature_kind", ptype);
  auto diag = diag_factory.create("PotentialTemperature",comm,params, grid);
  // Set the required fields for the diagnostic.
  std::map<std::string,Field> input_fields;
  {
    using namespace ShortFieldTagsNames;
    using namespace ekat::units;
    auto scalar3d = grid->get_3d_scalar_layout(LEV);
    for (const auto& fname : diag->get_input_fields_names()) {
      FieldIdentifier fid(fname, scalar3d, Pa, grid->name());
      Field f(fid);
      f.allocate_view();
      f.get_header().get_tracking().update_time_stamp(t0);
      diag->set_input_field(f);
      input_fields.emplace(fname, f);
    }
  }
  // Initialize the diagnostic
  diag->initialize();

  // Run tests
  {
    // Construct random data to use for test
    // Get views of input data and set to random values
    const auto& T_mid_f = input_fields["T_mid"];
    const auto& p_mid_f = input_fields["p_mid"];
    const auto& q_mid_f = input_fields["qc"];

    randomize_uniform(T_mid_f,seed++,200,400);
    randomize_uniform(p_mid_f,seed++,0,PC::P0.value);
    if (q_mid_f.is_allocated())
      randomize_uniform(q_mid_f,seed++,0,1e-2);

    // Run diagnostic and compare with manual calculation
    diag->compute_diagnostic(t0);
    const auto& diag_out = diag->get_diagnostic();
    Field theta_f = diag_out.clone();
    theta_f.deep_copy(0);
    const auto& theta_v = theta_f.get_view<Real**>();

    auto T_mid_v = T_mid_f.get_view<Real**>();
    auto p_mid_v = p_mid_f.get_view<Real**>();
    auto q_mid_v = q_mid_f.is_allocated() ? q_mid_f.get_view<Real**>() : decltype(T_mid_v){};

    auto manual = KOKKOS_LAMBDA(const int icol, const int ilev) {
      auto theta = PF::calculate_theta_from_T(T_mid_v(icol,ilev),p_mid_v(icol,ilev));
      if (int_ptype==1) {
        theta_v(icol,ilev) = theta - (theta/T_mid_v(icol,ilev)) * (PC::LatVap.value/PC::Cpair.value) * q_mid_v(icol,ilev);
      } else { theta_v(icol,ilev) = theta; }
    };
    MDRange policy({0,0},{ncols,nlevs});
    Kokkos::parallel_for("", policy, manual);
    Kokkos::fence();
    REQUIRE(views_are_equal(diag_out,theta_f));
  }
}

TEST_CASE("potential_temp_test", "potential_temp_test]"){
  using Device = scream::DefaultDevice;

  constexpr int num_runs = 5;

  auto engine = scream::setup_random_test();

  printf(" -> Number of randomized runs: %d\n\n", num_runs);

  printf(" -> Testing Pack<Real,%d> scalar type...",SCREAM_PACK_SIZE);
  for (int irun=0; irun<num_runs; ++irun) {
    for (const int int_ptype : {0,1}) {
    run<Device>(engine, int_ptype);
  }
  }
  printf("ok!\n");

  printf("\n");

} // TEST_CASE

} // namespace scream
