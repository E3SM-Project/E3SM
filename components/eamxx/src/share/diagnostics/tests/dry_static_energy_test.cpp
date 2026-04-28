#include "catch2/catch.hpp"

#include "share/grid/point_grid.hpp"
#include "share/diagnostics/register_diagnostics.hpp"

#include "share/physics/physics_constants.hpp"

#include "share/core/eamxx_setup_random_test.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"
#include "share/field/field_utils.hpp"

namespace scream {

template<typename DeviceT>
void run()
{
  using PF      = PhysicsFunctions<DeviceT>;
  using KT      = KokkosTypes<DeviceT>;
  using MDRange = Kokkos::MDRangePolicy<typename KT::ExeSpace,Kokkos::Rank<2>>;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // Create a grids manager - single column for these tests
  const int ncols = 1;
  const int nlevs = 33;
  auto grid = create_point_grid("physics",ncols,nlevs,comm);

  // For input randomization
  int seed = get_random_test_seed(&comm);

  // A time stamp
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Construct the Diagnostic
  ekat::ParameterList params;
  register_diagnostics();
  auto& diag_factory = DiagnosticFactory::instance();
  auto diag = diag_factory.create("DryStaticEnergy",comm,params, grid);
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
    create_and_add("height_mid", scalar3d, m);
    create_and_add("phis",  scalar2d, pow(m,2)/pow(s,2));
  }

  // Initialize the diagnostic
  diag->initialize();

  // Run tests
  {
    // Randomize inputs
    auto& T_mid = input_fields["T_mid"];
    auto& h_mid = input_fields["height_mid"];
    auto& phis  = input_fields["phis"];
    randomize_uniform(T_mid,seed++,200,400);
    randomize_uniform(h_mid,seed++,0,10000);
    randomize_uniform(phis, seed++,0,400);

    // Run diagnostic and compare with manual calculation
    diag->compute_diagnostic(t0);
    const auto& diag_out = diag->get_diagnostic();
    Field dse_f = diag_out.clone();
    dse_f.deep_copy(0);
    const auto& dse_v = dse_f.get_view<Real**>();

    MDRange policy({0,0},{ncols,nlevs});
    auto T_mid_v = T_mid.get_view<const Real**>();
    auto h_mid_v = h_mid.get_view<const Real**>();
    auto phis_v  = phis.get_view<const Real*>();
    auto manual = KOKKOS_LAMBDA(const int icol, const int ilev) {
      dse_v(icol,ilev) = PF::calculate_dse(T_mid_v(icol,ilev),h_mid_v(icol,ilev),phis_v(icol));
    };
    Kokkos::parallel_for("", policy, manual);
    Kokkos::fence();
    REQUIRE(views_are_equal(diag_out,dse_f));
  }
}

TEST_CASE("dry_static_energy_test", "dry_static_energy_test]"){
  using Device = scream::DefaultDevice;

  constexpr int num_runs = 5;

  printf(" -> Number of randomized runs: %d\n\n", num_runs);

  printf(" -> Testing Pack<Real,%d> scalar type...",SCREAM_PACK_SIZE);
  for (int irun=0; irun<num_runs; ++irun) {
    run<Device>();
  }
} // TEST_CASE

} // namespace scream
