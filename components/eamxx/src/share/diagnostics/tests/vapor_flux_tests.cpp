#include "catch2/catch.hpp"

#include "share/grid/point_grid.hpp"
#include "share/diagnostics/register_diagnostics.hpp"

#include "share/physics/physics_constants.hpp"

#include "share/core/eamxx_setup_random_test.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"
#include "share/field/field_utils.hpp"

#include <ekat_team_policy_utils.hpp>

namespace scream {

template<typename DeviceT>
void run(std::mt19937_64& engine)
{
  using PC         = scream::physics::Constants<Real>;
  using KT         = ekat::KokkosTypes<DeviceT>;
  using ExecSpace  = typename KT::ExeSpace;
  using TPF        = ekat::TeamPolicyFactory<ExecSpace>;
  using MemberType = typename KT::MemberType;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // For input randomization
  int seed = get_random_test_seed(&comm);

  // Create a grids manager - single column for these tests
  const int ncols = 1; //TODO should be set to the size of the communication group.
  const int nlevs = 33;
  auto grid = create_point_grid("physics",ncols,nlevs,comm);

  // Kokkos Policy
  auto policy = TPF::get_default_team_policy(ncols, nlevs);

  // A time stamp
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  ekat::ParameterList params;
  register_diagnostics();
  auto& diag_factory = DiagnosticFactory::instance();

  REQUIRE_THROWS (diag_factory.create("VaporFlux",comm,params, grid)); // No 'wind_component'
  params.set<std::string>("wind_component","foo");
  REQUIRE_THROWS (diag_factory.create("VaporFlux",comm,params, grid)); // Invalid 'wind_component'
  for (const std::string which_comp : {"Zonal", "Meridional"}) {
    // Construct the Diagnostic
    params.set<std::string>("wind_component",which_comp);
    auto diag = diag_factory.create("VaporFlux",comm,params, grid);
    // Set the required fields for the diagnostic.
    std::map<std::string,Field> input_fields;
    {
      using namespace ShortFieldTagsNames;
      using namespace ekat::units;
      auto scalar3d = grid->get_3d_scalar_layout(LEV);
      auto vec3d = grid->get_3d_vector_layout(LEV, 2, "cmp");
      for (const auto& fname : diag->get_input_fields_names()) {
        auto layout = (fname == "horiz_winds") ? vec3d : scalar3d;
        FieldIdentifier fid(fname, layout, Pa, grid->name());
        Field f(fid);
        f.allocate_view();
        f.get_header().get_tracking().update_time_stamp(t0);
        diag->set_input_field(f);
        input_fields[fname] = f;
      }
    }

    // Initialize the diagnostic
    diag->initialize();

    // Run tests
    {
      // Construct random data to use for test
      // Get views of input data and set to random values
      const auto& qv_f = input_fields["qv"];
      const auto& dp_f = input_fields["pseudo_density"];
      const auto& uv_f = input_fields["horiz_winds"];

      const auto& qv_v = qv_f.get_view<Real**>();
      const auto& dp_v = dp_f.get_view<Real**>();
      const auto& uv_v = uv_f.get_view<Real***>();

      randomize_uniform (qv_f,seed++,0,1e-3);
      randomize_uniform (dp_f,seed++,1,100);
      randomize_uniform (uv_f,seed++,0,200);

      // Run diagnostic and compare with manual calculation
      diag->compute_diagnostic(t0);
      const auto& diag_out = diag->get_diagnostic();
      Field qv_vert_integrated_flux_u_f = diag_out.clone();
      qv_vert_integrated_flux_u_f.deep_copy(0);
      const auto& qv_vert_integrated_flux_u_v = qv_vert_integrated_flux_u_f.get_view<Real*>();
      constexpr Real g = PC::gravit.value;
      int comp = which_comp=="Zonal" ? 0 : 1;

      Kokkos::parallel_for("", policy, KOKKOS_LAMBDA(const MemberType& team) {
        const int icol = team.league_rank();
        auto wind = ekat::subview(uv_v,icol,comp);
        Kokkos::parallel_reduce(Kokkos::TeamVectorRange(team, nlevs),
                                [&] (const int& ilev, Real& lsum) {
          lsum += wind(ilev) * qv_v(icol,ilev) * dp_v(icol,ilev) / g;
        },qv_vert_integrated_flux_u_v(icol));
      });
      Kokkos::fence();
      REQUIRE(views_are_equal(diag_out,qv_vert_integrated_flux_u_f));
    }
  }
}

TEST_CASE("zonal_vapor_flux_test", "zonal_vapor_flux_test]"){
  using Device = scream::DefaultDevice;

  constexpr int num_runs = 5;

  auto engine = scream::setup_random_test();

  printf(" -> Number of randomized runs: %d\n\n", num_runs);

  for (int irun=0; irun<num_runs; ++irun) {
    run<Device>(engine);
  }
}

} // namespace
