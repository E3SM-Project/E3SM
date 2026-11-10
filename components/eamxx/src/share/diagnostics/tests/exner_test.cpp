#include "catch2/catch.hpp"

#include "share/grid/point_grid.hpp"
#include "share/diagnostics/register_diagnostics.hpp"

#include "share/physics/physics_constants.hpp"

#include "share/core/eamxx_setup_random_test.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"
#include "share/field/field_utils.hpp"

#include <ekat_pack.hpp>
#include <ekat_team_policy_utils.hpp>
#include <ekat_view_utils.hpp>

#include <iomanip>

namespace scream {


//-----------------------------------------------------------------------------------------------//
template<typename DeviceT>
void run(std::mt19937_64& engine)
{
  using PF         = scream::PhysicsFunctions<DeviceT>;
  using PC         = scream::physics::Constants<Real>;
  using KT         = ekat::KokkosTypes<DeviceT>;
  using ExecSpace  = typename KT::ExeSpace;
  using MemberType = typename KT::MemberType;
  using TPF        = ekat::TeamPolicyFactory<ExecSpace>;

  const int packsize = SCREAM_PACK_SIZE;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // Create a grids manager - single column for these tests
  const int ncols = 1;
  const int nlevs = packsize*2 + 1; // Number of levels to use for tests, make sure the last pack can also have some empty slots (packsize>1).
  auto grid = create_point_grid("physics",ncols,nlevs,comm);

  // Kokkos Policy
  auto policy = TPF::get_default_team_policy(ncols, nlevs);

  // Construct random input data
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF pdf_pres(0.0,PC::P0.value);

  // A time stamp
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Construct the Diagnostic
  ekat::ParameterList params;
  register_diagnostics();
  auto& diag_factory = DiagnosticFactory::instance();
  auto diag = diag_factory.create("Exner",comm,params, grid);
  // Set the required fields for the diagnostic.
  std::map<std::string,Field> input_fields;
  {
    using namespace ShortFieldTagsNames;
    using namespace ekat::units;
    FieldIdentifier fid("p_mid", grid->get_3d_scalar_layout(LEV), Pa, grid->name());
    Field f(fid);
    f.allocate_view();
    f.get_header().get_tracking().update_time_stamp(t0);
    diag->set_input_field(f);
    input_fields.emplace("p_mid", f);
  }

  // Initialize the diagnostic
  diag->initialize();

  // Run tests
  {
    // Construct random data to use for test
    // Get views of input data and set to random values
    const auto& p_mid_f = input_fields["p_mid"];
    const auto& p_mid_v = p_mid_f.get_view<Real**>();
    for (int icol=0;icol<ncols;icol++) {
      const auto& p_sub = ekat::subview(p_mid_v,icol);
      ekat::genRandArray(p_sub, engine, pdf_pres);
    }

    // Run diagnostic and compare with manual calculation
    diag->compute(t0);
    const auto& diag_out = diag->get();
    Field exner_f = diag_out.clone();
    exner_f.deep_copy(0);
    const auto& exner_v = exner_f.get_view<Real**>();
    Kokkos::parallel_for("", policy, KOKKOS_LAMBDA(const MemberType& team) {
      const int icol = team.league_rank();
      Kokkos::parallel_for(Kokkos::TeamVectorRange(team,nlevs), [&] (const Int& ilev) {
        exner_v(icol,ilev) = PF::exner_function(p_mid_v(icol,ilev));
      });
      team.team_barrier();
    });
    Kokkos::fence();
    REQUIRE(views_are_equal(diag_out,exner_f));
  }
 
} // run()

TEST_CASE("exner_test", "exner_test]"){
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

} // namespace
