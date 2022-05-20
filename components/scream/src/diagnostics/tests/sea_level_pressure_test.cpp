#include "catch2/catch.hpp"

#include "diagnostics/tests/diagnostic_test_util.hpp"
#include "diagnostics/sea_level_pressure.hpp"
#include "diagnostics/register_diagnostics.hpp"

#include "physics/share/physics_constants.hpp"

#include "share/util/scream_setup_random_test.hpp"
#include "share/util/scream_common_physics_functions.hpp"
#include "share/field/field_utils.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/util/ekat_test_utils.hpp"

#include <iomanip>

namespace scream {

//-----------------------------------------------------------------------------------------------//
template<typename ScalarT, typename DeviceT>
void run(std::mt19937_64& engine)
{
  using STraits    = ekat::ScalarTraits<ScalarT>;
  using RealType   = typename STraits::scalar_type;

  using PF         = scream::PhysicsFunctions<DeviceT>;
  using PC         = scream::physics::Constants<RealType>;

  using KT         = ekat::KokkosTypes<DeviceT>;
  using ExecSpace  = typename KT::ExeSpace;
  using TeamPolicy = typename KT::TeamPolicy;
  using MemberType = typename KT::MemberType;
  using rview_1d   = typename KT::template view_1d<RealType>;

  constexpr int num_levs = 32; // Number of levels to use for tests.

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // Create a grids manager - single column for these tests
  const int ncols = 128;
  auto gm = create_gm(comm,ncols,num_levs);

  // Input (randomized) views
  rview_1d temperature("temperature",num_levs),
           pressure("pressure",num_levs);

  // Construct random input data
  using RPDF = std::uniform_real_distribution<RealType>;
  RPDF pdf_pres(0.0,PC::P0),
       pdf_temp(200.0,400.0),
       pdf_surface(100.0,400.0);

  // A time stamp
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Construct the Diagnostic
  ekat::ParameterList params;
  params.set<std::string>("Diagnostic Name", "Sea Level Pressure");
  params.set<std::string>("Grid", "Point Grid");
  register_diagnostics();
  auto& diag_factory = AtmosphereDiagnosticFactory::instance();
  auto diag = diag_factory.create("SeaLevelPressure",comm,params);
  diag->set_grids(gm);


  // Set the required fields for the diagnostic.
  std::map<std::string,Field> input_fields;
  for (const auto& req : diag->get_required_field_requests()) {
    Field f(req.fid);
    f.allocate_view();
    const auto name = f.name();
    f.get_header().get_tracking().update_time_stamp(t0);
    diag->set_required_field(f.get_const());
    REQUIRE_THROWS(diag->set_computed_field(f));
    input_fields.emplace(name,f);
  }

  // Initialize the diagnostic
  diag->initialize(t0,RunType::Initial);

  // Run tests
  // Get views of input data and set to random values
  const auto& T_mid_f       = input_fields["T_mid"];
  const auto& T_mid_v       = T_mid_f.get_view<Real**>();
  const auto& p_mid_f       = input_fields["p_mid"];
  const auto& p_mid_v       = p_mid_f.get_view<Real**>();
  const auto& phis_f        = input_fields["phis"];
  const auto& phis_v        = phis_f.get_view<Real*>();

  // The output from the diagnostic should match what would happen if we called "calculate_psl" directly
  {
  for (int icol = 0; icol<ncols;++icol) {
    const auto& T_sub      = ekat::subview(T_mid_v,icol);
    const auto& p_sub      = ekat::subview(p_mid_v,icol);
    ekat::genRandArray(temperature,   engine, pdf_temp);
    ekat::genRandArray(pressure,      engine, pdf_pres);
    Kokkos::deep_copy(T_sub,temperature);
    Kokkos::deep_copy(p_sub,pressure);
  } 
  ekat::genRandArray(phis_v, engine, pdf_surface);

  const auto& diag_out = diag->get_diagnostic(100.0);
  Field psl_f = phis_f.clone(); //diag_out.clone();
  const auto& psl_v = psl_f.get_view<Real*>();

  Kokkos::parallel_for("", ncols, KOKKOS_LAMBDA(const int i) {
    psl_v(i) = PF::calculate_psl(T_mid_v(i,num_levs-1),p_mid_v(i,num_levs-1),phis_v(i));
  });
  Kokkos::fence();
  REQUIRE(views_are_equal(diag_out,psl_f));
  }
 
  // Finalize the diagnostic
  diag->finalize(); 

} // run()

TEST_CASE("sea_level_pressure_test", "diagnostics"){
  // Run tests for both Real and Pack, and for (potentially) different pack sizes
  using scream::Real;
  using Device = scream::DefaultDevice;

  constexpr int num_runs = 5;

  auto engine = scream::setup_random_test();

  printf(" -> Number of randomized runs: %d\n\n", num_runs);

  printf(" -> Testing Real scalar type...");
  for (int irun=0; irun<num_runs; ++irun) {
    run<Real,Device>(engine);
  }
  printf("ok!\n");

  printf(" -> Testing Pack<Real,%d> scalar type...",SCREAM_SMALL_PACK_SIZE);
  for (int irun=0; irun<num_runs; ++irun) {
    run<ekat::Pack<Real,SCREAM_SMALL_PACK_SIZE>,Device>(engine);
  }
  printf("ok!\n");

  if (SCREAM_PACK_SIZE!=SCREAM_SMALL_PACK_SIZE) {
    printf(" -> Testing Pack<Real,%d> scalar type...",SCREAM_PACK_SIZE);
    for (int irun=0; irun<num_runs; ++irun) {
      run<ekat::Pack<Real,SCREAM_PACK_SIZE>,Device>(engine);
    }
    printf("ok!\n");
  }

  printf("\n");

} // TEST_CASE

} // namespace
