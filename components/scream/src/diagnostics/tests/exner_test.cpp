#include "catch2/catch.hpp"

#include "diagnostics/tests/diagnostic_test_util.hpp"
#include "diagnostics/exner.hpp"
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
  using view_1d    = typename KT::template view_1d<ScalarT>;
  using rview_1d   = typename KT::template view_1d<RealType>;


  constexpr int pack_size = sizeof(ScalarT) / sizeof(RealType);
  using pack_info = ekat::PackInfo<pack_size>;

  constexpr int num_levs = 32; // Number of levels to use for tests.
  const     int num_mid_packs = pack_info::num_packs(num_levs);

  using Check = ChecksHelpers<ScalarT,num_levs>;


  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // Create a grids manager - single column for these tests
  const int ncols = 1;
  auto gm = create_gm(comm,ncols,num_levs);

  // Kokkos Policy
  auto policy = ekat::ExeSpaceUtils<ExecSpace>::get_default_team_policy(ncols, num_mid_packs);

  // Input (randomized) views
  view_1d pressure("pressure",num_mid_packs);

  auto dview_as_real = [&] (const view_1d& v) -> rview_1d {
    return rview_1d(reinterpret_cast<RealType*>(v.data()),v.size()*pack_size);
  };

  // Construct random input data
  using RPDF = std::uniform_real_distribution<RealType>;
  RPDF pdf_pres(0.0,PC::P0);

  //contruct random integers
  using IPDF = std::uniform_int_distribution<int>;
  IPDF pdf_rand_int(1,100);

  ekat::genRandArray(dview_as_real(pressure),  engine,pdf_pres);

  // A time stamp
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Construct the Diagnostic
  ekat::ParameterList params;
  params.set<std::string>("Diagnostic Name", "Exner");
  params.set<std::string>("Grid", "Point Grid");
  register_diagnostics();
  auto& diag_factory = AtmosphereDiagnosticFactory::instance();
  auto diag = diag_factory.create("Exner",comm,params);
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
  const ScalarT p0   = PC::P0;
  const ScalarT zero = 0.0;

  // Get views of input data and set to random values
  const auto& p_mid_f = input_fields["p_mid"];
  const auto& p_mid_v = p_mid_f.get_view<ScalarT**>();

  // Test 1 - property tests 
  //  - exner(pmid=0) = 0
  {
  Field zero_f = p_mid_f.clone();  // Field with only zeros
  zero_f.deep_copy(0.0);
  for (int icol = 0; icol<ncols;++icol) {
    const auto& p_sub = ekat::subview(p_mid_v,icol);
    Kokkos::deep_copy(p_sub,zero);
  }
  const auto& diag_out = diag->get_diagnostic(100.0);
  REQUIRE(views_are_equal(diag_out,zero_f));
  }
  //  - exner=1 when p=p0
  {
  Field one_f = p_mid_f.clone();  // Field with only ones
  one_f.deep_copy(1.0);
  for (int icol = 0; icol<ncols;++icol) {
    const auto& p_sub = ekat::subview(p_mid_v,icol);
    Kokkos::deep_copy(p_sub,p0);
  } 
  const auto& diag_out = diag->get_diagnostic(100.0);
  REQUIRE(views_are_equal(diag_out,one_f));
  }
  // The output from the diagnostic should match what would happen if we called "exner" directly
  {
  for (int icol = 0; icol<ncols;++icol) {
    const auto& p_sub = ekat::subview(p_mid_v,icol);
    ekat::genRandArray(dview_as_real(pressure), engine, pdf_pres);
    Kokkos::deep_copy(p_sub,pressure);
  } 
  Field exner_f = p_mid_f;
  exner_f.deep_copy<double,Host>(0.0);
  const auto& exner_v = exner_f.get_view<ScalarT**>();
  Kokkos::parallel_for("", policy, KOKKOS_LAMBDA(const MemberType& team) {
    const int i = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team,num_mid_packs), [&] (const Int& k) {
      exner_v(i,k) = PF::exner_function(p_mid_v(i,k));
    });
    team.team_barrier();
  });
  Kokkos::fence();
  const auto& diag_out = diag->get_diagnostic(100.0);
  REQUIRE(views_are_equal(diag_out,exner_f));
  }
 
  // Finalize the diagnostic
  diag->finalize(); 

} // run()

TEST_CASE("exner_diagnostic_test", "diagnostics"){
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
