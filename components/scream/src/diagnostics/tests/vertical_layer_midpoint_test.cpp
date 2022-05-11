#include "catch2/catch.hpp"

#include "diagnostics/tests/diagnostic_test_util.hpp"
#include "diagnostics/vertical_layer_midpoint.hpp"
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
  const     int num_mid_packs_p1 = pack_info::num_packs(num_levs+1);

  using Check = ChecksHelpers<ScalarT,num_levs>;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // Create a grids manager - single column for these tests
  const int ncols = 1;
  auto gm = create_gm(comm,ncols,num_levs);

  // Kokkos Policy
  auto policy = ekat::ExeSpaceUtils<ExecSpace>::get_default_team_policy(ncols, num_mid_packs);

  // Input (randomized) views
  view_1d temperature("temperature",num_mid_packs),
          pseudodensity("pseudodensity",num_mid_packs),
          pressure("pressure",num_mid_packs),
          watervapor("watervapor",num_mid_packs);

  auto dview_as_real = [&] (const view_1d& v) -> rview_1d {
    return rview_1d(reinterpret_cast<RealType*>(v.data()),v.size()*pack_size);
  };

  // Construct random input data
  using RPDF = std::uniform_real_distribution<RealType>;
  RPDF pdf_qv(1e-6,1e-3),
       pdf_pseudodens(1.0,100.0),
       pdf_pres(0.0,PC::P0),
       pdf_temp(200.0,400.0);

  ekat::genRandArray(dview_as_real(temperature),     engine, pdf_temp);
  ekat::genRandArray(dview_as_real(pseudodensity),   engine, pdf_pseudodens);
  ekat::genRandArray(dview_as_real(pressure),        engine, pdf_pres);
  ekat::genRandArray(dview_as_real(watervapor),      engine, pdf_qv);

  // A time stamp
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Construct the Diagnostic
  ekat::ParameterList params;
  params.set<std::string>("Diagnostic Name", "Vertical Layer Midpoint");
  params.set<std::string>("Grid", "Point Grid");
  register_diagnostics();
  auto& diag_factory = AtmosphereDiagnosticFactory::instance();
  auto diag = diag_factory.create("VerticalLayerMidpoint",comm,params);
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
  const auto& T_mid_v       = T_mid_f.get_view<ScalarT**>();
  const auto& pseudo_dens_f = input_fields["pseudo_density"];
  const auto& pseudo_dens_v = pseudo_dens_f.get_view<ScalarT**>();
  const auto& p_mid_f       = input_fields["p_mid"];
  const auto& p_mid_v       = p_mid_f.get_view<ScalarT**>();
  const auto& qv_mid_f      = input_fields["qv"];
  const auto& qv_mid_v      = qv_mid_f.get_view<ScalarT**>();

  // The output from the diagnostic should match what would happen if we called "calculate_dz" directly
  {
  for (int icol = 0; icol<ncols;++icol) {
    const auto& T_sub      = ekat::subview(T_mid_v,icol);
    const auto& pseudo_sub = ekat::subview(pseudo_dens_v,icol);
    const auto& p_sub      = ekat::subview(p_mid_v,icol);
    const auto& qv_sub     = ekat::subview(qv_mid_v,icol);
    ekat::genRandArray(dview_as_real(temperature),   engine, pdf_temp);
    ekat::genRandArray(dview_as_real(pseudodensity), engine, pdf_pseudodens);
    ekat::genRandArray(dview_as_real(pressure),      engine, pdf_pres);
    ekat::genRandArray(dview_as_real(watervapor),    engine, pdf_qv);
    Kokkos::deep_copy(T_sub,temperature);
    Kokkos::deep_copy(pseudo_sub,pseudodensity);
    Kokkos::deep_copy(p_sub,pressure);
    Kokkos::deep_copy(qv_sub,watervapor);
  } 
  const auto& diag_out = diag->get_diagnostic(100.0);
  Field z_mid_f = diag_out.clone();
  const auto& z_mid_v = z_mid_f.get_view<ScalarT**>();

  const auto& z_int_v = view_1d("",num_mid_packs_p1);
  const auto& dz_v = view_1d("",num_mid_packs);
  Kokkos::parallel_for("", policy, KOKKOS_LAMBDA(const MemberType& team) {
    const int i = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team,num_mid_packs), [&] (const Int& k) {
      dz_v(k) = PF::calculate_dz(pseudo_dens_v(i,k),p_mid_v(i,k),T_mid_v(i,k),qv_mid_v(i,k));
    });
    team.team_barrier();
    const auto& z_mid_sub = ekat::subview(z_mid_v,i);
    PF::calculate_z_int(team,num_levs,dz_v,0.0,z_int_v);
    PF::calculate_z_mid(team,num_levs,z_int_v,z_mid_sub);
  });
  Kokkos::fence();
  REQUIRE(views_are_equal(diag_out,z_mid_f));
  }
 
  // Finalize the diagnostic
  diag->finalize(); 

} // run()

TEST_CASE("vertical_layer_midpoint_test", "diagnostics"){
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
