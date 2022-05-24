#include "catch2/catch.hpp"

#include "share/grid/mesh_free_grids_manager.hpp"
#include "diagnostics/potential_temperature.hpp"
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

std::shared_ptr<GridsManager>
create_gm (const ekat::Comm& comm, const int ncols, const int nlevs) {

  const int num_local_elems = 4;
  const int np = 4;
  const int num_global_cols = ncols*comm.size();

  ekat::ParameterList gm_params;
  gm_params.set<std::string>("Reference Grid", "Point Grid");
  gm_params.set<int>("Number of Global Columns", num_global_cols);
  gm_params.set<int>("Number of Local Elements", num_local_elems);
  gm_params.set<int>("Number of Vertical Levels", nlevs);
  gm_params.set<int>("Number of Gauss Points", np);

  auto gm = create_mesh_free_grids_manager(comm,gm_params);
  gm->build_all_grids();

  return gm;
}

template<typename ScalarT, int NumLevels>
struct ChecksHelpers {

  static bool is_non_negative (const ScalarT& s, const int k) {
    return not ( k<NumLevels && (s<0 || std::isnan(s)) );
  }
  static bool equal (const ScalarT& lhs, const ScalarT& rhs) {
    return lhs==rhs;
  }
  static bool approx_equal (const ScalarT lhs, const ScalarT rhs,
                            const int k, const ScalarT tol) {
    using std::abs;
    return not ( k<NumLevels && abs(lhs-rhs)>=tol );
  }
  static bool approx_equal (const ScalarT computed, const ScalarT expected, const ScalarT tol) {
    using std::abs;
    return abs(computed-expected)/abs(expected) < tol;
  }
};

template<typename T, int N, int NumLevels>
struct ChecksHelpers<ekat::Pack<T,N>,NumLevels> {
  using ScalarT = ekat::Pack<T,N>;

  static bool is_non_negative (const ScalarT& s, const int k) {
    const auto range = ekat::range<ScalarT>(k*N);
    const auto range_mask = range < NumLevels;
    return ( range_mask && (s<0 || isnan(s) ) ).none();
  }
  static bool equal (const ScalarT& lhs, const ScalarT& rhs) {
    return (lhs==rhs).all();
  }
  static bool approx_equal (const ScalarT& lhs, const ScalarT& rhs,
                            const int k, const T tol) {
    const auto range = ekat::range<ScalarT>(k*N);
    const auto range_mask = range < NumLevels;
    return (range_mask && abs(lhs-rhs)>=tol).none();
  }
  static bool approx_equal (const ScalarT& computed, const ScalarT& expected, const T tol) {
    return (abs(computed-expected)/abs(expected) < tol).all();
  }
};

// Helper function. Create Mirror View and Deep-Copy (CMVDC)
template<typename ViewT>
auto cmvdc (const ViewT& v_d) -> typename ViewT::HostMirror {
  auto v_h = Kokkos::create_mirror_view(v_d);
  Kokkos::deep_copy(v_h,v_d);
  return v_h;
}

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
  using MemberType = typename KT::MemberType;
  using view_1d    = typename KT::template view_1d<ScalarT>;
  using rview_1d   = typename KT::template view_1d<RealType>;


  constexpr int pack_size = sizeof(ScalarT) / sizeof(RealType);
  using pack_info = ekat::PackInfo<pack_size>;

  constexpr int num_levs = 32; // Number of levels to use for tests.
  const     int num_mid_packs = pack_info::num_packs(num_levs);

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // Create a grids manager - single column for these tests
  const int ncols = 1;
  auto gm = create_gm(comm,ncols,num_levs);

  // Kokkos Policy
  auto policy = ekat::ExeSpaceUtils<ExecSpace>::get_default_team_policy(ncols, num_mid_packs);

  // Input (randomized) views
  view_1d temperature("temperature",num_mid_packs),
          pressure("pressure",num_mid_packs);

  auto dview_as_real = [&] (const view_1d& v) -> rview_1d {
    return rview_1d(reinterpret_cast<RealType*>(v.data()),v.size()*pack_size);
  };

  // Construct random input data
  using RPDF = std::uniform_real_distribution<RealType>;
  RPDF pdf_pres(0.0,PC::P0),
       pdf_temp(200.0,400.0);

  //contruct random integers
  using IPDF = std::uniform_int_distribution<int>;
  IPDF pdf_rand_int(1,100);

  ekat::genRandArray(dview_as_real(temperature),     engine,pdf_temp);
  ekat::genRandArray(dview_as_real(pressure),        engine,pdf_pres);

  // A time stamp
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Construct the Diagnostic
  ekat::ParameterList params;
  params.set<std::string>("Diagnostic Name", "Potential Temperature");
  params.set<std::string>("Grid", "Point Grid");
  register_diagnostics();
  auto& diag_factory = AtmosphereDiagnosticFactory::instance();
  auto diag = diag_factory.create("PotentialTemperature",comm,params);
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
  const auto& T_mid_f = input_fields["T_mid"];
  const auto& T_mid_v = T_mid_f.get_view<ScalarT**>();
  const auto& p_mid_f = input_fields["p_mid"];
  const auto& p_mid_v = p_mid_f.get_view<ScalarT**>();

  // Test 1 - property tests 
  //  - theta(T=0) = 0
  {
    Field zero_f = T_mid_f;  // Field with only zeros
    zero_f.deep_copy(0.0);
    for (int icol = 0; icol<ncols;++icol) {
      const auto& T_sub = ekat::subview(T_mid_v,icol);
      const auto& p_sub = ekat::subview(p_mid_v,icol);
      Kokkos::deep_copy(T_sub,zero);
      Kokkos::deep_copy(p_sub,p0);
    }
    diag->run();
    const auto& diag_out = diag->get_diagnostic();
    REQUIRE(views_are_equal(diag_out,zero_f));
  }
  //  - theta=T when p=p0
  {
    for (int icol = 0; icol<ncols;++icol) {
      const auto& T_sub = ekat::subview(T_mid_v,icol);
      const auto& p_sub = ekat::subview(p_mid_v,icol);
      Kokkos::deep_copy(T_sub,temperature);
      Kokkos::deep_copy(p_sub,p0);
    } 
    diag->run();
    const auto& diag_out = diag->get_diagnostic();
    REQUIRE(views_are_equal(diag_out,T_mid_f));
  }
  // The output from the diagnostic should match what would happen if we called "calculate_theta_from_T" directly
  {
    Field theta_f = T_mid_f;
    theta_f.deep_copy<double,Host>(0.0);
    const auto& theta_v = theta_f.get_view<ScalarT**>();
    Kokkos::parallel_for("", policy, KOKKOS_LAMBDA(const MemberType& team) {
      const int i = team.league_rank();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team,num_mid_packs), [&] (const Int& k) {
        theta_v(i,k) = PF::calculate_theta_from_T(T_mid_v(i,k),p_mid_v(i,k));
      });
      team.team_barrier();
    });
    Kokkos::fence();
    diag->run();
    const auto& diag_out = diag->get_diagnostic();
    REQUIRE(views_are_equal(diag_out,theta_f));
  }
 
  // Finalize the diagnostic
  diag->finalize(); 

} // run()

TEST_CASE("potential_temp_test", "potential_temp_test]"){
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
