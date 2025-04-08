#include "catch2/catch.hpp"

#include "share/grid/mesh_free_grids_manager.hpp"
#include "diagnostics/sea_level_pressure.hpp"
#include "diagnostics/register_diagnostics.hpp"

#include "physics/share/physics_constants.hpp"

#include "share/util/eamxx_setup_random_test.hpp"
#include "share/util/eamxx_common_physics_functions.hpp"
#include "share/field/field_utils.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/util/ekat_test_utils.hpp"

#include <iomanip>

namespace scream {

std::shared_ptr<GridsManager>
create_gm (const ekat::Comm& comm, const int ncols, const int nlevs) {

  const int num_global_cols = ncols*comm.size();

  using vos_t = std::vector<std::string>;
  ekat::ParameterList gm_params;
  gm_params.set("grids_names",vos_t{"Point Grid"});
  auto& pl = gm_params.sublist("Point Grid");
  pl.set<std::string>("type","point_grid");
  pl.set("aliases",vos_t{"Physics"});
  pl.set<int>("number_of_global_columns", num_global_cols);
  pl.set<int>("number_of_vertical_levels", nlevs);

  auto gm = create_mesh_free_grids_manager(comm,gm_params);
  gm->build_grids();

  return gm;
}

//-----------------------------------------------------------------------------------------------//
template<typename DeviceT>
void run(std::mt19937_64& engine)
{
  using PF          = scream::PhysicsFunctions<DeviceT>;
  using PC          = scream::physics::Constants<Real>;
  using Pack        = ekat::Pack<Real,SCREAM_PACK_SIZE>;
  using KT          = ekat::KokkosTypes<DeviceT>;
  using RangePolicy = typename KT::RangePolicy;

  const     int packsize = SCREAM_PACK_SIZE;
  constexpr int num_levs = packsize*2 + 1; // Number of levels to use for tests, make sure the last pack can also have some empty slots (packsize>1).

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // Create a grids manager - single column for these tests
  const int ncols = 1;
  auto gm = create_gm(comm,ncols,num_levs);

  // Construct random input data
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF pdf_pres(0.0,PC::P0),
       pdf_temp(200.0,400.0),
       pdf_surface(100.0,400.0);

  // A time stamp
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Construct the Diagnostic
  ekat::ParameterList params;
  register_diagnostics();
  auto& diag_factory = AtmosphereDiagnosticFactory::instance();
  auto diag = diag_factory.create("SeaLevelPressure",comm,params);
  diag->set_grids(gm);


  // Set the required fields for the diagnostic.
  std::map<std::string,Field> input_fields;
  for (const auto& req : diag->get_required_field_requests()) {
    Field f(req.fid);
    auto & f_ap = f.get_header().get_alloc_properties();
    f_ap.request_allocation(packsize);
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
  {
    // Construct random data to use for test
    // Get views of input data and set to random values
    const auto& T_mid_f       = input_fields["T_mid"];
    const auto& T_mid_v       = T_mid_f.get_view<Pack**>();
    const auto& p_mid_f       = input_fields["p_mid"];
    const auto& p_mid_v       = p_mid_f.get_view<Pack**>();
    const auto& phis_f        = input_fields["phis"];
    const auto& phis_v        = phis_f.get_view<Real*>();
    for (int icol=0;icol<ncols;icol++) {
      const auto& T_sub      = ekat::subview(T_mid_v,icol);
      const auto& p_sub      = ekat::subview(p_mid_v,icol);
      ekat::genRandArray(T_sub, engine, pdf_temp);
      ekat::genRandArray(p_sub, engine, pdf_pres);
    }
    ekat::genRandArray(phis_v, engine, pdf_surface);

    // Run diagnostic and compare with manual calculation
    diag->compute_diagnostic();
    const auto& diag_out = diag->get_diagnostic();
    Field p_sealevel_f = diag_out.clone();
    p_sealevel_f.deep_copy(0);
    const int surf_lev = num_levs - 1;
    const auto& p_sealevel_v = p_sealevel_f.get_view<Real*>();
    auto T_mid_s = ekat::scalarize(T_mid_v);
    auto p_mid_s = ekat::scalarize(p_mid_v);
    auto policy = RangePolicy(0,ncols);
    Kokkos::parallel_for("", policy, KOKKOS_LAMBDA(const int icol) {
      p_sealevel_v(icol) = PF::calculate_psl(T_mid_s(icol,surf_lev),p_mid_s(icol,surf_lev),phis_v(icol));
    });
    Kokkos::fence();
    REQUIRE(views_are_equal(diag_out,p_sealevel_f));
  }
 
  // Finalize the diagnostic
  diag->finalize(); 

} // run()

TEST_CASE("sea_level_pressure_test", "sea_level_pressure_test]"){
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
