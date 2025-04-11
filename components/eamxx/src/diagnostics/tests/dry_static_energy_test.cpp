#include "catch2/catch.hpp"

#include "share/grid/mesh_free_grids_manager.hpp"
#include "diagnostics/dry_static_energy.hpp"
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
  using PF         = scream::PhysicsFunctions<DeviceT>;
  using PC         = scream::physics::Constants<Real>;
  using KT         = ekat::KokkosTypes<DeviceT>;
  using ExecSpace  = typename KT::ExeSpace;
  using MemberType = typename KT::MemberType;
  using view_1d    = typename KT::template view_1d<Real>;

  const     int packsize = SCREAM_PACK_SIZE;
  constexpr int num_levs = packsize*2 + 1; // Number of levels to use for tests, make sure the last pack can also have some empty slots (packsize>1).

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // Create a grids manager - single column for these tests
  const int ncols = 1;
  auto gm = create_gm(comm,ncols,num_levs);

  // Kokkos Policy
  auto policy = ekat::ExeSpaceUtils<ExecSpace>::get_default_team_policy(ncols, num_levs);

  // Construct random input data
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF pdf_qv(1e-6,1e-3),
       pdf_pseudodens(1.0,100.0),
       pdf_pres(0.0,PC::P0),
       pdf_temp(200.0,400.0),
       pdf_surface(100.0,400.0);

  // A time stamp
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Construct the Diagnostic
  ekat::ParameterList params;
  register_diagnostics();
  auto& diag_factory = AtmosphereDiagnosticFactory::instance();
  auto diag = diag_factory.create("DryStaticEnergy",comm,params);
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
    const auto& T_mid_v       = T_mid_f.get_view<Real**>();
    const auto& pseudo_dens_f = input_fields["pseudo_density"];
    const auto& pseudo_dens_v = pseudo_dens_f.get_view<Real**>();
    const auto& p_mid_f       = input_fields["p_mid"];
    const auto& p_mid_v       = p_mid_f.get_view<Real**>();
    const auto& qv_mid_f      = input_fields["qv"];
    const auto& qv_mid_v      = qv_mid_f.get_view<Real**>();
    const auto& phis_f        = input_fields["phis"];
    const auto& phis_v        = phis_f.get_view<Real*>();
    for (int icol=0;icol<ncols;icol++) {
      const auto& T_sub      = ekat::subview(T_mid_v,icol);
      const auto& pseudo_sub = ekat::subview(pseudo_dens_v,icol);
      const auto& p_sub      = ekat::subview(p_mid_v,icol);
      const auto& qv_sub     = ekat::subview(qv_mid_v,icol);
      ekat::genRandArray(T_sub,      engine, pdf_temp);
      ekat::genRandArray(pseudo_sub, engine, pdf_pseudodens);
      ekat::genRandArray(p_sub,      engine, pdf_pres);
      ekat::genRandArray(qv_sub,     engine, pdf_qv);
    }
    ekat::genRandArray(phis_v, engine, pdf_surface);

    // Run diagnostic and compare with manual calculation
    diag->compute_diagnostic();
    const auto& diag_out = diag->get_diagnostic();
    Field dse_f = diag_out.clone();
    dse_f.deep_copy(0);
    const auto& dse_v = dse_f.get_view<Real**>();
    // We need some temporary views
    const auto& z_int_v = view_1d("",num_levs+1);
    const auto& dz_v    = view_1d("",num_levs);
    const auto& z_mid_v = view_1d("",num_levs);
    Kokkos::parallel_for("", policy, KOKKOS_LAMBDA(const MemberType& team) {
      const int icol = team.league_rank();
      Kokkos::parallel_for(Kokkos::TeamVectorRange(team,num_levs), [&] (const Int& ilev) {
        dz_v(ilev) = PF::calculate_dz(pseudo_dens_v(icol,ilev),p_mid_v(icol,ilev),T_mid_v(icol,ilev),qv_mid_v(icol,ilev));
      });
      team.team_barrier();
      const auto& dse_sub = ekat::subview(dse_v,icol);
      PF::calculate_z_int(team,num_levs,dz_v,0.0,z_int_v);
      team.team_barrier();
      PF::calculate_z_mid(team,num_levs,z_int_v,z_mid_v);
      team.team_barrier();
      Kokkos::parallel_for(Kokkos::TeamVectorRange(team,num_levs), [&] (const Int& ilev) {
        dse_sub(ilev) = PF::calculate_dse(T_mid_v(icol,ilev),z_mid_v(ilev),phis_v(icol));
      });
      team.team_barrier();
    });
    Kokkos::fence();
    REQUIRE(views_are_equal(diag_out,dse_f));
  }
 
  // Finalize the diagnostic
  diag->finalize(); 

} // run()

TEST_CASE("dry_static_energy_test", "dry_static_energy_test]"){
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
