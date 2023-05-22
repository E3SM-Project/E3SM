#include "catch2/catch.hpp"

#include "share/grid/mesh_free_grids_manager.hpp"
#include "diagnostics/meridional_vapor_flux.hpp"
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
  gm_params.set<int>("number_of_global_columns", num_global_cols);
  gm_params.set<int>("number_of_local_elements", num_local_elems);
  gm_params.set<int>("number_of_vertical_levels", nlevs);
  gm_params.set<int>("number_of_gauss_points", np);

  auto gm = create_mesh_free_grids_manager(comm,gm_params);
  gm->build_grids();

  return gm;
}

//-----------------------------------------------------------------------------------------------//
template<typename DeviceT>
void run(std::mt19937_64& engine)
{
  using PC         = scream::physics::Constants<Real>;
  using Pack       = ekat::Pack<Real,SCREAM_PACK_SIZE>;
  using KT         = ekat::KokkosTypes<DeviceT>;
  using ExecSpace  = typename KT::ExeSpace;
  using MemberType = typename KT::MemberType;
  using view_1d    = typename KT::template view_1d<Pack>;
  using rview_1d   = typename KT::template view_1d<Real>;

  const     int packsize = SCREAM_PACK_SIZE;
  constexpr int num_levs = packsize*2 + 1; // Number of levels to use for tests, make sure the last pack can also have some empty slots (packsize>1).
  const     int num_mid_packs    = ekat::npack<Pack>(num_levs);

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // Create a grids manager - single column for these tests
  const int ncols = 1; //TODO should be set to the size of the communication group.
  auto gm = create_gm(comm,ncols,num_levs);

  // Kokkos Policy
  auto policy = ekat::ExeSpaceUtils<ExecSpace>::get_default_team_policy(ncols, num_mid_packs);

  // Input (randomized) views
  view_1d qv("qv",num_mid_packs),
          pseudo_density("pseudo_density",num_mid_packs),
          v("v",num_mid_packs);


  auto dview_as_real = [&] (const view_1d& vo) -> rview_1d {
    return rview_1d(reinterpret_cast<Real*>(vo.data()),vo.size()*packsize);
  };

  // Construct random input data
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF pdf_qv(0.0,1e-3),
  pdf_v(0.0,200),
  pdf_pseudo_density(1.0,100.0);

  // A time stamp
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Construct the Diagnostic
  ekat::ParameterList params;
  register_diagnostics();
  auto& diag_factory = AtmosphereDiagnosticFactory::instance();
  auto diag = diag_factory.create("MeridionalVapFlux",comm,params);
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
    input_fields.emplace(name,f);
  }

  // Initialize the diagnostic
  diag->initialize(t0,RunType::Initial);

  // Run tests
  {
    // Construct random data to use for test
    // Get views of input data and set to random values
    const auto& qv_f          = input_fields["qv"];
    const auto& qv_v          = qv_f.get_view<Pack**>();
    const auto& pseudo_density_f          = input_fields["pseudo_density"];
    const auto& pseudo_density_v          = pseudo_density_f.get_view<Pack**>();
    const auto& horiz_winds_f = input_fields["horiz_winds"];
    const auto& horiz_winds_v = horiz_winds_f.get_view<Pack***>();


    for (int icol=0;icol<ncols;icol++) {
      const auto& qv_sub      = ekat::subview(qv_v,icol);
      const auto& pseudo_density_sub      = ekat::subview(pseudo_density_v,icol);
      const auto& u_sub      = ekat::subview(horiz_winds_v,icol,0);
      const auto& v_sub      = ekat::subview(horiz_winds_v,icol,1);

      ekat::genRandArray(dview_as_real(qv),              engine, pdf_qv);
      ekat::genRandArray(dview_as_real(pseudo_density),              engine, pdf_pseudo_density);
      ekat::genRandArray(dview_as_real(v),              engine, pdf_v);
      Kokkos::deep_copy(qv_sub,qv);
      Kokkos::deep_copy(pseudo_density_sub,pseudo_density);
      // We shouldnt be using u, so make it a bad number.
      Kokkos::deep_copy(u_sub,std::numeric_limits<Real>::quiet_NaN());
      Kokkos::deep_copy(v_sub,v);

    }
   // ekat::genRandArray(phis_v, engine, pdf_surface);

    // Run diagnostic and compare with manual calculation
    diag->compute_diagnostic();
    const auto& diag_out = diag->get_diagnostic();
    Field qv_vert_integrated_flux_v_f = diag_out.clone();
    qv_vert_integrated_flux_v_f.deep_copy<double,Host>(0.0);
    qv_vert_integrated_flux_v_f.sync_to_dev();
    const auto& qv_vert_integrated_flux_v_v = qv_vert_integrated_flux_v_f.get_view<Real*>();
    constexpr Real gravit = PC::gravit;
    Kokkos::parallel_for("", policy, KOKKOS_LAMBDA(const MemberType& team) {
      const int icol = team.league_rank();
         Kokkos::parallel_reduce(Kokkos::TeamVectorRange(team, num_levs), [&] (const Int& idx, Real& lsum) {
      const int jpack = idx / Pack::n;
      const int klev  = idx % Pack::n;
      lsum += horiz_winds_v(icol,1,jpack)[klev] * qv_v(icol,jpack)[klev] * pseudo_density_v(icol,jpack)[klev]/gravit;
    },qv_vert_integrated_flux_v_v(icol));

    });
    Kokkos::fence();
    REQUIRE(views_are_equal(diag_out,qv_vert_integrated_flux_v_f));
  }
 
  // Finalize the diagnostic
  diag->finalize(); 

} // run()

TEST_CASE("meridional_vapor_flux_test", "meridional_vapor_flux_test]"){
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
