#include "catch2/catch.hpp"

#include "share/grid/mesh_free_grids_manager.hpp"
#include "diagnostics/vapor_flux.hpp"
#include "diagnostics/register_diagnostics.hpp"

#include "physics/share/physics_constants.hpp"

#include "share/util/eamxx_setup_random_test.hpp"
#include "share/util/eamxx_common_physics_functions.hpp"
#include "share/field/field_utils.hpp"

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
  using PC         = scream::physics::Constants<Real>;
  using KT         = ekat::KokkosTypes<DeviceT>;
  using ExecSpace  = typename KT::ExeSpace;
  using ESU        = ekat::ExeSpaceUtils<ExecSpace>;
  using MemberType = typename KT::MemberType;
  using view_1d    = typename KT::template view_1d<Real>;

  constexpr int num_levs = 33;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // Create a grids manager - single column for these tests
  const int ncols = 1; //TODO should be set to the size of the communication group.
  auto gm = create_gm(comm,ncols,num_levs);

  // Kokkos Policy
  auto policy = ESU::get_default_team_policy(ncols, num_levs);

  // Input (randomized) views
  view_1d qv("qv",num_levs),
          pseudo_density("pseudo_density",num_levs),
          u("u",num_levs),
          v("v",num_levs);


  // Construct random input data
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF pdf_qv(0.0,1e-3),
  pdf_uv(0.0,200),
  pdf_pseudo_density(1.0,100.0);

  // A time stamp
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  ekat::ParameterList params;
  register_diagnostics();
  auto& diag_factory = AtmosphereDiagnosticFactory::instance();

  REQUIRE_THROWS (diag_factory.create("VaporFlux",comm,params)); // No 'Wind Component'
  params.set<std::string>("Wind Component","foo");
  REQUIRE_THROWS (diag_factory.create("VaporFlux",comm,params)); // Invalid 'Wind Component'
  for (const std::string which_comp : {"Zonal", "Meridional"}) {
    // Construct the Diagnostic
    params.set<std::string>("Wind Component",which_comp);
    auto diag = diag_factory.create("VaporFlux",comm,params);
    diag->set_grids(gm);

    // Set the required fields for the diagnostic.
    std::map<std::string,Field> input_fields;
    for (const auto& req : diag->get_required_field_requests()) {
      Field f(req.fid);
      f.allocate_view();
      f.get_header().get_tracking().update_time_stamp(t0);
      diag->set_required_field(f.get_const());
      input_fields.emplace(f.name(),f);
    }

    // Initialize the diagnostic
    diag->initialize(t0,RunType::Initial);

    // Run tests
    {
      // Construct random data to use for test
      // Get views of input data and set to random values
      const auto& qv_f             = input_fields["qv"];
      const auto& pseudo_density_f = input_fields["pseudo_density"];
      const auto& horiz_winds_f    = input_fields["horiz_winds"];

      const auto& qv_v             = qv_f.get_view<Real**>();
      const auto& pseudo_density_v = pseudo_density_f.get_view<Real**>();
      const auto& horiz_winds_v    = horiz_winds_f.get_view<Real***>();

      for (int icol=0;icol<ncols;icol++) {
        const auto& qv_sub             = ekat::subview(qv_v,icol);
        const auto& pseudo_density_sub = ekat::subview(pseudo_density_v,icol);
        const auto& u_sub              = ekat::subview(horiz_winds_v,icol,0);
        const auto& v_sub              = ekat::subview(horiz_winds_v,icol,1);

        ekat::genRandArray(qv,             engine, pdf_qv);
        ekat::genRandArray(pseudo_density, engine, pdf_pseudo_density);
        ekat::genRandArray(u,              engine, pdf_uv);
        ekat::genRandArray(v,              engine, pdf_uv);
        Kokkos::deep_copy(qv_sub,qv);
        Kokkos::deep_copy(pseudo_density_sub,pseudo_density);
        Kokkos::deep_copy(u_sub,u);
        Kokkos::deep_copy(v_sub,v);
      }

      // Run diagnostic and compare with manual calculation
      diag->compute_diagnostic();
      const auto& diag_out = diag->get_diagnostic();
      Field qv_vert_integrated_flux_u_f = diag_out.clone();
      qv_vert_integrated_flux_u_f.deep_copy<double,Host>(0.0);
      qv_vert_integrated_flux_u_f.sync_to_dev();
      const auto& qv_vert_integrated_flux_u_v = qv_vert_integrated_flux_u_f.get_view<Real*>();
      constexpr Real g = PC::gravit;
      int comp = which_comp=="Zonal" ? 0 : 1;

      Kokkos::parallel_for("", policy, KOKKOS_LAMBDA(const MemberType& team) {
        const int icol = team.league_rank();
        auto wind = ekat::subview(horiz_winds_v,icol,comp);
        Kokkos::parallel_reduce(Kokkos::TeamVectorRange(team, num_levs),
                                [&] (const int& ilev, Real& lsum) {
          lsum += wind(ilev) * qv_v(icol,ilev) * pseudo_density_v(icol,ilev) / g;
        },qv_vert_integrated_flux_u_v(icol));

      });
      Kokkos::fence();
      REQUIRE(views_are_equal(diag_out,qv_vert_integrated_flux_u_f));
    }
   
    // Finalize the diagnostic
    diag->finalize(); 
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
  printf("ok!\n");
}

} // namespace
