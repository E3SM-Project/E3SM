#include "catch2/catch.hpp"

#include "share/grid/mesh_free_grids_manager.hpp"
#include "diagnostics/relative_humidity.hpp"
#include "diagnostics/register_diagnostics.hpp"

#include "physics/share/physics_constants.hpp"
#include "physics/share/physics_functions.hpp" 

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
  using PC         = scream::physics::Constants<Real>;
  using Pack       = ekat::Pack<Real,SCREAM_PACK_SIZE>;
  using KT         = ekat::KokkosTypes<DeviceT>;
  using ExecSpace  = typename KT::ExeSpace;
  using MemberType = typename KT::MemberType;
  using view_1d    = typename KT::template view_1d<Pack>;
  using rview_1d   = typename KT::template view_1d<Real>;

  const     int packsize = SCREAM_PACK_SIZE;
  constexpr int num_levs = packsize*2 + 1; // Number of levels to use for tests, make sure the last pack can also have some empty slots (packsize>1).
  const     int num_mid_packs = ekat::npack<Pack>(num_levs);

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // Create a grids manager - single column for these tests
  const int ncols = 1;
  auto gm = create_gm(comm,ncols,num_levs);

  // Kokkos Policy
  auto policy = ekat::ExeSpaceUtils<ExecSpace>::get_default_team_policy(ncols, num_mid_packs);

  // Input (randomized) views, device
  view_1d temperature("temperature",num_mid_packs),
          pressure("pressure",num_mid_packs),
          pseudo_density("pseudo_density",num_mid_packs),
          pseudo_density_dry("pseudo_density_dry",num_mid_packs),
          qv("qv",num_mid_packs);

  auto dview_as_real = [&] (const view_1d& v) -> rview_1d {
    return rview_1d(reinterpret_cast<Real*>(v.data()),v.size()*packsize);
  };

  // Construct random input data
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF pdf_pres(10.0,PC::P0),
       pdf_temp(200.0,400.0),
       pdf_qv(0.0,1e-2),
       pdf_pseudo_density(1.0,100.0);

  // A time stamp
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Construct the Diagnostic
  ekat::ParameterList params;
  register_diagnostics();
  auto& diag_factory = AtmosphereDiagnosticFactory::instance();
  auto diag = diag_factory.create("RelativeHumidity",comm,params);
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

    // Field
    const auto& T_mid_f     = input_fields["T_mid"];
    // its device view
    const auto& T_mid_v     = T_mid_f.get_view<Pack**>();

    const auto& p_dry_mid_f = input_fields["p_dry_mid"];
    const auto& p_dry_mid_v = p_dry_mid_f.get_view<Pack**>();

    const auto& dpwet_f = input_fields["pseudo_density"];
    const auto& dpwet_v = dpwet_f.get_view<Pack**>();

    const auto& dpdry_f = input_fields["pseudo_density_dry"];
    const auto& dpdry_v = dpdry_f.get_view<Pack**>();

    const auto& qv_f = input_fields["qv"];
    const auto& qv_v = qv_f.get_view<Pack**>();

    for (int icol=0;icol<ncols;icol++) {
      const auto& T_sub = ekat::subview(T_mid_v,    icol);
      const auto& p_sub = ekat::subview(p_dry_mid_v,icol);
      const auto& dpwet_sub = ekat::subview(dpwet_v,icol);
      const auto& dpdry_sub = ekat::subview(dpdry_v,icol);
      const auto& qv_sub = ekat::subview(qv_v,icol);

      // init device arrays as random except for dpdry
      ekat::genRandArray(dview_as_real(temperature),   engine, pdf_temp);
      ekat::genRandArray(dview_as_real(pressure),      engine, pdf_pres);
      ekat::genRandArray(dview_as_real(pseudo_density),engine, pdf_pseudo_density);
      ekat::genRandArray(dview_as_real(qv),            engine, pdf_qv);

      Kokkos::deep_copy(T_sub,temperature);
      Kokkos::deep_copy(p_sub,pressure);
      Kokkos::deep_copy(qv_sub,qv);
      Kokkos::deep_copy(dpwet_sub,pseudo_density);

    }

    // Run diagnostic and compare with manual calculation
    Field rh_f = T_mid_f.clone();
    rh_f.deep_copy(0);
    const auto& rh_v = rh_f.get_view<Pack**>();
    using physics = scream::physics::Functions<Real, DefaultDevice>;
    using Smask = ekat::Mask<Pack::n>;
    Smask range_mask(true);
    Kokkos::parallel_for("", policy, KOKKOS_LAMBDA(const MemberType& team) {
      const int icol = team.league_rank();

      const auto& dpwet_sub = ekat::subview(dpwet_v,icol);
      const auto& dpdry_sub = ekat::subview(dpdry_v,icol);
      const auto& qv_sub = ekat::subview(qv_v,icol);

      Kokkos::parallel_for(Kokkos::TeamVectorRange(team,num_mid_packs), [&] (const Int& jpack) {
        dpdry_sub(jpack) = dpwet_sub(jpack) - dpwet_sub(jpack)*qv_sub(jpack);
        auto qv_sat_l = physics::qv_sat_dry(T_mid_v(icol,jpack), p_dry_mid_v(icol,jpack), true, range_mask);
        qv_sat_l *=  dpdry_v(icol,jpack) ;
        qv_sat_l /=  dpwet_v(icol,jpack) ;
        rh_v(icol,jpack) = qv_v(icol,jpack)/qv_sat_l;
      });
      team.team_barrier();
    });
    Kokkos::fence();

    // Run diagnostic and compare with manual calculation
    diag->compute_diagnostic();
    const auto& diag_out = diag->get_diagnostic();

    //in case one needs to look at values 
    //print_field_hyperslab(rh_f);
    //print_field_hyperslab(diag_out);

    REQUIRE(views_are_equal(diag_out,rh_f));
  }
 
  // Finalize the diagnostic
  diag->finalize(); 

} // run()

TEST_CASE("relative_humidity_test", "relative_humidity_test]"){
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
