#include "catch2/catch.hpp"

#include "share/grid/mesh_free_grids_manager.hpp"
#include "diagnostics/vertical_layer.hpp"
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
void run(std::mt19937_64& engine, std::string diag_type, const bool from_sea_level = false)
{
  using PF         = scream::PhysicsFunctions<DeviceT>;
  using PC         = scream::physics::Constants<Real>;
  using Pack       = ekat::Pack<Real,SCREAM_PACK_SIZE>;
  using KT         = ekat::KokkosTypes<DeviceT>;
  using ExecSpace  = typename KT::ExeSpace;
  using MemberType = typename KT::MemberType;
  using view_1d    = typename KT::template view_1d<Pack>;
  using rview_1d   = typename KT::template view_1d<Real>;
  using view_2d    = typename KT::template view_2d<Pack>;

  const     int packsize = SCREAM_PACK_SIZE;
  constexpr int num_levs = packsize*2 + 1; // Number of levels to use for tests, make sure the last pack can also have some empty slots (packsize>1).
  const     int num_mid_packs    = ekat::npack<Pack>(num_levs);
  const     int num_mid_packs_p1 = ekat::npack<Pack>(num_levs+1);

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
    return rview_1d(reinterpret_cast<Real*>(v.data()),v.size()*packsize);
  };

  // Construct random input data
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF pdf_qv(1e-6,1e-3),
       pdf_pseudodens(1.0,100.0),
       pdf_pres(0.0,PC::P0),
       pdf_temp(200.0,400.0),
       pdf_phis(0.0,10000.0);

  // A time stamp
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Construct the Diagnostic
  ekat::ParameterList params;
  std::string diag_name;
  if (diag_type == "thickness") {
    diag_name = "dz";
  }
  else if (diag_type == "interface") {
    diag_name = from_sea_level ? "z_int" : "geopotential_int";
  } else if (diag_type == "midpoint") {
    diag_name = from_sea_level ? "z_mid" : "geopotential_mid";
  }
  params.set<std::string>("diag_name", diag_name);
  register_diagnostics();
  auto& diag_factory = AtmosphereDiagnosticFactory::instance();
  auto diag = diag_factory.create(diag_name,comm,params);
  diag->set_grids(gm);

  // Helpful bools
  const bool only_compute_dz = (diag_type == "thickness");
  const bool is_interface_layout = (diag_type == "interface");
  const bool generate_phis_data = (not only_compute_dz and not from_sea_level);

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
    const auto& pseudo_dens_f = input_fields["pseudo_density"];
    const auto& pseudo_dens_v = pseudo_dens_f.get_view<Pack**>();
    const auto& p_mid_f       = input_fields["p_mid"];
    const auto& p_mid_v       = p_mid_f.get_view<Pack**>();
    const auto& qv_mid_f      = input_fields["qv"];
    const auto& qv_mid_v      = qv_mid_f.get_view<Pack**>();
    Field phis_f;
    rview_1d phis_v;
    if (generate_phis_data) {
      phis_f = input_fields["phis"];
      phis_v = phis_f.get_view<Real*>();
    }

    for (int icol=0;icol<ncols;icol++) {
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
    if (generate_phis_data) {
      ekat::genRandArray(phis_v, engine, pdf_phis);
    }

    // Run diagnostic and compare with manual calculation
    diag->compute_diagnostic();
    const auto& diag_out = diag->get_diagnostic();

    // Need to generate temporary values for calculation
    const auto& dz_v   = view_2d("",ncols, num_mid_packs);
    const auto& zmid_v = view_2d("",ncols, num_mid_packs);
    const auto& zint_v = view_2d("",ncols, num_mid_packs_p1);

    Kokkos::parallel_for("", policy, KOKKOS_LAMBDA(const MemberType& team) {
      const int icol = team.league_rank();

      const auto& dz_s = ekat::subview(dz_v,icol);
      Kokkos::parallel_for(Kokkos::TeamVectorRange(team,num_mid_packs), [&] (const Int& jpack) {
        dz_s(jpack) = PF::calculate_dz(pseudo_dens_v(icol,jpack),p_mid_v(icol,jpack),T_mid_v(icol,jpack),qv_mid_v(icol,jpack));
      });
      team.team_barrier();

      if (not only_compute_dz) {
        const auto& zint_s = ekat::subview(zint_v,icol);
        const Real surf_geopotential = from_sea_level ? 0.0 : phis_v(icol);
        PF::calculate_z_int(team,num_levs,dz_s,surf_geopotential,zint_s);

        if (not is_interface_layout) {
          const auto& zmid_s = ekat::subview(zmid_v,icol);
          PF::calculate_z_mid(team,num_levs,zint_s,zmid_s);
        }
      }
    });
    Kokkos::fence();

    Field diag_calc = diag_out.clone();
    auto field_v = diag_calc.get_view<Pack**>();

    if (diag_type == "thickness")      Kokkos::deep_copy(field_v, dz_v);
    else if (diag_type == "interface") Kokkos::deep_copy(field_v, zint_v);
    else if (diag_type == "midpoint")  Kokkos::deep_copy(field_v, zmid_v);
    diag_calc.sync_to_host();
    REQUIRE(views_are_equal(diag_out,diag_calc));
  }

  // Finalize the diagnostic
  diag->finalize();

} // run()

TEST_CASE("vertical_layer_test", "vertical_layer_test]"){
  // Run tests for both Real and Pack, and for (potentially) different pack sizes
  using scream::Real;
  using Device = scream::DefaultDevice;

  constexpr int num_runs = 5;

  auto engine = scream::setup_random_test();

  printf(" -> Number of randomized runs: %d, Pack<Real,%d> scalar type\n\n", num_runs, SCREAM_PACK_SIZE);

  printf(" -> Testing dz...");
  for (int irun=0; irun<num_runs; ++irun) {
    run<Device>(engine, "thickness");
  }
  printf("ok!\n");
  printf(" -> Testing z_int/geopotential_int...");
  for (int irun=0; irun<num_runs; ++irun) {
    run<Device>(engine, "interface", irun%2==0); // alternate from_sea_level=true/false
  }
  printf("ok!\n");
  printf(" -> Testing z_mid/geopotential_mid...");
  for (int irun=0; irun<num_runs; ++irun) {
    run<Device>(engine, "midpoint", irun%2==0); // alternate from_sea_level=true/false
  }
  printf("ok!\n");

  printf("\n");

} // TEST_CASE

} // namespace
