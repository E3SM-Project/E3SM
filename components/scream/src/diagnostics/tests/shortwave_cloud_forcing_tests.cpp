#include "catch2/catch.hpp"

#include "share/grid/mesh_free_grids_manager.hpp"
#include "diagnostics/shortwave_cloud_forcing.hpp"
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

template<typename VT>
typename VT::HostMirror
cmvdc (const VT& v) {
  auto vh = Kokkos::create_mirror_view(v);
  Kokkos::deep_copy(vh,v);
  return vh;
}

//-----------------------------------------------------------------------------------------------//
template<typename DeviceT>
void run(std::mt19937_64& engine)
{
  using PF         = scream::PhysicsFunctions<DeviceT>;
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
  constexpr Real gravit = PC::gravit;
  constexpr Real macheps = PC::macheps;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // Create a grids manager - single column for these tests
  const int ncols = 1;
  auto gm = create_gm(comm,ncols,num_levs);

  // Kokkos Policy
  auto policy = ekat::ExeSpaceUtils<ExecSpace>::get_default_team_policy(ncols, num_mid_packs);

  // Input (randomized) views
  view_1d
    pseudo_density("pseudo_density",num_mid_packs),
    SW_flux_dn("SW_flux_dn",num_mid_packs),
    SW_flux_up("SW_flux_up",num_mid_packs),
    SW_clrsky_flux_dn("SW_clrsky_flux_dn",num_mid_packs),
    SW_clrsky_flux_up("SW_clrsky_flux_up",num_mid_packs);

  auto dview_as_real = [&] (const view_1d& v) -> rview_1d {
    return rview_1d(reinterpret_cast<Real*>(v.data()),v.size()*packsize);
  };

  // Construct random input data
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF 
    pdf_SW_flux_x(0.0,400),
    pdf_pseudodens(1.0,100.0),
    pdf_alpha(0.1,0.9);
  // Liquid

  // A time stamp
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Construct the Diagnostics
  std::map<std::string,std::shared_ptr<AtmosphereDiagnostic>> diags;
  auto& diag_factory = AtmosphereDiagnosticFactory::instance();
  register_diagnostics();
  ekat::ParameterList params;
  params.set<std::string>("Grid", "Point Grid");
  // Vapor
  params.set<std::string>("Diagnostic Name", "Shortwave Cloud Forcing");
  auto diag_vap = diag_factory.create("ShortwaveCloudForcing",comm,params);
  diag_vap->set_grids(gm);
//  diags.emplace("vwp",diag_vap);


  // Set the required fields for the diagnostic.
  std::map<std::string,Field> input_fields;
  for (const auto& dd : diags) {
    const auto& diag = dd.second;
    for (const auto& req : diag->get_required_field_requests()) {
      if (input_fields.find(req.fid.name())==input_fields.end()) {
        Field f(req.fid);
        auto & f_ap = f.get_header().get_alloc_properties();
        f_ap.request_allocation(packsize);
        f.allocate_view();
        const auto name = f.name();
        f.get_header().get_tracking().update_time_stamp(t0);
        diag->set_required_field(f.get_const());
        REQUIRE_THROWS(diag->set_computed_field(f));
        input_fields.emplace(name,f);
      } else {
        auto& f = input_fields[req.fid.name()];
        const auto name = f.name();
        f.get_header().get_tracking().update_time_stamp(t0);
        diag->set_required_field(f.get_const());
        REQUIRE_THROWS(diag->set_computed_field(f));
        input_fields.emplace(name,f);
      }
    }
    // Initialize the diagnostic
    diag->initialize(t0,RunType::Initial);
  }


  // Run tests
  {
    // Construct random data to use for test
    // Get views of input data and set to random values
    const auto& SW_flux_dn_f          = input_fields["SW_flux_dn"];
    const auto& SW_flux_dn_v          = SW_flux_dn_f.get_view<Pack**>();
    const auto& SW_flux_up_f          = input_fields["SW_flux_up"];
    const auto& SW_flux_up_v          = SW_flux_up_f.get_view<Pack**>();
    const auto& SW_clrsky_flux_dn_f          = input_fields["SW_clrsky_flux_dn"];
    const auto& SW_clrsky_flux_dn_v          = SW_clrsky_flux_dn_f.get_view<Pack**>();
    const auto& SW_clrsky_flux_up_f          = input_fields["SW_clrsky_flux_up"];
    const auto& SW_clrsky_flux_up_v          = SW_clrsky_flux_up_f.get_view<Pack**>();
    const auto& pseudo_dens_f = input_fields["pseudo_density"];
    const auto& pseudo_dens_v = pseudo_dens_f.get_view<Pack**>();
    for (int icol=0;icol<ncols;icol++) {
      const auto& SW_flux_dn_sub      = ekat::subview(SW_flux_dn_v,icol);
      const auto& SW_flux_up_sub      = ekat::subview(SW_flux_up_v,icol);
      const auto& SW_clrsky_flux_dn_sub      = ekat::subview(SW_clrsky_flux_dn_v,icol);
      const auto& SW_clrsky_flux_up_sub      = ekat::subview(SW_clrsky_flux_up_v,icol);
      const auto& dp_sub      = ekat::subview(pseudo_dens_v,icol);
      ekat::genRandArray(dview_as_real(pseudo_density),  engine, pdf_pseudodens);
      Kokkos::deep_copy(dp_sub,pseudo_density);

      ekat::genRandArray(dview_as_real(SW_flux_dn),              engine, pdf_SW_flux_x);
      ekat::genRandArray(dview_as_real(SW_flux_up),              engine, pdf_SW_flux_x);
      ekat::genRandArray(dview_as_real(SW_clrsky_flux_dn),              engine, pdf_SW_flux_x);
      ekat::genRandArray(dview_as_real(SW_clrsky_flux_up),              engine, pdf_SW_flux_x);
      Kokkos::deep_copy(SW_flux_dn_sub,SW_flux_dn);
      Kokkos::deep_copy(SW_flux_up_sub,SW_flux_up);
      Kokkos::deep_copy(SW_clrsky_flux_dn_sub,SW_clrsky_flux_dn);
      Kokkos::deep_copy(SW_clrsky_flux_up_sub,SW_clrsky_flux_up);
    }
    // Grab views for each of the water path diagnostics
    const auto& SWCF = diags["SWCF"]->get_diagnostic();
    const auto& SWCF_v = SWCF.get_view<const Real*>();
    const auto& SWCF_h = SWCF.get_view<const Real*,Host>();

    for (const auto& dd : diags) {
      dd.second->compute_diagnostic();
    }
    // Test 2: If the cell-wise mass is scaled by constant alpha then the water
    //         path should also be scaled by alpha.
    {
      Field SWCF_copy_f = diags["SWCF"]->get_diagnostic().clone();
      const auto& SWCF_copy_v = SWCF_copy_f.get_view<Real*>();

      const auto alpha_SW_flux_dn = pdf_alpha(engine);
      const auto alpha_SW_flux_up = pdf_alpha(engine);
      const auto alpha_SW_clrsky_flux_dn = pdf_alpha(engine);
      const auto alpha_SW_clrsky_flux_up = pdf_alpha(engine);
      REQUIRE(alpha_SW_flux_dn*alpha_SW_flux_up*alpha_SW_clrsky_flux_dn*alpha_SW_clrsky_flux_up != 1.0);

      Kokkos::parallel_for("",ncols*num_mid_packs,KOKKOS_LAMBDA(const int& idx) {
        const int icol  = idx / num_mid_packs;
        const int jpack = idx % num_mid_packs;

        SW_flux_dn_v(icol,jpack) *= alpha_SW_flux_dn;
        SW_flux_up_v(icol,jpack) *= alpha_SW_flux_up;
        SW_clrsky_flux_dn_v(icol,jpack) *= alpha_SW_clrsky_flux_dn;
        SW_clrsky_flux_up_v(icol,jpack) *= alpha_SW_clrsky_flux_up;

        if (jpack==0) {
          SW_flux_dn_copy_v(icol) *= alpha_SW_flux_dn;
          SW_flux_up_copy_v(icol) *= alpha_SW_flux_up;
          SW_clrsky_flux_dn_copy_v(icol) *= alpha_SW_clrsky_flux_dn;
          SW_clrsky_flux_up_copy_v(icol) *= alpha_SW_clrsky_flux_up;
        }
      });
      Kokkos::fence();
      for (const auto& dd : diags) {
        dd.second->compute_diagnostic();
      }
      SWCF_copy_f.sync_to_host();
      const auto& SWCF_copy_h = SWCF_copy_f.get_view<Real*,Host>();
      SWCF.sync_to_host();
      for (int icol=0;icol<ncols;icol++) {
        REQUIRE(std::abs(SWCF_copy_h(icol)-SWCF_h(icol))<macheps);
      }
    }
  // Finalize the diagnostic
  for (const auto& dd : diags) {
    const auto& diag = dd.second;
    diag->finalize();
  } 

} // run()

TEST_CASE("shortwave_cloud_forcing_test", "shortwave_cloud_forcing_test]"){
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
