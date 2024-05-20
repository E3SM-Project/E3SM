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
template<typename DeviceT, typename Engine>
void run(      Engine& engine,
               std::string  diag_name,
         const std::string& location)
{
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

  // A time stamp
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Kokkos Policy
  auto policy = ekat::ExeSpaceUtils<ExecSpace>::get_default_team_policy(ncols, num_mid_packs);

  // Input (randomized) views
  view_1d temperature("temperature",num_mid_packs),
          pseudodensity("pseudo_density",num_mid_packs),
          pressure("pressure",num_mid_packs),
          watervapor("watervapor",num_mid_packs);

  auto dview_as_real = [&] (const view_1d& v) -> rview_1d {
    return rview_1d(reinterpret_cast<Real*>(v.data()),v.size()*packsize);
  };

  register_diagnostics();
  auto& diag_factory = AtmosphereDiagnosticFactory::instance();
  
  // Construct the Diagnostic
  ekat::ParameterList params;
  if (location=="midpoints") {
    diag_name += "_mid";
  } else if (location=="interfaces") {
    diag_name += "_int";
  }
  params.set<std::string>("diag_name", diag_name);
  auto diag = diag_factory.create(diag_name,comm,params);
  diag->set_grids(gm);

  const bool needs_phis = diag_name=="z" or diag_name=="geopotential";

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

  // Note: we are not testing the calculate_dz utility. We are testing
  //       the diag class, so use some inputs that make checking results easier
  //       With these inputs, T_virt=T, and dz=8*rd/g
  const Real dz = PC::RD/PC::gravit;
  const Real phis = 3;
  input_fields["T_mid"].deep_copy(Real(4));
  input_fields["p_mid"].deep_copy(Real(2));
  input_fields["pseudo_density"].deep_copy(Real(4));
  input_fields["qv"].deep_copy(Real(0));
  if (needs_phis) {
    input_fields["phis"].deep_copy(phis);
  }

  // Initialize and run the diagnostic
  diag->initialize(t0,RunType::Initial);
  diag->compute_diagnostic();
  const auto& diag_out = diag->get_diagnostic();
  diag_out.sync_to_host();
  auto d_h = diag_out.get_view<Real**,Host>();

  // Compare against expecte value
  const auto last_lev = location=="interfaces" ? num_levs : num_levs-1;
  for (int icol=0; icol<ncols; ++icol) {
    for (int ilev=last_lev; ilev>=0; --ilev) {
      int num_mid_levs_below = ilev-last_lev;
      Real tgt_mid, tgt_int;
      if (diag_name=="dz") {
        tgt_mid = dz;
      } else if (diag_name=="altitude") {
        tgt_int = phis/PC::gravit + num_mid_levs_below*dz;
        tgt_mid = tgt_int + dz/2;
      } else if (diag_name=="geopotential") {
        tgt_int = phis + num_mid_levs_below*dz*PC::gravit;
        tgt_mid = tgt_int + PC::gravit*dz/2;
      } else {
        tgt_int = num_mid_levs_below*dz;
        tgt_mid = tgt_int + dz/2;
      }

      if (location=="interfaces") {
        REQUIRE (d_h(icol,ilev)==tgt_int);
      } else {
        REQUIRE (d_h(icol,ilev)==tgt_mid);
      }
    }
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

  printf("Test specs\n");
  printf(" - number of randomized runs: %d\n",num_runs);
  printf(" - scalar type: Pack<Real,%d>\n\n", SCREAM_PACK_SIZE);

  for (int irun=0; irun<num_runs; ++irun) {
    for (std::string loc : {"midpoints","interfaces"}) {
      for (std::string diag : {"geopotential","altitude","z"}) {
        printf(" -> Testing diag=%s at %s ...\n",diag.c_str(),loc.c_str());
        run<Device>(engine, loc, diag);
        printf(" -> Testing diag=%s at %s ... PASS!\n",diag.c_str(),loc.c_str());
      }
    }
    printf(" -> Testing diag=dz ...\n");
    run<Device>(engine, "", "dz");
    printf(" -> Testing diag=dz ... PASS!\n");
  }

} // TEST_CASE

} // namespace
