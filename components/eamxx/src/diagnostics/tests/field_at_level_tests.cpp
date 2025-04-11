#include "catch2/catch.hpp"

#include "diagnostics/field_at_level.hpp"

#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/field/field_utils.hpp"
#include "share/util/eamxx_setup_random_test.hpp"

namespace scream {

std::shared_ptr<GridsManager>
create_gm (const ekat::Comm& comm, const int ncols, const int nlevs) {

  const int num_global_cols = ncols*comm.size();

  auto gm = create_mesh_free_grids_manager(comm,0,0,nlevs,num_global_cols);
  gm->build_grids();

  return gm;
}

TEST_CASE("field_at_level")
{
  using namespace ShortFieldTagsNames;
  using FL = FieldLayout;

  constexpr int packsize = SCREAM_PACK_SIZE;

  ekat::Comm comm(MPI_COMM_WORLD);

  auto engine = scream::setup_random_test(&comm);

  // Create a grids manager
  const int ncols = 3;
  const int nlevs = packsize*2 + 1;
  auto gm = create_gm(comm,ncols,nlevs);
  auto grid = gm->get_grid("Point Grid");

  // A time stamp
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Create input fields
  const auto units = ekat::units::Units::invalid();

  FieldIdentifier fid_mid ("M",FL({COL,CMP,LEV},{ncols,2,nlevs}),units,grid->name());
  FieldIdentifier fid_int ("I",FL({COL,LEV},{ncols,nlevs}),units,grid->name());

  Field f_mid (fid_mid);
  Field f_int (fid_int);

  f_mid.get_header().get_alloc_properties().request_allocation(packsize);
  f_int.get_header().get_alloc_properties().request_allocation(packsize);

  f_mid.allocate_view();
  f_int.allocate_view();

  f_mid.get_header().get_tracking().update_time_stamp(t0);
  f_int.get_header().get_tracking().update_time_stamp(t0);

  // Construct random input data
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF pdf(-1.0,1.0);

  randomize(f_mid,engine,pdf);
  randomize(f_int,engine,pdf);

  auto f_mid_1 = f_mid.get_component(1);

  ekat::ParameterList params_mid, params_int;
  params_mid.set("field_name",f_mid_1.name());
  params_int.set("field_name",f_int.name());
  params_mid.set("grid_name",grid->name());
  params_int.set("grid_name",grid->name());

  using IPDF = std::uniform_int_distribution<int>;
  IPDF ipdf (1,nlevs-2);
  for (const std::string lev_loc : {"model_top", "model_bot", "rand"}) {
    int lev;
    std::string lev_str;
    if (lev_loc=="bot") {
      lev = nlevs-1;
      lev_str = lev_loc;
    } else if (lev_loc=="model_top") {
      lev = 0;
      lev_str = lev_loc;
    } else {
      lev = ipdf(engine);
      lev_str = "lev_" + std::to_string(lev);
    }
    printf (" -> testing extraction at level: %s\n",lev_str.c_str());

    // Create and setup diagnostics
    params_mid.set<std::string>("vertical_location",lev_str);
    params_int.set<std::string>("vertical_location",lev_str);
    auto diag_mid = std::make_shared<FieldAtLevel>(comm,params_mid);
    auto diag_int = std::make_shared<FieldAtLevel>(comm,params_int);

    diag_mid->set_grids(gm);
    diag_int->set_grids(gm);

    diag_mid->set_required_field(f_mid_1);
    diag_int->set_required_field(f_int);

    diag_mid->initialize(t0,RunType::Initial);
    diag_int->initialize(t0,RunType::Initial);

    // Run diagnostics
    diag_mid->compute_diagnostic();
    diag_int->compute_diagnostic();

    // Check output
    auto d_mid = diag_mid->get_diagnostic();
    auto d_int = diag_int->get_diagnostic();
    d_mid.sync_to_host();
    d_int.sync_to_host();

    auto f_mid_v = f_mid_1.get_view<const Real**,Host>();
    auto f_int_v = f_int.get_view<const Real**,Host>();
    auto d_mid_v = d_mid.get_view<const Real*,Host>();
    auto d_int_v = d_int.get_view<const Real*,Host>();
    for (int icol=0; icol<ncols; ++icol) {
      REQUIRE (d_mid_v(icol)==f_mid_v(icol,lev));
      REQUIRE (d_int_v(icol)==f_int_v(icol,lev));
    }
  }
}

} // namespace scream
