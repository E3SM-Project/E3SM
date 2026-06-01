#include "catch2/catch.hpp"

#include "share/diagnostics/field_at_level.hpp"
#include "share/grid/point_grid.hpp"
#include "share/field/field_utils.hpp"
#include "share/core/eamxx_setup_random_test.hpp"

namespace scream {


TEST_CASE("field_at_level")
{
  using namespace ShortFieldTagsNames;

  constexpr int packsize = SCREAM_PACK_SIZE;

  ekat::Comm comm(MPI_COMM_WORLD);

  int seed = get_random_test_seed(&comm);
  std::mt19937_64 engine(seed);

  // Create a grids manager
  const int ncols = 3;
  const int nlevs = packsize*2 + 1;
  auto grid = create_point_grid("physics",ncols,nlevs,comm);

  // A time stamp
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Create input fields
  auto scalar_3d = grid->get_3d_scalar_layout(LEV);
  auto vector_3d = grid->get_3d_vector_layout(LEV,2,"dim2");
  FieldIdentifier fid_mid ("M",vector_3d,ekat::units::none,grid->name());
  FieldIdentifier fid_int ("I",scalar_3d,ekat::units::none,grid->name());

  Field f_mid (fid_mid);
  Field f_int (fid_int);

  f_mid.get_header().get_alloc_properties().request_allocation(packsize);
  f_int.get_header().get_alloc_properties().request_allocation(packsize);

  f_mid.allocate_view();
  f_int.allocate_view();

  f_mid.get_header().get_tracking().update_time_stamp(t0);
  f_int.get_header().get_tracking().update_time_stamp(t0);

  // Construct random input data
  randomize_uniform(f_mid,seed++,-1,1);
  randomize_uniform(f_int,seed++,-1,1);

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
    auto diag_mid = std::make_shared<FieldAtLevel>(comm,params_mid,grid);
    auto diag_int = std::make_shared<FieldAtLevel>(comm,params_int,grid);

    diag_mid->set_input_field(f_mid_1);
    diag_int->set_input_field(f_int);

    diag_mid->initialize();
    diag_int->initialize();

    // Run diagnostics
    diag_mid->compute(t0);
    diag_int->compute(t0);

    // Check output
    auto d_mid = diag_mid->get();
    auto d_int = diag_int->get();
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
