#include <catch2/catch.hpp>

#include "share/grid/se_grid.hpp"
#include "share/grid/grid_utils.hpp"
#include "share/core/eamxx_types.hpp"

#include <algorithm>

namespace scream {


TEST_CASE("se_grid", "") {
  using namespace scream::ShortFieldTagsNames;
  ekat::Comm comm(MPI_COMM_WORLD);

  const int num_local_elems = 10;
  const int num_gp = 4;
  const int num_levels = 72;

  // SE grid
  auto se_grid = std::make_shared<SEGrid>("se_grid",num_local_elems,num_gp,num_levels,comm);
  {
    auto nldofs = se_grid->get_num_local_dofs();
    auto dof_offset = nldofs;
    comm.scan(&dof_offset,1,MPI_SUM);
    dof_offset -= nldofs; // comm.scan(..) is an INCLUSIVE prefix op
    auto dofs_gids = se_grid->get_dofs_gids();
    auto h_dofs_gids = dofs_gids.get_view<AbstractGrid::gid_type*,Host>();
    std::iota(h_dofs_gids.data(),h_dofs_gids.data()+nldofs,dof_offset);
    dofs_gids.sync_to_dev();
  }

  REQUIRE(se_grid->type() == GridType::SE);
  REQUIRE(se_grid->name() == "se_grid");
  REQUIRE(se_grid->get_num_vertical_levels() == num_levels);
  REQUIRE(se_grid->get_num_local_dofs() == num_local_elems*num_gp*num_gp);

  auto layout = se_grid->get_2d_scalar_layout();
  REQUIRE(layout.tags().size() == 3);
  REQUIRE(layout.tag(0) == EL);
  REQUIRE(layout.tag(1) == GP);
  REQUIRE(layout.tag(2) == GP);

  REQUIRE (se_grid->is_unique());

  const auto max_gid = se_grid->get_global_max_dof_gid();
  const auto min_gid = se_grid->get_global_min_dof_gid();
  REQUIRE( (max_gid-min_gid+1)==se_grid->get_num_global_dofs() );

  auto shallow_copy = se_grid->clone("shallow",true);
  auto deep_copy    = se_grid->clone("deep",false);

  using gid_type = AbstractGrid::gid_type;

  auto grid_gids = se_grid->get_dofs_gids().get_view<const gid_type*,Host>();
  auto scopy_gids = shallow_copy->get_dofs_gids().get_view<const gid_type*,Host>();
  auto dcopy_gids = deep_copy->get_dofs_gids().get_view<const gid_type*,Host>();
  REQUIRE (scopy_gids.data()==grid_gids.data());
  REQUIRE (dcopy_gids.data()!=grid_gids.data());
  for (int i=0; i<se_grid->get_num_local_dofs(); ++i) {
    REQUIRE (dcopy_gids[i]==grid_gids[i]);
  }
}

} // anonymous namespace
