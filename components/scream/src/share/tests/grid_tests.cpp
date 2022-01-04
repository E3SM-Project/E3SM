#include <catch2/catch.hpp>

#include "share/grid/point_grid.hpp"
#include "share/grid/se_grid.hpp"
#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/grid/grid_utils.hpp"
#include "share/scream_types.hpp"

#include "ekat/ekat_pack.hpp"

namespace {

using namespace scream;
using namespace scream::ShortFieldTagsNames;

TEST_CASE("point_grid", "") {

  ekat::Comm comm(MPI_COMM_WORLD);

  const int num_procs = comm.size();
  const int num_local_cols = 128;
  const int num_global_cols = num_local_cols*num_procs;
  const int num_levels = 72;

  auto grid = create_point_grid("my_grid", num_global_cols, num_levels, comm);
  REQUIRE(grid->type() == GridType::Point);
  REQUIRE(grid->name() == "my_grid");
  REQUIRE(grid->get_num_vertical_levels() == num_levels);
  REQUIRE(grid->get_num_local_dofs()  == num_local_cols);
  REQUIRE(grid->get_num_global_dofs() == num_global_cols);
  REQUIRE(grid->is_unique());

  // Point grids should have (global) gids spanning the interval [min_gid, min_gid+num_global_dofs)
  const auto max_gid = grid->get_global_max_dof_gid();
  const auto min_gid = grid->get_global_min_dof_gid();
  REQUIRE( (max_gid-min_gid+1)==grid->get_num_global_dofs() );

  auto lid_to_idx = grid->get_lid_to_idx_map();
  auto host_lid_to_idx = Kokkos::create_mirror_view(lid_to_idx);
  Kokkos::deep_copy(host_lid_to_idx, lid_to_idx);
  for (int i = 0; i < grid->get_num_local_dofs(); ++i) {
    REQUIRE(host_lid_to_idx.extent_int(1) == 1);
    REQUIRE(i == host_lid_to_idx(i, 0));
  }

  auto layout = grid->get_2d_scalar_layout();
  REQUIRE(layout.tags().size() == 1);
  REQUIRE(layout.tag(0) == COL);
}

TEST_CASE("se_grid", "") {
  ekat::Comm comm(MPI_COMM_WORLD);

  const int num_local_elems = 10;
  const int num_gp = 4;
  const int num_levels = 72;

  auto gm = create_mesh_free_grids_manager(comm,num_local_elems,num_gp,num_levels,0);
  gm->build_grids(std::set<std::string>{"SE Grid"});

  // SE grid
  auto se_grid = gm->get_grid("SE Grid");

  REQUIRE(se_grid->type() == GridType::SE);
  REQUIRE(se_grid->name() == "SE Grid");
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
}

} // anonymous namespace
