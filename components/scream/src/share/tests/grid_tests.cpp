#include <catch2/catch.hpp>

#include "share/grid/simple_grid.hpp"
#include "share/grid/se_grid.hpp"
#include "share/grid/grids_manager.hpp"
#include "share/grid/grid_utils.hpp"

#include "ekat/ekat_pack.hpp"

namespace {

TEST_CASE("simple_grid", "") {
  using namespace scream;

  ekat::Comm comm(MPI_COMM_WORLD);
  int num_procs = comm.size();
  int num_cols = 128, num_levels = 72;
  SimpleGrid grid("my_grid", num_cols, num_levels, comm);
  REQUIRE(grid.type() == GridType::MeshFree);
  REQUIRE(grid.name() == "my_grid");
  REQUIRE(grid.get_num_vertical_levels() == num_levels);
  REQUIRE(grid.get_num_local_dofs() == num_cols / num_procs);

  auto layout = grid.get_native_dof_layout();
  REQUIRE(layout.tags().size() == 1);
  REQUIRE(layout.tag(0) == COL);
}

TEST_CASE("se_cell_based_grid", "") {
  using namespace scream;

  // Make the grid and check its initial state.
  int num_elems = 64;
  int num_gp = 4, num_levels = 72;
  SEGrid cell_grid("se_cell_grid", GridType::SE_CellBased,
                   num_elems, num_gp, num_levels);
  REQUIRE(grid.type() == GridType::CellBased);
  REQUIRE(grid.name() == "cell_se_grid");
  REQUIRE(grid.get_num_vertical_levels() == num_levels);
  REQUIRE(grid.get_num_local_dofs() == 0);

  // Set the degrees of freedom.
  SEGrid::dofs_list_type dofs("", num_elems*num_gp*num_gp);
  SEGrid::lid_to_idx_map_type dofs_map("", num_elems*num_gp*num_gp, 3);
  grid.set_dofs(dofs, dofs_map);

  REQUIRE(grid.get_num_local_dofs() == num_elems*num_gp*num_gp);

  // Test element x gauss pt x gauss pt mappings.
  auto dofs_gids = grid.get_dofs_gids();
  auto lid_to_idx = grid.get_lid_to_idx_map();
  grid.create_elgpgp_to_lid_map();
  auto idx_to_lid = grid.get_idx_to_lid_map();
  Kokkos::parallel_for(kokkos_types::RangePolicy(0,dofs_gids.extent_int(0)),
                       KOKKOS_LAMBDA(const int i) {
    REQUIRE(i == idx_to_lid(lid_to_idx(i, 0), lid_to_idx(i, 1), lid_to_idx(i, 2)));
  }
}

TEST_CASE("se_node_based_grid", "") {
  using namespace scream;

  // Make the grid and check its initial state.
  int num_elems = 64;
  int num_gp = 4, num_levels = 72;
  int num_cols = 6*num_elems*num_elems*9 + 2;
  SEGrid node_grid("se_node_grid", GridType::SE_NodeBased,
                   num_elems, num_gp, num_levels);
  REQUIRE(grid.type() == GridType::CellBased);
  REQUIRE(grid.name() == "cell_se_grid");
  REQUIRE(grid.get_num_vertical_levels() == num_levels);
  REQUIRE(grid.get_num_local_dofs() == num_cols);

  // Set the degrees of freedom.
  SEGrid::dofs_list_type dofs("", num_elems*num_gp*num_gp);
  SEGrid::lid_to_idx_map_type dofs_map("", num_elems*num_gp*num_gp, 3);
  grid.set_dofs(dofs, dofs_map);

  REQUIRE(grid.get_num_local_dofs() == num_cols);
}

} // anonymous namespace
