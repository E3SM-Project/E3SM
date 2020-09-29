#include <catch2/catch.hpp>

#include "share/grid/simple_grid.hpp"
#include "share/grid/se_grid.hpp"
#include "share/grid/grids_manager.hpp"
#include "share/grid/grid_utils.hpp"
#include "share/scream_types.hpp"

#include "ekat/ekat_pack.hpp"

namespace {

using namespace scream;
using namespace scream::ShortFieldTagsNames;
using device_type = DefaultDevice;
using kokkos_types = KokkosTypes<device_type>;

TEST_CASE("simple_grid", "") {

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

  // Make the grid and check its initial state.
  int num_elems = 64;
  int num_gp = 4, num_levels = 72;
  SEGrid grid("se_cell_grid", GridType::SE_CellBased,
              num_elems, num_gp, num_levels);
  REQUIRE(grid.type() == GridType::SE_CellBased);
  REQUIRE(grid.name() == "se_cell_grid");
  REQUIRE(grid.get_num_vertical_levels() == num_levels);
  REQUIRE(grid.get_num_local_dofs() == 0);

  // Set up the degrees of freedom.
  SEGrid::dofs_list_type dofs("", num_elems*num_gp*num_gp);
  SEGrid::lid_to_idx_map_type dofs_map("", num_elems*num_gp*num_gp, 3);
  for (int ie = 0; ie < num_elems; ++ie) {
    for (int igp = 0; igp < num_gp; ++igp) {
      for (int jgp = 0; jgp < num_gp; ++jgp) {
        int idof = ie*num_gp*num_gp + igp*num_gp + jgp;
        dofs_map(idof, 0) = ie;
        dofs_map(idof, 1) = igp;
        dofs_map(idof, 2) = jgp;
      }
    }
  }
  grid.set_dofs(dofs, dofs_map);

  REQUIRE(grid.get_num_local_dofs() == num_elems*num_gp*num_gp);

  auto layout = grid.get_native_dof_layout();
  REQUIRE(layout.tags().size() == 3);
  REQUIRE(layout.tag(0) == EL);
  REQUIRE(layout.tag(1) == GP);
  REQUIRE(layout.tag(2) == GP);

  // Test element x gauss pt x gauss pt mappings.
  grid.create_elgpgp_to_lid_map();
  auto dofs_gids = grid.get_dofs_gids();
  auto lid_to_idx = grid.get_lid_to_idx_map();
  auto idx_to_lid = grid.get_idx_to_lid_map();
  Kokkos::parallel_for(kokkos_types::RangePolicy(0,dofs_gids.extent_int(0)),
                       KOKKOS_LAMBDA(const int i) {
    int ie  = lid_to_idx(i,0);
    int igp = lid_to_idx(i,1);
    int jgp = lid_to_idx(i,2);
    REQUIRE(i == idx_to_lid(ie, igp, jgp));
  });
}

TEST_CASE("se_node_based_grid", "") {

  // Make the grid and check its initial state.
  int num_elems = 64;
  int num_gp = 4, num_levels = 72;
  int num_cols = 6*num_elems*num_elems*9 + 2;
  SEGrid grid("se_node_grid", GridType::SE_NodeBased,
              num_elems, num_gp, num_levels);
  REQUIRE(grid.type() == GridType::SE_NodeBased);
  REQUIRE(grid.name() == "se_node_grid");
  REQUIRE(grid.get_num_vertical_levels() == num_levels);

  // Set the degrees of freedom.
  SEGrid::dofs_list_type dofs("", num_cols);
  SEGrid::lid_to_idx_map_type dofs_map("", num_cols, 1);
  grid.set_dofs(dofs, dofs_map);
  REQUIRE(grid.get_num_local_dofs() == num_cols);

  auto layout = grid.get_native_dof_layout();
  REQUIRE(layout.tags().size() == 1);
  REQUIRE(layout.tag(0) == COL);

  REQUIRE(grid.get_num_local_dofs() == num_cols);
}

} // anonymous namespace
