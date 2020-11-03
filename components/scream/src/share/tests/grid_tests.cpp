#include <catch2/catch.hpp>

#include "share/grid/point_grid.hpp"
#include "share/grid/se_grid.hpp"
#include "share/grid/grids_manager.hpp"
#include "share/grid/grid_utils.hpp"
#include "share/scream_types.hpp"

#include "ekat/ekat_pack.hpp"

namespace {

using namespace scream;
using namespace scream::ShortFieldTagsNames;

TEST_CASE("point_grid", "") {

  ekat::Comm comm(MPI_COMM_WORLD);
  int num_procs = comm.size();

  const int num_cols = 128*num_procs, num_levels = 72;

  auto grid = create_point_grid("my_grid", num_cols, num_levels, comm);
  REQUIRE(grid.type() == GridType::Point);
  REQUIRE(grid.name() == "my_grid");
  REQUIRE(grid.get_num_vertical_levels() == num_levels);
  REQUIRE(grid.get_num_local_dofs() == 128);

  auto lid_to_idx = grid.get_lid_to_idx_map();
  auto host_lid_to_idx = Kokkos::create_mirror_view(lid_to_idx);
  Kokkos::deep_copy(host_lid_to_idx, lid_to_idx);
  for (int i = 0; i < grid.get_num_local_dofs(); ++i) {
    REQUIRE(host_lid_to_idx.extent_int(1) == 1);
    REQUIRE(i == host_lid_to_idx(i, 0));
  }

  auto layout = grid.get_native_dof_layout();
  REQUIRE(layout.tags().size() == 1);
  REQUIRE(layout.tag(0) == COL);
}

TEST_CASE("se_grid", "") {

  // Make the grid and check its initial state.
  int num_elems = 64;
  int num_gp = 4, num_levels = 72;
  SEGrid grid("se_grid",num_elems, num_gp, num_levels);
  REQUIRE(grid.type() == GridType::SE);
  REQUIRE(grid.name() == "se_grid");
  REQUIRE(grid.get_num_vertical_levels() == num_levels);
  REQUIRE(grid.get_num_local_dofs() == num_elems*num_gp*num_gp);

  auto layout = grid.get_native_dof_layout();
  REQUIRE(layout.tags().size() == 3);
  REQUIRE(layout.tag(0) == EL);
  REQUIRE(layout.tag(1) == GP);
  REQUIRE(layout.tag(2) == GP);

  // Set up the degrees of freedom.
  SEGrid::dofs_list_type dofs("", num_elems*num_gp*num_gp);
  auto host_dofs = Kokkos::create_mirror_view(dofs);
  SEGrid::lid_to_idx_map_type dofs_map("", num_elems*num_gp*num_gp, 3);
  auto host_dofs_map = Kokkos::create_mirror_view(dofs_map);
  for (int ie = 0; ie < num_elems; ++ie) {
    for (int igp = 0; igp < num_gp; ++igp) {
      for (int jgp = 0; jgp < num_gp; ++jgp) {
        int idof = ie*num_gp*num_gp + igp*num_gp + jgp;
        host_dofs(idof) = idof;
        host_dofs_map(idof, 0) = ie;
        host_dofs_map(idof, 1) = igp;
        host_dofs_map(idof, 2) = jgp;
      }
    }
  }

  // Move the data to the device and set the DOFs.
  Kokkos::deep_copy(dofs, host_dofs);
  Kokkos::deep_copy(dofs_map, host_dofs_map);
  grid.set_dofs(dofs, dofs_map);
}

} // anonymous namespace
