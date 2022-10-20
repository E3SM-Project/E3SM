#include <catch2/catch.hpp>

#include "share/grid/point_grid.hpp"
#include "share/grid/se_grid.hpp"
#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/grid/grid_utils.hpp"
#include "share/util/scream_setup_random_test.hpp"
#include "share/scream_types.hpp"

#include "ekat/ekat_pack.hpp"

#include <algorithm>

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

  auto shallow_copy = grid->clone("shallow",true);
  auto deep_copy    = grid->clone("deep",false);

  REQUIRE (shallow_copy->get_dofs_gids().data()==grid->get_dofs_gids().data());
  REQUIRE (deep_copy->get_dofs_gids().data()!=grid->get_dofs_gids().data());
  for (int i=0; i<grid->get_num_local_dofs(); ++i) {
    REQUIRE (deep_copy->get_dofs_gids_host()[i]==grid->get_dofs_gids_host()[i]);
  }

  shallow_copy->reset_num_vertical_lev(4);
  REQUIRE (shallow_copy->get_num_vertical_levels()==4);
}

TEST_CASE("se_grid", "") {
  ekat::Comm comm(MPI_COMM_WORLD);

  const int num_local_elems = 10;
  const int num_gp = 4;
  const int num_levels = 72;

  auto gm = create_mesh_free_grids_manager(comm,num_local_elems,num_gp,num_levels,0);
  gm->build_grids();

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

  auto shallow_copy = se_grid->clone("shallow",true);
  auto deep_copy    = se_grid->clone("deep",false);

  REQUIRE (shallow_copy->get_dofs_gids().data()==se_grid->get_dofs_gids().data());
  REQUIRE (deep_copy->get_dofs_gids().data()!=se_grid->get_dofs_gids().data());
  for (int i=0; i<se_grid->get_num_local_dofs(); ++i) {
    REQUIRE (deep_copy->get_dofs_gids_host()[i]==se_grid->get_dofs_gids_host()[i]);
  }
}

TEST_CASE ("get_owners") {
  ekat::Comm comm(MPI_COMM_WORLD);

  auto engine = setup_random_test(&comm);

  const int num_local_dofs = 10;
  const int num_global_dofs = num_local_dofs*comm.size();;
  const int offset = num_local_dofs*comm.rank();
  auto grid = std::make_shared<PointGrid>("grid",num_local_dofs,2,comm);

  // Create dofs, shuffled them around across ranks.
  using dofs_dev_t = AbstractGrid::dofs_list_type;
  using dofs_host_t = AbstractGrid::dofs_list_h_type;
  using gid_t = AbstractGrid::gid_type;

  dofs_host_t all_dofs_h ("",num_global_dofs);
  if (comm.am_i_root()) {
    std::iota(all_dofs_h.data(),all_dofs_h.data()+all_dofs_h.size(),0);
    std::shuffle(all_dofs_h.data(),all_dofs_h.data()+num_global_dofs,engine);
  }
  comm.broadcast(all_dofs_h.data(),num_global_dofs,comm.root_rank());

  dofs_dev_t dofs_d ("",num_local_dofs);
  auto dofs_h = Kokkos::create_mirror_view(dofs_d);
  std::memcpy (dofs_h.data(),all_dofs_h.data()+offset,num_local_dofs*sizeof(gid_t));
  Kokkos::deep_copy(dofs_d,dofs_h);

  grid->set_dofs(dofs_d);

  // Now, ask each rank to retrieve owners, and verify
  auto dofs_owners = grid->get_owners(all_dofs_h);
  REQUIRE (dofs_owners.extent(0)==all_dofs_h.size());

  for (int i=0; i<num_global_dofs; ++i) {
    const int pid = dofs_owners(i);
    const int expected_pid = i / num_local_dofs;
    REQUIRE (pid==expected_pid);
  }
}

} // anonymous namespace
