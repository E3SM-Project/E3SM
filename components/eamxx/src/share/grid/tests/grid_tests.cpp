#include <catch2/catch.hpp>

#include "share/grid/point_grid.hpp"
#include "share/core/eamxx_setup_random_test.hpp"
#include "share/core/eamxx_types.hpp"

#include <algorithm>

namespace scream {

TEST_CASE ("get_owners") {
  ekat::Comm comm(MPI_COMM_WORLD);

  auto engine = setup_random_test(&comm);

  const int num_local_dofs = 10;
  const int num_global_dofs = num_local_dofs*comm.size();;
  const int offset = num_local_dofs*comm.rank();
  auto grid = std::make_shared<PointGrid>("grid",num_local_dofs,2,comm);

  // Create dofs, shuffled them around across ranks.
  using gid_type = AbstractGrid::gid_type;

  std::vector<gid_type> all_dofs (num_global_dofs);
  if (comm.am_i_root()) {
    std::iota(all_dofs.data(),all_dofs.data()+all_dofs.size(),0);
    std::shuffle(all_dofs.data(),all_dofs.data()+num_global_dofs,engine);
  }
  comm.broadcast(all_dofs.data(),num_global_dofs,comm.root_rank());

  auto dofs = grid->get_dofs_gids();
  auto dofs_h = dofs.get_view<gid_type*,Host>();
  std::memcpy (dofs_h.data(),all_dofs.data()+offset,num_local_dofs*sizeof(gid_type));
  dofs.sync_to_dev();

  // Now, ask each rank to retrieve owners, and verify
  auto dofs_owners = grid->get_owners(all_dofs);
  REQUIRE (dofs_owners.size()==all_dofs.size());

  for (int i=0; i<num_global_dofs; ++i) {
    const int pid = dofs_owners[i];
    const int expected_pid = i / num_local_dofs;
    REQUIRE (pid==expected_pid);
  }
}

TEST_CASE ("gid2lid_map") {
  using gid_type = AbstractGrid::gid_type;

  ekat::Comm comm(MPI_COMM_WORLD);

  auto engine = setup_random_test(&comm);

  const int num_local_dofs = 10;
  const int num_global_dofs = num_local_dofs*comm.size();;
  // Create dofs, shuffled them around across ranks.
  std::vector<gid_type> all_dofs (num_global_dofs);
  if (comm.am_i_root()) {
    std::iota(all_dofs.data(),all_dofs.data()+all_dofs.size(),0);
    std::shuffle(all_dofs.data(),all_dofs.data()+num_global_dofs,engine);
  }
  comm.broadcast(all_dofs.data(),num_global_dofs,comm.root_rank());

  // Create a grid, grabbing my portion of the all_dofs array
  const int offset = num_local_dofs*comm.rank();
  auto grid = std::make_shared<PointGrid>("grid",num_local_dofs,0,comm);
  auto dofs = grid->get_dofs_gids();
  auto dofs_h = dofs.get_view<gid_type*,Host>();
  std::memcpy (dofs_h.data(),all_dofs.data()+offset,num_local_dofs*sizeof(gid_type));
  dofs.sync_to_dev();

  auto gid2lid = grid->get_gid2lid_map();
  for (const auto& it : gid2lid) {
    REQUIRE (it.first==dofs_h[it.second]);
  }
}

TEST_CASE ("get_remote_pids_and_lids") {
  using gid_type = AbstractGrid::gid_type;

  ekat::Comm comm(MPI_COMM_WORLD);

  auto engine = setup_random_test(&comm);

  const int num_local_dofs = 10;
  const int num_global_dofs = num_local_dofs*comm.size();;
  // Create dofs, shuffled them around across ranks.
  std::vector<gid_type> all_dofs (num_global_dofs);
  if (comm.am_i_root()) {
    std::iota(all_dofs.data(),all_dofs.data()+all_dofs.size(),0);
    std::shuffle(all_dofs.data(),all_dofs.data()+num_global_dofs,engine);
  }
  comm.broadcast(all_dofs.data(),num_global_dofs,comm.root_rank());

  // Create a grid, grabbing my portion of the all_dofs array
  const int offset = num_local_dofs*comm.rank();
  auto grid = std::make_shared<PointGrid>("grid",num_local_dofs,0,comm);
  auto dofs = grid->get_dofs_gids();
  auto dofs_h = dofs.get_view<gid_type*,Host>();
  std::memcpy (dofs_h.data(),all_dofs.data()+offset,num_local_dofs*sizeof(gid_type));
  dofs.sync_to_dev();

  // Now, ask each rank to retrieve owners and local ids, and verify
  std::vector<int> pids, lids;
  grid->get_remote_pids_and_lids(all_dofs,pids,lids);
  REQUIRE (pids.size()==all_dofs.size());
  REQUIRE (lids.size()==all_dofs.size());

  for (int i=0; i<num_global_dofs; ++i) {
    const int expected_pid = i / num_local_dofs;
    const int expected_lid = i % num_local_dofs;

    REQUIRE (pids[i]==expected_pid);
    REQUIRE (lids[i]==expected_lid);
  }
}

} // namespace scream
