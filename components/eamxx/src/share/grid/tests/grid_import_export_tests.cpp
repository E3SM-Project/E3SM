// Access some internal catch2 types.
// #define CATCH_CONFIG_EXTERNAL_INTERFACES
#include <catch2/catch.hpp>

#include "share/grid/point_grid.hpp"
#include "share/grid/grid_import_export.hpp"
#include "share/util/eamxx_alltoallv.hpp"
#include "share/core/eamxx_setup_random_test.hpp"
#include "share/core/eamxx_types.hpp"

#include <algorithm>

namespace {

using namespace scream;
using namespace scream::ShortFieldTagsNames;

TEST_CASE ("grid_import_export") {
  using gid_type = AbstractGrid::gid_type;

  ekat::Comm comm(MPI_COMM_WORLD);
  MPI_Comm_set_errhandler(comm.mpi_comm(),MPI_ERRORS_RETURN);

  auto& catch_capture = Catch::getResultCapture();

  // return;
  auto engine = setup_random_test(&comm);

  bool ok;
  const int overlap = 5;
  const int nldofs  = 10;
  const int ngdofs  = nldofs*comm.size();;

  // Create the unique grid
  auto grid = create_point_grid("src",ngdofs,0,comm);
  auto gids = grid->get_dofs_gids().get_view<const gid_type*, Host>();

  // For the dst grid, shuffle dofs around randomly. Then,
  // have each rank grab a few extra dofs
  std::vector<gid_type> all_dofs (ngdofs);
  if (comm.am_i_root()) {
    std::iota(all_dofs.data(),all_dofs.data()+all_dofs.size(),0);
    std::shuffle(all_dofs.data(),all_dofs.data()+ngdofs,engine);
  }
  comm.broadcast(all_dofs.data(),ngdofs,comm.root_rank());

  const bool first = comm.rank()==0;
  const bool last  = comm.rank()==(comm.size()-1);
  auto start = all_dofs.data() + nldofs*comm.rank();
  auto end   = start + nldofs;
  end   += last ? 0 : overlap;
  start -= first ? 0 : overlap;
  const int nldofs_ov = nldofs + (first ? 0 : overlap) + (last ? 0 : overlap);

  auto ov_grid = std::make_shared<PointGrid>("grid",nldofs_ov,0,comm);
  auto ov_gids_field = ov_grid->get_dofs_gids();
  auto ov_gids = ov_gids_field.get_view<gid_type*,Host>();

  std::copy (start,end,ov_gids.data());
  ov_gids_field.sync_to_dev();

  GridImportExport imp_exp(grid,ov_grid);

  // Test import views
  if (comm.am_i_root()) {
    printf(" -> Testing import views ......\n");
  }
  ok = true;
  auto to_host = [](const auto& v) { return Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),v); };
  auto imp_pids = to_host(imp_exp.import_pids());
  auto imp_lids = to_host(imp_exp.import_lids());
  std::set<int> all_imp_lids;
  for (size_t i=0; i<imp_pids.size(); ++i) {
    auto lid = imp_lids[i];
    auto gid = ov_gids[lid];
    auto pid = gid / nldofs;
    CHECK (imp_pids[i]==pid);
    ok &= catch_capture.lastAssertionPassed();
    all_imp_lids.insert(lid);
  }
  CHECK (all_imp_lids.size()==imp_pids.size());
  ok &= catch_capture.lastAssertionPassed();
  if (comm.am_i_root()) {
    printf(" -> Testing import views ...... %s\n",ok ? "PASS" : "FAIL");
  }
  
  // Test export views
  if (comm.am_i_root()) {
    printf(" -> Testing export views ......\n");
  }
  ok = true;
  auto exp_pids = to_host(imp_exp.export_pids());
  auto exp_lids = imp_exp.export_lids_h();
  std::map<int,std::vector<int>> pid2lid;
  auto gid2lid = grid->get_gid2lid_map();
  for (int pid=0; pid<comm.size(); ++pid) {
    pid2lid[pid].resize(0);
    int pid_start = pid*nldofs;
    int pid_end   = pid_start+nldofs;
    if (pid!=0) {
      pid_start -= overlap;
    }
    if (pid!=comm.size()-1) {
      pid_end += overlap;
    }
    for (int i=pid_start;i<pid_end; ++i) {
      auto gid = all_dofs[i];
      if (gid2lid.count(gid)==1) {
        pid2lid[pid].push_back(gid2lid.at(gid));
      }
    }
    std::sort(pid2lid[pid].begin(),pid2lid[pid].end());
  }

  for (int pid=0,k=0; pid<comm.size(); ++pid) {
    const auto& lids = pid2lid[pid];
    const int n = lids.size();
    for (int i=0; i<n; ++i,++k) {
      CHECK (exp_pids[k]==pid);
      ok &= catch_capture.lastAssertionPassed();
      CHECK (exp_lids[k]==lids[i]);
      ok &= catch_capture.lastAssertionPassed();
    }
  }
  if (comm.am_i_root()) {
    printf(" -> Testing export views ...... %s\n",ok ? "PASS" : "FAIL");
  }

  // Test that the p2p messages are consistent: the gids that a rank imports from another rank
  // must be, in the same order, the ones that the other rank exports to it.
  if (comm.am_i_root()) {
    printf(" -> Testing import/export consistency ......\n");
  }
  ok = true;
  {
    const int nranks = comm.size();
    auto imp_offsets = imp_exp.import_pid_offsets_h();
    auto exp_offsets = imp_exp.export_pid_offsets_h();
    std::vector<int> send_counts(nranks), recv_counts(nranks);
    for (int p=0; p<nranks; ++p) {
      send_counts[p] = imp_offsets[p+1]-imp_offsets[p];
      recv_counts[p] = exp_offsets[p+1]-exp_offsets[p];
    }
    std::vector<gid_type> send_gids (imp_lids.size()), recv_gids;
    for (size_t k=0; k<imp_lids.size(); ++k) {
      send_gids[k] = ov_gids[imp_lids[k]];
    }
    std::vector<int> recv_counts_out = recv_counts;
    alltoallv(comm,send_gids,send_counts,recv_gids,recv_counts_out,true);

    CHECK (recv_gids.size()==exp_lids.size());
    ok &= catch_capture.lastAssertionPassed();
    for (size_t k=0; k<exp_lids.size() && k<recv_gids.size(); ++k) {
      CHECK (recv_gids[k]==gids[exp_lids[k]]);
      ok &= catch_capture.lastAssertionPassed();
    }
  }
  if (comm.am_i_root()) {
    printf(" -> Testing import/export consistency ...... %s\n",ok ? "PASS" : "FAIL");
  }
}

} // anonymous namespace
