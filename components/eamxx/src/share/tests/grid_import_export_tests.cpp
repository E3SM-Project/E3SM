// Access some internal catch2 types.
// #define CATCH_CONFIG_EXTERNAL_INTERFACES
#include <catch2/catch.hpp>

#include "share/grid/point_grid.hpp"
#include "share/grid/grid_import_export.hpp"
#include "share/util/eamxx_setup_random_test.hpp"
#include "share/eamxx_types.hpp"

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
  auto imp_pids = imp_exp.import_pids_h();
  auto imp_lids = imp_exp.import_lids_h();
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
  auto exp_pids = imp_exp.export_pids_h();
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

  // Test gather
  if (comm.am_i_root()) {
    printf(" -> Testing gather routine ....\n");
  }
  ok  = true;

  // We design the tests in such a way that, serially, the data stored for
  // gid N is [0,...,N]. We want this to also be the case in parallel, but the
  // same GID may be owned by 2+ processors in the ov grid. This would cause duplicate
  // entries upon gather completion, making checking the results harder.
  // To avoid this, we store a map gid->pids so that all ranks know how many ranks
  // have a give gid in the ov grid. When filling up the ov_data for gather,
  // rank P will only fill data for gid N if it is the smallest gid owning it.
  std::map<gid_type,std::vector<int>> gid2pids;
  for (int pid=0; pid<comm.size(); ++pid) {
    std::vector<gid_type> pid_gids;
    int n = ov_gids.size();
    comm.broadcast(&n,1,pid);
    pid_gids.resize(n);
    if (pid==comm.rank()) {
      for (int i=0; i<n; ++i) {
        pid_gids[i] = ov_gids[i];
      }
    }
    comm.broadcast(pid_gids.data(),n,pid);

    for (auto g : pid_gids) {
      gid2pids[g].push_back(pid);
    }
  }
  for (auto& [gid,pids] : gid2pids) {
    std::sort(pids.begin(),pids.end());
  }

  std::map<int,std::vector<Real>> ov_data;
  for (int i=0; i<ov_grid->get_num_local_dofs(); ++i) {
    auto& v = ov_data[i];
    auto gid = ov_gids[i];
    if (comm.rank()==gid2pids[gid][0]) {
      for (int k=0; k<=gid; ++k) {
        v.push_back(k);
      }
    }
  }

  std::map<int,std::vector<Real>> data;
  imp_exp.gather(ekat::get_mpi_type<Real>(),ov_data,data);

  for (auto [lid,v] : data) {
    auto gid = gids[lid];
    CHECK (static_cast<int>(v.size())==(gid+1));
    ok &= catch_capture.lastAssertionPassed();

    std::vector<Real> expected(gid+1);
    std::iota(expected.begin(),expected.end(),0);

    // Values may be gathered in a way that is just a permutation
    // of what we'd get serially, which is [0,...,gid-1]
    std::sort(v.begin(),v.end());
    CHECK (std::equal(expected.begin(),expected.end(),v.begin()));
    ok &= catch_capture.lastAssertionPassed();
  }
  if (comm.am_i_root()) {
    printf(" -> Testing gather routine .... %s\n",ok ? "PASS" : "FAIL");
  }

  // Test scatter
  if (comm.am_i_root()) {
    printf(" -> Testing scatter routine ...\n");
  }
  ok = true;
  data.clear();
  for (int i=0; i<grid->get_num_local_dofs(); ++i) {
    auto& v = data[i];
    auto gid = gids[i];
    v.resize(gid+1);
    std::iota(v.begin(),v.end(),0);
  }
  ov_data.clear();
  imp_exp.scatter(ekat::get_mpi_type<Real>(),data,ov_data);

  for (const auto& [lid,v] : ov_data) {
    auto gid = ov_gids[lid];
    CHECK (static_cast<int>(v.size())==(gid+1));
    ok &= catch_capture.lastAssertionPassed();

    std::vector<Real> expected(gid+1);
    std::iota(expected.begin(),expected.end(),0);
    CHECK (std::equal(expected.begin(),expected.end(),v.begin()));
    ok &= catch_capture.lastAssertionPassed();
  }
  if (comm.am_i_root()) {
    printf(" -> Testing scatter routine ... %s\n",ok ? "PASS" : "FAIL");
  }
}

} // anonymous namespace
