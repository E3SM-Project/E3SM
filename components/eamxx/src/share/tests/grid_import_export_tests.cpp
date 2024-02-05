// Access some internal catch2 types.
// #define CATCH_CONFIG_EXTERNAL_INTERFACES
#include <catch2/catch.hpp>

#include "share/grid/point_grid.hpp"
#include "share/grid/grid_import_export.hpp"
#include "share/util/scream_setup_random_test.hpp"
#include "share/scream_types.hpp"

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
  const int overlap    = 5;
  const int nldofs_src = 10;
  const int ngdofs     = nldofs_src*comm.size();;

  // Create the unique grid
  auto src_grid = create_point_grid("src",ngdofs,0,comm);
  auto src_gids = src_grid->get_dofs_gids().get_view<const gid_type*, Host>();

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
  auto start = all_dofs.data() + nldofs_src*comm.rank();
  auto end   = start + nldofs_src;
  end   += last ? 0 : overlap;
  start -= first ? 0 : overlap;
  const int ndofs_dst = nldofs_src + (first ? 0 : overlap) + (last ? 0 : overlap);

  auto dst_grid = std::make_shared<PointGrid>("grid",ndofs_dst,0,comm);
  auto dst_gids_field = dst_grid->get_dofs_gids();
  auto dst_gids = dst_gids_field.get_view<gid_type*,Host>();

  std::copy (start,end,dst_gids.data());
  dst_gids_field.sync_to_dev();

  GridImportExport imp_exp(src_grid,dst_grid);

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
    auto gid = dst_gids[lid];
    auto pid = gid / nldofs_src;
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
  auto src_gid2lid = src_grid->get_gid2lid_map();
  for (int pid=0; pid<comm.size(); ++pid) {
    pid2lid[pid].resize(0);
    int pid_start = pid*nldofs_src;
    int pid_end   = pid_start+nldofs_src;
    if (pid!=0) {
      pid_start -= overlap;
    }
    if (pid!=comm.size()-1) {
      pid_end += overlap;
    }
    for (int i=pid_start;i<pid_end; ++i) {
      auto gid = all_dofs[i];
      if (src_gid2lid.count(gid)==1) {
        pid2lid[pid].push_back(src_gid2lid.at(gid));
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
  std::map<int,std::vector<Real>> src_data;
  for (int i=0; i<dst_grid->get_num_local_dofs(); ++i) {
    auto& v = src_data[i];
    auto gid = dst_gids[i];
    v.resize(gid);
    std::iota(v.begin(),v.end(),0);
  }
  std::map<int,std::vector<Real>> dst_data;
  imp_exp.gather(ekat::get_mpi_type<Real>(),src_data,dst_data);

  std::vector<size_t> num_imp_per_lid (nldofs_src,0);
  for (size_t i=0; i<exp_lids.size(); ++i) {
    ++num_imp_per_lid[exp_lids[i]];
  }
  // std::cout << "num imp per lid: " << ekat::join(num_imp_per_lid," ") << "\n";
  for (const auto& it : dst_data) {
    auto lid = it.first;
    auto gid = src_gids[lid];
    const auto& v = it.second;
    CHECK (v.size()==gid*num_imp_per_lid[lid]);
    ok &= catch_capture.lastAssertionPassed();

    std::vector<Real> expected(gid);
    std::iota(expected.begin(),expected.end(),0);
    for (size_t imp=0; imp<num_imp_per_lid[lid]; ++imp) {
      CHECK (std::equal(expected.begin(),expected.end(),v.begin()+gid*imp));
      ok &= catch_capture.lastAssertionPassed();
    }
  }
  if (comm.am_i_root()) {
    printf(" -> Testing gather routine .... %s\n",ok ? "PASS" : "FAIL");
  }

  // Test scatter
  if (comm.am_i_root()) {
    printf(" -> Testing scatter routine ...\n");
  }
  ok = true;
  src_data.clear();
  // std::cout << "src_data:\n";
  for (int i=0; i<src_grid->get_num_local_dofs(); ++i) {
    auto& v = src_data[i];
    auto gid = src_gids[i];
    v.resize(gid);
    std::iota(v.begin(),v.end(),0);
    // std::cout << " " << gid << ":" << ekat::join(v," ") << "\n";
  }
  dst_data.clear();
  imp_exp.scatter(ekat::get_mpi_type<Real>(),src_data,dst_data);

  // std::cout << "dst_data:\n";
  for (const auto& it : dst_data) {
    auto lid = it.first;
    auto gid = dst_gids[lid];
    const auto& v = it.second;
    CHECK (static_cast<int>(v.size())==gid);
    ok &= catch_capture.lastAssertionPassed();
    // std::cout << " " << gid << ":" << ekat::join(v," ") << "\n";

    std::vector<Real> expected(gid);
    std::iota(expected.begin(),expected.end(),0);
    CHECK (std::equal(expected.begin(),expected.end(),v.begin()));
    ok &= catch_capture.lastAssertionPassed();
  }
  if (comm.am_i_root()) {
    printf(" -> Testing scatter routine ... %s\n",ok ? "PASS" : "FAIL");
  }
}

} // anonymous namespace
