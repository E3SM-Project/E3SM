#include <catch2/catch.hpp>

#include "share/grid/point_grid.hpp"
#include "share/grid/grid_import_export.hpp"

#include <ekat_assert.hpp>
#include <ekat_comm.hpp>

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <limits>
#include <map>
#include <numeric>
#include <set>

/*
 * Compare the construction of a GridImportExport (and AbstractGrid::is_unique) with
 * the "legacy" algorithm, where each rank loops over all other ranks, which broadcast
 * their gids (so that the cost grows with the number of ranks), against the current one,
 * which uses a distributed directory of gids.
 *
 * The test checks that the two algorithms give exactly the same import/export data, and
 * prints how long each takes (max over ranks).
 *
 * To compare the performance on a large number of ranks, run the executable directly
 * (it does not need to be launched by ctest), e.g.,
 *
 *   IMPEXP_NLOCAL=20000 mpirun -n 4096 ./grid_imp_exp_bench
 *
 * The following env vars can be used to change the size of the problem (all optional):
 *   IMPEXP_NLOCAL: the number of gids owned by each rank in the unique grid (default: 100)
 *   IMPEXP_HALO:   how many extra gids (in each direction) each rank needs in the overlapped grid (default: 5)
 *   IMPEXP_NRAND:  how many random gids each rank also needs, in the overlapped grid (default: 5)
 * The overlapped grid on each rank contains gids owned by that rank, by "neighboring" ranks (the halo),
 * and by random ranks, mimicking what happens with the grids in a horizontal remapper.
 */

namespace {

using namespace scream;
using gid_type = AbstractGrid::gid_type;
using clock_type = std::chrono::steady_clock;

int get_env_int (const char* name, const int def) {
  const char* v = std::getenv(name);
  return v ? std::atoi(v) : def;
}

double seconds_since (const clock_type::time_point& t0) {
  return std::chrono::duration<double>(clock_type::now()-t0).count();
}

// ----------------- The legacy algorithms (to be used as a reference) ------------------ //

struct LegacyImpExp {
  std::vector<int> imp_pids, imp_lids, exp_pids, exp_lids;
};

// This is what AbstractGrid::is_unique used to do
bool legacy_is_unique (const ekat::Comm& comm, std::vector<gid_type> dofs)
{
  const int nldofs = dofs.size();
  std::sort(dofs.begin(),dofs.end());
  auto unique_end = std::unique(dofs.begin(),dofs.end());

  int locally_unique = unique_end==dofs.end();
  int unique;
  comm.all_reduce(&locally_unique,&unique,1,MPI_PROD);
  if (unique==0) {
    return false;
  }

  int max_dofs;
  comm.all_reduce(&nldofs,&max_dofs,1,MPI_MAX);
  std::vector<gid_type> gids(max_dofs,-1);
  int unique_gids = 1;

  for (int pid=0; pid<comm.size(); ++pid) {
    // Rank pid broadcasts its gids, everyone else checks if there are duplicates
    if (pid==comm.rank()) {
      std::copy(dofs.begin(),dofs.end(),gids.data());
    }

    int ndofs = nldofs;
    comm.broadcast(&ndofs,1,pid);
    comm.broadcast(gids.data(),ndofs,pid);

    int my_unique_gids = 1;
    if (pid!=comm.rank()) {
      int i=0, j=0;
      while (i<nldofs && j<ndofs && my_unique_gids==1) {
        if (dofs[i]<gids[j]) {
          ++i;
        } else if (dofs[i]>gids[j]) {
          ++j;
        } else {
          my_unique_gids = 0;
          break;
        }
      }
    }
    comm.all_reduce(&my_unique_gids,&unique_gids,1,MPI_PROD);
    if (unique_gids==0) {
      break;
    }
  }
  return unique_gids==1;
}

// This is what the GridImportExport constructor used to do
LegacyImpExp legacy_import_export (const AbstractGrid& unique, const AbstractGrid& overlapped)
{
  const auto& comm = unique.get_comm();

  auto gids = unique.get_dofs_gids().clone(CloneFlags::CopyData);
  auto ov_gids = overlapped.get_dofs_gids().clone(CloneFlags::CopyData);
  auto gids_h = gids.get_view<gid_type*,Host>();
  auto ov_gids_h = ov_gids.get_view<gid_type*,Host>();
  const int num_ov_gids = ov_gids_h.size();

  LegacyImpExp result;

  gid_type* data;
  std::vector<gid_type> pid_gids;
  std::map<int,std::vector<int>> pid2lids;
  int num_imports = 0;

  // Imports
  auto ov_gid2lid = overlapped.get_gid2lid_map();
  for (int pid=0; pid<comm.size(); ++pid) {
    int num_gids_pid = gids_h.size();
    comm.broadcast(&num_gids_pid,1,pid);
    if (pid==comm.rank()) {
      data = gids_h.data();
    } else {
      pid_gids.resize(num_gids_pid);
      data = pid_gids.data();
    }
    comm.broadcast(data,num_gids_pid,pid);

    for (int i=0; i<num_gids_pid; ++i) {
      auto it = ov_gid2lid.find(data[i]);
      if (it!=ov_gid2lid.end()) {
        pid2lids[pid].push_back(it->second);
        ++num_imports;
      }
    }
  }
  EKAT_REQUIRE_MSG (num_ov_gids==num_imports,
      "Error! Could not locate the owner of one of the dst grid GIDs.\n");
  for (const auto& [pid,lids] : pid2lids) {
    for (auto lid : lids) {
      result.imp_lids.push_back(lid);
      result.imp_pids.push_back(pid);
    }
  }

  // Exports
  auto gid2lid = unique.get_gid2lid_map();
  pid2lids.clear();
  for (int pid=0; pid<comm.size(); ++pid) {
    int num_gids_pid = ov_gids_h.size();
    comm.broadcast(&num_gids_pid,1,pid);
    if (pid==comm.rank()) {
      data = ov_gids_h.data();
    } else {
      pid_gids.resize(num_gids_pid);
      data = pid_gids.data();
    }
    comm.broadcast(data,num_gids_pid,pid);

    for (int i=0; i<num_gids_pid; ++i) {
      auto it = gid2lid.find(data[i]);
      if (it!=gid2lid.end()) {
        pid2lids[pid].push_back(it->second);
      }
    }
    std::sort(pid2lids[pid].begin(),pid2lids[pid].end());
  }
  for (const auto& [pid,lids] : pid2lids) {
    for (auto lid : lids) {
      result.exp_lids.push_back(lid);
      result.exp_pids.push_back(pid);
    }
  }

  return result;
}

// ----------------- The test ------------------ //

template<typename ViewT>
std::vector<int> to_vector (const ViewT& v) {
  auto h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),v);
  return std::vector<int>(h.data(),h.data()+h.size());
}

TEST_CASE ("grid_import_export_bench") {
  ekat::Comm comm(MPI_COMM_WORLD);

  const int nranks = comm.size();
  const int rank   = comm.rank();

  const int nl    = get_env_int("IMPEXP_NLOCAL",100);
  const int halo  = get_env_int("IMPEXP_HALO",5);
  const int nrand = get_env_int("IMPEXP_NRAND",5);

  const long long N = static_cast<long long>(nl)*nranks;
  REQUIRE (N<static_cast<long long>(std::numeric_limits<gid_type>::max()));

  // A bijection of [0,N), so that the gids are not owned by ranks in the order they are numbered.
  // Pick a number that is coprime with N
  long long a = 7919;
  while (std::gcd(a,N)!=1) { ++a; }
  auto perm = [&](const long long j) -> gid_type {
    return static_cast<gid_type>(((a%N)*(((j%N)+N)%N) + 17) % N);
  };

  // Unique grid: rank r owns the (scrambled) gids in [r*nl,(r+1)*nl)
  auto unique = std::make_shared<PointGrid>("unique",nl,0,comm);
  {
    auto f = unique->get_dofs_gids();
    auto h = f.get_view<gid_type*,Host>();
    for (int i=0; i<nl; ++i) {
      h(i) = perm(static_cast<long long>(rank)*nl+i);
    }
    f.sync_to_dev();
  }

  // Overlapped grid: gids of this rank, plus a halo on both sides, plus a few random ones
  std::set<long long> ov_idx;
  for (long long j=static_cast<long long>(rank)*nl-halo; j<static_cast<long long>(rank+1)*nl+halo; ++j) {
    ov_idx.insert(((j%N)+N)%N);
  }
  unsigned long long state = 12345 + 7919ULL*rank;
  for (int k=0; k<nrand; ++k) {
    state = state*6364136223846793005ULL + 1442695040888963407ULL;
    ov_idx.insert((state>>33) % N);
  }
  auto overlapped = std::make_shared<PointGrid>("overlapped",ov_idx.size(),0,comm);
  {
    auto f = overlapped->get_dofs_gids();
    auto h = f.get_view<gid_type*,Host>();
    int i = 0;
    for (auto j : ov_idx) {
      h(i++) = perm(j);
    }
    f.sync_to_dev();
  }

  // Legacy algorithm
  comm.barrier();
  auto t0 = clock_type::now();
  auto unique_gids_h = unique->get_dofs_gids().get_view<const gid_type*,Host>();
  const bool legacy_unique = legacy_is_unique(comm,std::vector<gid_type>(unique_gids_h.data(),unique_gids_h.data()+nl));
  double t_legacy_unique = seconds_since(t0);

  comm.barrier();
  t0 = clock_type::now();
  auto legacy = legacy_import_export(*unique,*overlapped);
  double t_legacy = seconds_since(t0);

  // Current algorithm
  comm.barrier();
  t0 = clock_type::now();
  const bool new_unique = unique->is_unique();
  double t_new_unique = seconds_since(t0);

  comm.barrier();
  t0 = clock_type::now();
  GridImportExport imp_exp (unique,overlapped);
  double t_new = seconds_since(t0);

  // Check that they agree
  CHECK (legacy_unique);
  CHECK (new_unique);
  CHECK (legacy.imp_pids==to_vector(imp_exp.import_pids()));
  CHECK (legacy.imp_lids==to_vector(imp_exp.import_lids()));
  CHECK (legacy.exp_pids==to_vector(imp_exp.export_pids()));
  CHECK (legacy.exp_lids==to_vector(imp_exp.export_lids()));

  // Report the max time over ranks
  double times[4] = {t_legacy_unique, t_new_unique, t_legacy, t_new};
  double max_times[4];
  comm.all_reduce(times,max_times,4,MPI_MAX);
  int counts[2] = {imp_exp.num_imports(), imp_exp.num_exports()};
  int max_counts[2];
  comm.all_reduce(counts,max_counts,2,MPI_MAX);
  if (comm.am_i_root()) {
    printf("\n  GridImportExport construction: ranks=%d, gids per rank=%d (global: %lld), halo=%d, random=%d\n",
           nranks,nl,N,halo,nrand);
    printf("  (max imports/rank: %d, max exports/rank: %d. Times are the max over ranks, in seconds)\n",
           max_counts[0],max_counts[1]);
    printf("                       legacy      directory   speedup\n");
    printf("    is_unique        %9.3e   %9.3e   %6.1fx\n",max_times[0],max_times[1],max_times[0]/max_times[1]);
    printf("    import/export    %9.3e   %9.3e   %6.1fx\n",max_times[2],max_times[3],max_times[2]/max_times[3]);
    printf("    (the legacy import/export constructor also called is_unique first;\n"
           "     the new one checks uniqueness while building the directory)\n\n");
  }
}

} // anonymous namespace
