#include <catch2/catch.hpp>

#include "share/grid/gid_directory.hpp"
#include "share/core/eamxx_setup_random_test.hpp"

#include <algorithm>
#include <map>
#include <numeric>
#include <random>

namespace {

using namespace scream;
using gid_type = GidDirectory::gid_type;

// NOTE: we use CHECK (rather than REQUIRE), so that all ranks follow the same path in case of failure
//       on a subset of ranks (all the methods of GidDirectory are collective).

TEST_CASE ("gid_directory") {
  ekat::Comm comm(MPI_COMM_WORLD);
  auto engine = setup_random_test(&comm);

  const int nranks = comm.size();
  const int rank   = comm.rank();

  // Each rank owns a random set of gids (the number of gids can differ across ranks, and be 0)
  auto count_for = [](const int r) { return r%3==1 ? 0 : 5 + (r*7)%11; };
  std::vector<int> counts(nranks);
  int ngids = 0;
  for (int r=0; r<nranks; ++r) {
    counts[r] = count_for(r);
    ngids += counts[r];
  }

  std::vector<gid_type> all_gids (ngids);
  std::iota(all_gids.begin(),all_gids.end(),10);   // Do not start at 0, to not hide off-by-one errors
  if (comm.am_i_root()) {
    std::shuffle(all_gids.begin(),all_gids.end(),engine);
  }
  if (ngids>0) {
    comm.broadcast(all_gids.data(),ngids,comm.root_rank());
  }

  // Where each rank's gids start in the global list, and the global gid->(pid,lid) map
  std::map<gid_type,std::pair<int,int>> owners;
  int start = 0;
  for (int r=0; r<nranks; ++r) {
    for (int i=0; i<counts[r]; ++i) {
      owners[all_gids[start+i]] = {r,i};
    }
    start += counts[r];
  }
  start = std::accumulate(counts.begin(),counts.begin()+rank,0);
  std::vector<gid_type> my_gids (all_gids.begin()+start,all_gids.begin()+start+counts[rank]);

  SECTION ("unique_gids") {
    GidDirectory dir (comm,my_gids.data(),my_gids.size());
    CHECK (not dir.has_duplicates());

    // Ask for random gids. Some are not registered, and some are repeated
    std::uniform_int_distribution<int> pdf (0,ngids+ngids/2+3);
    std::vector<gid_type> queries (20);
    for (auto& q : queries) {
      q = pdf(engine);
    }
    queries.push_back(queries[0]);
    queries.push_back(queries[1]);

    std::vector<int> pids, lids;
    dir.lookup(queries.data(),queries.size(),pids,lids);
    CHECK (pids.size()==queries.size());
    CHECK (lids.size()==queries.size());
    for (size_t i=0; i<queries.size(); ++i) {
      auto it = owners.find(queries[i]);
      if (it==owners.end()) {
        CHECK (pids[i]==-1);
        CHECK (lids[i]==-1);
      } else {
        CHECK (pids[i]==it->second.first);
        CHECK (lids[i]==it->second.second);
      }
    }

    // Ranks can also ask nothing (and the call is still collective)
    dir.lookup(nullptr,0,pids,lids);
    CHECK (pids.size()==0);
  }

  SECTION ("duplicate_gids") {
    // Every rank also registers the gid 5, which makes it be registered several times
    // (by different ranks if nranks>1, or twice by the same rank otherwise).
    auto dup_gids = my_gids;
    dup_gids.push_back(5);
    if (nranks==1) {
      dup_gids.push_back(5);
    }
    GidDirectory dir (comm,dup_gids.data(),dup_gids.size());
    CHECK (dir.has_duplicates());

    std::vector<gid_type> queries = {5};
    if (not my_gids.empty()) {
      queries.push_back(my_gids[0]);
    }
    std::vector<int> pids, lids;
    dir.lookup(queries.data(),queries.size(),pids,lids);
    CHECK (pids[0]==-2);
    CHECK (lids[0]==-2);
    if (not my_gids.empty()) {
      auto it = owners.at(my_gids[0]);
      CHECK (pids[1]==it.first);
      CHECK (lids[1]==it.second);
    }
  }
}

} // anonymous namespace
