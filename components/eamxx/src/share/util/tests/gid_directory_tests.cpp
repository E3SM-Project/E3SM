#include <catch2/catch.hpp>

#include "share/util/gid_directory.hpp"
#include "share/core/eamxx_setup_random_test.hpp"

#include <algorithm>
#include <map>
#include <numeric>
#include <random>

namespace {

using namespace scream;

// Wrap a vector in a (unmanaged) host view
template<typename T>
typename GidDirectory<T>::gid_view_h view (const std::vector<T>& v) {
  return typename GidDirectory<T>::gid_view_h(v.data(),v.size());
}

// NOTE: we use CHECK (rather than REQUIRE), so that all ranks follow the same path in case of failure
//       on a subset of ranks (all the methods of GidDirectory are collective).

TEMPLATE_TEST_CASE ("gid_directory","",int,long long) {
  using gid_type = TestType;
  using Directory = GidDirectory<gid_type>;

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
    Directory dir (comm,view(my_gids));
    CHECK (not dir.has_duplicates());

    // Ask for random gids. Some are not registered, and some are repeated
    std::uniform_int_distribution<int> pdf (0,ngids+ngids/2+3);
    std::vector<gid_type> queries (20);
    for (auto& q : queries) {
      q = pdf(engine);
    }
    queries.push_back(queries[0]);
    queries.push_back(queries[1]);

    auto res = dir.lookup(view(queries));
    const auto& pids = res.pids;
    const auto& lids = res.lids;
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
    auto empty = dir.lookup(typename Directory::gid_view_h());
    CHECK (empty.pids.size()==0);
  }

  SECTION ("duplicate_gids") {
    // Every rank also registers the gid 5, which makes it be registered several times
    // (by different ranks if nranks>1, or twice by the same rank otherwise).
    auto dup_gids = my_gids;
    dup_gids.push_back(5);
    if (nranks==1) {
      dup_gids.push_back(5);
    }
    Directory dir (comm,view(dup_gids));
    CHECK (dir.has_duplicates());

    std::vector<gid_type> queries = {5};
    if (not my_gids.empty()) {
      queries.push_back(my_gids[0]);
    }
    auto res = dir.lookup(view(queries));
    CHECK (res.pids[0]==-2);
    CHECK (res.lids[0]==-2);
    if (not my_gids.empty()) {
      auto it = owners.at(my_gids[0]);
      CHECK (res.pids[1]==it.first);
      CHECK (res.lids[1]==it.second);
    }
  }

  SECTION ("negative_gids") {
    // Only one rank registers a negative gid, but all ranks must throw
    auto neg_gids = my_gids;
    if (rank==0) {
      neg_gids.push_back(-1);
    }
    CHECK_THROWS (Directory(comm,view(neg_gids)));

    // Negative gids can't be registered, so looking them up just reports "not found"
    Directory dir (comm,view(my_gids));
    std::vector<gid_type> queries = {-1,-7};
    auto res = dir.lookup(view(queries));
    for (size_t i=0; i<queries.size(); ++i) {
      CHECK (res.pids[i]==-1);
      CHECK (res.lids[i]==-1);
    }
  }
}

} // anonymous namespace
