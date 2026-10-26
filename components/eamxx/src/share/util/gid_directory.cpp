#include "share/util/gid_directory.hpp"

#include "share/util/eamxx_alltoallv.hpp"

#include <ekat_assert.hpp>

#include <algorithm>
#include <iterator>
#include <limits>
#include <numeric>

namespace scream
{

template<typename GidType>
int GidDirectory<GidType>::directory_rank (const gid_type gid) const
{
  // Gids are (nearly) contiguous in all our grids, so a simple modulo balances the load
  return gid % m_comm.size();
}

template<typename GidType>
GidDirectory<GidType>::
GidDirectory (const ekat::Comm& comm, const gid_view_h& gids)
 : m_comm (comm)
{
  const int nranks = comm.size();
  const int num_gids = gids.size();

  // Gids must be non-negative (directory_rank relies on it). Check it collectively,
  // so that either all ranks throw or none does (rather than some ranks hanging in MPI calls)
  int local_neg = 0;
  for (int i=0; i<num_gids; ++i) {
    if (gids[i]<0) {
      local_neg = 1;
      break;
    }
  }
  int global_neg;
  comm.all_reduce(&local_neg,&global_neg,1,MPI_MAX);
  EKAT_REQUIRE_MSG (global_neg==0,
      "Error! GidDirectory can only register non-negative gids.\n");

  // Group the registrations by directory rank
  std::vector<int> send_counts(nranks,0);
  std::vector<int> dest(num_gids);
  for (int i=0; i<num_gids; ++i) {
    dest[i] = directory_rank(gids[i]);
    ++send_counts[dest[i]];
  }
  std::vector<int> cursor(nranks,0);
  for (int p=1; p<nranks; ++p) {
    cursor[p] = cursor[p-1] + send_counts[p-1];
  }
  std::vector<gid_type> send_gids(num_gids);
  std::vector<int>      send_lids(num_gids);
  for (int i=0; i<num_gids; ++i) {
    const int pos = cursor[dest[i]]++;
    send_gids[pos] = gids[i];
    send_lids[pos] = i;
  }

  std::vector<gid_type> recv_gids;
  std::vector<int>      recv_lids, recv_counts;
  alltoallv(comm,send_gids,send_counts,recv_gids,recv_counts,false);
  alltoallv(comm,send_lids,send_counts,recv_lids,recv_counts,true);

  // Store what we received, sorted by gid. The owner of each entry is the rank that sent it.
  m_entries.reserve(recv_gids.size());
  int pos = 0;
  for (int src=0; src<nranks; ++src) {
    for (int k=0; k<recv_counts[src]; ++k, ++pos) {
      m_entries.push_back({recv_gids[pos],src,recv_lids[pos]});
    }
  }
  std::sort(m_entries.begin(),m_entries.end(),
            [](const Entry& a, const Entry& b) {
              return a.gid<b.gid or (a.gid==b.gid and (a.pid<b.pid or (a.pid==b.pid and a.lid<b.lid)));
            });

  int local_dups = 0;
  for (size_t i=1; i<m_entries.size(); ++i) {
    if (m_entries[i].gid==m_entries[i-1].gid) {
      local_dups = 1;
      break;
    }
  }
  int global_dups;
  comm.all_reduce(&local_dups,&global_dups,1,MPI_MAX);
  m_has_duplicates = global_dups==1;
}

template<typename GidType>
typename GidDirectory<GidType>::LookupResult
GidDirectory<GidType>::
lookup (const gid_view_h& gids) const
{
  const int nranks = m_comm.size();
  const int num_gids = gids.size();

  LookupResult result;
  result.pids.assign(num_gids,-1);
  result.lids.assign(num_gids,-1);

  // Ask only once for each gid. Negative gids can't be registered, so there's no need to ask.
  std::vector<gid_type> queries;
  queries.reserve(num_gids);
  for (int i=0; i<num_gids; ++i) {
    if (gids[i]>=0) {
      queries.push_back(gids[i]);
    }
  }
  std::sort(queries.begin(),queries.end());
  queries.erase(std::unique(queries.begin(),queries.end()),queries.end());
  const int nq = queries.size();

  // Group the queries by directory rank. Keep track of the order, to be able to map the answers back
  std::vector<int> send_counts(nranks,0);
  std::vector<int> dest(nq);
  for (int i=0; i<nq; ++i) {
    dest[i] = directory_rank(queries[i]);
    ++send_counts[dest[i]];
  }
  std::vector<int> cursor(nranks,0);
  for (int p=1; p<nranks; ++p) {
    cursor[p] = cursor[p-1] + send_counts[p-1];
  }
  std::vector<gid_type> send_queries(nq);
  std::vector<int>      query_idx(nq);   // For each query in the send buffer, its index in 'queries'
  for (int i=0; i<nq; ++i) {
    const int pos = cursor[dest[i]]++;
    send_queries[pos] = queries[i];
    query_idx[pos] = i;
  }

  std::vector<gid_type> recv_queries;
  std::vector<int>      recv_counts;
  alltoallv(m_comm,send_queries,send_counts,recv_queries,recv_counts,false);

  // Answer the queries we received
  const int nrecv = recv_queries.size();
  std::vector<int> ans_pids(nrecv), ans_lids(nrecv);
  for (int j=0; j<nrecv; ++j) {
    const auto gid = recv_queries[j];
    auto it = std::lower_bound(m_entries.begin(),m_entries.end(),gid,
                               [](const Entry& e, const gid_type g) { return e.gid<g; });
    if (it==m_entries.end() or it->gid!=gid) {
      ans_pids[j] = ans_lids[j] = -1;
    } else if ((it+1)!=m_entries.end() and (it+1)->gid==gid) {
      ans_pids[j] = ans_lids[j] = -2;
    } else {
      ans_pids[j] = it->pid;
      ans_lids[j] = it->lid;
    }
  }

  // Send the answers back (in the same order as we got the queries)
  std::vector<int> back_pids, back_lids;
  alltoallv(m_comm,ans_pids,recv_counts,back_pids,send_counts,true);
  alltoallv(m_comm,ans_lids,recv_counts,back_lids,send_counts,true);

  std::vector<int> query_pids(nq), query_lids(nq);
  for (int pos=0; pos<nq; ++pos) {
    query_pids[query_idx[pos]] = back_pids[pos];
    query_lids[query_idx[pos]] = back_lids[pos];
  }

  for (int i=0; i<num_gids; ++i) {
    if (gids[i]<0) {
      continue;
    }
    const int iq = std::lower_bound(queries.begin(),queries.end(),gids[i]) - queries.begin();
    result.pids[i] = query_pids[iq];
    result.lids[i] = query_lids[iq];
  }

  return result;
}

// Explicit instantiations, for the gid types we may use
template class GidDirectory<int>;
template class GidDirectory<long long>;

} // namespace scream
