#include "share/grid/gid_directory.hpp"

#include "share/grid/abstract_grid.hpp"
#include "share/util/eamxx_utils.hpp"

#include <ekat_assert.hpp>

#include <algorithm>
#include <limits>
#include <numeric>
#include <type_traits>

namespace scream
{

static_assert (std::is_same<GidDirectory::gid_type,AbstractGrid::gid_type>::value,
               "Error! GidDirectory::gid_type must match AbstractGrid::gid_type.\n");

namespace {

// Exchange data among all ranks. The data in 'send' must be grouped by destination rank, with
// send_counts[p] entries for rank p. If recv_counts_known=false, recv_counts is computed
// (with a MPI_Alltoall). Otherwise, it must already contain the right counts.
template<typename T>
void alltoallv (const ekat::Comm& comm,
                const std::vector<T>& send, const std::vector<int>& send_counts,
                std::vector<T>& recv, std::vector<int>& recv_counts,
                const bool recv_counts_known)
{
  const int nranks = comm.size();
  const auto mpi_comm = comm.mpi_comm();

  if (not recv_counts_known) {
    recv_counts.resize(nranks);
    check_mpi_call(MPI_Alltoall(send_counts.data(),1,MPI_INT,recv_counts.data(),1,MPI_INT,mpi_comm),
                   "[GidDirectory] exchanging counts");
  }

  std::vector<int> send_displs(nranks+1,0), recv_displs(nranks+1,0);
  long long total_send = 0, total_recv = 0;
  for (int p=0; p<nranks; ++p) {
    send_displs[p] = total_send;
    recv_displs[p] = total_recv;
    total_send += send_counts[p];
    total_recv += recv_counts[p];
  }
  EKAT_REQUIRE_MSG (total_send<=std::numeric_limits<int>::max() and
                    total_recv<=std::numeric_limits<int>::max(),
      "Error! [GidDirectory] Too many entries to exchange with MPI_Alltoallv.\n");

  recv.resize(total_recv);
  const auto mpi_t = ekat::get_mpi_type<T>();
  check_mpi_call(MPI_Alltoallv(send.data(),send_counts.data(),send_displs.data(),mpi_t,
                               recv.data(),recv_counts.data(),recv_displs.data(),mpi_t,mpi_comm),
                 "[GidDirectory] exchanging data");
}

} // anonymous namespace

int GidDirectory::directory_rank (const gid_type gid) const
{
  // Gids are (nearly) contiguous in all our grids, so a simple modulo balances the load
  return gid % m_comm.size();
}

GidDirectory::
GidDirectory (const ekat::Comm& comm,
              const gid_type* gids, const int num_gids)
 : m_comm (comm)
{
  const int nranks = comm.size();

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

void GidDirectory::
lookup (const gid_type* gids, const int num_gids,
        std::vector<int>& pids, std::vector<int>& lids) const
{
  const int nranks = m_comm.size();

  pids.assign(num_gids,-1);
  lids.assign(num_gids,-1);

  // Ask only once for each gid
  std::vector<gid_type> queries(gids,gids+num_gids);
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
    const int iq = std::lower_bound(queries.begin(),queries.end(),gids[i]) - queries.begin();
    pids[i] = query_pids[iq];
    lids[i] = query_lids[iq];
  }
}

} // namespace scream
