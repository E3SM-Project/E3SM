#ifndef EAMXX_ALLTOALLV_HPP
#define EAMXX_ALLTOALLV_HPP

#include "share/util/eamxx_utils.hpp"   // For check_mpi_call

#include <ekat_comm.hpp>

#include <mpi.h>

#include <limits>
#include <vector>

namespace scream
{

// Exchange data among all the ranks of a comm, with MPI_Alltoallv.
//  - send: the data to send, grouped by destination rank (send_counts[p] entries for rank p)
//  - recv: the data received, grouped by source rank (recv_counts[p] entries from rank p)
//  - recv_counts_known: if false, recv_counts is computed with an MPI_Alltoall. Otherwise,
//    the caller already knows how many entries to receive from each rank, and stored that in recv_counts
//  - mpi_t: the MPI datatype of T. For basic types, it can be omitted.
// This is a collective call.
template<typename T>
void alltoallv (const ekat::Comm& comm,
                const std::vector<T>& send, const std::vector<int>& send_counts,
                std::vector<T>& recv, std::vector<int>& recv_counts,
                const bool recv_counts_known,
                const MPI_Datatype mpi_t = ekat::get_mpi_type<T>())
{
  const int nranks = comm.size();
  const auto mpi_comm = comm.mpi_comm();

  if (not recv_counts_known) {
    recv_counts.resize(nranks);
    check_mpi_call(MPI_Alltoall(send_counts.data(),1,MPI_INT,recv_counts.data(),1,MPI_INT,mpi_comm),
                   "[alltoallv] exchanging counts");
  }

  std::vector<int> send_displs(nranks,0), recv_displs(nranks,0);
  long long total_send = 0, total_recv = 0;
  for (int p=0; p<nranks; ++p) {
    send_displs[p] = total_send;
    recv_displs[p] = total_recv;
    total_send += send_counts[p];
    total_recv += recv_counts[p];
    EKAT_REQUIRE_MSG (total_send<=std::numeric_limits<int>::max() and
                      total_recv<=std::numeric_limits<int>::max(),
        "Error! [alltoallv] Too many entries to exchange with MPI_Alltoallv.\n");
  }

  recv.resize(total_recv);
  check_mpi_call(MPI_Alltoallv(send.data(),send_counts.data(),send_displs.data(),mpi_t,
                               recv.data(),recv_counts.data(),recv_displs.data(),mpi_t,mpi_comm),
                 "[alltoallv] exchanging data");
}

} // namespace scream

#endif // EAMXX_ALLTOALLV_HPP
