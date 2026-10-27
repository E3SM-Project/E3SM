#ifndef EAMXX_GID_DIRECTORY_HPP
#define EAMXX_GID_DIRECTORY_HPP

#include <ekat_comm.hpp>
#include <ekat_kokkos_types.hpp>

#include <vector>

namespace scream
{

/*
 * A distributed directory of gids
 *
 * Each rank registers the gids it owns (a gid is identified by the rank that owns
 * it, and its position in that rank's list). Then, any rank can ask who owns a
 * given list of gids, and where (at which position in the owner's list).
 *
 * The directory is spread among ranks: gid%P is the "directory rank", which is the
 * one that stores the owner of that gid. (Gids in our grids are contiguous, so the load is balanced.) Therefore
 *  - each rank stores (and sends/receives) O(N/P) entries, where N is the global
 *    number of gids, and P the number of ranks;
 *  - both the construction and each lookup use a constant number of collective
 *    calls (MPI_Alltoall(v)), rather than P broadcasts.
 *
 * All the methods in this class are collective.
 */

template<typename GidType>
class GidDirectory {
public:
  using gid_type   = GidType;
  using gid_view_h = ekat::KokkosTypes<ekat::HostDevice>::view_1d<const gid_type>;

  // The owners found by a lookup. For the i-th input gid, the owner is (pids[i],lids[i]).
  //  - For gids that are not registered (this includes negative gids), pids[i]=lids[i]=-1.
  //  - For gids that were registered more than once, pids[i]=lids[i]=-2.
  struct LookupResult {
    std::vector<int> pids;
    std::vector<int> lids;
  };

  // Register the gids owned by this rank. The owner of gids[i] is (comm.rank(),i).
  // Gids must be non-negative (an error is thrown on all ranks otherwise).
  GidDirectory (const ekat::Comm& comm, const gid_view_h& gids);

  // Whether some gid was registered more than once (by the same or different ranks).
  // All ranks get the same answer.
  bool has_duplicates () const { return m_has_duplicates; }

  // For each input gid, find the rank that owns it, and its position in the owner's list.
  // The input gids may contain repetitions, and need not be owned by this rank.
  LookupResult lookup (const gid_view_h& gids) const;

private:

  // The rank where the owner of a gid is stored
  int directory_rank (const gid_type gid) const;

  struct Entry {
    gid_type gid;
    int      pid;
    int      lid;
  };

  ekat::Comm          m_comm;
  std::vector<Entry>  m_entries;   // Sorted by gid
  bool                m_has_duplicates = false;
};

} // namespace scream

#endif // EAMXX_GID_DIRECTORY_HPP
