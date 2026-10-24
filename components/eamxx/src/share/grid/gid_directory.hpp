#ifndef EAMXX_GID_DIRECTORY_HPP
#define EAMXX_GID_DIRECTORY_HPP

#include <ekat_comm.hpp>

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

class GidDirectory {
public:
  // NOTE: this must be the same as AbstractGrid::gid_type (checked in the cpp file).
  //       We don't include abstract_grid.hpp here, to be able to use the directory in that class.
  using gid_type = int;

  // Register the gids owned by this rank. The owner of gids[i] is (comm.rank(),i).
  GidDirectory (const ekat::Comm& comm,
                const gid_type* gids, const int num_gids);

  // Whether some gid was registered more than once (by the same or different ranks).
  // All ranks get the same answer.
  bool has_duplicates () const { return m_has_duplicates; }

  // For each input gid, return the rank that owns it, and its position in the owner's list.
  // The input gids may contain repetitions, and need not be owned by this rank.
  // For gids that are not registered, pids[i]=lids[i]=-1.
  // For gids that were registered more than once, pids[i]=lids[i]=-2.
  void lookup (const gid_type* gids, const int num_gids,
               std::vector<int>& pids, std::vector<int>& lids) const;

  // The rank where the owner of a gid is stored
  int directory_rank (const gid_type gid) const;

private:

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
