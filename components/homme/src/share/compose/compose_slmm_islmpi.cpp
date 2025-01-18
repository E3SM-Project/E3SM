#include "compose.hpp"
#include "compose_slmm_islmpi.hpp"
#include "compose_slmm_departure_point.hpp"

#include <set>
#include <map>
#include <memory>

namespace homme {
namespace mpi {
#ifdef COMPOSE_DEBUG_MPI
Request::Request () : unfreed(0) {}
Request::~Request () {
  if (unfreed) {
    std::stringstream ss;
    ss << "Request is being deleted with unfreed = " << unfreed;
    int fin;
    MPI_Finalized(&fin);
    if (fin) {
      ss << "\n";
      std::cerr << ss.str();
    } else {
      pr(ss.str());
    }
  }
}
#endif

int waitany (int count, Request* reqs, int* index, MPI_Status* stats) {
#ifdef COMPOSE_DEBUG_MPI
  std::vector<MPI_Request> vreqs(count);
  for (int i = 0; i < count; ++i) vreqs[i] = reqs[i].request;
  const auto out = MPI_Waitany(count, vreqs.data(), index,
                               stats ? stats : MPI_STATUS_IGNORE);
  for (int i = 0; i < count; ++i) reqs[i].request = vreqs[i];
  reqs[*index].unfreed--;
  return out;
#else
  return MPI_Waitany(count, reinterpret_cast<MPI_Request*>(reqs), index,
                     stats ? stats : MPI_STATUS_IGNORE);
#endif
}

int waitall (int count, Request* reqs, MPI_Status* stats) {
#ifdef COMPOSE_DEBUG_MPI
  std::vector<MPI_Request> vreqs(count);
  for (int i = 0; i < count; ++i) vreqs[i] = reqs[i].request;
  const auto out = MPI_Waitall(count, vreqs.data(),
                               stats ? stats : MPI_STATUS_IGNORE);
  for (int i = 0; i < count; ++i) {
    reqs[i].request = vreqs[i];
    reqs[i].unfreed--;
  }
  return out;
#else
  return MPI_Waitall(count, reinterpret_cast<MPI_Request*>(reqs),
                     stats ? stats : MPI_STATUS_IGNORE);
#endif
}

int wait (Request* req, MPI_Status* stat) {
#ifdef COMPOSE_DEBUG_MPI
  const auto out = MPI_Wait(&req->request, stat ? stat : MPI_STATUS_IGNORE);
  req->unfreed--;
  return out;
#else
  return MPI_Wait(reinterpret_cast<MPI_Request*>(req), stat ? stat : MPI_STATUS_IGNORE);
#endif
}
} // namespace mpi

namespace islmpi {
namespace extend_halo {
// Extend halo by one layer. This has two parts: finding neighbor (gid, rank) in
// collect_gid_rank, and extending the Advecter local mesh geometry in
// extend_local_meshes. The two parts have the same comm pattern: in round 1,
// request data for lists of GIDs; in round 2, fulfill these requests.

typedef Int Gid;
typedef Int Rank;
struct GidRankPair {
  const Gid gid;
  const Rank rank;
  GidRankPair(const Gid& gid_, const Rank& rank_) : gid(gid_), rank(rank_) {}
  bool operator< (const GidRankPair& o) const {
    if (gid < o.gid) return true;
    if (gid > o.gid) return false;
    return rank < o.rank;
  }
};
typedef std::vector<GidRankPair> GidRankPairs;
typedef std::map<Gid, GidRankPairs> Gid2Nbrs;
typedef std::vector<Int> IntBuf;
typedef std::vector<Real> RealBuf;
typedef std::map<Gid, int> Gid2Count;

template <typename MT>
GidRankPairs all_1halo_nbrs_but_me (const typename IslMpi<MT>::ElemDataH& ed) {
  GidRankPairs gs;
  gs.reserve(ed.nin1halo - 1);
  for (Int i = 0; i < ed.nin1halo; ++i) {
    const auto& n = ed.nbrs(i);
    if (&n != ed.me)
      gs.push_back(GidRankPair(n.gid, n.rank));
  }
  return gs;
}

template <typename MT>
void fill_gid2nbrs (const mpi::Parallel& p, const typename IslMpi<MT>::ElemDataListH& eds,
                    Gid2Nbrs& gid2nbrs, Gid2Count& gid2ninprevhalo) {
  static const Int tag = 6;
  const Rank my_rank = p.rank();
  const Int n_owned = eds.size();

  // Fill in the ones we know.
  for (const auto& ed : eds) {
    slmm_assert(ed.me->rank == my_rank);
    gid2nbrs[ed.me->gid] = all_1halo_nbrs_but_me<MT>(ed);
  }

  std::vector<Rank> ranks;
  Int nrank;
  std::vector<IntBuf> req_sends, req_recvs;
  std::vector<mpi::Request> req_recv_reqs;
  {
    // Find the ranks that know the rest.
    std::map<Rank,Int> rank2rankidx;
    std::map<Gid,Rank> needgid2rank;
    {
      std::set<Rank> unique_ranks;
      for (const auto& ed : eds) {
        // We only need information for GIDs in the current outermost halo.
        const auto& it = gid2ninprevhalo.find(ed.me->gid);
        const Int i0 = it == gid2ninprevhalo.end() ? 0 : it->second;
        for (Int i = i0; i < ed.nbrs.size(); ++i) {
          const auto& n = ed.nbrs(i);
          if (n.rank == my_rank) continue;
          slmm_assert(gid2nbrs.find(n.gid) == gid2nbrs.end());
          needgid2rank.insert(std::make_pair(n.gid, n.rank));
          unique_ranks.insert(n.rank);
        }
      }
      nrank = unique_ranks.size();
      ranks.insert(ranks.begin(), unique_ranks.begin(), unique_ranks.end());
      Int i = 0;
      for (const auto& rank : ranks)
        rank2rankidx[rank] = i++;
    }

    // Send and receive neighbor queries.
    slmm_assert(ranks.size() == static_cast<size_t>(nrank));
    req_sends.resize(nrank);
    req_recvs.resize(nrank); req_recv_reqs.resize(nrank);
    for (Int i = 0; i < nrank; ++i) {
      auto& r = req_recvs[i];
      r.resize(n_owned+1); // upper bound on #requests, plus 1 for size datum
      mpi::irecv(p, r.data(), r.size(), ranks[i], tag, &req_recv_reqs[i]);
      req_sends[i].push_back(0);
    }
    for (const auto& item : needgid2rank) {
      const auto& gid = item.first;
      const auto& rank = item.second;
      const auto& it = rank2rankidx.find(rank);
      slmm_assert(it != rank2rankidx.end());
      const auto& rank_idx = it->second;
      slmm_assert(gid2nbrs.find(gid) == gid2nbrs.end());
      req_sends[rank_idx].push_back(gid);
    }
    for (Int i = 0; i < nrank; ++i) {
      auto& s = req_sends[i];
      s[0] = s.size() - 1; // #gids in request
      mpi::isend(p, s.data(), s.size(), ranks[i], tag, nullptr);
    }
  }

  // Fullfill queries and receive answers to our queries.
  std::vector<IntBuf> nbr_sends(nrank), nbr_recvs(nrank);
  std::vector<mpi::Request> nbr_send_reqs(nrank), nbr_recv_reqs(nrank);
  for (Int i = 0; i < nrank; ++i) {
    auto& r = nbr_recvs[i];
    // 20 is from dimensions_mod::set_mesh_dimensions, the maximum size of the
    // 1-halo minus the 0-halo; factor of 2 is to get (gid,rank); 1 is for the
    // size datum.
    r.resize((20*2 + 1)*(req_sends[i].size() - 1));
    mpi::irecv(p, r.data(), r.size(), ranks[i], tag, &nbr_recv_reqs[i]);
  }
  for (Int k = 0; k < nrank; ++k) {
    Int i;
    mpi::waitany(nrank, req_recv_reqs.data(), &i);
    const auto& r = req_recvs[i];
    auto& s = nbr_sends[i];
    const Int ngid = r[0];
    for (Int j = 0; j < ngid; ++j) {
      const auto gid = r[j+1];
      const auto it = gid2nbrs.find(gid);
      slmm_assert(it != gid2nbrs.end());
      const auto& nbrs = it->second;
      s.push_back(nbrs.size());
      for (const auto& n : nbrs) {
        s.push_back(n.gid);
        s.push_back(n.rank);
      }
    }
    mpi::isend(p, s.data(), s.size(), ranks[i], tag, &nbr_send_reqs[i]);
  }
  for (Int j = 0; j < nrank; ++j) {
    Int i;
    mpi::waitany(nrank, nbr_recv_reqs.data(), &i);
    const auto& s = req_sends[i];
    const Int nrequested = s.size() - 1;
    const auto& r = nbr_recvs[i];
    for (Int si = 0, ri = 0; si < nrequested; ++si) {
      const Gid gid = s[si+1];
      const Int nnbr = r[ri++];
      GidRankPairs nbrs;
      for (Int ni = 0; ni < nnbr; ++ni, ri += 2)
        nbrs.push_back(GidRankPair(r[ri], r[ri+1]));
      slmm_assert(gid2nbrs.find(gid) == gid2nbrs.end());
      gid2nbrs.insert(std::make_pair(gid, nbrs));
    }
  }
  mpi::waitall(nbr_send_reqs.size(), nbr_send_reqs.data());
}

template <typename MT>
void extend_nbrs (const Gid2Nbrs& gid2nbrs, typename IslMpi<MT>::ElemDataListH& eds,
                  Gid2Count& gid2ninprevhalo) {
  for (auto& ed : eds) {
    // Get all <=(n+1)-halo neighbors, where we already have <=n-halo neighbors.
    std::set<GidRankPair> new_nbrs; {
      const auto& it = gid2ninprevhalo.find(ed.me->gid);
      const Int i0 = it == gid2ninprevhalo.end() ? 0 : it->second;
      for (Int i = i0; i < ed.nbrs.size(); ++i) {
        const auto& n = ed.nbrs(i);
        if (&n == ed.me) continue;
        const auto& it = gid2nbrs.find(n.gid);
        slmm_assert(it != gid2nbrs.end());
        const auto& gid_nbrs = it->second;
        for (const auto& gn : gid_nbrs)
          new_nbrs.insert(gn);
      }
    }
    // Remove the already known ones.
    for (const auto& n : ed.nbrs)
      new_nbrs.erase(GidRankPair(n.gid, n.rank));
    // Find me b/c of the pointer reset.
    Int me = -1;
    for (Int i = 0; i < ed.nbrs.size(); ++i)
      if (&ed.nbrs(i) == ed.me) {
        me = i;
        break;
      }
    slmm_assert(me >= 0);
    // Append the, now only new, (n+1)-halo ones.
    Int i = ed.nbrs.size();
    gid2ninprevhalo[ed.me->gid] = i;
    ed.nbrs.reset_capacity(i + new_nbrs.size(), true);
    ed.me = &ed.nbrs(me);
    for (const auto& n : new_nbrs) {
      auto& en = ed.nbrs(i++);
      en.gid = n.gid;
      en.rank = n.rank;
      en.rank_idx = -1;
      en.lid_on_rank = -1;
      en.lid_on_rank_idx = -1;
    }
#ifndef NDEBUG
    {
      std::set<Gid> ugid;
      for (Int i = 0; i < ed.nbrs.size(); ++i) ugid.insert(ed.nbrs(i).gid);
      slmm_assert(ugid.size() == size_t(ed.nbrs.size()));
    }
#endif
  }
}

template <typename MT>
void collect_gid_rank (const mpi::Parallel& p, typename IslMpi<MT>::ElemDataListH& eds,
                       Gid2Count& gid2ninprevhalo) {
  Gid2Nbrs gid2nbrs;
  fill_gid2nbrs<MT>(p, eds, gid2nbrs, gid2ninprevhalo);
  extend_nbrs<MT>(gid2nbrs, eds, gid2ninprevhalo);
}

template <typename MT>
void extend_local_meshes (const mpi::Parallel& p,
                          const typename IslMpi<MT>::ElemDataListH& eds,
                          typename IslMpi<MT>::Advecter& advecter) {
  using slmm::slice;
  using slmm::nslices;
  using slmm::szslice;
  using slmm::len;
  using RealArray3 = ko::View<Real***, typename MT::HES>;

  static const Int tag = 24;
  const Int my_rank = p.rank();
  const Int n_owned = eds.size();

  RealArray3 corners;
  std::map<Gid,Int> owngid2lid, rmtgid2idx;
  {
    // Find relevant ranks.
    std::vector<Int> ranks;
    Int nrank;
    std::vector<IntBuf> req_sends, req_recvs;
    std::vector<mpi::Request> req_recv_reqs;
    {
      std::map<Int,Int> rank2rankidx;
      {
        std::set<Int> uranks;
        for (const auto& ed : eds)
          for (const auto& n : ed.nbrs)
            uranks.insert(n.rank);
        uranks.erase(my_rank);
        ranks.insert(ranks.begin(), uranks.begin(), uranks.end());
        nrank = ranks.size();
        Int i = 0;
        for (const auto& rank : ranks)
          rank2rankidx[rank] = i++;
      }

      // Trade requests.
      req_sends.resize(nrank); req_recvs.resize(nrank);
      req_recv_reqs.resize(nrank);
      for (Int i = 0; i < nrank; ++i) { // Set up recvs.
        auto& r = req_recvs[i];
        r.resize(n_owned+1);
        mpi::irecv(p, r.data(), r.size(), ranks[i], tag, &req_recv_reqs[i]);
      }
      std::set<Gid> unique_rmtgids;
      { // Collect the remote GIDs by rank.
        std::vector<std::set<Gid> > req_gids(nrank);
        for (const auto& ed : eds)
          for (Int i = ed.nin1halo; i < ed.nbrs.size(); ++i) {
            const auto& n = ed.nbrs(i);
            if (n.rank != my_rank) {
              req_gids[rank2rankidx[n.rank]].insert(n.gid);
              unique_rmtgids.insert(n.gid);
            }
          }
        for (Int i = 0; i < nrank; ++i) { // Set up sends.
          auto& s = req_sends[i];
          s.push_back(0);
          s.insert(s.end(), req_gids[i].begin(), req_gids[i].end());
          s[0] = s.size() - 1;
          mpi::isend(p, s.data(), s.size(), ranks[i], tag, nullptr);
        }
        Int i = 0;
        for (const auto& gid : unique_rmtgids)
          rmtgid2idx[gid] = i++;
      }
    }

    for (Int i = 0; i < eds.size(); ++i) {
      const Gid gid = eds(i).me->gid;
      owngid2lid[gid] = i;
    }

    // Fulfill our requests and obtain the data in reply to ours.
    std::vector<RealBuf> geo_recvs(nrank), geo_sends(nrank);
    std::vector<mpi::Request> geo_recv_reqs(nrank), geo_send_reqs(nrank);
    for (Int i = 0; i < nrank; ++i) {
      auto& r = geo_recvs[i];
      r.resize(12*(req_sends[i].size() - 1)); // 4 vertices x 3 dimensions
      mpi::irecv(p, r.data(), r.size(), ranks[i], tag, &geo_recv_reqs[i]);
    }
    for (Int k = 0; k < nrank; ++k) { // Pack cell corner points for requesters.
      Int i;
      mpi::waitany(nrank, req_recv_reqs.data(), &i);
      const auto& r = req_recvs[i];
      auto& s = geo_sends[i];
      const Int ngid = r[0];
      s.resize(12*ngid);
      Int si = 0;
      for (Int j = 0; j < ngid; ++j) {
        const auto gid = r[j+1];
        const auto it = owngid2lid.find(gid);
        slmm_assert(it != owngid2lid.end());
        const auto lid = it->second;
        const auto& local_mesh = advecter.local_mesh_host(lid);
        const auto cell = slice(local_mesh.e, local_mesh.tgt_elem);
        slmm_assert(szslice(local_mesh.e) == 4 && cell[3] != -1);
        const auto& p = local_mesh.p;
        for (Int v = 0; v < 4; ++v)
          for (Int d = 0; d < 3; ++d)
            s[si++] = p(cell[v],d);
      }
      slmm_assert(si == static_cast<Int>(s.size()));
      mpi::isend(p, s.data(), s.size(), ranks[i], tag, &geo_send_reqs[i]);
    }
    corners = RealArray3("corners", rmtgid2idx.size(), 4, 3);
    for (Int j = 0; j < nrank; ++j) { // Pack cell corner points for me.
      Int i;
      mpi::waitany(nrank, geo_recv_reqs.data(), &i);
      const auto& s = req_sends[i];
      const Int nrequested = s.size() - 1;
      const auto& r = geo_recvs[i];
      for (Int si = 0, ri = 0; si < nrequested; ++si) {
        const Gid gid = s[si+1];
        const auto it = rmtgid2idx.find(gid);
        slmm_assert(it != rmtgid2idx.end());
        const auto idx = it->second;
        for (Int v = 0; v < 4; ++v)
          for (Int d = 0; d < 3; ++d)
            corners(idx,v,d) = r[ri++];
      }
    }
    // Do this here so we can release the memory.
    mpi::waitall(geo_send_reqs.size(), geo_send_reqs.data());
  }

  // Extend the local meshes.
  for (Int lid = 0; lid < eds.size(); ++lid) {
    const auto& ed = eds(lid);
    auto& local_mesh = advecter.local_mesh_host(lid);
    auto p0 = local_mesh.p;
    auto e0 = local_mesh.e;
    const Int ncell = ed.nbrs.size();
    slmm_assert(szslice(p0) == 3 && szslice(e0) == 4);
    local_mesh.p = typename slmm::LocalMesh<typename MT::HES>::RealArray("p", 4*ncell);
    auto& p = local_mesh.p;
    local_mesh.e = typename slmm::LocalMesh<typename MT::HES>::IntArray("e", ncell, 4);
    auto& e = local_mesh.e;
    // Copy in old data.
    for (Int pi = 0; pi < nslices(p0); ++pi)
      for (Int d = 0; d < 3; ++d)
        p(pi,d) = p0(pi,d);
    for (Int ei = 0; ei < nslices(e0); ++ei)
      for (Int d = 0; d < 4; ++d)
        e(ei,d) = e0(ei,d);
    // Fill in new data.
    for (Int ni = ed.nin1halo; ni < ed.nbrs.size(); ++ni) {
      const auto& n = ed.nbrs(ni);
      const Gid gid = n.gid;
      if (n.rank != my_rank) {
        const auto it = rmtgid2idx.find(gid);
        slmm_assert(it != rmtgid2idx.end());
        const Int idx = it->second;
        for (Int v = 0; v < 4; ++v) {
          const Int slot = 4*ni+v;
          e(ni,v) = slot;
          slmm_assert(slot < nslices(p));
          for (Int d = 0; d < 3; ++d) p(slot,d) = corners(idx,v,d);
        }
      } else {
        const auto it = owngid2lid.find(gid);
        slmm_assert(it != owngid2lid.end());
        const Int lid = it->second;
        const auto& local_mesh_other = advecter.local_mesh_host(lid);
        const auto cell_other = slice(local_mesh_other.e, local_mesh_other.tgt_elem);
        const auto& p_other = local_mesh_other.p;
        for (Int v = 0; v < 4; ++v) {
          const Int slot = 4*ni+v;
          e(ni,v) = slot;
          slmm_assert(slot < nslices(p));
          for (Int d = 0; d < 3; ++d) p(slot,d) = p_other(cell_other[v],d);
        }
      }
    }
    // Recompute all normals.
    slmm::fill_normals(local_mesh);
  }
}

template void
extend_local_meshes<ko::MachineTraits>(
  const mpi::Parallel& p,
  const typename IslMpi<ko::MachineTraits>::ElemDataListH& eds,
  IslMpi<ko::MachineTraits>::Advecter& advecter);

} // namespace extend_halo

// Fill in (gid, rank), the list of owning rank per gid.
template <typename MT>
void collect_gid_rank (IslMpi<MT>& cm, const Int* nbr_id_rank, const Int* nirptr) {
  cm.ed_h.reset_capacity(cm.nelemd, true);
  for (Int i = 0; i < cm.nelemd; ++i) {
    auto& ed = cm.ed_h(i);
    const Int* nir = nbr_id_rank + nirptr[i];
    const Int nnir = (nirptr[i+1] - nirptr[i]) / 2;
    const Int mygid = nir[0];
    ed.me = nullptr;
    ed.nin1halo = nnir-1;
    ed.nbrs.reset_capacity(ed.nin1halo, true);
#ifndef COMPOSE_PORT
    ed.own.reset_capacity(cm.nlev * cm.np2);
#endif
    ed.rmt.reset_capacity(cm.nlev * cm.np2);
    for (Int j = 1; j < nnir; ++j) {
      auto& n = ed.nbrs(j-1);
      n.gid = nir[2*j];
      if (n.gid == mygid) {
        slmm_assert( ! ed.me);
        ed.me = &n;
      }
      n.rank = nir[2*j+1];
      n.rank_idx = -1;
      n.lid_on_rank = -1;
      n.lid_on_rank_idx = -1;
    }
    slmm_assert(ed.me);
  }
  if (cm.halo > 1) {
    extend_halo::Gid2Count gid2ninprevhalo;
    for (int halo = 2; halo <= cm.halo; ++halo)
      extend_halo::collect_gid_rank<MT>(*cm.p, cm.ed_h, gid2ninprevhalo);
  }
#ifdef COMPOSE_PORT
  cm.own_dep_mask = typename IslMpi<MT>::DepMask("own_dep_mask",
                                                 cm.nelemd, cm.nlev, cm.np2);
  cm.own_dep_list = typename IslMpi<MT>::DepList("own_dep_list",
                                                 cm.nelemd*cm.nlev*cm.np2);
#endif
}

typedef std::map<Int, std::set<Int> > Rank2Gids;

template <typename MT>
void get_rank2gids (const IslMpi<MT>& cm, Rank2Gids& rank2rmtgids,
                    Rank2Gids& rank2owngids) {
  const Int myrank = cm.p->rank();
  for (Int i = 0; i < cm.nelemd; ++i) {
    const auto& ed = cm.ed_h(i);
    for (const auto& n: ed.nbrs) {
      if (n.rank == myrank) continue;
      // I need this rmt gid's lid.
      rank2rmtgids[n.rank].insert(n.gid);
      // This rank needs this gid's lid.
      rank2owngids[n.rank].insert(ed.me->gid);
    }
  }
}

// Fill in nbrs.lid_on_rank, the lid on the remote rank corresponding to the gid
// I share but do not own.
template <typename MT>
void comm_lid_on_rank (IslMpi<MT>& cm, const Rank2Gids& rank2rmtgids,
                       const Rank2Gids& rank2owngids,
                       std::map<Int, Int>& gid2rmt_owning_lid) {
  const Int myrank = cm.p->rank();

  std::map<Int, Int> gid2mylid;
  for (Int i = 0; i < cm.nelemd; ++i)
    gid2mylid[cm.ed_h(i).me->gid] = i;
  
  // Set up to recv remote (gid, lid) lists.
  Int rn = 0;
  for (const auto& e: rank2rmtgids)
    rn += e.second.size();
  const Int nrecv = rank2rmtgids.size();
  std::vector<Int> recv(rn), recvptr(nrecv+1), recvrank(nrecv);
  std::vector<mpi::Request> reqs(nrecv);
  recvptr[0] = 0;
  Int ir = 0;
  rn = 0;
  for (const auto& e: rank2rmtgids) {
    const Int ne = e.second.size();
    const Int rank = e.first;
    slmm_assert(rank != myrank);
    recvrank[ir] = rank;
    mpi::irecv(*cm.p, recv.data() + rn, ne, rank, 42, &reqs[ir++]);
    rn += ne;
    recvptr[ir] = rn;
  }

  // Send my (gid, lid) lists.
  Int sn = 0, i = 0;
  for (const auto& e: rank2owngids)
    sn += e.second.size();
  std::vector<Int> send(sn);
  std::vector<mpi::Request> sendreqs(rank2owngids.size());
  sn = 0;
  for (const auto& e: rank2owngids) {
    // Iteration through a set gives increasing GID, which means the LID list is
    // consistent between communicating ranks.
    Int slot = sn, pgid = -1;
    for (auto gid: e.second) {
      slmm_assert(gid > pgid);
      pgid = gid;
      send[slot++] = gid2mylid[pgid];
    }
    const Int ne = e.second.size();
    const Int rank = e.first;
    mpi::isend(*cm.p, send.data() + sn, ne, rank, 42, &sendreqs[i++]);
    sn += ne;
  }

  // Actually recv.
  const Int count = rank2owngids.size();
  for (Int c = 0; c < count; ++c) {
    int idx;
    mpi::waitany(count, reqs.data(), &idx);
    const Int* rmtlids = recv.data() + recvptr[idx];
    const auto& gids = rank2rmtgids.at(recvrank[idx]);
    slmm_assert(recvptr[idx+1] - recvptr[idx] == static_cast<int>(gids.size()));
    Int i = 0;
    for (auto gid: gids)
      gid2rmt_owning_lid[gid] = rmtlids[i++];
  }

  // Fill lid_on_rank and mylid_with_comm.
  std::vector<Int> mylid_with_comm;
  for (Int i = 0; i < cm.nelemd; ++i) {
    auto& ed = cm.ed_h(i);
    bool has_comm = false;
    for (auto& n: ed.nbrs)
      if (n.rank == myrank) {
        const auto it = gid2mylid.find(n.gid);
        slmm_throw_if(it == gid2mylid.end(),
                      "comm_lid_on_rank: On rank " << myrank << ", gid " << n.gid
                      << " is not in gid2mylid.");
        n.lid_on_rank = it->second;
      } else {
        has_comm = true;
        const auto it = gid2rmt_owning_lid.find(n.gid);
        slmm_throw_if(it == gid2rmt_owning_lid.end(),
                      "comm_lid_on_rank: On rank " << myrank << ", gid " << n.gid
                      << ", which I believe to be owned by rank " << n.rank
                      << ", is not in gid2rmt_owning_lid.");
        n.lid_on_rank = it->second;
      }
    if (has_comm)
      mylid_with_comm.push_back(i);
  }
  cm.mylid_with_comm_d.reset_capacity(mylid_with_comm.size(), true);
  cm.mylid_with_comm_h = cm.mylid_with_comm_d.mirror();
  for (Int i = 0; i < cm.mylid_with_comm_h.n(); ++i)
    cm.mylid_with_comm_h(i) = mylid_with_comm[i];
  deep_copy(cm.mylid_with_comm_d, cm.mylid_with_comm_h);

  mpi::waitall(sendreqs.size(), sendreqs.data());
}

// Useful maps between a linear index space 1:K to a set of K unique
// integers. These obviate sorts and use of hash or binary maps during time
// stepping.
template <typename MT>
void set_idx2_maps (IslMpi<MT>& cm, const Rank2Gids& rank2rmtgids,
                    const std::map<Int, Int>& gid2rmt_owning_lid) {
  const Int myrank = cm.p->rank();
  std::map<Int, Int> ranks;
  Int i = 0;
  for (const auto& e: rank2rmtgids)
    if (ranks.find(e.first) == ranks.end())
      ranks[e.first] = i++;
  ranks[myrank] = i;

  cm.ranks.reset_capacity(i+1, true);
  cm.ranks(i) = myrank;
  for (const auto& e: ranks)
    cm.ranks(e.second) = e.first;
  cm.sendcount.reset_capacity(i, true);
  cm.sendcount_h = cm.sendcount.mirror();
  cm.x_bulkdata_offset.reset_capacity(i, true);
  cm.x_bulkdata_offset_h = cm.x_bulkdata_offset.mirror();
  cm.sendreq.reset_capacity(i, true);
  cm.recvreq.reset_capacity(i, true);
  cm.recvreq_ri.reset_capacity(i, true);

  const Int nrmtrank = static_cast<Int>(cm.ranks.size()) - 1;
  std::vector<std::map<Int, Int> > lor2idx(nrmtrank);
  std::vector<Int> nlid_on_rank(nrmtrank);
  for (Int ri = 0; ri < nrmtrank; ++ri) {
    const auto& rmtgids = rank2rmtgids.at(cm.ranks(ri));
    nlid_on_rank[ri] = rmtgids.size();
  }
  cm.lid_on_rank.init(nrmtrank, nlid_on_rank.data());
  cm.lid_on_rank_h = cm.lid_on_rank.mirror();
  for (Int ri = 0; ri < nrmtrank; ++ri) {
    const auto& rmtgids = rank2rmtgids.at(cm.ranks(ri));
    auto& lor2idxri = lor2idx[ri];
    Int idx = 0;
    for (const auto& gid: rmtgids) {
      const auto lid = gid2rmt_owning_lid.at(gid);
      lor2idxri[lid] = idx;
      cm.lid_on_rank_h(ri,idx) = lid;
      ++idx;
    }
    slmm_assert(idx == nlid_on_rank[ri]);
  }
  deep_copy(cm.lid_on_rank, cm.lid_on_rank_h);

  for (i = 0; i < cm.nelemd; ++i) {
    auto& ed = cm.ed_h(i);
    ed.me->rank_idx = ranks.at(ed.me->rank);
    ed.me->lid_on_rank_idx = i;
    for (auto& n: ed.nbrs) {
      n.rank_idx = ranks.at(n.rank);
      n.lid_on_rank_idx = n.rank == myrank ?
        i :
        lor2idx[n.rank_idx].at(n.lid_on_rank);
    }
    ed.src = typename IslMpi<MT>::template ArrayH<Int**>("src", cm.nlev, cm.np2);
    ed.q_extrema = typename IslMpi<MT>::template ArrayH<Real**[2]>(
      "q_extrema", cm.qsize, cm.nlev);
  }
}

// In the original MPI pattern that has been in HOMME for years, each owned cell
// has a 1-halo patch of bulk data. For a 1-halo, allocations in this routine
// use essentially the same amount of memory, but not more. We could use less if
// we were willing to realloc space at each SL time step.
template <typename MT>
void size_mpi_buffers (IslMpi<MT>& cm, const Rank2Gids& rank2rmtgids,
                       const Rank2Gids& rank2owngids) {
  const auto myrank = cm.p->rank();
  // sizeof real, int, single int (b/c of alignment)
  const Int sor = sizeof(Real), soi = sizeof(Int), sosi = sor;
  static_assert(sizeof(Real) >= sizeof(Int),
                "For buffer packing, we require sizeof(Real) >= sizeof(Int)");
  const bool calc_trajectory = cm.traj_nsubstep > 0;
  const Int ndim = calc_trajectory ? cm.dep_points_ndim : 3;
  const Int qsize = calc_trajectory ? std::max(cm.dep_points_ndim, cm.qsize) : cm.qsize;
  const auto xbufcnt = [&] (const std::set<Int>& rmtgids,
                            const std::set<Int>& owngids,
                            const bool include_bulk = true) -> Int {
    return (sosi + (2*soi + (2*soi)*cm.nlev)*rmtgids.size() +               // meta data
            (include_bulk ? 1 : 0)*owngids.size()*cm.nlev*cm.np2*ndim*sor); // bulk data
  };
  const auto qbufcnt = [&] (const std::set<Int>& rmtgids,
                            const std::set<Int>& owngids) -> Int {
    return ((rmtgids.size()*2 +      // min/max q
             owngids.size()*cm.np2)* // q
            qsize*cm.nlev*sor);
  };
  const auto bytes2real = [&] (const Int& bytes) {
    return (bytes + sor - 1)/sor;
  };

  slmm_assert(cm.ranks.back() == myrank);
  const Int nrmtrank = static_cast<Int>(cm.ranks.size()) - 1;
  cm.nlid_per_rank.resize(nrmtrank);
  cm.sendsz.resize(nrmtrank);
  cm.recvsz.resize(nrmtrank);
  cm.sendmetasz.resize(nrmtrank);
  cm.recvmetasz.resize(nrmtrank);
  Int rmt_xs_sz = 0, rmt_qse_sz = 0;
  for (Int ri = 0; ri < nrmtrank; ++ri) {
    const auto& rmtgids = rank2rmtgids.at(cm.ranks(ri));
    const auto& owngids = rank2owngids.at(cm.ranks(ri));
    cm.nlid_per_rank[ri] = rmtgids.size();
    cm.sendsz[ri] = bytes2real(std::max(xbufcnt(rmtgids, owngids),
                                        qbufcnt(owngids, rmtgids)));
    cm.recvsz[ri] = bytes2real(std::max(xbufcnt(owngids, rmtgids),
                                        qbufcnt(rmtgids, owngids)));
    rmt_xs_sz  += 5*cm.np2*cm.nlev*rmtgids.size();
    rmt_qse_sz += 4       *cm.nlev*rmtgids.size();
#ifdef COMPOSE_PORT_SEPARATE_VIEWS
    cm.sendmetasz[ri] = bytes2real(xbufcnt(rmtgids, owngids));
    cm.recvmetasz[ri] = bytes2real(xbufcnt(owngids, rmtgids));
#endif
  }
  cm.rmt_xs.reset_capacity(rmt_xs_sz, true);
  cm.rmt_qs_extrema.reset_capacity(rmt_qse_sz, true);
  cm.rmt_xs_h = cm.rmt_xs.mirror();
  cm.rmt_qs_extrema_h = cm.rmt_qs_extrema.mirror();
}

template <typename MT>
void alloc_mpi_buffers (IslMpi<MT>& cm, Real* sendbuf, Real* recvbuf) {
  const Int nrmtrank = static_cast<Int>(cm.ranks.size()) - 1;
  cm.nx_in_rank.reset_capacity(nrmtrank, true);
  cm.nx_in_rank_h = cm.nx_in_rank.mirror();
  cm.nx_in_lid.init(nrmtrank, cm.nlid_per_rank.data());
  cm.nx_in_lid_h = cm.nx_in_lid.mirror();
  cm.bla.init(nrmtrank, cm.nlid_per_rank.data(), cm.nlev);
  cm.bla_h = cm.bla.mirror();
  cm.sendbuf.init(nrmtrank, cm.sendsz.data(), sendbuf);
  cm.recvbuf.init(nrmtrank, cm.recvsz.data(), recvbuf);
#ifdef COMPOSE_MPI_ON_HOST
  cm.sendbuf_h = cm.sendbuf.mirror();
  cm.recvbuf_h = cm.recvbuf.mirror();
#endif
  cm.nlid_per_rank.clear();
  cm.sendsz.clear();
  cm.recvsz.clear();
#ifdef COMPOSE_PORT_SEPARATE_VIEWS
  slmm_assert(sendbuf == nullptr && recvbuf == nullptr);
  cm.sendbuf_meta_h.init(nrmtrank, cm.sendmetasz.data());
  cm.recvbuf_meta_h.init(nrmtrank, cm.recvmetasz.data());
#else
  cm.sendbuf_meta_h = cm.sendbuf;
  cm.recvbuf_meta_h = cm.recvbuf;
#endif
#ifdef COMPOSE_HORIZ_OPENMP
  cm.ri_lidi_locks.init(nrmtrank, cm.nlid_per_rank.data());
  for (Int ri = 0; ri < nrmtrank; ++ri) {
    auto&& locks = cm.ri_lidi_locks(ri);
    for (auto& lock: locks)
      omp_init_lock(&lock);
  }
#endif  
}

// At simulation initialization, set up a bunch of stuff to make the work at
// each step as small as possible.
template <typename MT>
void setup_comm_pattern (IslMpi<MT>& cm, const Int* nbr_id_rank, const Int* nirptr) {
  collect_gid_rank(cm, nbr_id_rank, nirptr);
  Rank2Gids rank2rmtgids, rank2owngids;
  get_rank2gids(cm, rank2rmtgids, rank2owngids);
  {
    std::map<Int, Int> gid2rmt_owning_lid;
    comm_lid_on_rank(cm, rank2rmtgids, rank2owngids, gid2rmt_owning_lid);
    set_idx2_maps(cm, rank2rmtgids, gid2rmt_owning_lid);
  }
  size_mpi_buffers(cm, rank2rmtgids, rank2owngids);
}

template void
alloc_mpi_buffers(IslMpi<ko::MachineTraits>& cm, Real* sendbuf, Real* recvbuf);
template void
setup_comm_pattern(IslMpi<ko::MachineTraits>& cm, const Int* nbr_id_rank,
                   const Int* nirptr);

} // namespace islmpi
} // namespace homme
