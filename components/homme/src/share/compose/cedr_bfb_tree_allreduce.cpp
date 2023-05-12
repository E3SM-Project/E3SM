// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

#include "cedr_bfb_tree_allreduce.hpp"

namespace cedr {
using mpi::Parallel;

template <typename ES>
BfbTreeAllReducer<ES>
::BfbTreeAllReducer (const Parallel::Ptr& p, const tree::Node::Ptr& tree,
                     const Int nleaf, const Int nfield) {
  init(p, tree, nleaf, nfield);
}

template <typename ES>
void BfbTreeAllReducer<ES>
::init (const Parallel::Ptr& p, const tree::Node::Ptr& tree,
        const Int nleaf, const Int nfield) {
  p_ = p;
  ns_ = tree::analyze(p, nleaf, tree);
  nlocal_ = static_cast<Int>(ns_->levels[0].nodes.size());
  nfield_ = nfield;
}

template <typename ES>
void BfbTreeAllReducer<ES>
::get_host_buffers_sizes (size_t& buf1, size_t& buf2) {
  buf1 = ((impl::OnGpu<ES>::value ? nlocal_ : 0) + ns_->nslots)*nfield_;
  buf2 = 0;
}

template <typename ES>
void BfbTreeAllReducer<ES>
::set_host_buffers (Real* buf1, Real* buf2) {
  if ( ! buf1) return;
  size_t s1, s2;
  get_host_buffers_sizes(s1, s2);
  bd_ = RealListHost(buf1, s1);
}

template <typename ES>
void BfbTreeAllReducer<ES>
::finish_setup () {
  size_t s1, s2;
  get_host_buffers_sizes(s1, s2);
  if (bd_.size() > 0) {
    cedr_assert(bd_.size() == s1);
    return;
  }
  bd_ = RealListHost("bd_", s1);
}

template <typename ES>
const Real* BfbTreeAllReducer<ES>
::get_send_host (const ConstRealList& send) const {
  cedr_assert(send.extent_int(0) == nlocal_*nfield_);
  if (impl::OnGpu<ES>::value) {
    RealListHost m(bd_.data() + ns_->nslots * nfield_, nlocal_*nfield_);
    Kokkos::deep_copy(m, send);
    return m.data();
  } else {
    return send.data();
  }
}

template <typename ES>
void BfbTreeAllReducer<ES>
::fill_recv (const RealList& recv) const {
  Real* const buf = impl::OnGpu<ES>::value ? bd_.data() + ns_->nslots * nfield_ : recv.data();
  const Int idx = ns_->levels[0].nodes[0];
  const auto n = ns_->node_h(idx);
  Real* const redval = &bd_(n->offset*nfield_);
  for (Int i = 0; i < nfield_; ++i) buf[i] = redval[i];
  if (impl::OnGpu<ES>::value) Kokkos::deep_copy(recv, RealListHost(buf, nfield_));
}

template <typename ES>
void BfbTreeAllReducer<ES>
::allreduce (const ConstRealList& send, const RealList& recv, const bool transpose) const {
  const auto mpitag = tree::NodeSets::mpitag;
  const auto& ns = *ns_;
  const auto nf = nfield_;
  cedr_assert(ns.levels[0].nodes.size() == static_cast<size_t>(nlocal_));

  { // We want to be behaviorally const but still permit lazy finish_setup.
    //   Cuda 10.1.105 with GCC 8.5.0 incorrectly misses the const_cast; for
    // some reason, adding remove_const fixes that.
    typedef typename std::remove_const<Me>::type MeLcl;
    if ( ! bd_.size()) const_cast<MeLcl*>(this)->finish_setup();
  }

  const auto buf_host = get_send_host(send);

  // Fill leaf data.
  for (Int id = 0; id < nlocal_; ++id) {
    const auto& idx = ns.levels[0].nodes[id];
    const auto n = ns.node_h(idx);
    if (transpose) {
      auto d = &bd_[n->offset * nf];
      for (Int j = 0; j < nf; ++j) d[j] = buf_host[nlocal_*j + id];
    } else {
      const auto s = buf_host + id*nf;
      auto d = &bd_[n->offset * nf];
      for (Int j = 0; j < nf; ++j) d[j] = s[j];
    }
  }
  // Leaves to root.
  for (size_t il = 0; il < ns.levels.size(); ++il) {
    auto& lvl = ns.levels[il];
    // Set up receives.
    for (size_t i = 0; i < lvl.kids.size(); ++i) {
      const auto& mmd = lvl.kids[i];
      mpi::irecv(*p_, &bd_[mmd.offset * nf], mmd.size * nf, mmd.rank, mpitag,
                 &lvl.kids_req[i]);
    }
    mpi::waitall(lvl.kids_req.size(), lvl.kids_req.data());
    // Combine kids' data.
    for (const auto& idx : lvl.nodes) {
      const auto n = ns.node_h(idx);
      if ( ! n->nkids) continue;
      const auto d = &bd_[n->offset * nf];
      for (Int j = 0; j < nf; ++j) d[j] = 0;
      for (Int i = 0; i < n->nkids; ++i) {
        const auto s = &bd_[ns.node_h(n->kids[i])->offset * nf];
        for (Int j = 0; j < nf; ++j) d[j] += s[j];
      }
    }
    // Send to parents.
    for (size_t i = 0; i < lvl.me.size(); ++i) {
      const auto& mmd = lvl.me[i];
      mpi::isend(*p_, &bd_[mmd.offset * nf], mmd.size * nf, mmd.rank, mpitag);
    }
  }
  // Root to leaves.
  for (size_t il = ns.levels.size(); il > 0; --il) {
    auto& lvl = ns.levels[il-1];
    // Get the global sum from parent.
    for (size_t i = 0; i < lvl.me.size(); ++i) {
      const auto& mmd = lvl.me[i];
      mpi::irecv(*p_, &bd_[mmd.offset * nf], mmd.size * nf, mmd.rank, mpitag,
                 &lvl.me_recv_req[i]);
    }    
    mpi::waitall(lvl.me_recv_req.size(), lvl.me_recv_req.data());
    // Pass to kids.
    for (const auto& idx : lvl.nodes) {
      const auto n = ns.node_h(idx);
      if ( ! n->nkids) continue;
      const auto s = &bd_[n->offset * nf];
      for (Int i = 0; i < n->nkids; ++i) {
        auto d = &bd_[ns.node_h(n->kids[i])->offset * nf];
        for (Int j = 0; j < nf; ++j) d[j] = s[j];
      }
    }
    // Send.
    for (size_t i = 0; i < lvl.kids.size(); ++i) {
      const auto& mmd = lvl.kids[i];
      mpi::isend(*p_, &bd_[mmd.offset * nf], mmd.size * nf, mmd.rank, mpitag);
    }
  }

  fill_recv(recv);
}

template <typename ES>
Int BfbTreeAllReducer<ES>::unittest (const Parallel::Ptr& p) {
  using Mesh = tree::oned::Mesh;
  static const Int nfield = 3;
  const Int myrank = p->rank();
  const Int szs[] = {p->size(), 3*p->size()};
  const Mesh::ParallelDecomp::Enum dists[] = {Mesh::ParallelDecomp::pseudorandom,
                                              Mesh::ParallelDecomp::contiguous};
  Int nerr = 0;
  for (size_t is = 0; is < sizeof(szs)/sizeof(*szs); ++is)
    for (size_t id = 0; id < sizeof(dists)/sizeof(*dists); ++id)
      for (bool imbalanced : {false, true})
        for (bool transpose : {false, true}) {
          const Int ncell = szs[is];
          Mesh m(ncell, p, dists[id]);
          tree::Node::Ptr tree = make_tree(m, imbalanced);

          Int nlocal = 0;
          for (Int i = 0; i < ncell; ++i) if (myrank == m.rank(i)) ++nlocal;

          BfbTreeAllReducer<>::RealList send("send", nlocal*nfield), recv("recv", nfield);
          std::vector<Real> mpi_recv(nfield);
          const auto send_m = Kokkos::create_mirror_view(send);
          const auto recv_m = Kokkos::create_mirror_view(recv);
          if (transpose)
            for (Int j = 0; j < nfield; ++j)
              for (Int i = 0; i < nlocal; ++i)
                send_m(nlocal*j + i) = util::urand();
          else
            for (Int i = 0; i < nlocal; ++i)
              for (Int j = 0; j < nfield; ++j)
                send_m(nfield*i + j) = util::urand();
          Kokkos::deep_copy(send, send_m);

          BfbTreeAllReducer<> ar(p, tree, m.ncell(), nfield);
          ar.allreduce(send, recv, transpose);
          Kokkos::deep_copy(recv_m, recv);

          std::vector<Real> lcl_red(nfield, 0);
          if (transpose) {
            for (Int j = 0; j < nfield; ++j)
              for (Int i = 0; i < nlocal; ++i)
                lcl_red[j] += send_m[nlocal*j + i];
          } else {
            for (Int i = 0; i < nlocal; ++i)
              for (Int j = 0; j < nfield; ++j)
                lcl_red[j] += send_m[nfield*i + j];
          }
          mpi::all_reduce(*p, lcl_red.data(), mpi_recv.data(), nfield, MPI_SUM);

          Real max_re = 0;
          for (Int j = 0; j < nfield; ++j)
            max_re = std::max(max_re, util::reldif(mpi_recv[j], recv_m(j)));
          if (max_re > 2*std::log(ncell)*std::numeric_limits<Real>::epsilon()) {
            printf("FAIL BfbTreeAllReducer<>::unittest max_re %1.3e\n", max_re);
            ++nerr;
          }
        }
  return nerr;
}

} // namespace cedr

#ifdef KOKKOS_ENABLE_SERIAL
template class cedr::BfbTreeAllReducer<Kokkos::Serial>;
#endif
#ifdef KOKKOS_ENABLE_OPENMP
template class cedr::BfbTreeAllReducer<Kokkos::OpenMP>;
#endif
#ifdef CEDR_ENABLE_GPU
template class cedr::BfbTreeAllReducer<CedrGpuExeSpace>;
#endif
#ifdef KOKKOS_ENABLE_THREADS
template class cedr::BfbTreeAllReducer<Kokkos::Threads>;
#endif
