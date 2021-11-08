// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

#ifndef INCLUDE_CEDR_BFB_TREE_ALLREDUCE
#define INCLUDE_CEDR_BFB_TREE_ALLREDUCE

#include "cedr_tree.hpp"

namespace cedr {

// Use a tree and point-to-point communication to implement all-reduce. If the
// tree is independent of process deomposition, then
// BfbTreeAllReducer::allreduce is BFB-invariant to process decomposition.
template <typename ExeSpace = Kokkos::DefaultExecutionSpace>
struct BfbTreeAllReducer {
  typedef typename cedr::impl::DeviceType<ExeSpace>::type Device;
  typedef BfbTreeAllReducer<ExeSpace> Me;
  typedef std::shared_ptr<Me> Ptr;
  typedef Kokkos::View<Real*, Device> RealList;
  typedef typename Kokkos::View<Real*, Device>::HostMirror RealListHost;
  typedef Kokkos::View<const Real*, Device> ConstRealList;

  BfbTreeAllReducer(const mpi::Parallel::Ptr& p, const tree::Node::Ptr& tree,
                    // A leaf is a leaf node in the reduction tree. The global
                    // tree has nleaf leaves. Each leaf has nfield scalars to
                    // reduce.
                    const Int nleaf, const Int nfield);

  // All three are optional.
  void get_host_buffers_sizes(size_t& buf1, size_t& buf2);
  void set_host_buffers(Real* buf1, Real* buf2);
  // If you don't call this, it will be called on the first allreduce. So call
  // this if you're doing timing runs before the time stepping loop.
  void finish_setup();

  // In Fortran, these are formatted as recv(nfield) and -- with fastest index
  // last -- send(nfield, nlocal) if transpose or send(nlocal, nfield)
  // otherwise, where nlocal is the number of leaf nodes on this rank. send and
  // recv can point to the same memory.
  void allreduce(const ConstRealList& send, const RealList& recv,
                 const bool transpose = false) const;

  static Int unittest(const mpi::Parallel::Ptr& p);

private:
  mpi::Parallel::Ptr p_;
  Int nlocal_, nfield_;
  std::shared_ptr<const tree::NodeSets> ns_;
  mutable RealListHost bd_;

  void init(const mpi::Parallel::Ptr& p, const tree::Node::Ptr& tree,
            const Int nleaf, const Int nfield);
  const Real* get_send_host(const ConstRealList& send) const;
  void fill_recv(const RealList& recv) const;
};

} // namespace cedr

#endif
