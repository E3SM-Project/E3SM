// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

#ifndef INCLUDE_CEDR_TREE_HPP
#define INCLUDE_CEDR_TREE_HPP

#include "cedr_tree_caller.hpp"
#include "cedr_util.hpp"

#include <vector>

namespace cedr {
namespace tree {
using cedr::mpi::Parallel;

struct NodeSets {
  typedef std::shared_ptr<const NodeSets> ConstPtr;
  
  enum : int { mpitag = 42 };

  // A node in the tree that is relevant to this rank.
  struct Node {
    // Rank of the node. If the node is in a level, then its rank is my rank. If
    // it's not in a level, then it is a comm partner of a node on this rank.
    Int rank;
    // cellidx if leaf node, ie, if nkids == 0; otherwise, undefined.
    Int id;
    // This node's parent, a comm partner, if such a partner is required.
    Int parent;
    // This node's kids, comm partners, if such partners are required. Parent
    // and kid nodes are pruned relative to the full tree over the mesh to
    // contain just the nodes that matter to this rank.
    Int nkids;
    Int kids[2];
    // Offset factor into bulk data. An offset is a unit; actual buffer sizes
    // are multiples of this unit.
    Int offset;

    KOKKOS_FUNCTION Node () : rank(-1), id(-1), parent(-1), nkids(0), offset(-1) {}
  };

  // A level in the level schedule that is constructed to orchestrate
  // communication. A node in a level depends only on nodes in lower-numbered
  // levels (l2r) or higher-numbered (r2l).
  //
  // The communication patterns are as follows:
  //   > l2r
  //   MPI rcv into kids
  //   sum into node
  //   MPI send from node
  //   > r2l
  //   MPI rcv into node
  //   solve QP for kids
  //   MPI send from kids
  struct Level {
    struct MPIMetaData {
      Int rank;   // Rank of comm partner.
      Int offset; // Offset to start of buffer for this comm.
      Int size;   // Size of this buffer in units of offsets.
    };
    
    // The nodes in the level.
    std::vector<Int> nodes;
    // MPI information for this level.
    std::vector<MPIMetaData> me, kids;
    mutable std::vector<mpi::Request> me_send_req, me_recv_req, kids_req;
  };
  
  // Levels. levels[0] is level 0, the leaf level.
  std::vector<Level> levels;
  // Number of data slots this rank needs. Each node owned by this rank, plus
  // kids on other ranks, have an associated slot.
  Int nslots;
  
  // Allocate a node. The list node_mem_ is the mechanism for memory ownership;
  // node_mem_ isn't used for anything other than owning nodes.
  Int alloc () {
    const Int idx = node_mem_.size();
    node_mem_.push_back(Node());
    return idx;
  }

  Int nnode () const { return node_mem_.size(); }

  Node* node_h (const Int& idx) {
    cedr_assert(idx >= 0 && idx < static_cast<Int>(node_mem_.size()));
    return &node_mem_[idx];
  }
  const Node* node_h (const Int& idx) const {
    return const_cast<NodeSets*>(this)->node_h(idx);
  }

  void print(std::ostream& os) const;
  
private:
  std::vector<Node> node_mem_;
};

template <typename ExeSpace>
struct NodeSetsDeviceData {
  typedef typename cedr::impl::DeviceType<ExeSpace>::type Device;
  typedef Kokkos::View<Int*, Device> IntList;
  typedef Kokkos::View<NodeSets::Node*, Device> NodeList;

  NodeList node;
  // lvl(lvlptr(l):lvlptr(l+1)-1) is the list of node indices into node for
  // level l.
  IntList lvl, lvlptr;
};

typedef NodeSetsDeviceData<Kokkos::DefaultHostExecutionSpace> NodeSetsHostData;

// Analyze the tree to extract levels. Levels are run from 0 to #level - 1. Each
// level has nodes whose corresponding operations depend on only nodes in
// lower-indexed levels. This mechanism prevents deadlock in the general case of
// multiple cells per rank, with multiple ranks appearing in a subtree other
// than the root.
//   In addition, the set of nodes collected into levels are just those owned by
// this rank, and those with which owned nodes must communicate.
NodeSets::ConstPtr analyze(const Parallel::Ptr& p, const Int& ncells,
                           const tree::Node::Ptr& tree);

Int unittest(const Parallel::Ptr& p);

// Tree for a 1-D periodic domain, for unit testing.
namespace oned {
struct Mesh {
  struct ParallelDecomp {
    enum Enum {
      // The obvious distribution of ranks: 1 rank takes exactly 1 contiguous
      // set of cell indices.
      contiguous,
      // For heavy-duty testing of QLT comm pattern, use a ridiculous assignment
      // of ranks to cell indices. This forces the QLT tree to communicate,
      // pack, and unpack in silly ways.
      pseudorandom
    };
  };
  
  Mesh (const Int nc, const Parallel::Ptr& p,
        const ParallelDecomp::Enum& parallel_decomp = ParallelDecomp::contiguous) {
    init(nc, p, parallel_decomp);
  }
  
  void init(const Int nc, const Parallel::Ptr& p,
            const ParallelDecomp::Enum& parallel_decomp);

  Int ncell () const { return nc_; }

  const Parallel::Ptr& parallel () const { return p_; }

  Int rank(const Int& ci) const;

  static Int unittest(const Parallel::Ptr& p);

private:
  Int nc_, nranks_;
  Parallel::Ptr p_;
  ParallelDecomp::Enum pd_;
};

tree::Node::Ptr make_tree(const Parallel::Ptr& p, const Int& ncells,
                          const bool imbalanced);
tree::Node::Ptr make_tree(const Mesh& m, const bool imbalanced);

Int unittest(const Parallel::Ptr& p);
} // namespace oned
} // namespace tree
} // namespace cedr

#endif
