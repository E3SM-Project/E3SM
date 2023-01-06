// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

#ifndef INCLUDE_CEDR_TREE_CALLER_HPP
#define INCLUDE_CEDR_TREE_CALLER_HPP

#include "cedr_mpi.hpp"

namespace cedr {
namespace tree {
// The caller builds a tree of these nodes to pass to QLT.
struct Node {
  typedef std::shared_ptr<Node> Ptr;
  const Node* parent; // (Can't be a shared_ptr: would be a circular dependency.)
  Int rank;           // Owning rank.
  Long cellidx;       // If a leaf, the cell to which this node corresponds.
  Int nkids;          // 0 at leaf, 1 or 2 otherwise.
  Node::Ptr kids[2];
  Int reserved;       // For internal use.
  Int level;          // If providing only partial trees, set level to
                      // the level of this node, with a leaf node at
                      // level 0.
  Node () : parent(nullptr), rank(-1), cellidx(-1), nkids(0), reserved(-1), level(-1) {}
};

// Utility to make a tree over a 1D mesh. For testing, it can be useful to
// create an imbalanced tree.
Node::Ptr make_tree_over_1d_mesh(const mpi::Parallel::Ptr& p, const Int& ncells,
                                 const bool imbalanced = false);

} // namespace tree
} // namespace cedr

#endif
