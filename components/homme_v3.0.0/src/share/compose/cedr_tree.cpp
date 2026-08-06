#include "cedr_tree.hpp"

#include <set>

namespace cedr {
namespace tree {

void NodeSets::print (std::ostream& os) const {
  std::stringstream ss;
  if (levels.empty()) return;
  const Int myrank = node_h(levels[0].nodes[0])->rank;
  ss << "pid " << myrank << ":";
  ss << " #levels " << levels.size();
  for (size_t i = 0; i < levels.size(); ++i) {
    const auto& lvl = levels[i];
    ss << "\n  " << i << ": " << lvl.nodes.size();
    std::set<Int> ps, ks;
    for (size_t j = 0; j < lvl.nodes.size(); ++j) {
      const auto n = node_h(lvl.nodes[j]);
      for (Int k = 0; k < n->nkids; ++k)
        if (node_h(n->kids[k])->rank != myrank)
          ks.insert(node_h(n->kids[k])->rank);
      if (n->parent >= 0 && node_h(n->parent)->rank != myrank)
        ps.insert(node_h(n->parent)->rank);
    }
    ss << " |";
    for (const auto& e : ks) ss << " " << e;
    if ( ! lvl.kids.empty()) ss << " (" << lvl.kids.size() << ") |";
    for (const auto& e : ps) ss << " " << e;
    if ( ! lvl.me.empty()) ss << " (" << lvl.me.size() << ")";
  }
  ss << "\n";
  os << ss.str();
}

// Find tree depth, assign ranks to non-leaf nodes, and init 'reserved'.
Int init_tree (const Int& my_rank, const tree::Node::Ptr& node, Int& id) {
  node->reserved = -1;
  Int depth = 0;
  for (Int i = 0; i < node->nkids; ++i) {
    cedr_assert(node.get() == node->kids[i]->parent);
    depth = std::max(depth, init_tree(my_rank, node->kids[i], id));
  }
  if (node->nkids) {
    if (node->rank < 0) node->rank = node->kids[0]->rank;
    node->cellidx = id++;
  } else {
    cedr_throw_if(node->rank == my_rank && (node->cellidx < 0 || node->cellidx >= id),
                  "cellidx is " << node->cellidx << " but should be between " <<
                  0 << " and " << id);
  }
  return depth + 1;
}

void level_schedule_and_collect (
  NodeSets& ns, const Int& my_rank, const tree::Node::Ptr& node, Int& level,
  bool& need_parent_ns_node)
{
  cedr_assert(node->rank != -1);
  level = -1;
  bool make_ns_node = false;
  for (Int i = 0; i < node->nkids; ++i) {
    Int kid_level;
    bool kid_needs_ns_node;
    level_schedule_and_collect(ns, my_rank, node->kids[i], kid_level,
                               kid_needs_ns_node);
    level = std::max(level, kid_level);
    if (kid_needs_ns_node) make_ns_node = true;
  }
  ++level;
  if (node->level >= 0) {
    // The caller built only partial trees and so must provide the
    // level.
    level = node->level;
    cedr_assert(level < static_cast<Int>(ns.levels.size()));
  }
  // Is parent node needed for isend?
  const bool node_is_owned = node->rank == my_rank;
  need_parent_ns_node = node_is_owned;
  if (node_is_owned || make_ns_node) {
    cedr_assert(node->reserved == -1);
    const auto ns_node_idx = ns.alloc();
    // Levels hold only owned nodes.
    if (node_is_owned) ns.levels[level].nodes.push_back(ns_node_idx);
    node->reserved = ns_node_idx;
    NodeSets::Node* ns_node = ns.node_h(ns_node_idx);
    ns_node->rank = node->rank;
    if (node->nkids == 0) ns_node->id = node->cellidx;
    if (node_is_owned) {
      // If this node is owned, it needs to have information about all kids.
      ns_node->nkids = node->nkids;
      for (Int i = 0; i < node->nkids; ++i) {
        const auto& kid = node->kids[i];
        if (kid->reserved == -1) {
          // This kid isn't owned by this rank. But need it for irecv.
          const auto ns_kid_idx = ns.alloc();
          NodeSets::Node* ns_kid = ns.node_h(ns_kid_idx);
          kid->reserved = ns_kid_idx;
          ns_node = ns.node_h(ns_node_idx);
          ns_node->kids[i] = ns_kid_idx;
          cedr_assert(kid->rank != my_rank);
          ns_kid->rank = kid->rank;
          ns_kid->id = kid->cellidx;
          // The kid may have kids in the original tree, but in the tree pruned
          // according to rank, it does not.
          ns_kid->nkids = 0;
        } else {
          // This kid is owned by this rank, so fill in its parent pointer.
          NodeSets::Node* ns_kid = ns.node_h(kid->reserved);
          ns_node = ns.node_h(ns_node_idx);
          ns_node->kids[i] = kid->reserved;
          ns_kid->parent = ns_node_idx;
        }
      }
    } else {
      // This node is not owned. Update the owned kids with its parent.
      ns_node->nkids = 0;
      for (Int i = 0; i < node->nkids; ++i) {
        const auto& kid = node->kids[i];
        if (kid->reserved >= 0 && kid->rank == my_rank) {
          const auto ns_kid_idx = kid->reserved;
          ns_node->kids[ns_node->nkids++] = ns_kid_idx;
          NodeSets::Node* ns_kid = ns.node_h(ns_kid_idx);
          ns_kid->parent = ns_node_idx;
        }
      }
    }
  }
}

void level_schedule_and_collect (NodeSets& ns, const Int& my_rank,
                                 const tree::Node::Ptr& tree) {
  Int iunused;
  bool bunused;
  level_schedule_and_collect(ns, my_rank, tree, iunused, bunused);
}

void consolidate (NodeSets& ns) {
  auto levels = ns.levels;
  ns.levels.clear();
  for (const auto& level : levels)
    if ( ! level.nodes.empty())
      ns.levels.push_back(level);
}

typedef std::pair<Int, NodeSets::Node*> RankNode;

void init_offsets (const Int my_rank, std::vector<RankNode>& rns,
                   std::vector<NodeSets::Level::MPIMetaData>& mmds, Int& offset) {
  // Set nodes on my rank to have rank -1 so that they sort first.
  for (auto& rn : rns)
    if (rn.first == my_rank)
      rn.first = -1;

  // Sort so that all comms with a given rank are contiguous. Stable sort so
  // that rns retains its order, in particular in the leaf node level.
  std::stable_sort(rns.begin(), rns.end());

  // Collect nodes into groups by rank and set up comm metadata for each group.
  Int prev_rank = -1;
  for (auto& rn : rns) {
    const Int rank = rn.first;
    if (rank == -1) {
      if (rn.second->offset == -1)
        rn.second->offset = offset++;
      continue;
    }
    if (rank != prev_rank) {
      cedr_assert(rank > prev_rank);
      prev_rank = rank;
      mmds.push_back(NodeSets::Level::MPIMetaData());
      auto& mmd = mmds.back();
      mmd.rank = rank;
      mmd.offset = offset;
      mmd.size = 0;
    }
    ++mmds.back().size;
    rn.second->offset = offset++;
  }
}

// Set up comm data. Consolidate so that there is only one message between me
// and another rank per level. Determine an offset for each node, to be
// multiplied by data-size factors later, for use in data buffers.
void init_comm (const Int my_rank, NodeSets& ns) {
  ns.nslots = 0;
  for (auto& lvl : ns.levels) {
    Int nkids = 0;
    for (const auto& idx : lvl.nodes) {
      const auto n = ns.node_h(idx);
      nkids += n->nkids;
    }

    std::vector<RankNode> me(lvl.nodes.size()), kids(nkids);
    for (size_t i = 0, mi = 0, ki = 0; i < lvl.nodes.size(); ++i) {
      const auto n = ns.node_h(lvl.nodes[i]);
      me[mi].first = n->parent >= 0 ? ns.node_h(n->parent)->rank : my_rank;
      me[mi].second = const_cast<NodeSets::Node*>(n);
      ++mi;
      for (Int k = 0; k < n->nkids; ++k) {
        kids[ki].first = ns.node_h(n->kids[k])->rank;
        kids[ki].second = ns.node_h(n->kids[k]);
        ++ki;
      }
    }

    init_offsets(my_rank, me, lvl.me, ns.nslots);
    lvl.me_send_req.resize(lvl.me.size());
    lvl.me_recv_req.resize(lvl.me.size());
    init_offsets(my_rank, kids, lvl.kids, ns.nslots);
    lvl.kids_req.resize(lvl.kids.size());
  }
}

NodeSets::ConstPtr analyze (const Parallel::Ptr& p, const Int& ncells,
                            const tree::Node::Ptr& tree) {
  const auto nodesets = std::make_shared<NodeSets>();
  cedr_assert( ! tree->parent);
  Int id = ncells;
  Int depth = init_tree(p->rank(), tree, id);
  if (tree->level >= 0) {
    // If level is provided, don't trust depth from init_tree. Partial trees can
    // make depth too small.
    depth = tree->level + 1;
  }
  nodesets->levels.resize(depth);
  level_schedule_and_collect(*nodesets, p->rank(), tree);
  consolidate(*nodesets);
  init_comm(p->rank(), *nodesets);
  return nodesets;
}

// Check that the offsets are self consistent.
Int check_comm (const NodeSets& ns) {
  Int nerr = 0;
  std::vector<Int> offsets(ns.nslots, 0);
  for (const auto& lvl : ns.levels)
    for (const auto& idx : lvl.nodes) {
      const auto n = ns.node_h(idx);
      cedr_assert(n->offset < ns.nslots);
      ++offsets[n->offset];
      for (Int i = 0; i < n->nkids; ++i) {
        const auto kid = ns.node_h(n->kids[i]);
        if (kid->rank != n->rank)
          ++offsets[kid->offset];
      }
    }
  for (const auto& e : offsets)
    if (e != 1) ++nerr;
  return nerr;
}

// Check that there are the correct number of leaf nodes, and that their offsets
// all come first and are ordered the same as ns->levels[0]->nodes.
Int check_leaf_nodes (const Parallel::Ptr& p, const NodeSets& ns,
                      const Int ncells) {
  Int nerr = 0;
  cedr_assert( ! ns.levels.empty());
  cedr_assert( ! ns.levels[0].nodes.empty());
  Int my_nleaves = 0;
  for (const auto& idx : ns.levels[0].nodes) {
    const auto n = ns.node_h(idx);
    cedr_assert( ! n->nkids);
    ++my_nleaves;
  }
  for (const auto& idx : ns.levels[0].nodes) {
    const auto n = ns.node_h(idx);
    cedr_assert(n->offset < my_nleaves);
    cedr_assert(n->id < ncells);
  }
  Int glbl_nleaves = 0;
  mpi::all_reduce(*p, &my_nleaves, &glbl_nleaves, 1, MPI_SUM);
  if (glbl_nleaves != ncells)
    ++nerr;
  return nerr;
}

// Sum cellidx using the QLT comm pattern.
Int test_comm_pattern (const Parallel::Ptr& p, const NodeSets& ns,
                       const Int ncells) {
  Int nerr = 0;
  // Rank-wide data buffer.
  std::vector<Int> data(ns.nslots);
  // Sum this rank's cellidxs.
  for (const auto& idx : ns.levels[0].nodes) {
    const auto n = ns.node_h(idx);
    data[n->offset] = n->id;
  }
  // Leaves to root.
  for (size_t il = 0; il < ns.levels.size(); ++il) {
    auto& lvl = ns.levels[il];
    // Set up receives.
    for (size_t i = 0; i < lvl.kids.size(); ++i) {
      const auto& mmd = lvl.kids[i];
      mpi::irecv(*p, &data[mmd.offset], mmd.size, mmd.rank, NodeSets::mpitag,
                 &lvl.kids_req[i]);
    }
    mpi::waitall(lvl.kids_req.size(), lvl.kids_req.data());
    // Combine kids' data.
    for (const auto& idx : lvl.nodes) {
      const auto n = ns.node_h(idx);
      if ( ! n->nkids) continue;
      data[n->offset] = 0;
      for (Int i = 0; i < n->nkids; ++i)
        data[n->offset] += data[ns.node_h(n->kids[i])->offset];
    }
    // Send to parents.
    for (size_t i = 0; i < lvl.me.size(); ++i) {
      const auto& mmd = lvl.me[i];
      mpi::isend(*p, &data[mmd.offset], mmd.size, mmd.rank, NodeSets::mpitag);
    }
  }
  // Root to leaves.
  for (size_t il = ns.levels.size(); il > 0; --il) {
    auto& lvl = ns.levels[il-1];
    // Get the global sum from parent.
    for (size_t i = 0; i < lvl.me.size(); ++i) {
      const auto& mmd = lvl.me[i];
      mpi::irecv(*p, &data[mmd.offset], mmd.size, mmd.rank, NodeSets::mpitag,
                 &lvl.me_recv_req[i]);
    }    
    mpi::waitall(lvl.me_recv_req.size(), lvl.me_recv_req.data());
    // Pass to kids.
    for (const auto& idx : lvl.nodes) {
      const auto n = ns.node_h(idx);
      if ( ! n->nkids) continue;
      for (Int i = 0; i < n->nkids; ++i)
        data[ns.node_h(n->kids[i])->offset] = data[n->offset];
    }
    // Send.
    for (size_t i = 0; i < lvl.kids.size(); ++i) {
      const auto& mmd = lvl.kids[i];
      mpi::isend(*p, &data[mmd.offset], mmd.size, mmd.rank, NodeSets::mpitag);
    }
  }
  { // Check that all leaf nodes have the right number.
    const Int desired_sum = (ncells*(ncells - 1)) / 2;
    for (const auto& idx : ns.levels[0].nodes) {
      const auto n = ns.node_h(idx);
      if (data[n->offset] != desired_sum) ++nerr;
    }
    if (false && p->amroot()) {
      std::cout << " " << data[ns.node_h(ns.levels[0].nodes[0])->offset];
      std::cout.flush();
    }
  }
  return nerr;
}

Node::Ptr make_tree_over_1d_mesh (const Parallel::Ptr& p, const Int& ncells,
                                  const bool imbalanced) {
  return oned::make_tree(oned::Mesh(ncells, p), imbalanced);
}

// Tree for a 1-D periodic domain, for unit testing.
namespace oned {
void Mesh::init (const Int nc, const Parallel::Ptr& p,
                 const ParallelDecomp::Enum& parallel_decomp) {
  nc_ = nc;
  nranks_ = p->size();
  p_ = p;
  pd_ = parallel_decomp;
  cedr_throw_if(nranks_ > nc_, "#GIDs < #ranks is not supported.");
}

Int Mesh::rank (const Int& ci) const {
  switch (pd_) {
  case ParallelDecomp::contiguous:
    return std::min(nranks_ - 1, ci / (nc_ / nranks_));
  default: {
    const auto chunk = ci / nranks_;
    return (ci + chunk) % nranks_;
  }
  }
}

Int Mesh::unittest (const Parallel::Ptr& p) {
  const Mesh::ParallelDecomp::Enum dists[] = { Mesh::ParallelDecomp::pseudorandom,
                                               Mesh::ParallelDecomp::contiguous };
  Int ne = 0;
  for (size_t id = 0; id < sizeof(dists)/sizeof(*dists); ++id) {
    Mesh m(std::max(42, 3*p->size()), p, dists[id]);
    const Int nc = m.ncell();
    for (Int ci = 0; ci < nc; ++ci)
      if (m.rank(ci) < 0 || m.rank(ci) >= p->size())
        ++ne;
  }
  return ne;
}

tree::Node::Ptr make_tree (const Mesh& m, const Int cs, const Int ce,
                           const tree::Node* parent, const bool imbalanced) {
  const Int
    cn = ce - cs,
    cn0 = ( imbalanced && cn > 2 ?
            cn/3 :
            cn/2 );
  tree::Node::Ptr n = std::make_shared<tree::Node>();
  n->parent = parent;
  if (cn == 1) {
    n->nkids = 0;
    n->rank = m.rank(cs);
    n->cellidx = cs;
    return n;
  }
  n->nkids = 2;
  n->kids[0] = make_tree(m, cs, cs + cn0, n.get(), imbalanced);
  n->kids[1] = make_tree(m, cs + cn0, ce, n.get(), imbalanced);
  return n;
}

tree::Node::Ptr make_tree (const Mesh& m, const bool imbalanced) {
  return make_tree(m, 0, m.ncell(), nullptr, imbalanced);
}

tree::Node::Ptr make_tree (const Parallel::Ptr& p, const Int& ncells,
                           const bool imbalanced) {
  Mesh m(ncells, p);
  return make_tree(m, imbalanced);
}

void mark_cells (const tree::Node::Ptr& node, std::vector<Int>& cells) {
  if ( ! node->nkids) {
    ++cells[node->cellidx];
    return;
  }
  for (Int i = 0; i < node->nkids; ++i)
    mark_cells(node->kids[i], cells);
}

Int unittest (const Parallel::Ptr& p) {
  const Mesh::ParallelDecomp::Enum dists[] = { Mesh::ParallelDecomp::pseudorandom,
                                               Mesh::ParallelDecomp::contiguous };
  Int ne = 0;
  for (size_t id = 0; id < sizeof(dists)/sizeof(*dists); ++id)
    for (bool imbalanced: {false, true}) {
      Mesh m(std::max(42, 3*p->size()), p, Mesh::ParallelDecomp::pseudorandom);
      tree::Node::Ptr tree = make_tree(m, imbalanced);
      std::vector<Int> cells(m.ncell(), 0);
      mark_cells(tree, cells);
      for (Int i = 0; i < m.ncell(); ++i)
        if (cells[i] != 1) ++ne;
    }
  return ne;
}
} // namespace oned

// Unit tests for NodeSets.
Int unittest_NodeSets (const Parallel::Ptr& p, const NodeSets::ConstPtr& ns,
                       const Int ncells) {
  Int nerr = 0, ne;
  ne = check_comm(*ns);
  if (ne && p->amroot()) pr("check_comm failed");
  nerr += ne;
  ne = check_leaf_nodes(p, *ns, ncells);
  if (ne && p->amroot()) pr("check_leaf_nodes failed");
  nerr += ne;
  for (Int trial = 0; trial < 11; ++trial) {
    ne = test_comm_pattern(p, *ns, ncells);
    if (ne && p->amroot()) pr("test_comm_pattern failed for trial" pu(trial));
    nerr += ne;
  }
  return nerr;
}

Int unittest (const Parallel::Ptr& p) {
  using Mesh = oned::Mesh;
  const Int szs[] = { p->size(), 3*p->size() };
  const Mesh::ParallelDecomp::Enum dists[] = { Mesh::ParallelDecomp::pseudorandom,
                                               Mesh::ParallelDecomp::contiguous };
  Int nerr = 0;
  for (size_t is = 0; is < sizeof(szs)/sizeof(*szs); ++is)
    for (size_t id = 0; id < sizeof(dists)/sizeof(*dists); ++id)
      for (bool imbalanced: {false, true}) {
        Mesh m(szs[is], p, dists[id]);
        tree::Node::Ptr tree = make_tree(m, imbalanced);
        tree::NodeSets::ConstPtr nodesets = analyze(p, m.ncell(), tree);
        tree = nullptr;
        nerr += unittest_NodeSets(p, nodesets, m.ncell());
      }
  return nerr;
}

} // namespace tree
} // namespace cedr
