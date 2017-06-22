#include <sys/time.h>
#include <mpi.h>
#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>
#include <map>
#include <set>
#include <list>
#include <vector>
#include <memory>
#include <limits>
#include <algorithm>

#include <Kokkos_Core.hpp>

// QLT: Quasi-local tree-based non-iterative tracer density reconstructor for
//      mass conservation, shape preservation, and tracer consistency.
//
// Implementation for use in Homme, as wrapped by compose_mod.F90.
//
// Build standalone with -DQLT_MAIN to make a program that runs correctness and
// performance tests.
//
// This implementation is intended to produce bit-for-bit output values
// independent of parallel decomposition. Compile with -DQLT_FASTEST to drop
// this constraint.
namespace qlt {
typedef int Int;
typedef size_t Size;
typedef double Real;

class Parallel {
  MPI_Comm comm_;
public:
  typedef std::shared_ptr<Parallel> Ptr;
  Parallel(MPI_Comm comm) : comm_(comm) {}
  MPI_Comm comm () const { return comm_; }
  Int size() const;
  Int rank() const;
  Int root () const { return 0; }
  bool amroot () const { return rank() == root(); }
};

namespace tree {
// The caller builds a tree of these nodes to pass to QLT.
struct Node {
  typedef std::shared_ptr<Node> Ptr;
  const Node* parent; // (Can't be a shared_ptr: would be a circular dependency.)
  Int rank;           // Owning rank.
  Int cellidx;        // If a leaf, the cell to which this node corresponds.
  Int nkids;          // 0 at leaf, 1 or 2 otherwise.
  Node::Ptr kids[2];
  void* reserved;     // For internal use.
  Node () : parent(nullptr), rank(-1), cellidx(-1), nkids(0), reserved(nullptr) {}
};
} // namespace tree

namespace mpi {
template <typename T> MPI_Datatype get_type();
template <> MPI_Datatype get_type<int>() { return MPI_INT; }
template <> MPI_Datatype get_type<double>() { return MPI_DOUBLE; }

template <typename T>
int reduce (const Parallel& p, const T* sendbuf, T* rcvbuf, int count, MPI_Op op,
            int root) {
  MPI_Datatype dt = get_type<T>();
  return MPI_Reduce(const_cast<T*>(sendbuf), rcvbuf, count, dt, op, root, p.comm());
}

template <typename T>
int all_reduce (const Parallel& p, const T* sendbuf, T* rcvbuf, int count, MPI_Op op) {
  MPI_Datatype dt = get_type<T>();
  return MPI_Allreduce(const_cast<T*>(sendbuf), rcvbuf, count, dt, op, p.comm());
}

template <typename T>
int isend (const Parallel& p, const T* buf, int count, int dest, int tag,
           MPI_Request* ireq) {
  MPI_Datatype dt = get_type<T>();
  MPI_Request ureq;
  MPI_Request* req = ireq ? ireq : &ureq;
  int ret = MPI_Isend(const_cast<T*>(buf), count, dt, dest, tag, p.comm(), req);
  if ( ! ireq) MPI_Request_free(req);
  return ret;
}

template <typename T>
int irecv (const Parallel& p, T* buf, int count, int src, int tag, MPI_Request* ireq) {
  MPI_Datatype dt = get_type<T>();
  MPI_Request ureq;
  MPI_Request* req = ireq ? ireq : &ureq;
  int ret = MPI_Irecv(buf, count, dt, src, tag, p.comm(), req);
  if ( ! ireq) MPI_Request_free(req);
  return ret;
}

int waitany (int count, MPI_Request* reqs, int* index, MPI_Status* stats = nullptr) {
  return MPI_Waitany(count, reqs, index, stats ? stats : MPI_STATUS_IGNORE);
}

int waitall (int count, MPI_Request* reqs, MPI_Status* stats = nullptr) {
  return MPI_Waitall(count, reqs, stats ? stats : MPI_STATUS_IGNORE);
}

template<typename T>
int gather (const Parallel& p, const T* sendbuf, int sendcount,
            T* recvbuf, int recvcount, int root) {
  MPI_Datatype dt = get_type<T>();
  return MPI_Gather(sendbuf, sendcount, dt, recvbuf, recvcount, dt, root, p.comm());
}

template <typename T>
int gatherv (const Parallel& p, const T* sendbuf, int sendcount,
             T* recvbuf, const int* recvcounts, const int* displs, int root) {
  MPI_Datatype dt = get_type<T>();
  return MPI_Gatherv(sendbuf, sendcount, dt, recvbuf, recvcounts, displs, dt, root,
                     p.comm());
}

bool all_ok (const Parallel& p, bool im_ok) {
  int ok = im_ok, msg;
  all_reduce<int>(p, &ok, &msg, 1, MPI_LAND);
  return static_cast<bool>(msg);
}
} // namespace mpi

Parallel::Ptr make_parallel (MPI_Comm comm) {
  return std::make_shared<Parallel>(comm);
}

Int Parallel::size () const {
  int sz = 0;
  MPI_Comm_size(comm_, &sz);
  return sz;
}

Int Parallel::rank () const {
  int pid = 0;
  MPI_Comm_rank(comm_, &pid);
  return pid;
}

namespace impl {
#define pr(m) do {                                  \
    int _pid_ = 0;                                      \
    MPI_Comm_rank(MPI_COMM_WORLD, &_pid_);              \
    std::stringstream _ss_;                             \
    _ss_.precision(15);                                 \
    _ss_ << "pid " << _pid_ << " " << m << std::endl;   \
    std::cerr << _ss_.str();                            \
  } while (0)
#define pr0(m) do {                                 \
    int _pid_; MPI_Comm_rank(MPI_COMM_WORLD, &_pid_);   \
    if (_pid_ != 0) break;                              \
    std::stringstream _ss_;                             \
    _ss_ << "pid " << _pid_ << " " << m << std::endl;   \
    std::cerr << _ss_.str();                            \
  } while (0)
#define prc(m) pr(#m << " | " << (m))
#define pr0c(m) pr0(#m << " | " << (m))
#define puf(m) "(" << #m << " " << (m) << ")"
#define pu(m) << " " << puf(m)
template <typename T>
void prarr (const std::string& name, const T* const v, const size_t n) {
  std::stringstream ss;
  ss.precision(15);
  ss << name << " = [";
  for (size_t i = 0; i < n; ++i) ss << " " << v[i];
  ss << "];";
  pr(ss.str());
}
#define mprarr(m) qlt::impl::prarr(#m, m.data(), m.size())

#define qlt_assert(condition) do {                                      \
    if ( ! (condition)) {                                               \
      std::stringstream _ss_;                                           \
      _ss_ << __FILE__ << ":" << __LINE__ << ": FAIL:\n" << #condition  \
        << "\n";                                                        \
        throw std::logic_error(_ss_.str());                             \
    }                                                                   \
  } while (0)
#define qlt_throw_if(condition, message) do {                           \
    if (condition) {                                                    \
      std::stringstream _ss_;                                           \
      _ss_ << __FILE__ << ":" << __LINE__ << ": The condition:\n"       \
           << #condition "\nled to the exception\n" << message << "\n"; \
        throw std::logic_error(_ss_.str());                             \
    }                                                                   \
  } while (0)
#define qlt_kernel_assert(condition) do {       \
    if ( ! (condition))                         \
      Kokkos::abort(#condition);                \
  } while (0)
#define qlt_kernel_throw_if(condition, message) do {                    \
    if (condition)                                                      \
      Kokkos::abort(#condition " led to the exception\n" message);      \
  } while (0)

// GPU-friendly replacements for std::min/max.
template <typename T> KOKKOS_INLINE_FUNCTION
const T& min (const T& a, const T& b) { return a < b ? a : b; }
template <typename T> KOKKOS_INLINE_FUNCTION
const T& max (const T& a, const T& b) { return a > b ? a : b; }

inline Real reldif (const Real a, const Real b)
{ return std::abs(b - a)/std::max(std::abs(a), std::abs(b)); }

struct FILECloser { void operator() (FILE* fh) { fclose(fh); } };

template <typename MemoryTraitsType, Kokkos::MemoryTraitsFlags flag>
using MemoryTraits = Kokkos::MemoryTraits<
  MemoryTraitsType::Unmanaged | MemoryTraitsType::RandomAccess |
  MemoryTraitsType::Atomic | flag>;

template <typename View>
using Unmanaged = Kokkos::View<
  typename View::data_type, typename View::array_layout,
  typename View::device_type, MemoryTraits<typename View::memory_traits,
                                           Kokkos::Unmanaged> >;
template <typename View>
using Const = Kokkos::View<
  typename View::const_data_type, typename View::array_layout,
  typename View::device_type, typename View::memory_traits>;
template <typename View>
using ConstUnmanaged = Const<Unmanaged<View> >;

template <typename ExeSpace>
struct DeviceType {
  typedef Kokkos::Device<typename ExeSpace::execution_space,
                         typename ExeSpace::memory_space> type;
};

#ifdef KOKKOS_HAVE_CUDA
// Because Tpetra and others in Trilinos use UVM and we do not want to, we need
// to be very precise about the (execution space, memory space) pair we want.
typedef Kokkos::Device<Kokkos::CudaSpace::execution_space,
                       Kokkos::CudaSpace::memory_space> DefaultDeviceType;

template <> struct DeviceType<Kokkos::Cuda> {
  typedef DefaultDeviceType type;
};
#else
typedef Kokkos::Device<Kokkos::DefaultExecutionSpace::execution_space,
                       Kokkos::DefaultExecutionSpace::memory_space> DefaultDeviceType;
#endif
} // namespace impl

class Timer {
public:
  enum Op { tree, analyze, trcrinit, trcrgen, trcrcheck,
            qltrun, qltrunl2r, qltrunr2l, snp, waitall,
            total, NTIMERS };
  static inline void init () {
#ifdef QLT_TIME
    for (int i = 0; i < NTIMERS; ++i) {
      et_[i] = 0;
      cnt_[i] = 0;
    }
#endif
  }
  static inline void reset (const Op op) {
#ifdef QLT_TIME
    et_[op] = 0;
    cnt_[op] = 0;
#endif
  }
  static inline void start (const Op op) {
#ifdef QLT_TIME
    gettimeofday(&t_start_[op], 0);
    ++cnt_[op];
#endif
  }
  static inline void stop (const Op op) {
#ifdef QLT_TIME
    timeval t2;
    gettimeofday(&t2, 0);
    const timeval& t1 = t_start_[op];
    static const double us = 1.0e6;
    et_[op] += (t2.tv_sec*us + t2.tv_usec - t1.tv_sec*us - t1.tv_usec)/us;
#endif
  }
# define tpr(op) do {                                                   \
    printf("%-20s %10.3e %10.1f (%4d %10.3e)\n",                        \
           #op, et_[op], 100*et_[op]/tot, cnt_[op], et_[op]/cnt_[op]);  \
  } while (0)
  static void print () {
#ifdef QLT_TIME
    const double tot = et_[total];
    tpr(tree); tpr(analyze);
    tpr(trcrinit); tpr(trcrgen); tpr(trcrcheck);
    tpr(qltrun); tpr(qltrunl2r); tpr(qltrunr2l); tpr(snp); tpr(waitall);
    printf("%-20s %10.3e %10.1f\n", "total", tot, 100.0);
#endif
  }
#undef tpr
private:
#ifdef QLT_TIME
  static timeval t_start_[NTIMERS];
  static double et_[NTIMERS];
  static int cnt_[NTIMERS];
#endif
};
#ifdef QLT_TIME
timeval Timer::t_start_[Timer::NTIMERS];
double Timer::et_[Timer::NTIMERS];
int Timer::cnt_[Timer::NTIMERS];
#endif

namespace impl {
struct NodeSets {
  typedef std::shared_ptr<const NodeSets> ConstPtr;
  
  enum : int { mpitag = 42 };

  // A node in the tree that is relevant to this rank.
  struct Node {
    // Rank of the node. If the node is in a level, then its rank is my rank. If
    // it's not in a level, then it is a comm partner of a node on this rank.
    Int rank;
    // Globally unique identifier; cellidx if leaf node, ie, if nkids == 0.
    Int id;
    // This node's parent, a comm partner, if such a partner is required.
    const Node* parent;
    // This node's kids, comm partners, if such partners are required. Parent
    // and kid nodes are pruned relative to the full tree over the mesh to
    // contain just the nodes that matter to this rank.
    Int nkids;
    const Node* kids[2];
    // Offset factor into bulk data. An offset is a unit; actual buffer sizes
    // are multiples of this unit.
    Int offset;

    Node () : rank(-1), id(-1), parent(nullptr), nkids(0), offset(-1) {}
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
    std::vector<Node*> nodes;
    // MPI information for this level.
    std::vector<MPIMetaData> me, kids;
    // Have to keep requests separate so we can call waitall if we want to.
    mutable std::vector<MPI_Request> me_req, kids_req;
  };
  
  // Levels. nodes[0] is level 0, the leaf level.
  std::vector<Level> levels;
  // Number of data slots this rank needs. Each node owned by this rank, plus
  // kids on other ranks, have an associated slot.
  Int nslots;
  
  // Allocate a node. The list node_mem_ is the mechanism for memory ownership;
  // node_mem_ isn't used for anything other than owning nodes.
  Node* alloc () {
    node_mem_.push_front(Node());
    return &node_mem_.front();
  }

  void print(std::ostream& os) const;
  
private:
  std::list<Node> node_mem_;
};

void NodeSets::print (std::ostream& os) const {
  std::stringstream ss;
  if (levels.empty()) return;
  const Int myrank = levels[0].nodes[0]->rank;
  ss << "pid " << myrank << ":";
  ss << " #levels " << levels.size();
  for (size_t i = 0; i < levels.size(); ++i) {
    const auto& lvl = levels[i];
    ss << "\n  " << i << ": " << lvl.nodes.size();
    std::set<Int> ps, ks;
    for (size_t j = 0; j < lvl.nodes.size(); ++j) {
      const auto n = lvl.nodes[j];
      for (Int k = 0; k < n->nkids; ++k)
        if (n->kids[k]->rank != myrank)
          ks.insert(n->kids[k]->rank);
      if (n->parent && n->parent->rank != myrank)
        ps.insert(n->parent->rank);
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
Int init_tree (const tree::Node::Ptr& node, Int& id) {
  node->reserved = nullptr;
  Int depth = 0;
  for (Int i = 0; i < node->nkids; ++i) {
    qlt_assert(node.get() == node->kids[i]->parent);
    depth = std::max(depth, init_tree(node->kids[i], id));
  }
  if (node->nkids) {
    node->rank = node->kids[0]->rank;
    node->cellidx = id++;
  } else {
    qlt_throw_if(node->cellidx < 0 || node->cellidx >= id,
                 "cellidx is " << node->cellidx << " but should be between " <<
                 0 << " and " << id);
  }
  return depth + 1;
}

void level_schedule_and_collect (
  NodeSets& ns, const Int& my_rank, const tree::Node::Ptr& node, Int& level,
  bool& need_parent_ns_node)
{
  qlt_assert(node->rank != -1);
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
  // Is parent node needed for isend?
  const bool node_is_owned = node->rank == my_rank;
  need_parent_ns_node = node_is_owned;
  if (node_is_owned || make_ns_node) {
    qlt_assert( ! node->reserved);
    NodeSets::Node* ns_node = ns.alloc();
    // Levels hold only owned nodes.
    if (node_is_owned) ns.levels[level].nodes.push_back(ns_node);
    node->reserved = ns_node;
    ns_node->rank = node->rank;
    ns_node->id = node->cellidx;
    ns_node->parent = nullptr;
    if (node_is_owned) {
      // If this node is owned, it needs to have information about all kids.
      ns_node->nkids = node->nkids;
      for (Int i = 0; i < node->nkids; ++i) {
        const auto& kid = node->kids[i];
        if ( ! kid->reserved) {
          // This kid isn't owned by this rank. But need it for irecv.
          NodeSets::Node* ns_kid;
          kid->reserved = ns_kid = ns.alloc();
          ns_node->kids[i] = ns_kid;
          qlt_assert(kid->rank != my_rank);
          ns_kid->rank = kid->rank;
          ns_kid->id = kid->cellidx;
          ns_kid->parent = nullptr; // Not needed.
          // The kid may have kids in the original tree, but in the tree pruned
          // according to rank, it does not.
          ns_kid->nkids = 0;
        } else {
          // This kid is owned by this rank, so fill in its parent pointer.
          NodeSets::Node* ns_kid = static_cast<NodeSets::Node*>(kid->reserved);
          ns_node->kids[i] = ns_kid;
          ns_kid->parent = ns_node;
        }
      }
    } else {
      // This node is not owned. Update the owned kids with its parent.
      ns_node->nkids = 0;
      for (Int i = 0; i < node->nkids; ++i) {
        const auto& kid = node->kids[i];
        if (kid->reserved && kid->rank == my_rank) {
          NodeSets::Node* ns_kid = static_cast<NodeSets::Node*>(kid->reserved);
          ns_node->kids[ns_node->nkids++] = ns_kid;
          ns_kid->parent = ns_node;
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
      qlt_assert(rank > prev_rank);
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
    for (const auto& n : lvl.nodes)
      nkids += n->nkids;

    std::vector<RankNode> me(lvl.nodes.size()), kids(nkids);
    for (size_t i = 0, mi = 0, ki = 0; i < lvl.nodes.size(); ++i) {
      const auto& n = lvl.nodes[i];
      me[mi].first = n->parent ? n->parent->rank : my_rank;
      me[mi].second = const_cast<NodeSets::Node*>(n);
      ++mi;
      for (Int k = 0; k < n->nkids; ++k) {
        kids[ki].first = n->kids[k]->rank;
        kids[ki].second = const_cast<NodeSets::Node*>(n->kids[k]);
        ++ki;
      }
    }

    init_offsets(my_rank, me, lvl.me, ns.nslots);
    lvl.me_req.resize(lvl.me.size());
    init_offsets(my_rank, kids, lvl.kids, ns.nslots);
    lvl.kids_req.resize(lvl.kids.size());
  }
}

// Analyze the tree to extract levels. Levels are run from 0 to #level - 1. Each
// level has nodes whose corresponding operations depend on only nodes in
// lower-indexed levels. This mechanism prevents deadlock in the general case of
// multiple cells per rank, with multiple ranks appearing in a subtree other
// than the root.
//   In addition, the set of nodes collected into levels are just those owned by
// this rank, and those with which owned nodes must communicate.
//   Once this function is done, the tree can be deleted.
NodeSets::ConstPtr analyze (const Parallel::Ptr& p, const Int& ncells,
                            const tree::Node::Ptr& tree) {
  const auto nodesets = std::make_shared<NodeSets>();
  qlt_assert( ! tree->parent);
  Int id = ncells;
  const Int depth = init_tree(tree, id);
  nodesets->levels.resize(depth);
  level_schedule_and_collect(*nodesets, p->rank(), tree);
  consolidate(*nodesets);
  init_comm(p->rank(), *nodesets);
  return nodesets;
}

// Check that the offsets are self consistent.
Int check_comm (const NodeSets::ConstPtr& ns) {
  Int nerr = 0;
  std::vector<Int> offsets(ns->nslots, 0);
  for (const auto& lvl : ns->levels)
    for (const auto& n : lvl.nodes) {
      qlt_assert(n->offset < ns->nslots);
      ++offsets[n->offset];
      for (Int i = 0; i < n->nkids; ++i)
        if (n->kids[i]->rank != n->rank)
          ++offsets[n->kids[i]->offset];
    }
  for (const auto& e : offsets)
    if (e != 1) ++nerr;
  return nerr;
}

// Check that there are the correct number of leaf nodes, and that their offsets
// all come first and are ordered the same as ns->levels[0]->nodes.
Int check_leaf_nodes (const Parallel::Ptr& p, const NodeSets::ConstPtr& ns,
                      const Int ncells) {
  Int nerr = 0;
  qlt_assert( ! ns->levels.empty());
  qlt_assert( ! ns->levels[0].nodes.empty());
  Int my_nleaves = 0;
  for (const auto& n : ns->levels[0].nodes) {
    qlt_assert( ! n->nkids);
    ++my_nleaves;
  }
  for (const auto& n : ns->levels[0].nodes) {
    qlt_assert(n->offset < my_nleaves);
    qlt_assert(n->id < ncells);
  }
  Int glbl_nleaves = 0;
  mpi::all_reduce(*p, &my_nleaves, &glbl_nleaves, 1, MPI_SUM);
  if (glbl_nleaves != ncells)
    ++nerr;
  return nerr;
}

// Sum cellidx using the QLT comm pattern.
Int test_comm_pattern (const Parallel::Ptr& p, const NodeSets::ConstPtr& ns,
                       const Int ncells) {
  Int nerr = 0;
  // Rank-wide data buffer.
  std::vector<Int> data(ns->nslots);
  // Sum this rank's cellidxs.
  for (auto& n : ns->levels[0].nodes)
    data[n->offset] = n->id;
  // Leaves to root.
  for (size_t il = 0; il < ns->levels.size(); ++il) {
    auto& lvl = ns->levels[il];
    // Set up receives.
    for (size_t i = 0; i < lvl.kids.size(); ++i) {
      const auto& mmd = lvl.kids[i];
      mpi::irecv(*p, &data[mmd.offset], mmd.size, mmd.rank, NodeSets::mpitag,
                 &lvl.kids_req[i]);
    }
    //todo Replace with simultaneous waitany and isend.
    mpi::waitall(lvl.kids_req.size(), lvl.kids_req.data());
    // Combine kids' data.
    for (auto& n : lvl.nodes) {
      if ( ! n->nkids) continue;
      data[n->offset] = 0;
      for (Int i = 0; i < n->nkids; ++i)
        data[n->offset] += data[n->kids[i]->offset];
    }
    // Send to parents.
    for (size_t i = 0; i < lvl.me.size(); ++i) {
      const auto& mmd = lvl.me[i];
      mpi::isend(*p, &data[mmd.offset], mmd.size, mmd.rank, NodeSets::mpitag,
                 &lvl.me_req[i]);
    }
    if (il+1 == ns->levels.size())
      mpi::waitall(lvl.me_req.size(), lvl.me_req.data());
  }
  // Root to leaves.
  for (size_t il = ns->levels.size(); il > 0; --il) {
    auto& lvl = ns->levels[il-1];
    // Get the global sum from parent.
    for (size_t i = 0; i < lvl.me.size(); ++i) {
      const auto& mmd = lvl.me[i];
      mpi::irecv(*p, &data[mmd.offset], mmd.size, mmd.rank, NodeSets::mpitag,
                 &lvl.me_req[i]);
    }    
    //todo Replace with simultaneous waitany and isend.
    mpi::waitall(lvl.me_req.size(), lvl.me_req.data());
    // Pass to kids.
    for (auto& n : lvl.nodes) {
      if ( ! n->nkids) continue;
      for (Int i = 0; i < n->nkids; ++i)
        data[n->kids[i]->offset] = data[n->offset];
    }
    // Send.
    for (size_t i = 0; i < lvl.kids.size(); ++i) {
      const auto& mmd = lvl.kids[i];
      mpi::isend(*p, &data[mmd.offset], mmd.size, mmd.rank, NodeSets::mpitag,
                 &lvl.kids_req[i]);
    }
  }
  // Wait on sends to clean up.
  for (size_t il = 0; il < ns->levels.size(); ++il) {
    auto& lvl = ns->levels[il];
    if (il+1 < ns->levels.size())
      mpi::waitall(lvl.me_req.size(), lvl.me_req.data());
    mpi::waitall(lvl.kids_req.size(), lvl.kids_req.data());
  }
  { // Check that all leaf nodes have the right number.
    const Int desired_sum = (ncells*(ncells - 1)) / 2;
    for (const auto& n : ns->levels[0].nodes)
      if (data[n->offset] != desired_sum) ++nerr;
    if (p->amroot()) {
      std::cout << " " << data[ns->levels[0].nodes[0]->offset];
      std::cout.flush();
    }
  }
  return nerr;
}

// Unit tests for NodeSets.
Int unittest (const Parallel::Ptr& p, const NodeSets::ConstPtr& ns,
              const Int ncells) {
  Int nerr = 0;
  nerr += check_comm(ns);
  if (nerr) return nerr;
  nerr += check_leaf_nodes(p, ns, ncells);
  if (nerr) return nerr;
  nerr += test_comm_pattern(p, ns, ncells);
  if (nerr) return nerr;
  return nerr;
}
} // namespace impl

namespace slv {
KOKKOS_INLINE_FUNCTION
Real get_xbd (const Real* xbd, const Int i, const bool xbds_scalar)
{ return xbds_scalar ? *xbd : xbd[i]; }

KOKKOS_INLINE_FUNCTION
bool is_inside (const Real xi, const Real* xlo, const Real* xhi, const Int i,
                const bool xbds_scalar) {
  return (xi > get_xbd(xlo, i, xbds_scalar) &&
          xi < get_xbd(xhi, i, xbds_scalar));
}

KOKKOS_INLINE_FUNCTION
bool is_outside (const Real xi, const Real* xlo, const Real* xhi, const Int i,
                 const bool xbds_scalar) {
  return (xi < get_xbd(xlo, i, xbds_scalar) ||
          xi > get_xbd(xhi, i, xbds_scalar));
}

KOKKOS_INLINE_FUNCTION
Real calc_r_tol (const Real b, const Real* a, const Real* y, const Int n) {
  Real ab = std::abs(b);
  for (Int i = 0; i < n; ++i) ab = std::max(ab, std::abs(a[i]*y[i]));
  return 1e1*std::numeric_limits<Real>::epsilon()*std::abs(ab);
}

KOKKOS_INLINE_FUNCTION
void calc_r (const Int n, const Real* w, const Real* a, const Real b,
             const Real* xlo, const Real* xhi, const bool xbds_scalar,
             const Real* y, const Real& lambda, Real* x, Real& r, Real& r_lambda) {
  r = 0;
  r_lambda = 0;
  for (Int i = 0; i < n; ++i) {
    const Real q = a[i]/w[i];
    const Real x_trial = y[i] + lambda*q;
    Real xtmp;
    if (x_trial < (xtmp = get_xbd(xlo, i, xbds_scalar)))
      x[i] = xtmp;
    else if (x_trial > (xtmp = get_xbd(xhi, i, xbds_scalar)))
      x[i] = xtmp;
    else {
      x[i] = x_trial;
      r_lambda += a[i]*q;
    }
    r += a[i]*x[i];
  }
  r -= b;
}

// Solve
//     min_x sum_i w(i) (x(i) - y(i))^2
//      st   a' x = b
//           xlo <= x <= xhi.
// This function assumes w > 0 to save a few operations. Return 0 on success and
// x == y, 1 on success and x != y, -1 if infeasible, -2 if max_its hit with no
// solution. See Section 3 of Bochev, Ridzal, Shashkov, Fast optimization-based
// conservative remap of scalar fields through aggregate mass transfer. lambda
// is used in check_1eq_bc_qp_foc.
//todo 2D version of this function that takes advantage of 2D.
KOKKOS_FUNCTION
Int solve_1eq_bc_qp (const Int n, const Real* w, const Real* a, const Real b,
                     const Real* xlo, const Real* xhi, const bool xbds_scalar,
                     const Real* y, Real* x, const Int max_its = 100) {
  const Real r_tol = calc_r_tol(b, a, y, n);

  { // Check for a quick exit.
    bool all_in = true;
    Real r = 0;
    for (Int i = 0; i < n; ++i) {
      if (is_outside(x[i], xlo, xhi, i, xbds_scalar)) {
        all_in = false;
        break;
      }
      r += a[i]*x[i];
    }
    if (all_in) {
      r -= b;
      if (std::abs(r) <= r_tol)
        return 0;
    }
  }

  { // Eval r at end points to check for feasibility, and also possibly a quick
    // exit on a common case.
    Real r = -b;
    for (Int i = 0; i < n; ++i) {
      x[i] = get_xbd(xlo, i, xbds_scalar);
      r += a[i]*x[i];
    }
    if (std::abs(r) <= r_tol) return 1;
    if (r > 0) return -1;
    r = -b;
    for (Int i = 0; i < n; ++i) {
      x[i] = get_xbd(xhi, i, xbds_scalar);
      r += a[i]*x[i];
    }
    if (std::abs(r) <= r_tol) return 1;
    if (r < 0) return -1;
  }

  { // Check for a quick exit: the bounds are so tight that the midpoint of the
    // box satisfies r_tol.
    Real r = -b;
    for (Int i = 0; i < n; ++i) {
      x[i] = 0.5*(get_xbd(xlo, i, xbds_scalar) + get_xbd(xhi, i, xbds_scalar));
      r += a[i]*x[i];
    }
    if (std::abs(r) <= r_tol) return 1;
  }

  const Real wall_dist = 1e-3;

  // Get lambda endpoints.
  Real lamlo = 0, lamhi = 0;
  for (Int i = 0; i < n; ++i) {
    const Real rq = w[i]/a[i];
    const Real lamlo_i = rq*(get_xbd(xlo, i, xbds_scalar) - y[i]);
    const Real lamhi_i = rq*(get_xbd(xhi, i, xbds_scalar) - y[i]);
    if (i == 0) {
      lamlo = lamlo_i;
      lamhi = lamhi_i;
    } else {
      lamlo = std::min(lamlo, lamlo_i);
      lamhi = std::max(lamhi, lamhi_i);
    }
  }
  const Real lamlo_feas = lamlo, lamhi_feas = lamhi;
  Real lambda = lamlo <= 0 && lamhi >= 0 ? 0 : lamlo;

  Int info = -2;

  // Bisection-safeguarded Newton iteration for r(lambda) = 0.
  bool prev_step_bisect = false;
  Int nbisect = 0;
  for (Int iteration = 0; iteration < max_its; ++iteration) {
    // Compute x, r, r_lambda.
    Real r, r_lambda;
    calc_r(n, w, a, b, xlo, xhi, xbds_scalar, y, lambda, x, r, r_lambda);
    // Is r(lambda) - b sufficiently == 0?
    if (std::abs(r) <= r_tol) {
      info = 1;
      break;
    }
    // Check if the lambda bounds are too close.
    if (nbisect > 64) {
      if (lamhi == lamhi_feas || lamlo == lamlo_feas) {
        // r isn't small enough and one lambda bound is on the feasibility
        // limit. The QP must not be feasible.
        info = -1;
        break;
      }
      info = 1;
      break;
    }
    // Adjust lambda bounds.
    if (r > 0)
      lamhi = lambda;
    else
      lamlo = lambda;
    if (r_lambda != 0) {
      // Newton step.
      lambda -= r/r_lambda;
    } else {
      // Force bisection.
      lambda = lamlo;
    }
    // Safeguard. The wall distance check assures progress, but use it only
    // every other potential bisection.
    const Real D = prev_step_bisect ? 0 : wall_dist*(lamhi - lamlo);
    if (lambda - lamlo < D || lamhi - lambda < D) {
      lambda = 0.5*(lamlo + lamhi);
      ++nbisect;
      prev_step_bisect = true;
    } else {
      prev_step_bisect = false;
    }
  }

  return info;
}

KOKKOS_FUNCTION
void r2l_nl_adjust_bounds (Real Qm_bnd[2], const Real rhom[2], Real Qm_extra) {
  Real q[2];
  for (Int i = 0; i < 2; ++i) q[i] = Qm_bnd[i] / rhom[i];
  if (Qm_extra < 0) {
    Int i0, i1;
    if (q[0] >= q[1]) { i0 = 0; i1 = 1; } else { i0 = 1; i1 = 0; }
    const Real Qm_gap = (q[i1] - q[i0])*rhom[i0];
    if (Qm_gap <= Qm_extra) {
      Qm_bnd[i0] += Qm_extra;
      return;
    }
  } else {
    Int i0, i1;
    if (q[0] <= q[1]) { i0 = 0; i1 = 1; } else { i0 = 1; i1 = 0; }
    const Real Qm_gap = (q[i1] - q[i0])*rhom[i0];
    if (Qm_gap >= Qm_extra) {
      Qm_bnd[i0] += Qm_extra;
      return;
    }
  }
  { // Have to adjust both. Adjust so that the q bounds are the same. This
    // procedure assures that as long as rhom is conservative, then the
    // adjustment never pushes q_{min,max} out of the safety bounds.
    const Real Qm_tot = Qm_bnd[0] + Qm_bnd[1] + Qm_extra;
    const Real rhom_tot = rhom[0] + rhom[1];
    const Real q_tot = Qm_tot / rhom_tot;
    for (Int i = 0; i < 2; ++i)
      Qm_bnd[i] = q_tot*rhom[i];
  }
}

KOKKOS_FUNCTION
void r2l_l_adjust_bounds (const Int np, Real* q_min, Real* q_max, const Real* rhom,
                          Real Qm_extra) {
  static constexpr int max_np = 16;
  Real* const q_bnd = Qm_extra < 0 ? q_min : q_max;
  // Try solving a QP that adjusts a q bound.
  Real Qm = Qm_extra;
  Real w[max_np], q_bnd_min[max_np], q_bnd_max[max_np], q_bnd_orig[max_np];
  q_bnd_min[0] = q_min[0];
  q_bnd_max[0] = q_max[0];
  for (Int i = 0; i < np; ++i) {
    const Real rhomi = rhom[i];
    Qm += q_bnd[i]*rhomi;
    q_bnd_orig[i] = q_bnd[i];
    w[i] = rhomi;
    if (Qm_extra < 0) {
      q_bnd_min[0] = std::min(q_bnd_min[0], q_min[i]);
      q_bnd_max[i] = q_max[i];
    } else {
      q_bnd_min[i] = q_min[i];
      q_bnd_max[0] = std::max(q_bnd_max[0], q_max[i]);
    }
  }
  if (Qm_extra < 0)
    for (Int i = 1; i < np; ++i) q_bnd_min[i] = q_bnd_min[0];
  else
    for (Int i = 1; i < np; ++i) q_bnd_max[i] = q_bnd_max[0];
  // Check for feasibility.
  bool feasible; {
    Real Qm_lo = 0, Qm_hi = 0;
    for (Int i = 0; i < np; ++i) {
      Qm_lo += q_bnd_min[i]*w[i];
      Qm_hi += q_bnd_max[i]*w[i];
    }
    feasible = Qm_lo <= Qm && Qm <= Qm_hi;
  }
  if (feasible) {
    solve_1eq_bc_qp(np, w, w, Qm, q_bnd_min, q_bnd_max, false, q_bnd_orig, q_bnd);
  } else {
    // The QP isn't feasible, so set the bound to a constant.
    Real rhom_tot = 0, Qm_tot = Qm_extra;
    for (Int i = 0; i < np; ++i) {
      const Real rhomi = rhom[i];
      rhom_tot += rhomi;
      Qm_tot += q_bnd_orig[i]*rhomi;
    }
    const Real q_tot = Qm_tot / rhom_tot;
    for (Int i = 0; i < np; ++i)
      q_bnd[i] = q_tot;
#pragma message "CHECK THIS 1"
    //return;
    // Assert that this constant is outside of all previous bound values. That's
    // why the QP wasn't feasible.
    if (Qm_extra < 0)
      for (Int i = 0; i < np; ++i)
        assert(q_tot <= q_bnd_orig[i]);
    else
      for (Int i = 0; i < np; ++i)
        assert(q_tot >= q_bnd_orig[i]);
  }
}

KOKKOS_FUNCTION
void solve_node_problem (const Real& rhom, const Real* pd, const Real& Qm,
                         const Real& rhom0, const Real* k0d, Real& Qm0,
                         const Real& rhom1, const Real* k1d, Real& Qm1) {
  Real Qm_min_kids [] = {k0d[0], k1d[0]};
  Real Qm_orig_kids[] = {k0d[1], k1d[1]};
  Real Qm_max_kids [] = {k0d[2], k1d[2]};
  { // Set the target values so that mass gets redistributed in a relative sense
    // rather than absolute. If a kid doesn't have much mass, don't give it too
    // much.
    const Real Qm_orig = pd[1], Qm_extra = Qm - Qm_orig;
    if (Qm_orig != 0)
      for (Int i = 0; i < 2; ++i)
        Qm_orig_kids[i] += (Qm_orig_kids[i] / Qm_orig) * Qm_extra;
  }
  { // The ideal problem is not assuredly feasible. Test for feasibility. If not
    // feasible, adjust bounds to solve the safety problem, which is assuredly
    // feasible if the total density field rho is mass conserving (Q doesn't
    // have to be mass conserving, of course; achieving mass conservation is one
    // use for QLT).
    const Real Qm_min = pd[0], Qm_max = pd[2];
    const bool lo = Qm < Qm_min, hi = Qm > Qm_max;
    if (lo || hi) {
      const Real rhom_kids[] = {rhom0, rhom1};
      r2l_nl_adjust_bounds(lo ? Qm_min_kids : Qm_max_kids,
                           rhom_kids,
                           Qm - (lo ? Qm_min : Qm_max));
    }
  }
  { // Solve the node's QP.
    static const Real ones[] = {1, 1};
    Real Qm_kids[2] = {k0d[1], k1d[1]};
    solve_1eq_bc_qp(2, ones, ones, Qm, Qm_min_kids, Qm_max_kids, false, Qm_orig_kids,
                    Qm_kids);
    Qm0 = Qm_kids[0];
    Qm1 = Qm_kids[1];
  }
}
} // namespace slv

template <typename ExeSpace>
class QLT {
public:
  typedef typename impl::DeviceType<ExeSpace>::type Device;
  typedef QLT<ExeSpace> Me;
  typedef std::shared_ptr<Me> Ptr;
  
  struct ProblemType {
    enum : Int { conserve = 1, shapepreserve = 1 << 1, consistent = 1 << 2 };
  };

private:
  typedef Kokkos::View<Int*, Kokkos::LayoutLeft, Device> IntList;
  typedef impl::Const<IntList> ConstIntList;
  typedef impl::ConstUnmanaged<IntList> ConstUnmanagedIntList;

  static void init (const std::string& name, IntList& d,
                    typename IntList::HostMirror& h, size_t n) {
    d = IntList(name, n);
    h = Kokkos::create_mirror_view(d);
  }

  struct MetaDataBuilder {
    typedef std::shared_ptr<MetaDataBuilder> Ptr;
    std::vector<int> trcr2prob;
  };

  struct MetaData {
    enum : Int { nprobtypes = 4 };

    template <typename IntListT>
    struct Arrays {
      // trcr2prob(i) is the ProblemType of tracer i.
      IntListT trcr2prob;
      // bidx2trcr(prob2trcrptr(i) : prob2trcrptr(i+1)-1) is the list of
      // tracers having ProblemType index i. bidx2trcr is the permutation
      // from the user's tracer index to the bulk data's ordering (bidx).
      Int prob2trcrptr[nprobtypes+1];
      IntListT bidx2trcr;
      // Inverse of bidx2trcr.
      IntListT trcr2bidx;
      // Points to the start of l2r bulk data for each problem type, within a
      // slot.
      Int prob2bl2r[nprobtypes + 1];
      // Point to the start of l2r bulk data for each tracer, within a slot.
      IntListT trcr2bl2r;
      // Same for r2l bulk data.
      Int prob2br2l[nprobtypes + 1];
      IntListT trcr2br2l;
    };

    static int get_problem_type (const int& idx) { return problem_type_[idx]; }
    
    // icpc doesn't let us use problem_type_ here, even though it's constexpr.
    KOKKOS_INLINE_FUNCTION
    static int get_problem_type_idx (const int& mask) {
      switch (mask) {
      case CPT::s:  case CPT::st:  return 0;
      case CPT::cs: case CPT::cst: return 1;
      case CPT::t:  return 2;
      case CPT::ct: return 3;
      default: qlt_kernel_throw_if(true, "Invalid problem type."); return -1;
      }
    }

    KOKKOS_INLINE_FUNCTION
    static int get_problem_type_l2r_bulk_size (const int& mask) {
      if (mask & ProblemType::conserve) return 4;
      return 3;
    }

    KOKKOS_INLINE_FUNCTION
    static int get_problem_type_r2l_bulk_size (const int& mask) {
      if (mask & ProblemType::shapepreserve) return 1;
      return 3;
    }

    struct CPT {
      // We could make the l2r buffer smaller by one entry, Qm. However, the
      // l2r comm is more efficient if it's done with one buffer. Similarly,
      // we separate the r2l data into a separate buffer for packing and MPI
      // efficiency.
      //   There are 7 possible problems.
      //   The only problem not supported is conservation alone. It makes very
      // little sense to use QLT for conservation alone.
      //   The remaining 6 fall into 4 categories of details. These 4 categories
      // are traceked by QLT; which of the original 6 problems being solved is
      // not important.
      enum {
        // l2r: rhom, (Qm_min, Qm, Qm_max)*; l2r, r2l: Qm*
        s  = ProblemType::shapepreserve,
        st = ProblemType::shapepreserve | ProblemType::consistent,
        // l2r: rhom, (Qm_min, Qm, Qm_max, Qm_prev)*; l2r, r2l: Qm*
        cs  = ProblemType::conserve | s,
        cst = ProblemType::conserve | st,
        // l2r: rhom, (q_min, Qm, q_max)*; l2r, r2l: Qm*
        t = ProblemType::consistent,
        // l2r: rhom, (q_min, Qm, q_max, Qm_prev)*; l2r, r2l: Qm*
        ct = ProblemType::conserve | t
      };
    };

    Arrays<typename ConstUnmanagedIntList::HostMirror> a_h;
    Arrays<ConstUnmanagedIntList> a_d;

    void init (const MetaDataBuilder& mdb) {
      const Int ntracers = mdb.trcr2prob.size();

      Me::init("trcr2prob", a_d_.trcr2prob, a_h_.trcr2prob, ntracers);
      std::copy(mdb.trcr2prob.begin(), mdb.trcr2prob.end(), a_h_.trcr2prob.data());
      Kokkos::deep_copy(a_d_.trcr2prob, a_h_.trcr2prob);

      Me::init("bidx2trcr", a_d_.bidx2trcr, a_h_.bidx2trcr, ntracers);
      Me::init("trcr2bl2r", a_d_.trcr2bl2r, a_h_.trcr2bl2r, ntracers);
      Me::init("trcr2br2l", a_d_.trcr2br2l, a_h_.trcr2br2l, ntracers);
      a_h_.prob2trcrptr[0] = 0;
      a_h_.prob2bl2r[0] = 1; // rho is at 0.
      a_h_.prob2br2l[0] = 0;
      for (Int pi = 0; pi < nprobtypes; ++pi) {
        a_h_.prob2trcrptr[pi+1] = a_h_.prob2trcrptr[pi];
        const Int l2rbulksz = get_problem_type_l2r_bulk_size(get_problem_type(pi));
        const Int r2lbulksz = get_problem_type_r2l_bulk_size(get_problem_type(pi));
        for (Int ti = 0; ti < ntracers; ++ti) {
          const auto problem_type = a_h_.trcr2prob[ti];
          if (problem_type != problem_type_[pi]) continue;
          const auto tcnt = a_h_.prob2trcrptr[pi+1] - a_h_.prob2trcrptr[pi];
          a_h_.trcr2bl2r[ti] = a_h_.prob2bl2r[pi] + tcnt*l2rbulksz;
          a_h_.trcr2br2l[ti] = a_h_.prob2br2l[pi] + tcnt*r2lbulksz;
          a_h_.bidx2trcr[a_h_.prob2trcrptr[pi+1]++] = ti;
        }
        Int ni = a_h_.prob2trcrptr[pi+1] - a_h_.prob2trcrptr[pi];
        a_h_.prob2bl2r[pi+1] = a_h_.prob2bl2r[pi] + ni*l2rbulksz;
        a_h_.prob2br2l[pi+1] = a_h_.prob2br2l[pi] + ni*r2lbulksz;
      }
      Kokkos::deep_copy(a_d_.bidx2trcr, a_h_.bidx2trcr);
      Kokkos::deep_copy(a_d_.trcr2bl2r, a_h_.trcr2bl2r);
      Kokkos::deep_copy(a_d_.trcr2br2l, a_h_.trcr2br2l);

      Me::init("trcr2bidx", a_d_.trcr2bidx, a_h_.trcr2bidx, ntracers);
      for (Int ti = 0; ti < ntracers; ++ti)
        a_h_.trcr2bidx(a_h_.bidx2trcr(ti)) = ti;
      Kokkos::deep_copy(a_d_.trcr2bidx, a_h_.trcr2bidx);
            
      a_h = a_h_;

      // Won't default construct Unmanaged, so have to do pointer stuff and raw
      // array copy explicitly.
      a_d.trcr2prob = a_d_.trcr2prob;
      a_d.bidx2trcr = a_d_.bidx2trcr;
      a_d.trcr2bidx = a_d_.trcr2bidx;
      a_d.trcr2bl2r = a_d_.trcr2bl2r;
      a_d.trcr2br2l = a_d_.trcr2br2l;
      std::copy(a_h_.prob2trcrptr, a_h_.prob2trcrptr + nprobtypes + 1,
                a_d.prob2trcrptr);
      std::copy(a_h_.prob2bl2r, a_h_.prob2bl2r + nprobtypes + 1, a_d.prob2bl2r);
      std::copy(a_h_.prob2br2l, a_h_.prob2br2l + nprobtypes + 1, a_d.prob2br2l);
      qlt_assert(a_d.prob2trcrptr[nprobtypes] == ntracers);
    }

  private:
    static constexpr Int problem_type_[] = { CPT::st, CPT::cst, CPT::t, CPT::ct };
    Arrays<typename IntList::HostMirror> a_h_;
    Arrays<IntList> a_d_;
  };

  struct BulkData {
    typedef Kokkos::View<Real*, Kokkos::LayoutLeft, Device> RealList;
    typedef impl::Unmanaged<RealList> UnmanagedRealList;

    UnmanagedRealList l2r_data, r2l_data;

    void init (const MetaData& md, const Int& nslots) {
      l2r_data_ = RealList("l2r_data", md.a_h.prob2bl2r[md.nprobtypes]*nslots);
      r2l_data_ = RealList("r2l_data", md.a_h.prob2br2l[md.nprobtypes]*nslots);
      l2r_data = l2r_data_;
      r2l_data = r2l_data_;
    }

  private:
    RealList l2r_data_, r2l_data_;
  };

private:
  Parallel::Ptr p_;
  // Tree and communication topology.
  impl::NodeSets::ConstPtr ns_;
  // Globally unique cellidx -> rank-local index.
  std::map<Int,Int> gci2lci_;
  // Temporary to collect caller's tracer information prior to calling
  // end_tracer_declarations().
  typename MetaDataBuilder::Ptr mdb_;
  // Constructed in end_tracer_declarations().
  MetaData md_;
  BulkData bd_;

private:
  void init (const Parallel::Ptr& p, const Int& ncells, const tree::Node::Ptr& tree) {
    p_ = p;
    Timer::start(Timer::analyze);
    ns_ = impl::analyze(p, ncells, tree);
    init_ordinals();
    Timer::stop(Timer::analyze);
    mdb_ = std::make_shared<MetaDataBuilder>();
  }

  void init_ordinals () {
    for (const auto& n : ns_->levels[0].nodes)
      gci2lci_[n->id] = n->offset;
  }

  KOKKOS_INLINE_FUNCTION
  static void solve_node_problem (const Int problem_type,
                                  const Real& rhom, const Real* pd, const Real& Qm,
                                  const Real& rhom0, const Real* k0d, Real& Qm0,
                                  const Real& rhom1, const Real* k1d, Real& Qm1) {
    if ( ! (problem_type & ProblemType::shapepreserve)) {      
      Real mpd[3], mk0d[3], mk1d[3];
      mpd[0]  = pd [0]*rhom ; mpd [1] = pd[1] ; mpd [2] = pd [2]*rhom ;
      mk0d[0] = k0d[0]*rhom0; mk0d[1] = k0d[1]; mk0d[2] = k0d[2]*rhom0;
      mk1d[0] = k1d[0]*rhom1; mk1d[1] = k1d[1]; mk1d[2] = k1d[2]*rhom1;
      slv::solve_node_problem(rhom, mpd, Qm, rhom0, mk0d, Qm0, rhom1, mk1d, Qm1);
      return;
    }
    slv::solve_node_problem(rhom, pd, Qm, rhom0, k0d, Qm0, rhom1, k1d, Qm1);
  }

public:
  // Set up QLT topology and communication data structures based on a tree.
  QLT (const Parallel::Ptr& p, const Int& ncells, const tree::Node::Ptr& tree) {
    init(p, ncells, tree);
  }

  void print (std::ostream& os) const {
    ns_->print(os);
  }

  // Number of cells owned by this rank.
  Int nlclcells () const { return ns_->levels[0].nodes.size(); }

  // Cells owned by this rank, in order of local numbering. Thus,
  // gci2lci(gcis[i]) == i. Ideally, the caller never actually calls gci2lci(),
  // and instead uses the information from get_owned_glblcells to determine
  // local cell indices.
  void get_owned_glblcells (std::vector<Int>& gcis) const {
    gcis.resize(ns_->levels[0].nodes.size());
    for (const auto& n : ns_->levels[0].nodes)
      gcis[n->offset] = n->id;
  }

  // For global cell index cellidx, i.e., the globally unique ordinal associated
  // with a cell in the caller's tree, return this rank's local index for
  // it. This is not an efficient operation.
  Int gci2lci (const Int& gci) const {
    const auto it = gci2lci_.find(gci);
    if (it == gci2lci_.end()) {
      pr(puf(gci));
      std::vector<Int> gcis;
      get_owned_glblcells(gcis);
      mprarr(gcis);
    }
    qlt_throw_if(it == gci2lci_.end(), "gci " << gci << " not in gci2lci map.");
    return it->second;
  }

  // Set up QLT tracer metadata. Once end_tracer_declarations is called, it is
  // an error to call declare_tracer again. Call declare_tracer in order of the
  // tracer index in the caller's numbering.
  void declare_tracer (int problem_type) {
    qlt_throw_if( ! mdb_, "end_tracer_declarations was already called; "
                  "it is an error to call declare_tracer now.");
    // For its exception side effect, and to get canonical problem type, since
    // some possible problem types map to the same canonical one:
    problem_type = md_.get_problem_type(md_.get_problem_type_idx(problem_type));
    mdb_->trcr2prob.push_back(problem_type);
  }

  void end_tracer_declarations () {
    md_.init(*mdb_);
    mdb_ = nullptr;
    bd_.init(md_, ns_->nslots);
  }

  int get_problem_type (const Int& tracer_idx) const {
    qlt_throw_if(tracer_idx < 0 || tracer_idx > md_.a_h.trcr2prob.extent_int(0),
                 "tracer_idx is out of bounds: " << tracer_idx);
    return md_.a_h.trcr2prob[tracer_idx];
  }

  Int get_num_tracers () const {
    return md_.a_h.trcr2prob.size();
  }

  // set_{rho,Q}: Set cell values prior to running the QLT algorithm.
  //   set_rho must be called before set_Q.
  //   lclcellidx is gci2lci(cellidx).
  //   Notation:
  //     rho: Total density.
  //       Q: Tracer density.
  //       q: Tracer mixing ratio = Q/rho.
  //      *m: Mass corresponding to the density; results from an integral over a
  //          region, such as a cell.
  KOKKOS_INLINE_FUNCTION
  void set_rho (const Int& lclcellidx,
                // Current total mass in this cell.
                const Real& rhom) {
    qlt_throw_if(mdb_, "Call end_tracer_declarations before this function.");
    const Int ndps = md_.a_d.prob2bl2r[md_.nprobtypes];
    bd_.l2r_data(ndps*lclcellidx) = rhom;
  }

  KOKKOS_INLINE_FUNCTION
  void set_Q (const Int& lclcellidx, const Int& tracer_idx,
              // Current tracer mass in this cell.
              const Real& Qm,
              // Minimum and maximum permitted tracer mass in this cell.
              const Real& Qm_min, const Real& Qm_max,
              // If mass conservation is requested, provide the previous Qm,
              // which will be summed to give the desired global mass.
              const Real Qm_prev = -1) {
    qlt_throw_if(mdb_, "Call end_tracer_declarations before this function.");
    const Int ndps = md_.a_d.prob2bl2r[md_.nprobtypes];
    Real* bd; {
      const Int bdi = md_.a_d.trcr2bl2r(tracer_idx);
      bd = &bd_.l2r_data(ndps*lclcellidx + bdi);
    }
    bd[1] = Qm;
    {
      const Int problem_type = md_.a_d.trcr2prob(tracer_idx);
      if (problem_type & ProblemType::shapepreserve) {
        bd[0] = Qm_min;
        bd[2] = Qm_max;
      } else if (problem_type & ProblemType::consistent) {
        const Real rhom = bd_.l2r_data(ndps*lclcellidx);
        bd[0] = Qm_min / rhom;
        bd[2] = Qm_max / rhom;
      } else {
        qlt_throw_if(true, "set_Q: invalid problem_type.");
      }
      if (problem_type & ProblemType::conserve) {
        qlt_kernel_throw_if(Qm_prev < 0, "Qm_prev was not provided to set_Q.");
        bd[3] = Qm_prev;
      }
    }
  }

  // Run the QLT algorithm with the values set by set_{rho,Q}.
  void run () {
    Timer::start(Timer::qltrunl2r);
    using namespace impl;
    // Number of data per slot.
    const Int l2rndps = md_.a_d.prob2bl2r[md_.nprobtypes];
    const Int r2lndps = md_.a_d.prob2br2l[md_.nprobtypes];
    // Leaves to root.
    for (size_t il = 0; il < ns_->levels.size(); ++il) {
      auto& lvl = ns_->levels[il];
      // Set up receives.
      for (size_t i = 0; i < lvl.kids.size(); ++i) {
        const auto& mmd = lvl.kids[i];
        mpi::irecv(*p_, &bd_.l2r_data(mmd.offset*l2rndps), mmd.size*l2rndps, mmd.rank,
                   NodeSets::mpitag, &lvl.kids_req[i]);
      }
      //todo Replace with simultaneous waitany and isend.
      Timer::start(Timer::waitall);
      mpi::waitall(lvl.kids_req.size(), lvl.kids_req.data());
      Timer::stop(Timer::waitall);
      // Combine kids' data.
      //todo Kernelize, interacting with waitany todo above.
      Timer::start(Timer::snp);
      for (const auto& n : lvl.nodes) {
        if ( ! n->nkids) continue;
        qlt_kernel_assert(n->nkids == 2);
        // Total density.
        bd_.l2r_data(n->offset*l2rndps) = (bd_.l2r_data(n->kids[0]->offset*l2rndps) +
                                           bd_.l2r_data(n->kids[1]->offset*l2rndps));
        // Tracers.
        for (Int pti = 0; pti < md_.nprobtypes; ++pti) {
          const Int problem_type = md_.get_problem_type(pti);
          const bool sum_only = problem_type & ProblemType::shapepreserve;
          const Int bsz = md_.get_problem_type_l2r_bulk_size(problem_type);
          const Int bis = md_.a_d.prob2trcrptr[pti], bie = md_.a_d.prob2trcrptr[pti+1];
          for (Int bi = bis; bi < bie; ++bi) {
            const Int bdi = md_.a_d.trcr2bl2r(md_.a_d.bidx2trcr(bi));
            Real* const me = &bd_.l2r_data(n->offset*l2rndps + bdi);
            const Real* const k0 = &bd_.l2r_data(n->kids[0]->offset*l2rndps + bdi);
            const Real* const k1 = &bd_.l2r_data(n->kids[1]->offset*l2rndps + bdi);
            me[0] = sum_only ? k0[0] + k1[0] : impl::min(k0[0], k1[0]);
            me[1] =            k0[1] + k1[1] ;
            me[2] = sum_only ? k0[2] + k1[2] : impl::max(k0[2], k1[2]);
            if (bsz == 4)
              me[3] =          k0[3] + k1[3] ;
          }
        }
      }
      Timer::stop(Timer::snp);
      // Send to parents.
      for (size_t i = 0; i < lvl.me.size(); ++i) {
        const auto& mmd = lvl.me[i];
        mpi::isend(*p_, &bd_.l2r_data(mmd.offset*l2rndps), mmd.size*l2rndps, mmd.rank,
                   NodeSets::mpitag, &lvl.me_req[i]);
      }
      if (il+1 == ns_->levels.size()) {
        Timer::start(Timer::waitall);
        mpi::waitall(lvl.me_req.size(), lvl.me_req.data());
        Timer::stop(Timer::waitall);
      }
    }
    Timer::stop(Timer::qltrunl2r); Timer::start(Timer::qltrunr2l);
    // Root.
    if ( ! ns_->levels.empty() && ns_->levels.back().nodes.size() == 1 &&
         ! ns_->levels.back().nodes[0]->parent) {
      const auto& n = ns_->levels.back().nodes[0];
      for (Int pti = 0; pti < md_.nprobtypes; ++pti) {
        const Int problem_type = md_.get_problem_type(pti);
        const Int bis = md_.a_d.prob2trcrptr[pti], bie = md_.a_d.prob2trcrptr[pti+1];
        for (Int bi = bis; bi < bie; ++bi) {
          const Int l2rbdi = md_.a_d.trcr2bl2r(md_.a_d.bidx2trcr(bi));
          const Int r2lbdi = md_.a_d.trcr2br2l(md_.a_d.bidx2trcr(bi));
          // If QLT is enforcing global mass conservation, set the root's r2l Qm
          // value to the l2r Qm_prev's sum; otherwise, copy the l2r Qm value to
          // the r2l one.
          const Int os = problem_type & ProblemType::conserve ? 3 : 1;
          bd_.r2l_data(n->offset*r2lndps + r2lbdi) =
            bd_.l2r_data(n->offset*l2rndps + l2rbdi + os);
          if ( ! (problem_type & ProblemType::shapepreserve)) {
            // We now know the global q_{min,max}. Start propagating it
            // leafward.
            bd_.r2l_data(n->offset*r2lndps + r2lbdi + 1) =
              bd_.l2r_data(n->offset*l2rndps + l2rbdi + 0);
            bd_.r2l_data(n->offset*r2lndps + r2lbdi + 2) =
              bd_.l2r_data(n->offset*l2rndps + l2rbdi + 2);
          }
        }
      }
    }
    // Root to leaves.
    for (size_t il = ns_->levels.size(); il > 0; --il) {
      auto& lvl = ns_->levels[il-1];
      for (size_t i = 0; i < lvl.me.size(); ++i) {
        const auto& mmd = lvl.me[i];
        mpi::irecv(*p_, &bd_.r2l_data(mmd.offset*r2lndps), mmd.size*r2lndps, mmd.rank,
                   NodeSets::mpitag, &lvl.me_req[i]);
      }
      //todo Replace with simultaneous waitany and isend.
      Timer::start(Timer::waitall);
      mpi::waitall(lvl.me_req.size(), lvl.me_req.data());
      Timer::stop(Timer::waitall);
      // Solve QP for kids' values.
      //todo Kernelize, interacting with waitany todo above.
      Timer::start(Timer::snp);
      for (const auto& n : lvl.nodes) {
        if ( ! n->nkids) continue;
        for (Int pti = 0; pti < md_.nprobtypes; ++pti) {
          const Int problem_type = md_.get_problem_type(pti);
          const Int bis = md_.a_d.prob2trcrptr[pti], bie = md_.a_d.prob2trcrptr[pti+1];
          for (Int bi = bis; bi < bie; ++bi) {
            const Int l2rbdi = md_.a_d.trcr2bl2r(md_.a_d.bidx2trcr(bi));
            const Int r2lbdi = md_.a_d.trcr2br2l(md_.a_d.bidx2trcr(bi));
            qlt_assert(n->nkids == 2);
            if ( ! (problem_type & ProblemType::shapepreserve)) {
              // Pass q_{min,max} info along. l2r data are updated for use in
              // solve_node_problem. r2l data are updated for use in isend.
              const Real q_min = bd_.r2l_data(n->offset*r2lndps + r2lbdi + 1);
              const Real q_max = bd_.r2l_data(n->offset*r2lndps + r2lbdi + 2);
              bd_.l2r_data(n->offset*l2rndps + l2rbdi + 0) = q_min;
              bd_.l2r_data(n->offset*l2rndps + l2rbdi + 2) = q_max;
              for (Int k = 0; k < 2; ++k) {
                bd_.l2r_data(n->kids[k]->offset*l2rndps + l2rbdi + 0) = q_min;
                bd_.l2r_data(n->kids[k]->offset*l2rndps + l2rbdi + 2) = q_max;
                bd_.r2l_data(n->kids[k]->offset*r2lndps + r2lbdi + 1) = q_min;
                bd_.r2l_data(n->kids[k]->offset*r2lndps + r2lbdi + 2) = q_max;
              }
            }
            const auto& k0 = n->kids[0];
            const auto& k1 = n->kids[1];
            solve_node_problem(
              problem_type,
               bd_.l2r_data( n->offset*l2rndps),
              &bd_.l2r_data( n->offset*l2rndps + l2rbdi),
               bd_.r2l_data( n->offset*r2lndps + r2lbdi),
               bd_.l2r_data(k0->offset*l2rndps),
              &bd_.l2r_data(k0->offset*l2rndps + l2rbdi),
               bd_.r2l_data(k0->offset*r2lndps + r2lbdi),
               bd_.l2r_data(k1->offset*l2rndps),
              &bd_.l2r_data(k1->offset*l2rndps + l2rbdi),
               bd_.r2l_data(k1->offset*r2lndps + r2lbdi));
          }
        }
      }
      Timer::stop(Timer::snp);
      // Send.
      for (size_t i = 0; i < lvl.kids.size(); ++i) {
        const auto& mmd = lvl.kids[i];
        mpi::isend(*p_, &bd_.r2l_data(mmd.offset*r2lndps), mmd.size*r2lndps, mmd.rank,
                   NodeSets::mpitag, &lvl.kids_req[i]);
      }
    }
    // Wait on sends to clean up.
    for (size_t il = 0; il < ns_->levels.size(); ++il) {
      auto& lvl = ns_->levels[il];
      if (il+1 < ns_->levels.size())
        mpi::waitall(lvl.me_req.size(), lvl.me_req.data());
      mpi::waitall(lvl.kids_req.size(), lvl.kids_req.data());
    }
    Timer::stop(Timer::qltrunr2l);
  }
  // Get a cell's tracer mass Qm after the QLT algorithm has run.
  KOKKOS_INLINE_FUNCTION
  Real get_Q (const Int& lclcellidx, const Int& tracer_idx) {
    const Int ndps = md_.a_d.prob2br2l[md_.nprobtypes];
    const Int bdi = md_.a_d.trcr2br2l(tracer_idx);
    return bd_.r2l_data(ndps*lclcellidx + bdi);
  }
}; // class QLT

template <typename ExeSpace>
constexpr Int QLT<ExeSpace>::MetaData::problem_type_[];

template class QLT<Kokkos::DefaultExecutionSpace>;

namespace test {
using namespace impl;

class TestQLT {
  typedef QLT<Kokkos::DefaultExecutionSpace> QLTT;
  typedef Kokkos::View<Real**, QLTT::Device> R2D;

  struct Tracer {
    typedef QLTT::ProblemType PT;
    
    Int idx;
    Int problem_type;
    Int perturbation_type;
    bool no_change_should_hold, safe_should_hold, local_should_hold;
    bool write;

    std::string str () const {
      std::stringstream ss;
      ss << "(ti " << idx;
      if (problem_type & PT::conserve) ss << " c";
      if (problem_type & PT::shapepreserve) ss << " s";
      if (problem_type & PT::consistent) ss << " t";
      ss << " pt " << perturbation_type << " ssh " << safe_should_hold
         << " lsh " << local_should_hold << ")";
      return ss.str();
    }

    Tracer ()
      : idx(-1), problem_type(-1), perturbation_type(-1), no_change_should_hold(false),
        safe_should_hold(true), local_should_hold(true), write(false)
    {}
  };

  struct Values {
    Values (const Int ntracers, const Int ncells)
      : ncells_(ncells), v_((4*ntracers + 1)*ncells)
    {}
    Int ncells () const { return ncells_; }
    Real* rhom () { return v_.data(); }
    Real* Qm_min  (const Int& ti) { return v_.data() + ncells_*(1 + 4*ti    ); }
    Real* Qm      (const Int& ti) { return v_.data() + ncells_*(1 + 4*ti + 1); }
    Real* Qm_max  (const Int& ti) { return v_.data() + ncells_*(1 + 4*ti + 2); }
    Real* Qm_prev (const Int& ti) { return v_.data() + ncells_*(1 + 4*ti + 3); }
    const Real* rhom () const { return const_cast<Values*>(this)->rhom(); }
    const Real* Qm_min  (const Int& ti) const
    { return const_cast<Values*>(this)->Qm_min (ti); }
    const Real* Qm      (const Int& ti) const
    { return const_cast<Values*>(this)->Qm     (ti); }
    const Real* Qm_max  (const Int& ti) const
    { return const_cast<Values*>(this)->Qm_max (ti); }
    const Real* Qm_prev (const Int& ti) const
    { return const_cast<Values*>(this)->Qm_prev(ti); }
  private:
    Int ncells_;
    std::vector<Real> v_;
  };

  // For solution output, if requested.
  struct Writer {
    std::unique_ptr<FILE, impl::FILECloser> fh;
    std::vector<Int> ngcis;  // Number of i'th rank's gcis_ array.
    std::vector<int> displs; // Cumsum of above.
    std::vector<Int> gcis;   // Global cell indices packed by rank's gcis_ vector.
    ~Writer () {
      if ( ! fh) return;
      fprintf(fh.get(), "  return s\n");
    }
  };

private:
  const Parallel::Ptr p_;
  const Int ncells_;
  QLTT qlt_;
  // Caller index (local cell index in the app code) -> QLT lclcellidx.
  std::vector<Int> gcis_, i2lci_;
  std::vector<Tracer> tracers_;
  // For optional output.
  bool write_inited_;
  std::shared_ptr<Writer> w_; // Only on root.

private:
  void init_numbering (const tree::Node::Ptr& node) {
    // TestQLT doesn't actually care about a particular ordering, as there is no
    // geometry to the test problem. However, use *some* ordering to model what
    // a real problem must do.
    if ( ! node->nkids) {
      if (node->rank == p_->rank()) {
        gcis_.push_back(node->cellidx);
        i2lci_.push_back(qlt_.gci2lci(gcis_.back()));
      }
      return;
    }
    for (Int i = 0; i < node->nkids; ++i)
      init_numbering(node->kids[i]);
  }

  void init_tracers () {
    Timer::start(Timer::trcrinit);
    typedef Tracer::PT PT;
    static const Int pts[] = {
      PT::conserve | PT::shapepreserve | PT::consistent,
      PT::shapepreserve, // Test a noncanonical problem type.
      PT::conserve | PT::consistent,
      PT::consistent
    };
    Int tracer_idx = 0;
    for (Int perturb = 0; perturb < 6; ++perturb)
      for (Int ti = 0; ti < 4; ++ti) {
        Tracer t;
        t.problem_type = pts[ti];
        const bool shapepreserve = t.problem_type & PT::shapepreserve;
        t.idx = tracer_idx++;
        t.perturbation_type = perturb;
        t.safe_should_hold = true;
        t.no_change_should_hold = perturb == 0;
        t.local_should_hold = perturb < 4 && shapepreserve;
        t.write = perturb == 2 && ti == 2;
        tracers_.push_back(t);
        qlt_.declare_tracer(t.problem_type);
      }
    qlt_.end_tracer_declarations();
    qlt_assert(qlt_.get_num_tracers() == static_cast<Int>(tracers_.size()));
    for (size_t i = 0; i < tracers_.size(); ++i)
      qlt_assert(qlt_.get_problem_type(i) == (tracers_[i].problem_type | PT::consistent));
    Timer::stop(Timer::trcrinit);
  }

  static Real urand () { return rand() / ((Real) RAND_MAX + 1.0); }

  static void generate_rho (Values& v) {
    auto r = v.rhom();
    const Int n = v.ncells();
    for (Int i = 0; i < n; ++i)
      r[i] = 0.5 + 1.5*urand();
  }

  static void generate_Q (const Tracer& t, Values& v) {
    Real* rhom = v.rhom(), * Qm_min = v.Qm_min(t.idx), * Qm = v.Qm(t.idx),
      * Qm_max = v.Qm_max(t.idx), * Qm_prev = v.Qm_prev(t.idx);
    const Int n = v.ncells();
    for (Int i = 0; i < n; ++i) {
      const Real
        q_min = 0.1 + 0.8*urand(),
        q_max = std::min<Real>(1, q_min + (0.9 - q_min)*urand()),
        q = q_min + (q_max - q_min)*urand();
      // Check correctness up to FP.
      assert(q_min >= 0 &&
             q_max <= 1 + 10*std::numeric_limits<Real>::epsilon() &&
             q_min <= q && q <= q_max);
      Qm_min[i] = q_min*rhom[i];
      Qm_max[i] = q_max*rhom[i];
      // Protect against FP error.
      Qm[i] = std::max<Real>(Qm_min[i], std::min(Qm_max[i], q*rhom[i]));
      // Set previous Qm to the current unperturbed value.
      Qm_prev[i] = Qm[i];
    }
  }

  static void gen_rand_perm (const size_t n, std::vector<Int>& p) {
    p.resize(n);
    for (size_t i = 0; i < n; ++i)
      p[i] = i;
    for (size_t i = 0; i < n; ++i) {
      const int j = urand()*n, k = urand()*n;
      std::swap(p[j], p[k]);
    }
  }

  // Permuting the Qm array, even just on a rank as long as there is > 1 cell,
  // produces a problem likely requiring considerable reconstruction, which
  // reconstruction assuredly satisfies the properties. But because this is a
  // local operation only, it doesn't test the 1 cell/rank case.
  static void permute_Q (const Tracer& t, Values& v) {
    Real* const Qm = v.Qm(t.idx);
    const Int N = v.ncells();
    std::vector<Int> p;
    gen_rand_perm(N, p);
    std::vector<Real> Qm_orig(N);
    std::copy(Qm, Qm + N, Qm_orig.begin());
    for (Int i = 0; i < N; ++i)
      Qm[i] = Qm_orig[p[i]];
  }

  void add_const_to_Q (const Tracer& t, Values& v,
                       // Move 0 < alpha <= 1 of the way to the QLT or safety
                       // feasibility bound.
                       const Real& alpha,
                       // Whether the modification should be done in a
                       // mass-conserving way.
                       const bool conserve_mass,
                       // Only safety problem is feasible.
                       const bool safety_problem) {
    // Some of these reductions aren't used at present. Might add more test
    // options later that use them.
    Real rhom, Qm, Qm_max; {
      Real Qm_sum_lcl[3] = {0};
      for (Int i = 0; i < v.ncells(); ++i) {
        Qm_sum_lcl[0] += v.rhom()[i];
        Qm_sum_lcl[1] += v.Qm(t.idx)[i];
        Qm_sum_lcl[2] += v.Qm_max(t.idx)[i];
      }
      Real Qm_sum_gbl[3] = {0};
      mpi::all_reduce(*p_, Qm_sum_lcl, Qm_sum_gbl, 3, MPI_SUM);
      rhom = Qm_sum_gbl[0]; Qm = Qm_sum_gbl[1]; Qm_max = Qm_sum_gbl[2];
    }
    Real Qm_max_safety = 0;
    if (safety_problem) {
      Real q_safety_lcl = v.Qm_max(t.idx)[0] / v.rhom()[0];
      for (Int i = 1; i < v.ncells(); ++i)
        q_safety_lcl = std::max(q_safety_lcl, v.Qm_max(t.idx)[i] / v.rhom()[i]);
      Real q_safety_gbl = 0;
      mpi::all_reduce(*p_, &q_safety_lcl, &q_safety_gbl, 1, MPI_MAX);
      Qm_max_safety = q_safety_gbl*rhom;
    }
    const Real dQm = safety_problem ?
      ((Qm_max - Qm) + alpha * (Qm_max_safety - Qm_max)) / ncells_ :
      alpha * (Qm_max - Qm) / ncells_;
    for (Int i = 0; i < v.ncells(); ++i)
      v.Qm(t.idx)[i] += dQm;
    // Now permute Qm so that it's a little more interesting.
    permute_Q(t, v);
    // Adjust Qm_prev. Qm_prev is used to test the PT::conserve case, and also
    // simply to record the correct total mass. The modification above modified
    // Q's total mass. If conserve_mass, then Qm_prev needs to be made to sum to
    // the same new mass. If ! conserve_mass, we want Qm_prev to be modified in
    // an interesting way, so that PT::conserve doesn't trivially undo the mod
    // that was made above when the root fixes the mass discrepancy.
    const Real
      relax = 0.9,
      dQm_prev = (conserve_mass ? dQm :
                  (safety_problem ?
                   ((Qm_max - Qm) + relax*alpha * (Qm_max_safety - Qm_max)) / ncells_ :
                   relax*alpha * (Qm_max - Qm) / ncells_));
    for (Int i = 0; i < v.ncells(); ++i)
      v.Qm_prev(t.idx)[i] += dQm_prev;
  }

  void perturb_Q (const Tracer& t, Values& v) {
    // QLT is naturally mass conserving. But if QLT isn't being asked to impose
    // mass conservation, then the caller better have a conservative
    // method. Here, we model that by saying that Qm_prev and Qm should sum to
    // the same mass.
    const bool cm = ! (t.problem_type & Tracer::PT::conserve);
    // For the edge cases, we cannot be exactly on the edge and still expect the
    // q-limit checks to pass to machine precision. Thus, back away from the
    // edge by an amount that bounds the error in the global mass due to FP,
    // assuming each cell's mass is O(1).
    const Real edg = 1 - ncells_*std::numeric_limits<Real>::epsilon();
    switch (t.perturbation_type) {
    case 0:
      // Do nothing, to test that QLT doesn't make any changes if none is
      // needed.
      break;
    case 1: permute_Q(t, v); break;
    case 2: add_const_to_Q(t, v, 0.5, cm, false); break;
    case 3: add_const_to_Q(t, v, edg, cm, false); break;
    case 4: add_const_to_Q(t, v, 0.5, cm, true ); break;
    case 5: add_const_to_Q(t, v, edg, cm, true ); break;
    }
  }

  static std::string get_tracer_name (const Tracer& t) {
    std::stringstream ss;
    ss << "t" << t.idx;
    return ss.str();
  }

  void init_writer () {
    if (p_->amroot()) {
      w_ = std::make_shared<Writer>();
      w_->fh = std::unique_ptr<FILE, impl::FILECloser>(fopen("QLT.py", "w"));
      int n = gcis_.size();
      w_->ngcis.resize(p_->size());
      mpi::gather(*p_, &n, 1, w_->ngcis.data(), 1, p_->root());
      w_->displs.resize(p_->size() + 1);
      w_->displs[0] = 0;
      for (size_t i = 0; i < w_->ngcis.size(); ++i)
        w_->displs[i+1] = w_->displs[i] + w_->ngcis[i];
      qlt_assert(w_->displs.back() == ncells_);
      w_->gcis.resize(ncells_);
      mpi::gatherv(*p_, gcis_.data(), gcis_.size(), w_->gcis.data(), w_->ngcis.data(),
                   w_->displs.data(), p_->root());
    } else {
      int n = gcis_.size();
      mpi::gather(*p_, &n, 1, static_cast<int*>(nullptr), 0, p_->root());
      Int* Inull = nullptr;
      const int* inull = nullptr;
      mpi::gatherv(*p_, gcis_.data(), gcis_.size(), Inull, inull, inull, p_->root());
    }
    write_inited_ = true;
  }

  void gather_field (const Real* Qm_lcl, std::vector<Real>& Qm_gbl,
                     std::vector<Real>& wrk) {
    if (p_->amroot()) {
      Qm_gbl.resize(ncells_);
      wrk.resize(ncells_);
      mpi::gatherv(*p_, Qm_lcl, gcis_.size(), wrk.data(), w_->ngcis.data(),
                   w_->displs.data(), p_->root());
      for (Int i = 0; i < ncells_; ++i)
        Qm_gbl[w_->gcis[i]] = wrk[i];
    } else {
      Real* rnull = nullptr;
      const int* inull = nullptr;
      mpi::gatherv(*p_, Qm_lcl, gcis_.size(), rnull, inull, inull, p_->root());
    }
  }

  void write_field (const std::string& tracer_name, const std::string& field_name,
                    const std::vector<Real>& Qm) {
    if ( ! p_->amroot()) return;
    fprintf(w_->fh.get(), "  s.%s.%s = [", tracer_name.c_str(), field_name.c_str());
    for (const auto& e : Qm)
      fprintf(w_->fh.get(), "%1.15e, ", e);
    fprintf(w_->fh.get(), "]\n");
  }

  void write_pre (const Tracer& t, Values& v) {
    if ( ! t.write) return;
    std::vector<Real> f, wrk;
    if ( ! write_inited_) {
      init_writer();
      if (w_)
        fprintf(w_->fh.get(),
                "def getsolns():\n"
                "  class Struct:\n"
                "    pass\n"
                "  s = Struct()\n"
                "  s.all = Struct()\n");
      gather_field(v.rhom(), f, wrk);
      write_field("all", "rhom", f);
    }
    const auto name = get_tracer_name(t);
    if (w_)
      fprintf(w_->fh.get(), "  s.%s = Struct()\n", name.c_str());
    gather_field(v.Qm_min(t.idx), f, wrk);
    write_field(name, "Qm_min", f);
    gather_field(v.Qm_prev(t.idx), f, wrk);
    write_field(name, "Qm_orig", f);
    gather_field(v.Qm(t.idx), f, wrk);
    write_field(name, "Qm_pre", f);
    gather_field(v.Qm_max(t.idx), f, wrk);
    write_field(name, "Qm_max", f);
  }

  void write_post (const Tracer& t, Values& v) {
    if ( ! t.write) return;
    const auto name = get_tracer_name(t);
    std::vector<Real> Qm, wrk;
    gather_field(v.Qm(t.idx), Qm, wrk);
    write_field(name, "Qm_qlt", Qm);
  }

  static void check (const QLTT& qlt) {
    const Int n = qlt.nlclcells();
    std::vector<Int> gcis;
    qlt.get_owned_glblcells(gcis);
    qlt_assert(static_cast<Int>(gcis.size()) == n);
    for (Int i = 0; i < n; ++i)
      qlt_assert(qlt.gci2lci(gcis[i]) == i);
  }

  static Int check (const Parallel& p, const std::vector<Tracer>& ts, const Values& v) {
    static const bool details = false;
    static const Real ulp2 = 2*std::numeric_limits<Real>::epsilon();
    Int nerr = 0;
    std::vector<Real> lcl_mass(2*ts.size()), q_min_lcl(ts.size()), q_max_lcl(ts.size());
    std::vector<Int> t_ok(ts.size(), 1), local_violated(ts.size(), 0);
    for (size_t ti = 0; ti < ts.size(); ++ti) {
      const auto& t = ts[ti];

      qlt_assert(t.safe_should_hold);
      const bool safe_only = ! t.local_should_hold;
      const Int n = v.ncells();
      const Real* rhom = v.rhom(), * Qm_min = v.Qm_min(t.idx), * Qm = v.Qm(t.idx),
        * Qm_max = v.Qm_max(t.idx), * Qm_prev = v.Qm_prev(t.idx);

      q_min_lcl[ti] = 1;
      q_max_lcl[ti] = 0;
      for (Int i = 0; i < n; ++i) {
        const bool lv = Qm[i] < Qm_min[i] || Qm[i] > Qm_max[i];
        if (lv) local_violated[ti] = 1;
        if ( ! safe_only && lv) {
          if (details)
            pr("check q " << t.str() << ": " << Qm[i] << " " <<
               (Qm[i] < Qm_min[i] ? Qm[i] - Qm_min[i] : Qm[i] - Qm_max[i]));
          t_ok[ti] = false;
          ++nerr;
        }
        if (t.no_change_should_hold && Qm[i] != Qm_prev[i]) {
          if (details)
            pr("Q should be unchanged but is not: " << Qm_prev[i] << " changed to " <<
               Qm[i] << " in " << t.str());
          t_ok[ti] = false;
          ++nerr;
        }
        lcl_mass[2*ti    ] += Qm_prev[i];
        lcl_mass[2*ti + 1] += Qm[i];
        q_min_lcl[ti] = std::min(q_min_lcl[ti], Qm_min[i]/rhom[i]);
        q_max_lcl[ti] = std::max(q_max_lcl[ti], Qm_max[i]/rhom[i]);
      }
    }

    std::vector<Real> q_min_gbl(ts.size(), 0), q_max_gbl(ts.size(), 0);
    mpi::all_reduce(p, q_min_lcl.data(), q_min_gbl.data(), q_min_lcl.size(), MPI_MIN);
    mpi::all_reduce(p, q_max_lcl.data(), q_max_gbl.data(), q_max_lcl.size(), MPI_MAX);

    for (size_t ti = 0; ti < ts.size(); ++ti) {
      // Check safety problem. If local_should_hold and it does, then the safety
      // problem is by construction also solved (since it's a relaxation of the
      // local problem).
      const auto& t = ts[ti];
      const bool safe_only = ! t.local_should_hold;
      if (safe_only) {
        const Int n = v.ncells();
        const Real* rhom = v.rhom(), * Qm_min = v.Qm_min(t.idx), * Qm = v.Qm(t.idx),
          * Qm_max = v.Qm_max(t.idx);
        const Real q_min = q_min_gbl[ti], q_max = q_max_gbl[ti];
        for (Int i = 0; i < n; ++i) {
          if (Qm[i] < q_min*rhom[i]*(1 - ulp2) ||
              Qm[i] > q_max*rhom[i]*(1 + ulp2)) {
            if (details)
              pr("check q " << t.str() << ": " << q_min*rhom[i] << " " << Qm_min[i] <<
                 " " << Qm[i] << " " << Qm_max[i] << " " << q_max*rhom[i] << " | " <<
                 (Qm[i] < q_min*rhom[i] ?
                  Qm[i] - q_min*rhom[i] :
                  Qm[i] - q_max*rhom[i]));
            t_ok[ti] = false;
            ++nerr;
          }
        }        
      }
    }

    std::vector<Real> glbl_mass(2*ts.size(), 0);
    mpi::reduce(p, lcl_mass.data(), glbl_mass.data(), lcl_mass.size(), MPI_SUM,
                p.root());
    std::vector<Int> t_ok_gbl(ts.size(), 0);
    mpi::reduce(p, t_ok.data(), t_ok_gbl.data(), t_ok.size(), MPI_MIN, p.root());
    // Right now we're not using these:
    std::vector<Int> local_violated_gbl(ts.size(), 0);
    mpi::reduce(p, local_violated.data(), local_violated_gbl.data(),
                local_violated.size(), MPI_MAX, p.root());

    if (p.amroot()) {
      const Real tol = 1e3*std::numeric_limits<Real>::epsilon();
      for (size_t ti = 0; ti < ts.size(); ++ti) {
        // Check mass conservation.
        const Real desired_mass = glbl_mass[2*ti], actual_mass = glbl_mass[2*ti+1],
          rd = reldif(desired_mass, actual_mass);
        const bool mass_failed = rd > tol;
        if (mass_failed) {
          ++nerr;
          t_ok_gbl[ti] = false;
        }
        if ( ! t_ok_gbl[ti]) {
          std::cout << "FAIL  " << ts[ti].str();
          if (mass_failed) std::cout << " mass re " << rd;
          std::cout << "\n";
        }
      }
    }

    return nerr;
  }
  
public:
  TestQLT (const Parallel::Ptr& p, const tree::Node::Ptr& tree,
           const Int& ncells, const bool verbose = false)
    : p_(p), ncells_(ncells), qlt_(p_, ncells, tree), write_inited_(false)
  {
    check(qlt_);
    init_numbering(tree);
    init_tracers();
    if (verbose) qlt_.print(std::cout);
  }

  Int run (const Int nrepeat = 1, const bool write=false) {
    Timer::start(Timer::trcrgen);
    const Int nt = qlt_.get_num_tracers(), nlclcells = qlt_.nlclcells();
    Values v(nt, nlclcells);
    generate_rho(v);
    {
      Real* rhom = v.rhom();
      for (Int i = 0; i < nlclcells; ++i)
        qlt_.set_rho(i2lci_[i], rhom[i]);
    }
    for (Int ti = 0; ti < nt; ++ti) {
      generate_Q(tracers_[ti], v);
      perturb_Q(tracers_[ti], v);
      if (write) write_pre(tracers_[ti], v);
    }
    Timer::stop(Timer::trcrgen);
    for (Int trial = 0; trial <= nrepeat; ++trial) {
      for (Int ti = 0; ti < nt; ++ti) {
        Real* Qm_min = v.Qm_min(ti), * Qm = v.Qm(ti), * Qm_max = v.Qm_max(ti),
          * Qm_prev = v.Qm_prev(ti);
        for (Int i = 0; i < nlclcells; ++i)
          qlt_.set_Q(i2lci_[i], ti, Qm[i], Qm_min[i], Qm_max[i], Qm_prev[i]);
      }
      MPI_Barrier(p_->comm());
      Timer::start(Timer::qltrun);
      qlt_.run();
      MPI_Barrier(p_->comm());
      Timer::stop(Timer::qltrun);
      if (trial == 0) {
        Timer::reset(Timer::qltrun);
        Timer::reset(Timer::qltrunl2r);
        Timer::reset(Timer::qltrunr2l);
        Timer::reset(Timer::waitall);
        Timer::reset(Timer::snp);
      }
    }
    Timer::start(Timer::trcrcheck);
    Int nerr = 0;
    for (Int ti = 0; ti < nt; ++ti) {
      Real* Qm = v.Qm(ti);
      for (Int i = 0; i < nlclcells; ++i)
        Qm[i] = qlt_.get_Q(i2lci_[i], ti);
      if (write) write_post(tracers_[ti], v);
    }
    nerr += check(*p_, tracers_, v);
    Timer::stop(Timer::trcrcheck);
    return nerr;
  }
};

// Test all QLT variations and situations.
Int test_qlt (const Parallel::Ptr& p, const tree::Node::Ptr& tree, const Int& ncells,
              const int nrepeat = 1,
              // Diagnostic output for dev and illustration purposes. To be
              // clear, no QLT unit test requires output to be checked; each
              // checks in-memory data and returns a failure count.
              const bool write = false,
              const bool verbose = false) {
  return TestQLT(p, tree, ncells, verbose).run(nrepeat, write);
}
} // namespace test

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
  
  void init (const Int nc, const Parallel::Ptr& p,
             const ParallelDecomp::Enum& parallel_decomp) {
    nc_ = nc;
    nranks_ = p->size();
    p_ = p;
    pd_ = parallel_decomp;
    qlt_assert(nranks_ <= nc_);
  }

  Int ncell () const { return nc_; }

  const Parallel::Ptr& parallel () const { return p_; }

  Int rank (const Int& ci) const {
    switch (pd_) {
    case ParallelDecomp::contiguous:
      return std::min(nranks_ - 1, ci / (nc_ / nranks_));
    default: {
      const auto chunk = ci / nranks_;
      return (ci + chunk) % nranks_;
    }
    }
  }

  static Int unittest (const Parallel::Ptr& p) {
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

private:
  Int nc_, nranks_;
  Parallel::Ptr p_;
  ParallelDecomp::Enum pd_;
};

tree::Node::Ptr make_tree (const Mesh& m, const Int cs, const Int ce,
                           const tree::Node* parent) {
  const Int cn = ce - cs, cn0 = cn/2;
  tree::Node::Ptr n = std::make_shared<tree::Node>();
  n->parent = parent;
  if (cn == 1) {
    n->nkids = 0;
    n->rank = m.rank(cs);
    n->cellidx = cs;
    return n;
  }
  n->nkids = 2;
  n->kids[0] = make_tree(m, cs, cs + cn0, n.get());
  n->kids[1] = make_tree(m, cs + cn0, ce, n.get());
  return n;
}

tree::Node::Ptr make_tree (const Mesh& m) {
  return make_tree(m, 0, m.ncell(), nullptr);
}

tree::Node::Ptr make_tree (const Parallel::Ptr& p, const Int& ncells) {
  Mesh m(ncells, p);
  return make_tree(m);
}

namespace test {
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
  for (size_t id = 0; id < sizeof(dists)/sizeof(*dists); ++id) {
    Mesh m(std::max(42, 3*p->size()), p, Mesh::ParallelDecomp::pseudorandom);
    tree::Node::Ptr tree = make_tree(m);
    std::vector<Int> cells(m.ncell(), 0);
    mark_cells(tree, cells);
    for (Int i = 0; i < m.ncell(); ++i)
      if (cells[i] != 1) ++ne;
  }
  return ne;
}
} // namespace test
} // namespace oned

namespace test {
Int unittest_NodeSets (const Parallel::Ptr& p) {
  using Mesh = oned::Mesh;
  const Int szs[] = { p->size(), 3*p->size() };
  const Mesh::ParallelDecomp::Enum dists[] = { Mesh::ParallelDecomp::pseudorandom,
                                               Mesh::ParallelDecomp::contiguous };
  Int nerr = 0;
  for (size_t is = 0; is < sizeof(szs)/sizeof(*szs); ++is)
    for (size_t id = 0; id < sizeof(dists)/sizeof(*dists); ++id) {
      Mesh m(szs[is], p, dists[id]);
      tree::Node::Ptr tree = make_tree(m);
      impl::NodeSets::ConstPtr nodesets = impl::analyze(p, m.ncell(), tree);
      tree = nullptr;
      nerr += impl::unittest(p, nodesets, m.ncell());
    }
  return nerr;
}

Int unittest_QLT (const Parallel::Ptr& p, const bool write_requested=false) {
  using Mesh = oned::Mesh;
  const Int szs[] = { p->size(), 2*p->size(), 7*p->size(), 21*p->size() };
  const Mesh::ParallelDecomp::Enum dists[] = { Mesh::ParallelDecomp::contiguous,
                                               Mesh::ParallelDecomp::pseudorandom };
  Int nerr = 0;
  for (size_t is = 0, islim = sizeof(szs)/sizeof(*szs); is < islim; ++is)
    for (size_t id = 0, idlim = sizeof(dists)/sizeof(*dists); id < idlim; ++id) {
      if (p->amroot()) {
        std::cout << " (" << szs[is] << ", " << id << ")";
        std::cout.flush();
      }
      Mesh m(szs[is], p, dists[id]);
      tree::Node::Ptr tree = make_tree(m);
      const bool write = (write_requested && m.ncell() < 3000 &&
                          is == islim-1 && id == idlim-1);
      nerr += test::test_qlt(p, tree, m.ncell(), 1, write);
    }
  return nerr;
}

inline bool eq (const std::string& a, const char* const b1, const char* const b2 = 0) {
  return (a == std::string(b1) || (b2 && a == std::string(b2)) ||
          a == std::string("-") + std::string(b1));
}

struct Input {
  bool quiet, unittest, write;
  Int ncells, ntracers, tracer_type, nrepeat;
  bool pseudorandom, verbose;

  class ArgAdvancer {
    const Int argc_;
    char const* const* argv_;
    Int i_;
  public:
    ArgAdvancer (int argc, char** argv) : argc_(argc), argv_(argv), i_(1) {}
    const char* advance () {
      qlt_assert(i_+1 < argc_);
      return argv_[++i_];
    }
    const char* token () const { return argv_[i_]; }
    void incr () { ++i_; }
    bool more () const { return i_ < argc_; }
  };
  
  Input (Int argc, char** argv, const Parallel::Ptr& p) {
    quiet = false;
    unittest = false;
    write = false;
    ncells = 0;
    ntracers = 1;
    tracer_type = 0;
    nrepeat = 1;
    pseudorandom = false;
    verbose = false;
    for (ArgAdvancer aa(argc, argv); aa.more(); aa.incr()) {
      const char* token = aa.token();
      if (eq(token, "-q", "--quiet")) quiet = true;
      else if (eq(token, "-t", "--unittest")) unittest = true;
      else if (eq(token, "-w", "--write")) write = true;
      else if (eq(token, "-nc", "--ncells")) ncells = std::atoi(aa.advance());
      else if (eq(token, "-nt", "--ntracers")) ntracers = std::atoi(aa.advance());
      else if (eq(token, "-tt", "--tracertype")) tracer_type = std::atoi(aa.advance());
      else if (eq(token, "-nr", "--nrepeat")) nrepeat = std::atoi(aa.advance());
      else if (eq(token, "--random")) pseudorandom = true;
      else if (eq(token, "-v", "--verbose")) verbose = true;
      else qlt_throw_if(true, "Invalid token " << token);
    }
    qlt_assert(tracer_type >= 0 && tracer_type < 4);
    qlt_assert(ntracers >= 1);
  }

  void print (std::ostream& os) const {
    os << "ncells " << ncells
       << " nrepeat " << nrepeat;
    //todo
    //<< " ntracers " << ntracers << " tracer_type " << tracer_type;
    if (pseudorandom) os << " random";
    os << "\n";
  }
};

Int run (const Parallel::Ptr& p, const Input& in) {
  Int nerr = 0;
  if (in.unittest) {
    Int ne;
    ne = oned::Mesh::unittest(p);
    if (ne && p->amroot()) std::cerr << "FAIL: Mesh::unittest()\n";
    nerr += ne;
    ne = oned::test::unittest(p);
    if (ne && p->amroot()) std::cerr << "FAIL: oned::unittest_tree()\n";
    nerr += ne;
    ne = unittest_NodeSets(p);
    if (ne && p->amroot()) std::cerr << "FAIL: oned::unittest_NodeSets()\n";
    nerr += ne;
    ne = unittest_QLT(p, in.write);
    if (ne && p->amroot()) std::cerr << "FAIL: oned::unittest_QLT()\n";
    nerr += ne;
    if (p->amroot()) std::cout << "\n";
  }
  if (nerr)
    return nerr;
  // Performance test.
  if (in.ncells > 0) {
    if (p->amroot()) in.print(std::cout);
    oned::Mesh m(in.ncells, p,
                 (in.pseudorandom ?
                  oned::Mesh::ParallelDecomp::pseudorandom :
                  oned::Mesh::ParallelDecomp::contiguous));
    Timer::init();
    Timer::start(Timer::total); Timer::start(Timer::tree);
    tree::Node::Ptr tree = make_tree(m);
    Timer::stop(Timer::tree);
    test::test_qlt(p, tree, in.ncells, in.nrepeat, false, in.verbose);
    Timer::stop(Timer::total);
    if (p->amroot()) Timer::print();
  }
  return nerr;
}

} // namespace test
} // namespace qlt

#ifdef QLT_MAIN
int main (int argc, char** argv) {
  int ret = 0;
  MPI_Init(&argc, &argv);
  auto p = qlt::make_parallel(MPI_COMM_WORLD);
  srand(p->rank());
  Kokkos::initialize(argc, argv);
  try {
    qlt::test::Input in(argc, argv, p);
    ret = qlt::test::run(p, in);
  } catch (const std::exception& e) {
    if (p->amroot())
      std::cerr << e.what();
  }
  Kokkos::finalize_all();
  MPI_Finalize();
  return ret;
}
#endif

// Homme-specific impl details.
//todo Move to a separate file, qlt_homme.cpp.
namespace homme {
typedef qlt::Int Int;
typedef qlt::Real Real;
  
void renumber (const Int* sc2gci, const Int* sc2rank, const qlt::tree::Node::Ptr& node) {
  if (node->nkids)
    for (Int k = 0; k < node->nkids; ++k)
      renumber(sc2gci, sc2rank, node->kids[k]);
  else {
    const Int ci = node->cellidx;
    node->cellidx = sc2gci[ci];
    node->rank = sc2rank[ci];
  }
}

qlt::tree::Node::Ptr
make_tree (const qlt::Parallel::Ptr& p, const Int nelem, const Int* sc2gci, const Int* sc2rank) {
  // Partition 0:nelem-1, the space-filling curve space.
  auto tree = qlt::oned::make_tree(p, nelem);
  // Renumber so that node->cellidx records the global element number, and
  // associate the correct rank with the element.
  renumber(sc2gci, sc2rank, tree);
  return tree;
}

struct QLT {
  typedef std::shared_ptr<QLT> Ptr;
  typedef qlt::QLT<Kokkos::DefaultExecutionSpace> QLTT;
  
  const Int ncell;
  const qlt::Parallel::Ptr p;
  qlt::tree::Node::Ptr tree; // Don't need this except for unit testing.
  qlt::QLT<Kokkos::DefaultExecutionSpace>::Ptr qlt;
  std::vector<Int> ie2gci; // Map Homme ie to Home global cell index.
  std::vector<Int> ie2lci; // Map Homme ie to QLT local cell index (lclcellidx).

  QLT (Int ncell_, const Int* sc2gci, const Int* sc2rank, const qlt::Parallel::Ptr& p_)
    : ncell(ncell_), p(p_), inited_tracers_(false)
  {
    tree = make_tree(p, ncell, sc2gci, sc2rank);
    qlt = std::make_shared<QLTT>(p, ncell, tree);
    ie2gci.resize(qlt->nlclcells());
  }

  void init_tracers (const Int nlev, const Int qsize, const bool need_conservation) {
    typedef QLTT::ProblemType PT;
    for (Int ti = 0, nt = nlev*qsize; ti < nt; ++ti)
      qlt->declare_tracer(PT::shapepreserve |
                          (need_conservation ? PT::conserve : 0));
    qlt->end_tracer_declarations();
  }

private:
  bool inited_tracers_;
};

void set_ie2gci (QLT& q, const Int ie, const Int gci) { q.ie2gci[ie] = gci; }

void init_ie2lci (QLT& q) {
  q.ie2lci.resize(q.ie2gci.size());
  for (size_t ie = 0; ie < q.ie2lci.size(); ++ie)
    q.ie2lci[ie] = q.qlt->gci2lci(q.ie2gci[ie]);
}

void init_tracers (QLT& q, const Int nlev, const Int qsize,
                   const bool need_conservation) {
  q.init_tracers(nlev, qsize, need_conservation);
}

namespace sl { // For sl_advection.F90
// Fortran array wrappers.
template <typename T> using FA2 =
  Kokkos::View<T**,    Kokkos::LayoutLeft, Kokkos::HostSpace>;
template <typename T> using FA4 =
  Kokkos::View<T****,  Kokkos::LayoutLeft, Kokkos::HostSpace>;
template <typename T> using FA5 =
  Kokkos::View<T*****, Kokkos::LayoutLeft, Kokkos::HostSpace>;

// Following are naming conventions in element_state and sl_advection:
//     elem(ie)%state%Q(:,:,k,q) is tracer mixing ratio.
//     elem(ie)%state%dp3d(:,:,k,tl%np1) is essentially total density.
//     elem(ie)%state%Qdp(:,:,k,q,n0_qdp) is Q*dp3d.
//     rho(:,:,k,ie) is spheremp*dp3d, essentially total mass at a GLL point.
//     Hence Q*rho = Q*spheremp*dp3d is tracer mass at a GLL point.
// We need to get pointers to some of these; elem can't be given the bind(C)
// attribute, so we can't take the elem array directly. We get these quantities
// at a mix of previous and current time steps.
//   In the code that follows, _p is previous and _c is current time step. Q is
// renamed to q, and Q is tracer mass at a GLL point.
struct Data {
  typedef std::shared_ptr<Data> Ptr;
  const Int np, nlev, qsize, qsize_d, timelevels;
  Int n0_qdp, n1_qdp, tl_np1;
  std::vector<const Real*> spheremp, dp3d_c;
  std::vector<Real*> q_c, qdp_pc;

  Data (Int lcl_ncell, Int np_, Int nlev_, Int qsize_, Int qsize_d_, Int timelevels_)
    : np(np_), nlev(nlev_), qsize(qsize_), qsize_d(qsize_d_), timelevels(timelevels_),
      spheremp(lcl_ncell, nullptr), dp3d_c(lcl_ncell, nullptr), q_c(lcl_ncell, nullptr),
      qdp_pc(lcl_ncell, nullptr)
  {}
};

static void check (const QLT& q, const Data& d) {
  qlt_assert(q.qlt->nlclcells() == static_cast<Int>(d.spheremp.size()));
}

template <typename T>
void insert (std::vector<T*>& r, const Int i, T* v) {
  qlt_assert(i >= 0 && i < static_cast<int>(r.size()));
  r[i] = v;
}

void insert (const Data::Ptr& d, const Int ie, const Int ptridx, Real* array,
             const Int i0 = 0, const Int i1 = 0) {
  qlt_assert(d);
  switch (ptridx) {
  case 0: insert<const double>(d->spheremp, ie, array); break;
  case 1: insert<      double>(d->qdp_pc,   ie, array); d->n0_qdp = i0; d->n1_qdp = i1; break;
  case 2: insert<const double>(d->dp3d_c,   ie, array); d->tl_np1 = i0; break;
  case 3: insert<      double>(d->q_c,      ie, array); break;
  default: qlt_throw_if(true, "Invalid pointer index " << ptridx);
  }
}

static void run_qlt (QLT& q) {
#ifdef HORIZ_OPENMP
# prama omp barrier
#endif
  q.qlt->run();
#ifdef HORIZ_OPENMP
# prama omp barrier
#endif
}

void run (QLT& qlt, const Data& d, const Real* q_min_r, const Real* q_max_r) {
  static constexpr Int max_np = 4;
  const Int np = d.np, nlev = d.nlev, qsize = d.qsize, ncell = d.spheremp.size();
  qlt_assert(np <= max_np);

  FA5<const Real>
    q_min(q_min_r, np, np, nlev, qsize, ncell),
    q_max(q_max_r, np, np, nlev, qsize, ncell);

  for (Int ie = 0; ie < ncell; ++ie) {
    FA2<const Real> spheremp(d.spheremp[ie], np, np);
    FA5<const Real> qdp_p(d.qdp_pc[ie], np, np, nlev, d.qsize_d, 2);
    FA4<const Real> dp3d_c(d.dp3d_c[ie], np, np, nlev, d.timelevels);
    FA4<const Real> q_c(d.q_c[ie], np, np, nlev, d.qsize_d);
    const Int lci = qlt.ie2lci[ie];

    for (Int k = 0; k < nlev; ++k)
      for (Int q = 0; q < qsize; ++q) {
        const Int ti = k*qsize + q;

        Real Qm = 0, Qm_min = 0, Qm_max = 0, Qm_prev = 0, rhom = 0;
        for (Int j = 0; j < np; ++j)
          for (Int i = 0; i < np; ++i) {
            const Real rhomij = dp3d_c(i,j,k,d.tl_np1) * spheremp(i,j);
            rhom += rhomij;
            Qm += q_c(i,j,k,q) * rhomij;
            Qm_min += q_min(i,j,k,q,ie) * rhomij;
            Qm_max += q_max(i,j,k,q,ie) * rhomij;
            Qm_prev += qdp_p(i,j,k,q,d.n0_qdp) * spheremp(i,j);
          }

        if (ti == 0) qlt.qlt->set_rho(lci, rhom);
        qlt.qlt->set_Q(lci, ti, Qm, Qm_min, Qm_max, Qm_prev);
      }
  }

  run_qlt(qlt);
}

void run_local (QLT& qlt, const Data& d, const Real* q_min_r, const Real* q_max_r,
                const bool scalar_bounds) {
  static constexpr Int max_np = 4, max_np2 = max_np*max_np;
  const Int np = d.np, np2 = np*np, nlev = d.nlev, qsize = d.qsize,
    ncell = d.spheremp.size();
  qlt_assert(np <= max_np);

  FA5<const Real>
    q_min(q_min_r, np, np, nlev, qsize, ncell),
    q_max(q_max_r, np, np, nlev, qsize, ncell);

  for (Int ie = 0; ie < ncell; ++ie) {
    FA2<const Real> spheremp(d.spheremp[ie], np, np);
    FA5<      Real> qdp_c(d.qdp_pc[ie], np, np, nlev, d.qsize_d, 2);
    FA4<const Real> dp3d_c(d.dp3d_c[ie], np, np, nlev, d.timelevels);
    FA4<      Real> q_c(d.q_c[ie], np, np, nlev, d.qsize_d);
    const Int lci = qlt.ie2lci[ie];

    for (Int k = 0; k < nlev; ++k)
      for (Int q = 0; q < qsize; ++q) {
        const Int ti = k*qsize + q;

        Real wa[max_np2], qlo[max_np2], qhi[max_np2], y[max_np2], x[max_np2];
        Real rhom = 0;
        for (Int j = 0, cnt = 0; j < np; ++j)
          for (Int i = 0; i < np; ++i, ++cnt) {
            const Real rhomij = dp3d_c(i,j,k,d.tl_np1) * spheremp(i,j);
            rhom += rhomij;
            wa[cnt] = rhomij;
            qlo[cnt] = q_min(i,j,k,q,ie);
            qhi[cnt] = q_max(i,j,k,q,ie);
            y[cnt] = q_c(i,j,k,q);
            x[cnt] = y[cnt];
          }
        const Real Qm = qlt.qlt->get_Q(lci, ti);

        int info = qlt::slv::solve_1eq_bc_qp(
          np2, wa, wa, Qm, qlo, qhi, scalar_bounds, y, x);
        if (info < 0) {
          if (scalar_bounds) {
            const Real q = Qm / rhom;
            *qlo = std::min(*qlo, q);
            *qhi = std::max(*qhi, q);
            info = qlt::slv::solve_1eq_bc_qp(
              np2, wa, wa, Qm, qlo, qhi, true, y, x);
          } else {
            if (false) {
              //todo Need to think about this more. The sign of Qm_extra is not
              // enough to determine which bound is the problem.
              Real Qm_before_qlt = 0;
              for (Int j = 0; j < np; ++j)
                for (Int i = 0; i < np; ++i) {
                  const Real rhomij = dp3d_c(i,j,k,d.tl_np1) * spheremp(i,j);
                  Qm_before_qlt += q_c(i,j,k,q) * rhomij;
                }
              qlt::slv::r2l_l_adjust_bounds(np2, qlo, qhi, wa, Qm - Qm_before_qlt);
              info = qlt::slv::solve_1eq_bc_qp(np2, wa, wa, Qm, qlo, qhi, false, y, x);

#pragma message "CHECK THIS 2"
              if (info < 0) {
                Real Qm_min = 0, Qm_max = 0;
                Real qlo_orig[max_np2], qhi_orig[max_np2];
                for (Int j = 0, os = 0; j < np; ++j)
                  for (Int i = 0; i < np; ++i, ++os) {
                    const Real rhomij = dp3d_c(i,j,k,d.tl_np1) * spheremp(i,j);
                    Qm_min += qlo[os] * rhomij;
                    Qm_max += qhi[os] * rhomij;
                    qlo_orig[os] = q_min(i,j,k,q,ie);
                    qhi_orig[os] = q_max(i,j,k,q,ie);
                  }
                using qlt::impl::prarr;
                const Real Qm_extra = Qm - Qm_before_qlt;
                pr(puf(ie) pu(k) pu(q) pu(ti) pu(info)
                   pu(Qm_min) pu(Qm-Qm_min) pu(Qm) pu(Qm_extra) pu(Qm_max) pu(Qm_max-Qm));
                prarr("wa", wa, np2);
                prarr("qlo_orig", qlo_orig, np2);
                prarr("qlo", qlo, np2);
                prarr("qhi_orig", qhi_orig, np2);
                prarr("qhi", qhi, np2);
                prarr("y", y, np2);
                prarr("x", x, np2);
                exit(-1);
              }
            } else {
              //todo For now, just make sure we can accommodate the cell
              // mean. QLT has guaranteed that this process does not violate the
              // safety problem; it just isn't as nice as precise as the
              // r2l_l_adjust_bounds method.
              const Real q = Qm / rhom;
              for (Int i = 0; i < np2; ++i) qlo[i] = std::min(qlo[i], q);
              for (Int i = 0; i < np2; ++i) qhi[i] = std::max(qhi[i], q);
              info = qlt::slv::solve_1eq_bc_qp(
                np2, wa, wa, Qm, qlo, qhi, true, y, x);
            }
          }
        }
        
        for (Int j = 0, cnt = 0; j < np; ++j)
          for (Int i = 0; i < np; ++i, ++cnt) {
            q_c(i,j,k,q) = x[cnt];
            qdp_c(i,j,k,q,d.n1_qdp) = q_c(i,j,k,q) * dp3d_c(i,j,k,d.tl_np1);
          }
      }
  }
}

void check (QLT& qlt, const Data& d, const Real* q_min_r, const Real* q_max_r) {
  using qlt::mpi::reduce;

  const Int np = d.np, nlev = d.nlev, qsize = d.qsize, ncell = d.spheremp.size();

  Kokkos::View<Real**, Kokkos::HostSpace>
    mass_p_g("mass_p_g", nlev, qsize), mass_c_g("mass_c_g", nlev, qsize),
    mass_lo_g("mass_lo_g", nlev, qsize), mass_hi_g("mass_hi_g", nlev, qsize),
    q_lo_g("q_lo_g", nlev, qsize), q_hi_g("q_hi_g", nlev, qsize),
    q_min_g("q_min_g", nlev, qsize), q_max_g("q_max_g", nlev, qsize),
    qd_lo_g("qd_lo_g", nlev, qsize), qd_hi_g("qd_hi_g", nlev, qsize);
  {
    Kokkos::View<Real**, Kokkos::HostSpace>
      mass_p("mass_p", nlev, qsize), mass_c("mass_c", nlev, qsize),
      mass_lo("mass_lo", nlev, qsize), mass_hi("mass_hi", nlev, qsize),
      q_lo("q_lo", nlev, qsize), q_hi("q_hi", nlev, qsize),
      q_min_l("q_min_l", nlev, qsize), q_max_l("q_max_l", nlev, qsize),
      qd_lo("qd_lo", nlev, qsize), qd_hi("qd_hi", nlev, qsize);
    FA5<const Real>
      q_min(q_min_r, np, np, nlev, qsize, ncell),
      q_max(q_max_r, np, np, nlev, qsize, ncell);
    Kokkos::deep_copy(q_lo,  1e200);
    Kokkos::deep_copy(q_hi, -1e200);
    Kokkos::deep_copy(q_min_l,  1e200);
    Kokkos::deep_copy(q_max_l, -1e200);
    Kokkos::deep_copy(qd_lo, 0);
    Kokkos::deep_copy(qd_hi, 0);

    for (Int ie = 0; ie < ncell; ++ie) {
      FA2<const Real> spheremp(d.spheremp[ie], np, np);
      FA5<const Real> qdp_pc(d.qdp_pc[ie], np, np, nlev, d.qsize_d, 2);
      FA4<const Real> dp3d_c(d.dp3d_c[ie], np, np, nlev, d.timelevels);
      FA4<const Real> q_c(d.q_c[ie], np, np, nlev, d.qsize_d);
      for (Int k = 0; k < nlev; ++k)
        for (Int q = 0; q < qsize; ++q)
          for (Int j = 0; j < np; ++j)
            for (Int i = 0; i < np; ++i) {
              mass_p(k,q) += qdp_pc(i,j,k,q,d.n0_qdp) * spheremp(i,j);
              mass_c(k,q) += qdp_pc(i,j,k,q,d.n1_qdp) * spheremp(i,j);
              if (q_c(i,j,k,q) < q_min(i,j,k,q,ie))
                qd_lo(k,q) = std::max(qd_lo(k,q), q_min(i,j,k,q,ie) - q_c(i,j,k,q));
              if (q_c(i,j,k,q) > q_max(i,j,k,q,ie))
                qd_hi(k,q) = std::max(qd_hi(k,q), q_c(i,j,k,q) - q_max(i,j,k,q,ie));
              // Safety problem.
              mass_lo(k,q) += q_min(i,j,k,q,ie) * dp3d_c(i,j,k,d.tl_np1) * spheremp(i,j);
              mass_hi(k,q) += q_max(i,j,k,q,ie) * dp3d_c(i,j,k,d.tl_np1) * spheremp(i,j);
              q_lo(k,q) = std::min(q_lo(k,q), q_min(i,j,k,q,ie));
              q_hi(k,q) = std::max(q_hi(k,q), q_max(i,j,k,q,ie));
              q_min_l(k,q) = std::min(q_min_l(k,q), q_min(i,j,k,q,ie));
              q_max_l(k,q) = std::max(q_max_l(k,q), q_max(i,j,k,q,ie));
            }
    }

    reduce(*qlt.p, mass_p.data(), mass_p_g.data(), nlev*qsize, MPI_SUM, qlt.p->root());
    reduce(*qlt.p, mass_c.data(), mass_c_g.data(), nlev*qsize, MPI_SUM, qlt.p->root());
    reduce(*qlt.p, qd_lo.data(), qd_lo_g.data(), nlev*qsize, MPI_MAX, qlt.p->root());
    reduce(*qlt.p, qd_hi.data(), qd_hi_g.data(), nlev*qsize, MPI_MAX, qlt.p->root());
    // Safety problem.
    reduce(*qlt.p, mass_lo.data(), mass_lo_g.data(), nlev*qsize, MPI_SUM, qlt.p->root());
    reduce(*qlt.p, mass_hi.data(), mass_hi_g.data(), nlev*qsize, MPI_SUM, qlt.p->root());
    reduce(*qlt.p, q_lo.data(), q_lo_g.data(), nlev*qsize, MPI_MIN, qlt.p->root());
    reduce(*qlt.p, q_hi.data(), q_hi_g.data(), nlev*qsize, MPI_MAX, qlt.p->root());
    reduce(*qlt.p, q_min_l.data(), q_min_g.data(), nlev*qsize, MPI_MIN, qlt.p->root());
    reduce(*qlt.p, q_max_l.data(), q_max_g.data(), nlev*qsize, MPI_MAX, qlt.p->root());
  }

  if (qlt.p->amroot()) {
    const Real tol = 1e4*std::numeric_limits<Real>::epsilon();
    for (Int k = 0; k < nlev; ++k)
      for (Int q = 0; q < qsize; ++q) {
        const Real rd = qlt::impl::reldif(mass_p_g(k,q), mass_c_g(k,q));
        if (rd > tol)
          pr(puf(k) pu(q) pu(mass_p_g(k,q)) pu(mass_c_g(k,q)) pu(rd));
        if (mass_lo_g(k,q) <= mass_c_g(k,q) && mass_c_g(k,q) <= mass_hi_g(k,q)) {
          // Local problems should be feasible.
          if (qd_lo_g(k,q) > 0)
            pr(puf(k) pu(q) pu(qd_lo_g(k,q)));
          if (qd_hi_g(k,q) > 0)
            pr(puf(k) pu(q) pu(qd_hi_g(k,q)));
        } else {
          // Safety problem must hold.
          if (q_lo_g(k,q) < q_min_g(k,q))
            pr(puf(k) pu(q) pu(q_lo_g(k,q) - q_min_g(k,q)) pu(q_min_g(k,q)));
          if (q_hi_g(k,q) > q_max_g(k,q))
            pr(puf(k) pu(q) pu(q_max_g(k,q) - q_hi_g(k,q)) pu(q_max_g(k,q)));
        }
      }
  }
}
} // namespace sl
} // namespace homme

// Interface for Homme, through compose_mod.F90.
extern "C" void kokkos_init_ () { Kokkos::initialize(); }

extern "C" void kokkos_finalize_ () { Kokkos::finalize_all(); }

static homme::QLT::Ptr g_qlt;

extern "C" void
qlt_init_impl_ (const homme::Int* fcomm, const homme::Int** sc2gci,
                const homme::Int** sc2rank, const homme::Int* gbl_ncell) {
  const auto p = qlt::make_parallel(MPI_Comm_f2c(*fcomm));
  g_qlt = std::make_shared<homme::QLT>(*gbl_ncell, *sc2gci, *sc2rank, p);
}

extern "C" void
qlt_unittest_ (const homme::Int* fcomm, homme::Int* nerrp) {
  qlt_assert(g_qlt);
  qlt_assert(g_qlt->tree);
  auto p = qlt::make_parallel(MPI_Comm_f2c(*fcomm));
  *nerrp = qlt::test::TestQLT(p, g_qlt->tree, g_qlt->ncell).run();
}

extern "C" void
qlt_set_ie2gci_ (const homme::Int* ie, const homme::Int* gci) {
  qlt_assert(g_qlt);
  // Now is a good time to drop the tree, whose persistence was used for unit
  // testing if at all.
  g_qlt->tree = nullptr;
  homme::set_ie2gci(*g_qlt, *ie - 1, *gci - 1);
}

static homme::sl::Data::Ptr g_sl;

//todo Need to be tid aware if HORIZ_OPENMP.
extern "C" homme::Int qlt_sl_set_pointers_begin_ (
  homme::Int* nets, homme::Int* nete, const homme::Int* np, const homme::Int* nlev,
  const homme::Int* qsize, const homme::Int* qsized, const homme::Int* timelevels,
  const homme::Int* need_conservation)
{
  if (g_sl) return 0;
  qlt_assert(g_qlt);
  g_sl = std::make_shared<homme::sl::Data>(g_qlt->qlt->nlclcells(), *np, *nlev, *qsize,
                                           *qsized, *timelevels);
  homme::init_ie2lci(*g_qlt);
  homme::init_tracers(*g_qlt, *nlev, *qsize, *need_conservation);
  homme::sl::check(*g_qlt, *g_sl);
  return 1;
}

extern "C" void qlt_sl_set_spheremp_ (homme::Int* ie, homme::Real* v)
{ homme::sl::insert(g_sl, *ie - 1, 0, v); }
extern "C" void qlt_sl_set_qdp_ (homme::Int* ie, homme::Real* v, homme::Int* n0_qdp,
                                 homme::Int* n1_qdp)
{ homme::sl::insert(g_sl, *ie - 1, 1, v, *n0_qdp - 1, *n1_qdp - 1); }
extern "C" void qlt_sl_set_dp3d_ (homme::Int* ie, homme::Real* v, homme::Int* tl_np1)
{ homme::sl::insert(g_sl, *ie - 1, 2, v, *tl_np1 - 1); }
extern "C" void qlt_sl_set_q_ (homme::Int* ie, homme::Real* v)
{ homme::sl::insert(g_sl, *ie - 1, 3, v); }

extern "C" void qlt_sl_set_pointers_end_ () {}

// Run QLT.
extern "C" void qlt_sl_run_ (const homme::Real* minq, const homme::Real* maxq,
                             homme::Int*, homme::Int*) {
  qlt_assert(minq != maxq);
  qlt_assert(g_qlt);
  qlt_assert(g_sl);
  homme::sl::run(*g_qlt, *g_sl, minq, maxq);
}

// Run the cell-local limiter problem.
extern "C" void qlt_sl_run_local_ (const homme::Real* minq, const homme::Real* maxq,
                                   homme::Int*, homme::Int*, homme::Int* use_ir) {
  qlt_assert(minq != maxq);
  qlt_assert(g_qlt);
  qlt_assert(g_sl);
  homme::sl::run_local(*g_qlt, *g_sl, minq, maxq, *use_ir);
}

// Check properties for this transport step.
extern "C" void qlt_sl_check_ (const homme::Real* minq, const homme::Real* maxq,
                               homme::Int*, homme::Int*) {
  qlt_assert(g_qlt);
  qlt_assert(g_sl);
  homme::sl::check(*g_qlt, *g_sl, minq, maxq);
}
