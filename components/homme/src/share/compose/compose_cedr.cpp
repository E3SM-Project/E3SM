// Uncomment this to look for MPI-related memory leaks.
//#define COMPOSE_DEBUG_MPI

//>> cedr_kokkos.hpp
// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

#ifndef INCLUDE_CEDR_KOKKOS_HPP
#define INCLUDE_CEDR_KOKKOS_HPP

#include <Kokkos_Core.hpp>

#define KIF KOKKOS_INLINE_FUNCTION

// Clarify that a class member type is meant to be private but is
// marked public for Cuda visibility.
#define PRIVATE_CUDA public
#define PROTECTED_CUDA public

#if defined KOKKOS_COMPILER_GNU
// See https://github.com/kokkos/kokkos-kernels/issues/129 
# define ConstExceptGnu
#else
# define ConstExceptGnu const
#endif

namespace cedr {
namespace impl {

// Turn a View's MemoryTraits (traits::memory_traits) into the equivalent
// unsigned int mask.
template <typename View>
struct MemoryTraitsMask {
  enum : unsigned int {
#ifndef KOKKOS_VERSION
    // pre-v3
    value = ((View::traits::memory_traits::RandomAccess ? Kokkos::RandomAccess : 0) |
             (View::traits::memory_traits::Atomic ? Kokkos::Atomic : 0) |
             (View::traits::memory_traits::Restrict ? Kokkos::Restrict : 0) |
             (View::traits::memory_traits::Aligned ? Kokkos::Aligned : 0) |
             (View::traits::memory_traits::Unmanaged ? Kokkos::Unmanaged : 0))
#else
    // >= v3
    value = ((View::traits::memory_traits::is_random_access ? Kokkos::RandomAccess : 0) |
             (View::traits::memory_traits::is_atomic ? Kokkos::Atomic : 0) |
             (View::traits::memory_traits::is_restrict ? Kokkos::Restrict : 0) |
             (View::traits::memory_traits::is_aligned ? Kokkos::Aligned : 0) |
             (View::traits::memory_traits::is_unmanaged ? Kokkos::Unmanaged : 0))
#endif
      };
};

// Make the input View Unmanaged, whether or not it already is. One might
// imagine that View::unmanaged_type would provide this.
//   Use: Unmanaged<ViewType>
template <typename View>
using Unmanaged =
  // Provide a full View type specification, augmented with Unmanaged.
  Kokkos::View<typename View::traits::scalar_array_type,
               typename View::traits::array_layout,
               typename View::traits::device_type,
               Kokkos::MemoryTraits<
                 // All the current values...
                 MemoryTraitsMask<View>::value |
                 // ... |ed with the one we want, whether or not it's
                 // already there.
                 Kokkos::Unmanaged> >;

template <typename View>
using Const = typename View::const_type;

template <typename View>
using ConstUnmanaged = Const<Unmanaged<View> >;

template <typename ExeSpace>
struct DeviceType {
  typedef Kokkos::Device<typename ExeSpace::execution_space,
                         typename ExeSpace::memory_space> type;
};

#ifdef KOKKOS_HAVE_CUDA
typedef Kokkos::Device<Kokkos::CudaSpace::execution_space,
                       Kokkos::CudaSpace::memory_space> DefaultDeviceType;

template <> struct DeviceType<Kokkos::Cuda> {
  typedef DefaultDeviceType type;
};
#else
typedef Kokkos::Device<Kokkos::DefaultExecutionSpace::execution_space,
                       Kokkos::DefaultExecutionSpace::memory_space> DefaultDeviceType;
#endif

template <typename ES> struct OnGpu {
  enum : bool { value =
#ifdef COMPOSE_MIMIC_GPU
                true
#else
                false
#endif
  };
};
#ifdef KOKKOS_ENABLE_CUDA
template <> struct OnGpu<Kokkos::Cuda> { enum : bool { value = true }; };
#endif

template <typename ExeSpace = Kokkos::DefaultExecutionSpace>
struct ExeSpaceUtils {
  using TeamPolicy = Kokkos::TeamPolicy<ExeSpace>;
  using Member = typename TeamPolicy::member_type;
  static TeamPolicy get_default_team_policy (int outer, int inner) {
#ifdef COMPOSE_MIMIC_GPU
    const int max_threads =
#ifdef KOKKOS_ENABLE_OPENMP
      ExeSpace::concurrency()
#else
      1
#endif
      ;
    const int team_size = max_threads < 7 ? max_threads : 7;
    return TeamPolicy(outer, team_size, 1);
#else
    return TeamPolicy(outer, 1, 1);
#endif
}
};

#ifdef KOKKOS_ENABLE_CUDA
template <>
struct ExeSpaceUtils<Kokkos::Cuda> {
  using TeamPolicy = Kokkos::TeamPolicy<Kokkos::Cuda>;
  using Member = typename TeamPolicy::member_type;
  static TeamPolicy get_default_team_policy (int outer, int inner) {
    return TeamPolicy(outer, std::min(128, 32*((inner + 31)/32)), 1);
  }
};
#endif

// GPU-friendly replacements for std::*.
template <typename T> KOKKOS_INLINE_FUNCTION
const T& min (const T& a, const T& b) { return a < b ? a : b; }
template <typename T> KOKKOS_INLINE_FUNCTION
const T& max (const T& a, const T& b) { return a > b ? a : b; }
template <typename T> KOKKOS_INLINE_FUNCTION
void swap (T& a, T& b) { const T tmp = a; a = b; b = tmp; }

} // namespace impl
} // namespace cedr

#endif

//>> cedr.hpp
// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

#ifndef INCLUDE_CEDR_HPP
#define INCLUDE_CEDR_HPP

//#include "cedr_kokkos.hpp"

// Communication-Efficient Constrained Density Reconstructors
namespace cedr {
typedef int Int;
typedef long int Long;
typedef std::size_t Size;
typedef double Real;

// CDRs in general implement
// * tracer mass, Qm, conservation;
// * mixing ratio, q, shape preservation: any of local bound preservation,
//   dynamic range preservation, or simply non-negativity; and
// * tracer consistency, which follows from dynamic range preservation or
//   stronger (including local bound preservation) with rhom coming from the
//   dynamics.
//
// One can solve a subset of these.
//   If !conserve, then the CDR does not alter the tracer mass, but it does not
// correct for any failure in mass conservation in the field given to it.
//   If consistent but !shapepreserve, then the CDR solves the dynamic range
// preservation problem rather than the local bound preservation problem.
struct ProblemType {
  enum : Int {
    conserve = 1, shapepreserve = 1 << 1, consistent = 1 << 2,
    // The 'nonnegative' problem type can be combined only with 'conserve'. The
    // caller can implement nonnegativity when running with 'shapepreserve' or
    // 'consistent' simply by setting Qm_min = 0. The 'nonnegativity' type is
    // reserved for a particularly efficient type of problem in which
    // Qm_{min,max} are not specified.
    nonnegative = 1 << 3
  };
};
}

#endif

//>> cedr_mpi.hpp
// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

#ifndef INCLUDE_CEDR_MPI_HPP
#define INCLUDE_CEDR_MPI_HPP

#include <memory>

#include <mpi.h>

//#include "compose_config.hpp"
//#include "cedr.hpp"

namespace cedr {
namespace mpi {

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

struct Request {
  MPI_Request request;

#ifdef COMPOSE_DEBUG_MPI
  int unfreed;
  Request();
  ~Request();
#endif
};

Parallel::Ptr make_parallel(MPI_Comm comm);

template <typename T> MPI_Datatype get_type();

template <typename T>
int reduce(const Parallel& p, const T* sendbuf, T* rcvbuf, int count, MPI_Op op,
           int root);

template <typename T>
int all_reduce(const Parallel& p, const T* sendbuf, T* rcvbuf, int count, MPI_Op op);

template <typename T>
int isend(const Parallel& p, const T* buf, int count, int dest, int tag,
          Request* ireq = nullptr);

template <typename T>
int irecv(const Parallel& p, T* buf, int count, int src, int tag,
          Request* ireq = nullptr);

int waitany(int count, Request* reqs, int* index, MPI_Status* stats = nullptr);

int waitall(int count, Request* reqs, MPI_Status* stats = nullptr);

template<typename T>
int gather(const Parallel& p, const T* sendbuf, int sendcount,
           T* recvbuf, int recvcount, int root);

template <typename T>
int gatherv(const Parallel& p, const T* sendbuf, int sendcount,
            T* recvbuf, const int* recvcounts, const int* displs, int root);

bool all_ok(const Parallel& p, bool im_ok);

struct Op {
  typedef std::shared_ptr<Op> Ptr;

  Op (MPI_User_function* function, bool commute) {
    MPI_Op_create(function, static_cast<int>(commute), &op_);
  }

  ~Op () { MPI_Op_free(&op_); }

  const MPI_Op& get () const { return op_; }

private:
  MPI_Op op_;
};

} // namespace mpi
} // namespace cedr

//#include "cedr_mpi_inl.hpp"

#endif

//>> cedr_util.hpp
// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

#ifndef INCLUDE_CEDR_UTIL_HPP
#define INCLUDE_CEDR_UTIL_HPP

#include <sstream>

//#include "cedr_kokkos.hpp"
//#include "cedr_mpi.hpp"

namespace cedr {
namespace util {

template <typename T> KOKKOS_INLINE_FUNCTION constexpr
T square (const T& x) { return x*x; }

bool eq(const std::string& a, const char* const b1, const char* const b2 = 0);

// Uniform rand in [0, 1).
Real urand();

#define pr(m) do {                                      \
    int _pid_ = 0;                                      \
    MPI_Comm_rank(MPI_COMM_WORLD, &_pid_);              \
    std::stringstream _ss_;                             \
    _ss_.precision(15);                                 \
    _ss_ << "pid " << _pid_ << " " << m << std::endl;   \
    std::cerr << _ss_.str();                            \
  } while (0)
#define pr0(m) do {                                     \
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
#define mprarr(m) cedr::util::prarr(#m, m.data(), m.size())

#ifndef NDEBUG
# define cedr_assert(condition) do {                                    \
    if ( ! (condition)) {                                               \
      std::stringstream _ss_;                                           \
      _ss_ << __FILE__ << ":" << __LINE__ << ": FAIL:\n" << #condition  \
        << "\n";                                                        \
        throw std::logic_error(_ss_.str());                             \
    }                                                                   \
  } while (0)
# define cedr_kernel_assert(condition) do {     \
    if ( ! (condition))                         \
      Kokkos::abort(#condition);                \
  } while (0)
#else
# define cedr_assert(condition)
# define cedr_kernel_assert(condition)
#endif
#define cedr_throw_if(condition, message) do {                          \
    if (condition) {                                                    \
      std::stringstream _ss_;                                           \
      _ss_ << __FILE__ << ":" << __LINE__ << ": The condition:\n"       \
           << #condition "\nled to the exception\n" << message << "\n"; \
        throw std::logic_error(_ss_.str());                             \
    }                                                                   \
  } while (0)
#define cedr_kernel_throw_if(condition, message) do {                   \
    if (condition)                                                      \
      Kokkos::abort(#condition " led to the exception\n" message);      \
  } while (0)

inline Real reldif (const Real a, const Real b)
{ return std::abs(b - a)/std::max(std::abs(a), std::abs(b)); }

Real reldif(const Real* a, const Real* b, const Int n);

struct FILECloser { void operator() (FILE* fh) { fclose(fh); } };

template <typename T, typename ExeSpace>
struct RawArrayRaft {
  typedef typename cedr::impl::DeviceType<ExeSpace>::type Device;
  typedef Kokkos::View<T*, Device> List;
  
  RawArrayRaft (T* a, const Int n)
    : a_(a), n_(n)
  {
    a_d_ = List("RawArrayRaft::a_d_", n_);
    a_h_ = typename List::HostMirror(a_, n_);
  }

  const List& sync_device () {
    Kokkos::deep_copy(a_d_, a_h_);
    return a_d_;
  }

  T* sync_host () {
    Kokkos::deep_copy(a_h_, a_d_);
    return a_;
  }

  T* device_ptr () { return a_d_.data(); }

private:
  T* a_;
  Int n_;
  List a_d_;
  typename List::HostMirror a_h_;
};

}
}

#endif

//>> cedr_cdr.hpp
// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

#ifndef INCLUDE_CEDR_CDR_HPP
#define INCLUDE_CEDR_CDR_HPP

//#include "cedr_mpi.hpp"

namespace cedr {
// Constrained Density Reconstructor interface.
//   Find a point Qm in the set
//      { Qm: ( i) e'Qm = Qm_global
//            (ii) Qm_min <= Qm <= Qm_max },
// where e is the vector of 1s. Each algorithm in CEDR has its own optimality
// principle to decide on the point.
struct CDR {
  typedef std::shared_ptr<CDR> Ptr;

  struct Options {
    // In exact arithmetic, mass is conserved exactly and bounds of a feasible
    // problem are not violated. In inexact arithmetic, these properties are
    // inexact, and approximation quality of one can be preferred to that of the
    // other. Use this option to prefer numerical mass conservation over
    // numerical bounds nonviolation. In practice, this option governs errrs at
    // the level of 10 to 1000 times numeric_limits<Real>::epsilon().
    bool prefer_numerical_mass_conservation_to_numerical_bounds;

    Options ()
      : prefer_numerical_mass_conservation_to_numerical_bounds(false)
    {}
  };

  CDR (const Options options = Options()) : options_(options) {}

  virtual void print(std::ostream& os) const {}

  const Options& get_options () const { return options_; }

  // Set up QLT tracer metadata. Call declare_tracer in order of the tracer
  // index in the caller's numbering. Once end_tracer_declarations is called, it
  // is an error to call declare_tracer again.
  //   Associate the tracer with a rhom index. In many problems, there will be
  // only one rhom, so rhomidx is always 0.
  //   It is an error to call this function from a parallel region.
  virtual void declare_tracer(int problem_type, const Int& rhomidx) = 0;

  // It is an error to call this function from a parallel region.
  virtual void end_tracer_declarations() = 0;

  // Optionally query the sizes of the two primary memory buffers. This
  // operation is valid after end_tracer_declarations was called.
  virtual void get_buffers_sizes(size_t& buf1, size_t& buf2) = 0;
  // Optionally provide the two primary memory buffers. These pointers must
  // point to memory having at least the sizes returned by get_buffer_sizes.
  virtual void set_buffers(Real* buf1, Real* buf2) = 0;

  // Call this method to finish the setup phase. The subsequent methods are
  // valid only after this method has been called.
  virtual void finish_setup() = 0;

  virtual int get_problem_type(const Int& tracer_idx) const = 0;

  virtual Int get_num_tracers() const = 0;

  // set_{rhom,Qm}: Set cell values prior to running the QLT algorithm.
  //
  //   Notation:
  //     rho: Total density.
  //       Q: Tracer density.
  //       q: Tracer mixing ratio = Q/rho.
  //      *m: Mass corresponding to the density; results from an integral over a
  //          region, such as a cell.
  //   Some CDRs have a nontrivial local <-> global cell index map. For these
  // CDRs, lclcellidx may be nontrivial. For others, the caller should provide
  // the index into the local cell.
  //
  //   set_rhom must be called before set_Qm.
  KOKKOS_FUNCTION
  virtual void set_rhom(
    const Int& lclcellidx, const Int& rhomidx,
    // Current total mass in this cell.
    const Real& rhom) const = 0;

  KOKKOS_FUNCTION
  virtual void set_Qm(
    const Int& lclcellidx, const Int& tracer_idx,
    // Current tracer mass in this cell.
    const Real& Qm,
    // Minimum and maximum permitted tracer mass in this cell. Ignored if
    // ProblemType is 'nonnegative'.
    const Real& Qm_min, const Real& Qm_max,
    // If mass conservation is requested, provide the previous Qm, which will be
    // summed to give the desired global mass.
    const Real Qm_prev = std::numeric_limits<Real>::infinity()) const = 0;

  // Run the QLT algorithm with the values set by set_{rho,Q}. It is an error to
  // call this function from a parallel region.
  virtual void run() = 0;

  // Get a cell's tracer mass Qm after the QLT algorithm has run.
  KOKKOS_FUNCTION
  virtual Real get_Qm(const Int& lclcellidx, const Int& tracer_idx) const = 0;

protected:
  Options options_;
};
} // namespace cedr

#endif

//>> cedr_qlt.hpp
// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

#ifndef INCLUDE_CEDR_QLT_HPP
#define INCLUDE_CEDR_QLT_HPP

#include <mpi.h>

#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <list>

//#include "cedr_cdr.hpp"
//#include "cedr_util.hpp"

namespace cedr {
// QLT: Quasi-local tree-based non-iterative tracer density reconstructor for
//      mass conservation, shape preservation, and tracer consistency.
namespace qlt {
using cedr::mpi::Parallel;

namespace impl {
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
  
  // Levels. nodes[0] is level 0, the leaf level.
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

typedef impl::NodeSetsDeviceData<Kokkos::DefaultHostExecutionSpace> NodeSetsHostData;
} // namespace impl

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
Node::Ptr make_tree_over_1d_mesh(const Parallel::Ptr& p, const Int& ncells,
                                 const bool imbalanced = false);
} // namespace tree

template <typename ExeSpace = Kokkos::DefaultExecutionSpace>
class QLT : public cedr::CDR {
public:
  typedef typename cedr::impl::DeviceType<ExeSpace>::type Device;
  typedef QLT<ExeSpace> Me;
  typedef std::shared_ptr<Me> Ptr;
  typedef Kokkos::View<Real*, Device> RealList;
  
  // Set up QLT topology and communication data structures based on a tree. Both
  // ncells and tree refer to the global mesh, not just this processor's
  // part. The tree must be identical across ranks.
  QLT(const Parallel::Ptr& p, const Int& ncells, const tree::Node::Ptr& tree,
      CDR::Options options = Options());

  void print(std::ostream& os) const override;

  // Number of cells owned by this rank.
  Int nlclcells() const;

  // Cells owned by this rank, in order of local numbering. Thus,
  // gci2lci(gcis[i]) == i. Ideally, the caller never actually calls gci2lci(),
  // and instead uses the information from get_owned_glblcells to determine
  // local cell indices.
  void get_owned_glblcells(std::vector<Long>& gcis) const;

  // For global cell index cellidx, i.e., the globally unique ordinal associated
  // with a cell in the caller's tree, return this rank's local index for
  // it. This is not an efficient operation.
  Int gci2lci(const Int& gci) const;

  void declare_tracer(int problem_type, const Int& rhomidx) override;

  void end_tracer_declarations() override;

  void get_buffers_sizes(size_t& buf1, size_t& buf2) override;

  void set_buffers(Real* buf1, Real* buf2) override;

  void finish_setup() override;

  int get_problem_type(const Int& tracer_idx) const override;

  Int get_num_tracers() const override;

  // lclcellidx is gci2lci(cellidx).
  KOKKOS_INLINE_FUNCTION
  void set_rhom(const Int& lclcellidx, const Int& rhomidx, const Real& rhom) const override;

  // lclcellidx is gci2lci(cellidx).
  KOKKOS_INLINE_FUNCTION
  void set_Qm(const Int& lclcellidx, const Int& tracer_idx,
              const Real& Qm, const Real& Qm_min, const Real& Qm_max,
              const Real Qm_prev = std::numeric_limits<Real>::infinity()) const override;

  void run() override;

  KOKKOS_INLINE_FUNCTION
  Real get_Qm(const Int& lclcellidx, const Int& tracer_idx) const override;

protected:
  typedef Kokkos::View<Int*, Device> IntList;
  typedef cedr::impl::Const<IntList> ConstIntList;
  typedef cedr::impl::ConstUnmanaged<IntList> ConstUnmanagedIntList;

  static void init(const std::string& name, IntList& d,
                   typename IntList::HostMirror& h, size_t n);

  struct MetaDataBuilder {
    typedef std::shared_ptr<MetaDataBuilder> Ptr;
    std::vector<int> trcr2prob;
  };

PROTECTED_CUDA:
  struct MetaData {
    enum : Int { nprobtypes = 6 };

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

    KOKKOS_INLINE_FUNCTION static int get_problem_type(const int& idx);
    
    // icpc doesn't let us use problem_type_ here, even though it's constexpr.
    static int get_problem_type_idx(const int& mask);

    KOKKOS_INLINE_FUNCTION static int get_problem_type_l2r_bulk_size(const int& mask);

    static int get_problem_type_r2l_bulk_size(const int& mask);

    struct CPT {
      //   The only problem not supported is conservation alone. It makes very
      // little sense to use QLT for conservation alone.
      //   The remaining problems fall into 6 categories of details. These
      // categories are tracked by QLT; which of the original problems being
      // solved is not important.
      enum {
        // l2r: rhom, (Qm_min, Qm, Qm_max)*; r2l: Qm*
        s  = ProblemType::shapepreserve,
        st = ProblemType::shapepreserve | ProblemType::consistent,
        // l2r: rhom, (Qm_min, Qm, Qm_max, Qm_prev)*; r2l: Qm*
        cs  = ProblemType::conserve | s,
        cst = ProblemType::conserve | st,
        // l2r: rhom, (q_min, Qm, q_max)*; r2l: (Qm, q_min, q_max)*
        t = ProblemType::consistent,
        // l2r: rhom, (q_min, Qm, q_max, Qm_prev)*; r2l: (Qm, q_min, q_max)*
        ct = ProblemType::conserve | t,
        // l2r: rhom, Qm*; r2l: Qm*
        nn = ProblemType::nonnegative,
        // l2r: rhom, (Qm, Qm_prev)*; r2l: Qm*
        cnn = ProblemType::conserve | nn
      };
    };

    Arrays<typename ConstUnmanagedIntList::HostMirror> a_h;
    Arrays<ConstUnmanagedIntList> a_d;

    void init(const MetaDataBuilder& mdb);

  private:
    Arrays<typename IntList::HostMirror> a_h_;
    Arrays<IntList> a_d_;
  };

  struct BulkData {
    typedef cedr::impl::Unmanaged<RealList> UnmanagedRealList;

    UnmanagedRealList l2r_data, r2l_data;

    BulkData () : inited_(false) {}

    bool inited () const { return inited_; }

    void init(const size_t& l2r_sz, const size_t& r2l_sz);
    void init(Real* l2r_buf, const size_t& l2r_sz,
              Real* r2l_buf, const size_t& r2l_sz);

  private:
    bool inited_;
    RealList l2r_data_, r2l_data_;
  };

protected:
  void init(const Parallel::Ptr& p, const Int& ncells, const tree::Node::Ptr& tree);

  void init_ordinals();

  /// Pointer data for initialization and host computation.
  Parallel::Ptr p_;
  // Tree and communication topology.
  std::shared_ptr<const impl::NodeSets> ns_;
  // Data extracted from ns_ for use in run() on device.
  std::shared_ptr<impl::NodeSetsDeviceData<ExeSpace> > nsdd_;
  std::shared_ptr<impl::NodeSetsHostData> nshd_;
  // Globally unique cellidx -> rank-local index.
  typedef std::map<Int,Int> Gci2LciMap;
  std::shared_ptr<Gci2LciMap> gci2lci_;
  // Temporary to collect caller's tracer information prior to calling
  // end_tracer_declarations().
  typename MetaDataBuilder::Ptr mdb_;
  /// View data for host and device computation.
  // Constructed in end_tracer_declarations().
  MetaData md_;
  BulkData bd_;

PRIVATE_CUDA:
  void l2r_recv(const impl::NodeSets::Level& lvl, const Int& l2rndps) const;
  void l2r_combine_kid_data(const Int& lvlidx, const Int& l2rndps) const;
  void l2r_send_to_parents(const impl::NodeSets::Level& lvl, const Int& l2rndps) const;
  void root_compute(const Int& l2rndps, const Int& r2lndps) const;
  void r2l_recv(const impl::NodeSets::Level& lvl, const Int& r2lndps) const;
  void r2l_solve_qp(const Int& lvlidx, const Int& l2rndps, const Int& r2lndps) const;
  void r2l_send_to_kids(const impl::NodeSets::Level& lvl, const Int& r2lndps) const;
};

namespace test {
struct Input {
  bool unittest, perftest, write;
  Int ncells, ntracers, tracer_type, nrepeat;
  bool pseudorandom, verbose;
};

Int run_unit_and_randomized_tests(const Parallel::Ptr& p, const Input& in);

Int test_qlt(const Parallel::Ptr& p, const tree::Node::Ptr& tree, const Int& ncells,
             const Int nrepeat,
             // Diagnostic output for dev and illustration purposes. To be
             // clear, no QLT unit test requires output to be checked; each
             // checks in-memory data and returns a failure count.
             const bool write,
             // Provide memory to QLT for its buffers.
             const bool external_memory,
             // Set CDR::Options.prefer_numerical_mass_conservation_to_numerical_bounds.
             const bool prefer_mass_con_to_bounds,
             const bool verbose);
} // namespace test
} // namespace qlt
} // namespace cedr

// These are the definitions that must be visible in the calling translation
// unit, unless Cuda relocatable device code is enabled.
//#include "cedr_qlt_inl.hpp"

#endif

//>> cedr_caas.hpp
// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

#ifndef INCLUDE_CEDR_CAAS_HPP
#define INCLUDE_CEDR_CAAS_HPP

//#include "cedr_cdr.hpp"

namespace cedr {
// ClipAndAssuredSum.
namespace caas {

template <typename ExeSpace = Kokkos::DefaultExecutionSpace>
class CAAS : public CDR {
public:
  typedef typename cedr::impl::DeviceType<ExeSpace>::type Device;
  typedef CAAS<ExeSpace> Me;
  typedef std::shared_ptr<Me> Ptr;
  typedef Kokkos::View<Real*, Kokkos::LayoutLeft, Device> RealList;

public:
  struct UserAllReducer {
    typedef std::shared_ptr<const UserAllReducer> Ptr;
    virtual int n_accum_in_place () const { return 1; }
    virtual int operator()(const mpi::Parallel& p,
                           // In Fortran, these are formatted as
                           //   sendbuf(nlocal, nfld)
                           //   rcvbuf(nfld)
                           // The implementation is permitted to modify sendbuf.
                           Real* sendbuf, Real* rcvbuf,
                           // nlocal is number of values to reduce in this rank.
                           // nfld is number of fields.
                           int nlocal, int nfld,
                           MPI_Op op) const = 0;
  };

  CAAS(const mpi::Parallel::Ptr& p, const Int nlclcells,
       const typename UserAllReducer::Ptr& r = nullptr);

  void declare_tracer(int problem_type, const Int& rhomidx) override;

  void end_tracer_declarations() override;

  void get_buffers_sizes(size_t& buf1, size_t& buf2) override;

  void set_buffers(Real* buf1, Real* buf2) override;

  void finish_setup() override;

  int get_problem_type(const Int& tracer_idx) const override;

  Int get_num_tracers() const override;

  // lclcellidx is trivial; it is the user's index for the cell.
  KOKKOS_INLINE_FUNCTION
  void set_rhom(const Int& lclcellidx, const Int& rhomidx, const Real& rhom) const override;

  KOKKOS_INLINE_FUNCTION
  void set_Qm(const Int& lclcellidx, const Int& tracer_idx,
              const Real& Qm, const Real& Qm_min, const Real& Qm_max,
              const Real Qm_prev = std::numeric_limits<Real>::infinity()) const override;

  void run() override;

  KOKKOS_INLINE_FUNCTION
  Real get_Qm(const Int& lclcellidx, const Int& tracer_idx) const override;

protected:
  typedef cedr::impl::Unmanaged<RealList> UnmanagedRealList;
  typedef Kokkos::View<Int*, Kokkos::LayoutLeft, Device> IntList;

  struct Decl {
    int probtype;
    Int rhomidx;
    Decl (const int probtype_, const Int rhomidx_)
      : probtype(probtype_), rhomidx(rhomidx_) {}
  };

  mpi::Parallel::Ptr p_;
  typename UserAllReducer::Ptr user_reducer_;

  Int nlclcells_, nrhomidxs_;
  std::shared_ptr<std::vector<Decl> > tracer_decls_;
  bool need_conserve_;
  IntList probs_, t2r_;
  typename IntList::HostMirror probs_h_;
  RealList d_, send_, recv_;
  bool finished_setup_;

  void reduce_globally();

PRIVATE_CUDA:
  void reduce_locally();
  void finish_locally();

private:
  void get_buffers_sizes(size_t& buf1, size_t& buf2, size_t& buf3);
};

namespace test {
Int unittest(const mpi::Parallel::Ptr& p);
} // namespace test
} // namespace caas
} // namespace cedr

//#include "cedr_caas_inl.hpp"

#endif

//>> cedr_caas_inl.hpp
// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

#ifndef INCLUDE_CEDR_CAAS_INL_HPP
#define INCLUDE_CEDR_CAAS_INL_HPP

//#include "cedr_util.hpp"

namespace cedr {
// ClipAndAssuredSum.
namespace caas {

template <typename ES> KOKKOS_INLINE_FUNCTION
void CAAS<ES>::set_rhom (const Int& lclcellidx, const Int& rhomidx,
                         const Real& rhom) const {
  cedr_kernel_assert(lclcellidx >= 0 && lclcellidx < nlclcells_);
  cedr_kernel_assert(rhomidx >= 0 && rhomidx < nrhomidxs_);
  d_(lclcellidx) = rhom;
}

template <typename ES> KOKKOS_INLINE_FUNCTION
void CAAS<ES>
::set_Qm (const Int& lclcellidx, const Int& tracer_idx,
          const Real& Qm, const Real& Qm_min, const Real& Qm_max,
          const Real Qm_prev) const {
  cedr_kernel_assert(lclcellidx >= 0 && lclcellidx < nlclcells_);
  cedr_kernel_assert(tracer_idx >= 0 && tracer_idx < probs_.extent_int(0));
  const Int nt = probs_.size();
  d_((1 +               tracer_idx)*nlclcells_ + lclcellidx) = Qm;
  d_((1 +   nt + tracer_idx)*nlclcells_ + lclcellidx) = Qm_min;
  d_((1 + 2*nt + tracer_idx)*nlclcells_ + lclcellidx) = Qm_max;
  if (need_conserve_)
    d_((1 + 3*nt + tracer_idx)*nlclcells_ + lclcellidx) = Qm_prev;
}

template <typename ES> KOKKOS_INLINE_FUNCTION
Real CAAS<ES>::get_Qm (const Int& lclcellidx, const Int& tracer_idx) const {
  cedr_kernel_assert(lclcellidx >= 0 && lclcellidx < nlclcells_);
  cedr_kernel_assert(tracer_idx >= 0 && tracer_idx < probs_.extent_int(0));
  return d_((1 + tracer_idx)*nlclcells_ + lclcellidx);
}

} // namespace caas
} // namespace cedr

#endif

//>> cedr_local.hpp
// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

#ifndef INCLUDE_CEDR_LOCAL_HPP
#define INCLUDE_CEDR_LOCAL_HPP

//#include "cedr.hpp"
//#include "cedr_kokkos.hpp"

namespace cedr {
namespace local {

// The following routines solve
//     min_x norm(x - y; w)
//      st   a'x = b
//           xlo <= x <= xhi,
// a > 0, w > 0.

// Minimize the weighted 2-norm. Return 0 on success and x == y, 1 on success
// and x != y, -1 if infeasible, -2 if max_its hit with no solution. See section
// 3 of Bochev, Ridzal, Shashkov, Fast optimization-based conservative remap of
// scalar fields through aggregate mass transfer.
KOKKOS_INLINE_FUNCTION
Int solve_1eq_bc_qp(const Int n, const Real* w, const Real* a, const Real b,
                    const Real* xlo, const Real* xhi,
                    const Real* y, Real* x, const Int max_its = 100);

KOKKOS_INLINE_FUNCTION
Int solve_1eq_bc_qp_2d(const Real* w, const Real* a, const Real b,
                       const Real* xlo, const Real* xhi,
                       const Real* y, Real* x,
                       const bool clip = true,
                       // The 2D algorithm doesn't need to terminate based on a
                       // tolerance, but by default it will exit early if the ND
                       // algorithm would. Set this to false to run the full 2D
                       // algorithm w/o checking the tolerance at start. If this
                       // is false, the feasibility check is also disabled.
                       const bool early_exit_on_tol = true);

// ClipAndAssuredSum. Minimize the 1-norm with w = 1s. Does not check for
// feasibility.
KOKKOS_INLINE_FUNCTION
void caas(const Int n, const Real* a, const Real b,
          const Real* xlo, const Real* xhi,
          const Real* y, Real* x,
          const bool clip = true);

struct Method { enum Enum { least_squares, caas }; };

// Solve
//     min_x norm(x - y; w)
//      st   a'x = b
//           x >= 0,
// a, w > 0. Return 0 on success and x == y, 1 on success and x != y, -1 if
// infeasible. w is used only if lcl_method = least_squares.
KOKKOS_INLINE_FUNCTION
Int solve_1eq_nonneg(const Int n, const Real* a, const Real b, const Real* y, Real* x,
                     const Real* w, const Method::Enum lcl_method);

Int unittest();

} // namespace local
} // namespace cedr

//#include "cedr_local_inl.hpp"

#endif

//>> cedr_mpi_inl.hpp
// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

#ifndef INCLUDE_CEDR_MPI_INL_HPP
#define INCLUDE_CEDR_MPI_INL_HPP

namespace cedr {
namespace mpi {

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
           Request* ireq) {
  MPI_Datatype dt = get_type<T>();
  MPI_Request ureq;
  MPI_Request* req = ireq ? &ireq->request : &ureq;
  int ret = MPI_Isend(const_cast<T*>(buf), count, dt, dest, tag, p.comm(), req);
  if ( ! ireq) MPI_Request_free(req);
#ifdef COMPOSE_DEBUG_MPI
  else ireq->unfreed++;
#endif
  return ret;
}

template <typename T>
int irecv (const Parallel& p, T* buf, int count, int src, int tag, Request* ireq) {
  MPI_Datatype dt = get_type<T>();
  MPI_Request ureq;
  MPI_Request* req = ireq ? &ireq->request : &ureq;
  int ret = MPI_Irecv(buf, count, dt, src, tag, p.comm(), req);
  if ( ! ireq) MPI_Request_free(req);
#ifdef COMPOSE_DEBUG_MPI
  else ireq->unfreed++;
#endif
  return ret;
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

} // namespace mpi
} // namespace cedr

#endif

//>> cedr_local_inl.hpp
// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

#ifndef INCLUDE_CEDR_LOCAL_INL_HPP
#define INCLUDE_CEDR_LOCAL_INL_HPP

//#include "cedr_util.hpp"

namespace cedr {
namespace local {

namespace impl {
KOKKOS_INLINE_FUNCTION
Real calc_r_tol (const Real b, const Real* a, const Real* y, const Int n) {
  Real ab = std::abs(b);
  for (Int i = 0; i < n; ++i) ab = cedr::impl::max(ab, std::abs(a[i]*y[i]));
  return 1e1*std::numeric_limits<Real>::epsilon()*std::abs(ab);
}

// Eval r at end points to check for feasibility, and also possibly a quick exit
// on a common case. Return -1 if infeasible, 1 if a corner is a solution, 0 if
// feasible and a corner is not.
KOKKOS_INLINE_FUNCTION
Int check_lu (const Int n, const Real* a, const Real& b, const Real* xlo,
              const Real* xhi, const Real& r_tol, Real* x) {
  Real r = -b;
  for (Int i = 0; i < n; ++i) {
    x[i] = xlo[i];
    r += a[i]*x[i];
  }
  if (std::abs(r) <= r_tol) return 1;
  if (r > 0) return -1;
  r = -b;
  for (Int i = 0; i < n; ++i) {
    x[i] = xhi[i];
    r += a[i]*x[i];
  }
  if (std::abs(r) <= r_tol) return 1;
  if (r < 0) return -1;
  return 0;
}

KOKKOS_INLINE_FUNCTION
void calc_r (const Int n, const Real* w, const Real* a, const Real b,
             const Real* xlo, const Real* xhi,  const Real* y, const Real& lambda,
             Real* x, Real& r, Real& r_lambda) {
  r = 0;
  r_lambda = 0;
  for (Int i = 0; i < n; ++i) {
    const Real q = a[i]/w[i];
    const Real x_trial = y[i] + lambda*q;
    Real xtmp;
    if (x_trial < (xtmp = xlo[i]))
      x[i] = xtmp;
    else if (x_trial > (xtmp = xhi[i]))
      x[i] = xtmp;
    else {
      x[i] = x_trial;
      r_lambda += a[i]*q;
    }
    r += a[i]*x[i];
  }
  r -= b;
}
} // namespace impl

// 2D special case for efficiency.
KOKKOS_INLINE_FUNCTION
Int solve_1eq_bc_qp_2d (const Real* w, const Real* a, const Real b,
                        const Real* xlo, const Real* xhi, 
                        const Real* y, Real* x,
                        const bool clip, const bool early_exit_on_tol) {
  Int info;
  if (early_exit_on_tol) {
    const Real r_tol = impl::calc_r_tol(b, a, y, 2);
    Int info = impl::check_lu(2, a, b, xlo, xhi, r_tol, x);
    if (info == -1) return info;
  }

  { // Check if the optimal point ignoring bound constraints is in bounds.
    Real qmass = 0, dm = b;
    for (int i = 0; i < 2; ++i) {
      const Real qi = a[i]/w[i];
      qmass += a[i]*qi;
      dm -= a[i]*y[i];
    }
    const Real lambda = dm/qmass;
    bool ok = true;
    for (int i = 0; i < 2; ++i) {
      x[i] = y[i] + lambda*(a[i]/w[i]);
      if (x[i] < xlo[i] || x[i] > xhi[i]) {
        ok = false;
        break;
      }
    }
    if (ok) return 1;
  }

  // Solve for intersection of a'x = b, given by the parameterized line
  //     p(alpa) = x_base + alpha x_dir,
  // with a bounding line.

  // Get parameterized line.
  Real x_base[2];
  for (int i = 0; i < 2; ++i)
    x_base[i] = 0.5*b/a[i];
  Real x_dir[] = {-a[1], a[0]};

  // Get the 4 alpha values.
  Real alphas[4];
  alphas[0] = (xlo[1] - x_base[1])/x_dir[1]; // bottom
  alphas[1] = (xhi[0] - x_base[0])/x_dir[0]; // right
  alphas[2] = (xhi[1] - x_base[1])/x_dir[1]; // top
  alphas[3] = (xlo[0] - x_base[0])/x_dir[0]; // left

  // Find the middle two in the sorted alphas.
  Real min = alphas[0], max = min;
  Int imin = 0, imax = 0;
  for (Int i = 1; i < 4; ++i) {
    const Real alpha = alphas[i];
    if (alpha < min) { min = alpha; imin = i; }
    if (alpha > max) { max = alpha; imax = i; }
  }
  Int ais[2];
  Int cnt = 0;
  for (Int i = 0; i < 4; ++i)
    if (i != imin && i != imax) {
      ais[cnt++] = i;
      if (cnt == 2) break;
    }

  Real objs[2];
  Real alpha_mid = 0;
  for (Int j = 0; j < 2; ++j) {
    const Real alpha = alphas[ais[j]];
    alpha_mid += alpha;
    Real obj = 0;
    for (Int i = 0; i < 2; ++i) {
      x[i] = x_base[i] + alpha*x_dir[i];
      obj += w[i]*cedr::util::square(y[i] - x[i]);
    }
    objs[j] = obj;
  }

  const Int ai = ais[objs[0] <= objs[1] ? 0 : 1];

  info = 1;
  Int i0 = 0;
  switch (ai) {
  case 0: case 2:
    x[1] = ai == 0 ? xlo[1] : xhi[1];
    i0 = 1;
    break;
  case 1: case 3:
    x[0] = ai == 1 ? xhi[0] : xlo[0];
    i0 = 0;
    break;
  default: cedr_kernel_assert(0); info = -2;
  }
  const Int i1 = (i0 + 1) % 2;
  x[i1] = (b - a[i0]*x[i0])/a[i1];
  if (clip)
    x[i1] = cedr::impl::min(xhi[i1], cedr::impl::max(xlo[i1], x[i1]));
  return info;
}

KOKKOS_INLINE_FUNCTION
Int solve_1eq_bc_qp (const Int n, const Real* w, const Real* a, const Real b,
                     const Real* xlo, const Real* xhi, const Real* y, Real* x,
                     const Int max_its) {
  const Real r_tol = impl::calc_r_tol(b, a, y, n);
  Int info = impl::check_lu(n, a, b, xlo, xhi, r_tol, x);
  if (info != 0) return info;

  for (int i = 0; i < n; ++i)
    if (x[i] != y[i]) {
      info = 1;
      x[i] = y[i];
    }

  // In our use case, the caller has already checked (more cheaply) for a quick
  // exit.
#if 0
  { // Check for a quick exit.
    bool all_in = true;
    Real r = 0;
    for (Int i = 0; i < n; ++i) {
      if (x[i] < xlo[i] || x[i] > xhi[i]) {
        all_in = false;
        break;
      }
      r += a[i]*x[i];
    }
    if (all_in) {
      r -= b;
      if (std::abs(r) <= r_tol)
        return info;
    }
  }
#endif

  const Real wall_dist = 1e-3;

  // Get lambda endpoints.
  Real lamlo = 0, lamhi = 0;
  for (Int i = 0; i < n; ++i) {
    const Real rq = w[i]/a[i];
    const Real lamlo_i = rq*(xlo[i] - y[i]);
    const Real lamhi_i = rq*(xhi[i] - y[i]);
    if (i == 0) {
      lamlo = lamlo_i;
      lamhi = lamhi_i;
    } else {
      lamlo = cedr::impl::min(lamlo, lamlo_i);
      lamhi = cedr::impl::max(lamhi, lamhi_i);
    }
  }
  const Real lamlo_feas = lamlo, lamhi_feas = lamhi;
  Real lambda = lamlo <= 0 && lamhi >= 0 ? 0 : lamlo;

  // Bisection-safeguarded Newton iteration for r(lambda) = 0.
  bool prev_step_bisect = false;
  Int nbisect = 0;
  info = -2;
  for (Int iteration = 0; iteration < max_its; ++iteration) {
    // Compute x, r, r_lambda.
    Real r, r_lambda;
    impl::calc_r(n, w, a, b, xlo, xhi, y, lambda, x, r, r_lambda);
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

KOKKOS_INLINE_FUNCTION
void caas (const Int n, const Real* a, const Real b,
           const Real* xlo, const Real* xhi,
           const Real* y, Real* x,
           const bool clip) {
  Real dm = b;
  for (Int i = 0; i < n; ++i) {
    x[i] = cedr::impl::max(xlo[i], cedr::impl::min(xhi[i], y[i]));
    dm -= a[i]*x[i];
  }
  if (dm == 0) return;
  if (dm > 0) {
    Real fac = 0;
    for (Int i = 0; i < n; ++i)
      fac += a[i]*(xhi[i] - x[i]);
    if (fac > 0) {
      fac = dm/fac;
      for (Int i = 0; i < n; ++i)
        x[i] += fac*(xhi[i] - x[i]);
    }
  } else if (dm < 0) {
    Real fac = 0;
    for (Int i = 0; i < n; ++i)
      fac += a[i]*(x[i] - xlo[i]);
    if (fac > 0) {
      fac = dm/fac;
      for (Int i = 0; i < n; ++i)
        x[i] += fac*(x[i] - xlo[i]);
    }
  }
  // Clip again for numerics.
  if (clip)
    for (Int i = 0; i < n; ++i)
      x[i] = cedr::impl::max(xlo[i], cedr::impl::min(xhi[i], x[i]));
}

KOKKOS_INLINE_FUNCTION
Int solve_1eq_nonneg (const Int n, const Real* a, const Real b, const Real* y, Real* x,
                      const Real* w,  const Method::Enum method) {
  cedr_kernel_assert(n <= 16);
  if (b < 0) return -1;

  const Real zero[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  // Set the upper bound to the value that implies that just one slot gets all
  // of the mass.
  Real xhi[16];
  for (int i = 0; i < n; ++i)
    xhi[i] = b/a[i];

  if (method == Method::caas) {
    caas(n, a, b, zero, xhi, y, x);
    return 1;
  } else {
    if (n == 2)
      return solve_1eq_bc_qp_2d(w, a, b, zero, xhi, y, x);
    else
      return solve_1eq_bc_qp(n, w, a, b, zero, xhi, y, x);
  }
}

} // namespace local
} // namespace cedr

#endif

//>> cedr_qlt_inl.hpp
// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

#ifndef INCLUDE_CEDR_QLT_INL_HPP
#define INCLUDE_CEDR_QLT_INL_HPP

#include <cassert>

//#include "cedr_local.hpp"

namespace cedr {
namespace qlt {

template <typename ES> KOKKOS_INLINE_FUNCTION
void QLT<ES>::set_rhom (const Int& lclcellidx, const Int& rhomidx,
                        const Real& rhom) const {
  const Int ndps = md_.a_d.prob2bl2r[md_.nprobtypes];
  bd_.l2r_data(ndps*lclcellidx) = rhom;  
}

template <typename ES> KOKKOS_INLINE_FUNCTION
void QLT<ES>::set_Qm (const Int& lclcellidx, const Int& tracer_idx,
                      const Real& Qm,
                      const Real& Qm_min, const Real& Qm_max,
                      const Real Qm_prev) const {
  const Int ndps = md_.a_d.prob2bl2r[md_.nprobtypes];
  Real* bd; {
    const Int bdi = md_.a_d.trcr2bl2r(tracer_idx);
    bd = &bd_.l2r_data(ndps*lclcellidx + bdi);
  }
  {
    const Int problem_type = md_.a_d.trcr2prob(tracer_idx);
    Int next = 0;
    if (problem_type & ProblemType::shapepreserve) {
      bd[0] = Qm_min;
      bd[1] = Qm;
      bd[2] = Qm_max;
      next = 3;
    } else if (problem_type & ProblemType::consistent) {
      const Real rhom = bd_.l2r_data(ndps*lclcellidx);
      bd[0] = Qm_min / rhom;
      bd[1] = Qm;
      bd[2] = Qm_max / rhom;
      next = 3;
    } else if (problem_type & ProblemType::nonnegative) {
      bd[0] = Qm;
      next = 1;
    } else {
      cedr_kernel_throw_if(true, "set_Q: invalid problem_type.");
    }
    if (problem_type & ProblemType::conserve) {
      cedr_kernel_throw_if(Qm_prev == std::numeric_limits<Real>::infinity(),
                           "Qm_prev was not provided to set_Q.");
      bd[next] = Qm_prev;
    }
  }
}

template <typename ES> KOKKOS_INLINE_FUNCTION
Real QLT<ES>::get_Qm (const Int& lclcellidx, const Int& tracer_idx) const {
  const Int ndps = md_.a_d.prob2br2l[md_.nprobtypes];
  const Int bdi = md_.a_d.trcr2br2l(tracer_idx);
  return bd_.r2l_data(ndps*lclcellidx + bdi);
}

//todo Replace this and the calling code with ReconstructSafely.
KOKKOS_INLINE_FUNCTION
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

template <typename ES> KOKKOS_INLINE_FUNCTION
int QLT<ES>::MetaData::get_problem_type (const int& idx) {
  static const Int problem_type[] = {
    CPT::st, CPT::cst, CPT::t, CPT::ct, CPT::nn, CPT::cnn
  };
  return problem_type[idx];
}
    
template <typename ES> KOKKOS_INLINE_FUNCTION
int QLT<ES>::MetaData::get_problem_type_l2r_bulk_size (const int& mask) {
  if (mask & ProblemType::nonnegative) {
    if (mask & ProblemType::conserve) return 2;
    return 1;
  }
  if (mask & ProblemType::conserve) return 4;
  return 3;
}

namespace impl {
KOKKOS_INLINE_FUNCTION
void solve_node_problem (const Real& rhom, const Real* pd, const Real& Qm,
                         const Real& rhom0, const Real* k0d, Real& Qm0,
                         const Real& rhom1, const Real* k1d, Real& Qm1,
                         const bool prefer_mass_con_to_bounds) {
  Real Qm_min_kids [] = {k0d[0], k1d[0]};
  Real Qm_orig_kids[] = {k0d[1], k1d[1]};
  Real Qm_max_kids [] = {k0d[2], k1d[2]};
  { // The ideal problem is not assuredly feasible. Test for feasibility. If not
    // feasible, adjust bounds to solve the safety problem, which is assuredly
    // feasible if the total density field rho is mass conserving (Q doesn't
    // have to be mass conserving, of course; achieving mass conservation is one
    // use for QLT).
    const Real Qm_min = pd[0], Qm_max = pd[2];
    const bool lo = Qm < Qm_min, hi = Qm > Qm_max;
    if (lo || hi) {
      // If the discrepancy is numerical noise, don't act on it.
      const Real tol = 10*std::numeric_limits<Real>::epsilon();
      const Real discrepancy = lo ? Qm_min - Qm : Qm - Qm_max;
      if (discrepancy > tol*(Qm_max - Qm_min)) {
        const Real rhom_kids[] = {rhom0, rhom1};
        r2l_nl_adjust_bounds(lo ? Qm_min_kids : Qm_max_kids,
                             rhom_kids,
                             Qm - (lo ? Qm_min : Qm_max));
      }
    } else {
      // Quick exit if everything is OK as is. This is a speedup, and it also
      // lets the subnode solver make ~1 ulp changes instead of having to keep x
      // = y if y satisfies the conditions. Without this block, the
      // no_change_should_hold tests can fail.
      if (Qm == pd[1] && // Was our total tracer mass adjusted?
          // Are the kids' problems feasible?
          Qm_orig_kids[0] >= Qm_min_kids[0] && Qm_orig_kids[0] <= Qm_max_kids[0] &&
          Qm_orig_kids[1] >= Qm_min_kids[1] && Qm_orig_kids[1] <= Qm_max_kids[1]) {
        // Don't need to do anything, so skip even the math-based quick exits in
        // solve_node_problem.
        Qm0 = Qm_orig_kids[0];
        Qm1 = Qm_orig_kids[1];
        return;
      }
    }
  }
  { // Solve the node's QP.
    static const Real ones[] = {1, 1};
    const Real w[] = {1/rhom0, 1/rhom1};
    Real Qm_kids[] = {k0d[1], k1d[1]};
    local::solve_1eq_bc_qp_2d(w, ones, Qm, Qm_min_kids, Qm_max_kids,
                              Qm_orig_kids, Qm_kids,
                              ! prefer_mass_con_to_bounds /* clip */,
                              ! prefer_mass_con_to_bounds /* early_exit_on_tol */);
    Qm0 = Qm_kids[0];
    Qm1 = Qm_kids[1];
  }
}

KOKKOS_INLINE_FUNCTION
void solve_node_problem (const Int problem_type,
                         const Real& rhom, const Real* pd, const Real& Qm,
                         const Real& rhom0, const Real* k0d, Real& Qm0,
                         const Real& rhom1, const Real* k1d, Real& Qm1,
                         const bool prefer_mass_con_to_bounds) {
  if ((problem_type & ProblemType::consistent) &&
      ! (problem_type & ProblemType::shapepreserve)) {      
    Real mpd[3], mk0d[3], mk1d[3];
    mpd[0]  = pd [0]*rhom ; mpd [1] = pd[1] ; mpd [2] = pd [2]*rhom ;
    mk0d[0] = k0d[0]*rhom0; mk0d[1] = k0d[1]; mk0d[2] = k0d[2]*rhom0;
    mk1d[0] = k1d[0]*rhom1; mk1d[1] = k1d[1]; mk1d[2] = k1d[2]*rhom1;
    solve_node_problem(rhom, mpd, Qm, rhom0, mk0d, Qm0, rhom1, mk1d, Qm1,
                       prefer_mass_con_to_bounds);
    return;
  } else if (problem_type & ProblemType::nonnegative) {
    static const Real ones[] = {1, 1};
    const Real w[] = {1/rhom0, 1/rhom1};
    Real Qm_orig_kids[] = {k0d[0], k1d[0]};
    Real Qm_kids[2] = {k0d[0], k1d[0]};
    local::solve_1eq_nonneg(2, ones, Qm, Qm_orig_kids, Qm_kids, w,
                            local::Method::least_squares);
    Qm0 = Qm_kids[0];
    Qm1 = Qm_kids[1];
  } else {
    solve_node_problem(rhom, pd, Qm, rhom0, k0d, Qm0, rhom1, k1d, Qm1,
                       prefer_mass_con_to_bounds);
  }
}

} // namespace impl
} // namespace qlt
} // namespace cedr

#endif

//>> cedr_test_randomized.hpp
// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

#ifndef INCLUDE_CEDR_TEST_RANDOMIZED_HPP
#define INCLUDE_CEDR_TEST_RANDOMIZED_HPP

//#include "cedr_cdr.hpp"
//#include "cedr_mpi.hpp"
//#include "cedr_util.hpp"

namespace cedr {
namespace test {

class TestRandomized {
public:
  TestRandomized(const std::string& cdr_name, const mpi::Parallel::Ptr& p,
                 const Int& ncells, const bool verbose = false,
                 const CDR::Options options = CDR::Options());

  // The subclass should call this, probably in its constructor.
  void init();

  template <typename CDRT, typename ExeSpace = Kokkos::DefaultExecutionSpace>
  Int run(const Int nrepeat = 1, const bool write=false);

private:
  const std::string cdr_name_;
  const CDR::Options options_;

protected:
  struct Tracer {
    typedef ProblemType PT;
    
    Int idx;
    Int problem_type;
    Int perturbation_type;
    bool no_change_should_hold, safe_should_hold, local_should_hold;
    bool write;

    std::string str() const;

    Tracer ()
      : idx(-1), problem_type(-1), perturbation_type(-1), no_change_should_hold(false),
        safe_should_hold(true), local_should_hold(true), write(false)
    {}
  };

  struct ValuesPartition {
    Int ncells () const { return ncells_; }
    KIF Real* rhom () const { return v_; }
    KIF Real* Qm_min  (const Int& ti) const { return v_ + ncells_*(1 + 4*ti    ); }
    KIF Real* Qm      (const Int& ti) const { return v_ + ncells_*(1 + 4*ti + 1); }
    KIF Real* Qm_max  (const Int& ti) const { return v_ + ncells_*(1 + 4*ti + 2); }
    KIF Real* Qm_prev (const Int& ti) const { return v_ + ncells_*(1 + 4*ti + 3); }
  protected:
    void init (const Int ncells, Real* v) {
      ncells_ = ncells;
      v_ = v;
    }
  private:
    Int ncells_;
    Real* v_;
  };

  struct Values : public ValuesPartition {
    Values (const Int ntracers, const Int ncells)
      : v_((4*ntracers + 1)*ncells)
    { init(ncells, v_.data()); }
    Real* data () { return v_.data(); }
    size_t size () const { return v_.size(); }
  private:
    std::vector<Real> v_;
  };

PRIVATE_CUDA:
  template <typename ExeSpace>
  struct ValuesDevice : public ValuesPartition {
    // This Values object is the source of data and gets updated by sync_host.
    ValuesDevice (Values& v)
      : rar_(v.data(), v.size())
    { init(v.ncells(), rar_.device_ptr()); }
    // Values -> device.
    void sync_device () { rar_.sync_device(); }
    // Update Values from device.
    void sync_host () { rar_.sync_host(); }
  private:
    util::RawArrayRaft<Real, ExeSpace> rar_;
  };

protected:
  // For solution output, if requested.
  struct Writer {
    std::unique_ptr<FILE, cedr::util::FILECloser> fh;
    std::vector<Int> ngcis;  // Number of i'th rank's gcis_ array.
    std::vector<Long> gcis;  // Global cell indices packed by rank's gcis_ vector.
    std::vector<int> displs; // Cumsum of above.
    ~Writer();
  };

  const mpi::Parallel::Ptr p_;
  const Int ncells_;
  // Global mesh entity IDs, 1-1 with reduction array index or QLT leaf node.
  std::vector<Long> gcis_;
  std::vector<Tracer> tracers_;

  // Tell this class the CDR.
  virtual CDR& get_cdr() = 0;

  // Fill gcis_.
  virtual void init_numbering() = 0;

  // Using tracers_, the vector of Tracers, initialize the CDR's tracers.
  virtual void init_tracers() = 0;

  virtual void run_impl(const Int trial) = 0;

private:
  // For optional output.
  bool write_inited_;
  std::shared_ptr<Writer> w_; // Only on root.

  void init_tracers_vector();

  void add_const_to_Q(
    const Tracer& t, Values& v,
    // Move 0 < alpha <= 1 of the way to the QLT or safety feasibility bound.
    const Real& alpha,
    // Whether the modification should be done in a mass-conserving way.
    const bool conserve_mass,
    // Only safety problem is feasible.
    const bool safety_problem);

  void perturb_Q(const Tracer& t, Values& v);
  void init_writer();
  void gather_field(const Real* Qm_lcl, std::vector<Real>& Qm_gbl,
                    std::vector<Real>& wrk);
  void write_field(const std::string& tracer_name, const std::string& field_name,
                   const std::vector<Real>& Qm);
  void write_pre(const Tracer& t, Values& v);
  void write_post(const Tracer& t, Values& v);

  static void generate_rho(Values& v);
  static void generate_Q(const Tracer& t, Values& v);
  static void permute_Q(const Tracer& t, Values& v);
  static std::string get_tracer_name(const Tracer& t);
  Int check(const std::string& cdr_name, const mpi::Parallel& p,
            const std::vector<Tracer>& ts, const Values& v);
};

} // namespace test
} // namespace cedr

//#include "cedr_test_randomized_inl.hpp"

#endif

//>> cedr_test_randomized_inl.hpp
// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

#ifndef INCLUDE_CEDR_TEST_RANDOMIZED_INL_HPP
#define INCLUDE_CEDR_TEST_RANDOMIZED_INL_HPP

//#include "cedr_test_randomized.hpp"

namespace cedr {
namespace test {

template <typename CDRT, typename ES>
Int TestRandomized::run (const Int nrepeat, const bool write) {
  const Int nt = tracers_.size(), nlclcells = gcis_.size();

  Values v(nt, nlclcells);
  generate_rho(v);
  for (const auto& t : tracers_) {
    generate_Q(t, v);
    perturb_Q(t, v);
  }

  if (write)
    for (const auto& t : tracers_)
      write_pre(t, v);

  CDRT cdr = static_cast<CDRT&>(get_cdr());
  ValuesDevice<ES> vd(v);
  vd.sync_device();

  {
    const auto rhom = vd.rhom();
    const auto set_rhom = KOKKOS_LAMBDA (const Int& i) {
      cdr.set_rhom(i, 0, rhom[i]);
    };
    Kokkos::parallel_for(Kokkos::RangePolicy<ES>(0, nlclcells), set_rhom);
  }
  // repeat > 1 runs the same values repeatedly for performance
  // meaurement.
  for (Int trial = 0; trial <= nrepeat; ++trial) {
    const auto set_Qm = KOKKOS_LAMBDA (const Int& j) {
      const auto ti = j / nlclcells;
      const auto i = j % nlclcells;
      cdr.set_Qm(i, ti, vd.Qm(ti)[i], vd.Qm_min(ti)[i], vd.Qm_max(ti)[i],
                 vd.Qm_prev(ti)[i]);
    };
    Kokkos::parallel_for(Kokkos::RangePolicy<ES>(0, nt*nlclcells), set_Qm);
    run_impl(trial);
  }
  {
    const auto get_Qm = KOKKOS_LAMBDA (const Int& j) {
      const auto ti = j / nlclcells;
      const auto i = j % nlclcells;
      vd.Qm(ti)[i] = cdr.get_Qm(i, ti);
    };
    Kokkos::parallel_for(Kokkos::RangePolicy<ES>(0, nt*nlclcells), get_Qm);
  }
  vd.sync_host(); // => v contains computed values

  if (write)
    for (const auto& t : tracers_)
      write_post(t, v);
  return check(cdr_name_, *p_, tracers_, v);
}

} // namespace test
} // namespace cedr

#endif

//>> cedr_test.hpp
// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

#ifndef INCLUDE_CEDR_TEST_HPP
#define INCLUDE_CEDR_TEST_HPP

//#include "cedr.hpp"
//#include "cedr_mpi.hpp"

namespace cedr {
namespace test {
namespace transport1d {

struct Input {
  Int ncells;
  bool verbose;
};

Int run(const mpi::Parallel::Ptr& p, const Input& in);

} // namespace transport1d
} // namespace test
} // namespace cedr

#endif

//>> cedr_util.cpp
// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

//#include "cedr_util.hpp"

namespace cedr {
namespace util {

bool eq (const std::string& a, const char* const b1, const char* const b2) {
  return (a == std::string(b1) || (b2 && a == std::string(b2)) ||
          a == std::string("-") + std::string(b1));
}

Real urand () { return std::rand() / ((Real) RAND_MAX + 1.0); }

Real reldif (const Real* a, const Real* b, const Int n) {
  Real num = 0, den = 0;
  for (Int i = 0; i < n; ++i) {
    num += std::abs(a[i] - b[i]);
    den += std::abs(a[i]);
  }
  return num/den;
}

}
}

//>> cedr_local.cpp
// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

//#include "cedr_local.hpp"
//#include "cedr_local_inl.hpp"

namespace cedr {
namespace local {
namespace test {
// Check the first-order optimality conditions. Return true if OK, false
// otherwise. If quiet, don't print anything.
bool check_1eq_bc_qp_foc (
  const char* label, const Int n, const Real* w, const Real* a, const Real b,
  const Real* xlo, const Real* xhi, const Real* y, const Real* x, const bool verbose)
{
  auto& os = std::cout;
  bool ok = true;
  Real xtmp;
  // Check the bound constraints.
  for (Int i = 0; i < n; ++i)
    if (x[i] < (xtmp = xlo[i])) {
      if (verbose)
        os << "x[" << i << "] = " << x[i]
           << " but x[i] - xlo[i] = " << (x[i] - xtmp) << "\n";
      ok = false;
    }
  for (Int i = 0; i < n; ++i)
    if (x[i] > (xtmp = xhi[i])) {
      if (verbose)
        os << "x[" << i << "] = " << x[i]
           << " but xhi[i] - x[i] = " << (xtmp - x[i]) << "\n";
      ok = false;
    }
  // Check the equality constraint.
  Real r = 0;
  for (Int i = 0; i < n; ++i)
    r += a[i]*x[i];
  r -= b;
  if (std::abs(r) > impl::calc_r_tol(b, a, y, n)) {
    if (verbose)
      os << "r = " << r << "\n";
    ok = false;
  }
  // Check the gradient is 0 when projected into the constraints. Compute
  //     g = W (x - y)
  //     g_reduced = g - C ((C'C) \ (C'g))
  // where
  //     IA = I(:,A)
  //     C = [IA a],
  // and A is the active set.
  const Real padtol = 1e5*std::numeric_limits<Real>::epsilon();
  Real lambda = 0, den = 0;
  for (Int i = 0; i < n; ++i) {
    const Real pad = padtol*(xhi[i] - xlo[i]);
    if (xlo[i] + pad <= x[i] && x[i] <= xhi[i] - pad) {
      const Real gi = w[i]*(x[i] - y[i]);
      lambda += a[i]*gi;
      den += a[i]*a[i];
    }
  }
  lambda /= den;
  Real normg = 0, normy = 0;
  for (Int i = 0; i < n; ++i) {
    normy += cedr::util::square(y[i]);
    const Real pad = padtol*(xhi[i] - xlo[i]);
    if (xlo[i] + pad <= x[i] && x[i] <= xhi[i] - pad)
      normg += cedr::util::square(w[i]*(x[i] - y[i]) - a[i]*lambda);
  }
  normy = std::sqrt(normy);
  normg = std::sqrt(normg);
  const Real gtol = 1e4*std::numeric_limits<Real>::epsilon()*normy;
  if (normg > gtol) {
    if (verbose)
      os << "norm(g) = " << normg << " gtol = " << gtol << "\n";
    ok = false;
  }
  // Check the gradient at the active boundaries.
  for (Int i = 0; i < n; ++i) {
    const bool onlo = x[i] == xlo[i];
    const bool onhi = onlo ? false : x[i] == xhi[i];
    if (onlo || onhi) {
      const Real rg = w[i]*(x[i] - y[i]) - a[i]*lambda;
      if (onlo && rg < -gtol) {
        if (verbose)
          os << "onlo but rg = " << rg << "\n";
        ok = false;
      } else if (onhi && rg > gtol) {
        if (verbose)
          os << "onhi but rg = " << rg << "\n";
        ok = false;
      }
    }
  }
  if ( ! ok && verbose)
    os << "label: " << label << "\n";
  return ok;
}

Int test_1eq_bc_qp () {
  bool verbose = true;
  Int nerr = 0;

  Int n;
  static const Int N = 16;
  Real w[N], a[N], b, xlo[N], xhi[N], y[N], x[N], al, au;

  auto run = [&] () {
    const Int info = solve_1eq_bc_qp(n, w, a, b, xlo, xhi, y, x);
    const bool ok = test::check_1eq_bc_qp_foc(
      "unittest", n, w, a, b, xlo, xhi, y, x, verbose);
    if ( ! ok) ++nerr;

    if (n == 2) {
      // This version never returns 0.
      Real x2[2];
      const Int info2 = solve_1eq_bc_qp_2d(w, a, b, xlo, xhi, y, x2);
      if (info2 != 1 && (info == 0 || info == 1)) {
        if (verbose) pr(puf(info) pu(info2));
        ++nerr;
      }
      const Real rd = cedr::util::reldif(x, x2, 2);
      if (rd > 1e4*std::numeric_limits<Real>::epsilon()) {
        if (verbose)
          printf("%1.1e | y %1.15e %1.15e | x %1.15e %1.15e | "
                 "x2 %1.15e %1.15e | l %1.15e %1.15e | u %1.15e %1.15e\n",
                 rd, y[0], y[1], x[0], x[1], x2[0], x2[1],
                 xlo[0], xlo[1], xhi[0], xhi[1]);
        ++nerr;
      }
    }

    caas(n, a, b, xlo, xhi, y, x);
    Real m = 0, den = 0;
    for (Int i = 0; i < n; ++i) {
      m += a[i]*x[i];
      den += std::abs(a[i]*x[i]);
      if (x[i] < xlo[i]) ++nerr;
      else if (x[i] > xhi[i]) ++nerr;
    }
    const Real rd = std::abs(b - m)/den;
    if (rd > 1e3*std::numeric_limits<Real>::epsilon()) {
      if (verbose) pr(puf(rd) pu(n) pu(b) pu(m));
      ++nerr;
    }
  };

  auto gena = [&] () {
    for (Int i = 0; i < n; ++i)
      a[i] = 0.1 + cedr::util::urand();
  };
  auto genw = [&] () {
    for (Int i = 0; i < n; ++i)
      w[i] = 0.1 + cedr::util::urand();
  };
  auto genbnds = [&] () {
    al = au = 0;
    for (Int i = 0; i < n; ++i) {
      xlo[i] = cedr::util::urand() - 0.5;
      al += a[i]*xlo[i];
      xhi[i] = xlo[i] + cedr::util::urand();
      au += a[i]*xhi[i];
    }
  };
  auto genb = [&] (const bool in) {
    if (in) {
      const Real alpha = cedr::util::urand();
      b = alpha*al + (1 - alpha)*au;
    } else {
      if (cedr::util::urand() > 0.5)
        b = au + 0.01 + cedr::util::urand();
      else
        b = al - 0.01 - cedr::util::urand();
    }
  };
  auto geny = [&] (const bool in) {
    if (in) {
      for (Int i = 0; i < n; ++i) {
        const Real alpha = cedr::util::urand();
        y[i] = alpha*xlo[i] + (1 - alpha)*xhi[i];
      }
    } else if (cedr::util::urand() > 0.2) {
      for (Int i = 1; i < n; i += 2) {
        const Real alpha = cedr::util::urand();
        y[i] = alpha*xlo[i] + (1 - alpha)*xhi[i];
        cedr_assert(y[i] >= xlo[i] && y[i] <= xhi[i]);
      }      
      for (Int i = 0; i < n; i += 4)
        y[i] = xlo[i] - cedr::util::urand();
      for (Int i = 2; i < n; i += 4)
        y[i] = xhi[i] + cedr::util::urand();
    } else {
      for (Int i = 0; i < n; i += 2)
        y[i] = xlo[i] - cedr::util::urand();
      for (Int i = 1; i < n; i += 2)
        y[i] = xhi[i] + cedr::util::urand();
    }
  };
  auto b4y = [&] () {
    b = 0;
    for (Int i = 0; i < n; ++i)
      b += a[i]*y[i];
  };

  for (n = 2; n <= N; ++n) {
    const Int count = n == 2 ? 100 : 10;
    for (Int i = 0; i < count; ++i) {
      gena();
      genw();
      genbnds();
      genb(true);
      geny(true);
      run();
      b4y();
      run();
      genb(true);
      geny(false);
      run();
    }
  }

  return  nerr;
}

Int test_1eq_nonneg () {
  using cedr::util::urand;
  using cedr::util::reldif;

  bool verbose = true;
  Int nerr = 0;

  Int n;
  static const Int N = 16;
  Real w[N], a[N], b, xlo[N], xhi[N], y[N], x_ls[N], x_caas[N], x1_ls[N], x1_caas[N];

  for (n = 2; n <= 2; ++n) {
    const Int count = 20;
    for (Int trial = 0; trial < count; ++trial) {
      b = 0.5*n*urand();
      for (Int i = 0; i < n; ++i) {
        w[i] = 0.1 + urand();
        a[i] = 0.1 + urand();
        xlo[i] = 0;
        xhi[i] = b/a[i];
        y[i] = urand();
        if (urand() > 0.8) y[i] *= -1;
        x1_caas[i] = urand() > 0.5 ? y[i] : -1;
        x1_ls[i] = urand() > 0.5 ? y[i] : -1;
      }
      solve_1eq_nonneg(n, a, b, y, x1_caas, w, Method::caas);
      caas(n, a, b, xlo, xhi, y, x_caas);
      solve_1eq_nonneg(n, a, b, y, x1_ls, w, Method::least_squares);
      solve_1eq_bc_qp(n, w, a, b, xlo, xhi, y, x_ls);
      const Real rd_caas = reldif(x_caas, x1_caas, 2);
      const Real rd_ls = reldif(x_ls, x1_ls, 2);
      if (rd_ls > 5e1*std::numeric_limits<Real>::epsilon() ||
          rd_caas > 1e1*std::numeric_limits<Real>::epsilon()) {
        pr(puf(rd_ls) pu(rd_caas));
        if (verbose) {
          using cedr::util::prarr;
          prarr("w", w, n);
          prarr("a", a, n);
          prarr("xhi", xhi, n);
          prarr("y", y, n);
          prarr("x_ls", x_ls, n);
          prarr("x1_ls", x1_ls, n);
          prarr("x_caas", x_caas, n);
          prarr("x1_caas", x1_caas, n);
          prc(b);
          Real mass = 0;
          for (Int i = 0; i < n; ++i) mass += a[i]*x_ls[i];
          prc(mass);
          mass = 0;
          for (Int i = 0; i < n; ++i) mass += a[i]*x_ls[i];
          prc(mass);
          mass = 0;
          for (Int i = 0; i < n; ++i) mass += a[i]*x_caas[i];
          prc(mass);
        }
        ++nerr;
      }
    }
  }

  return nerr;
}

} // namespace test

Int unittest () {
  Int nerr = 0;
  nerr += test::test_1eq_bc_qp();
  nerr += test::test_1eq_nonneg();
  return nerr;
}

} // namespace local
} // namespace cedr

//>> cedr_mpi.cpp
// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

//#include "cedr_mpi.hpp"
//#include "cedr_util.hpp"

namespace cedr {
namespace mpi {

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

template <> MPI_Datatype get_type<int>() { return MPI_INT; }
template <> MPI_Datatype get_type<double>() { return MPI_DOUBLE; }
template <> MPI_Datatype get_type<long>() { return MPI_LONG_INT; }

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

bool all_ok (const Parallel& p, bool im_ok) {
  int ok = im_ok, msg;
  all_reduce<int>(p, &ok, &msg, 1, MPI_LAND);
  return static_cast<bool>(msg);
}

} // namespace mpi
} // namespace cedr

//>> cedr_qlt.cpp
// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

//#include "cedr_qlt.hpp"
//#include "cedr_test_randomized.hpp"

#include <sys/time.h>

#include <cassert>
#include <cmath>

#include <set>
#include <limits>
#include <algorithm>

namespace cedr {
namespace qlt {

class Timer {
public:
  enum Op { tree, analyze, qltrun, qltrunl2r, qltrunr2l, snp, waitall,
            total, NTIMERS };
  static inline void init () {
#ifdef COMPOSE_QLT_TIME
    for (int i = 0; i < NTIMERS; ++i) {
      et_[i] = 0;
      cnt_[i] = 0;
    }
#endif
  }
  static inline void reset (const Op op) {
#ifdef COMPOSE_QLT_TIME
    et_[op] = 0;
    cnt_[op] = 0;
#endif
  }
  static inline void start (const Op op) {
#ifdef COMPOSE_QLT_TIME
    gettimeofday(&t_start_[op], 0);
    ++cnt_[op];
#endif
  }
  static inline void stop (const Op op) {
#ifdef COMPOSE_QLT_TIME
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
#ifdef COMPOSE_QLT_TIME
    const double tot = et_[total];
    tpr(tree); tpr(analyze);
    tpr(qltrun); tpr(qltrunl2r); tpr(qltrunr2l); tpr(snp); tpr(waitall);
    printf("%-20s %10.3e %10.1f\n", "total", tot, 100.0);
#endif
  }
#undef tpr
private:
#ifdef COMPOSE_QLT_TIME
  static timeval t_start_[NTIMERS];
  static double et_[NTIMERS];
  static int cnt_[NTIMERS];
#endif
};
#ifdef COMPOSE_QLT_TIME
timeval Timer::t_start_[Timer::NTIMERS];
double Timer::et_[Timer::NTIMERS];
int Timer::cnt_[Timer::NTIMERS];
#endif

namespace impl {
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

// Unit tests for NodeSets.
Int unittest (const Parallel::Ptr& p, const NodeSets::ConstPtr& ns,
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
} // namespace impl

template <typename ES>
void QLT<ES>::init (const std::string& name, IntList& d,
                    typename IntList::HostMirror& h, size_t n) {
  d = IntList("QLT " + name, n);
  h = Kokkos::create_mirror_view(d);
}

template <typename ES>
int QLT<ES>::MetaData::get_problem_type_idx (const int& mask) {
  switch (mask) {
  case CPT::s:   case CPT::st:  return 0;
  case CPT::cs:  case CPT::cst: return 1;
  case CPT::t:   return 2;
  case CPT::ct:  return 3;
  case CPT::nn:  return 4;
  case CPT::cnn: return 5;
  default: cedr_kernel_throw_if(true, "Invalid problem type."); return -1;
  }
}

template <typename ES>
int QLT<ES>::MetaData::get_problem_type_r2l_bulk_size (const int& mask) {
  if (mask & ProblemType::shapepreserve) return 1;
  return 3;
}

template <typename ES>
void QLT<ES>::MetaData::init (const MetaDataBuilder& mdb) {
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
      if (problem_type != get_problem_type(pi)) continue;
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
  cedr_assert(a_d.prob2trcrptr[nprobtypes] == ntracers);
}

template <typename ES>
void QLT<ES>::BulkData::init (const size_t& l2r_sz, const size_t& r2l_sz) {
  l2r_data_ = RealList("QLT l2r_data", l2r_sz);
  r2l_data_ = RealList("QLT r2l_data", r2l_sz);
  l2r_data = l2r_data_;
  r2l_data = r2l_data_;
  inited_ = true;
}

template <typename ES>
void QLT<ES>::BulkData::init (Real* buf1, const size_t& l2r_sz,
                              Real* buf2, const size_t& r2l_sz) {
  l2r_data_ = RealList(buf1, l2r_sz);
  r2l_data_ = RealList(buf2, r2l_sz);
  l2r_data = l2r_data_;
  r2l_data = r2l_data_;
  inited_ = true;
}

template <typename ES>
void init_device_data (const impl::NodeSets& ns, impl::NodeSetsHostData& h,
                       impl::NodeSetsDeviceData<ES>& d) {
  typedef impl::NodeSetsDeviceData<ES> NSDD;
  d.node = typename NSDD::NodeList("NSDD::node", ns.nnode());
  h.node = Kokkos::create_mirror_view(d.node);
  d.lvlptr = typename NSDD::IntList("NSDD::lvlptr", ns.levels.size() + 1);
  h.lvlptr = Kokkos::create_mirror_view(d.lvlptr);
  Int nnode = 0;
  for (const auto& lvl : ns.levels)
    nnode += lvl.nodes.size();
  d.lvl = typename NSDD::IntList("NSDD::lvl", nnode);
  h.lvl = Kokkos::create_mirror_view(d.lvl);
  h.lvlptr(0) = 0;
  for (size_t il = 0; il < ns.levels.size(); ++il) {
    const auto& level = ns.levels[il];
    h.lvlptr(il+1) = h.lvlptr(il) + level.nodes.size();
    for (Int os = h.lvlptr(il), i = 0; i < h.lvlptr(il+1) - os; ++i)
      h.lvl(os+i) = level.nodes[i];
    nnode += level.nodes.size();
  }
  for (Int i = 0; i < ns.nnode(); ++i)
    h.node(i) = *ns.node_h(i);
  Kokkos::deep_copy(d.node, h.node);
  Kokkos::deep_copy(d.lvl, h.lvl);
  Kokkos::deep_copy(d.lvlptr, h.lvlptr);
}

template <typename ES>
void QLT<ES>::init (const Parallel::Ptr& p, const Int& ncells,
                    const tree::Node::Ptr& tree) {
  p_ = p;
  Timer::start(Timer::analyze);
  ns_ = impl::analyze(p, ncells, tree);
  nshd_ = std::make_shared<impl::NodeSetsHostData>();
  nsdd_ = std::make_shared<impl::NodeSetsDeviceData<ES> >();
  init_device_data(*ns_, *nshd_, *nsdd_);
  init_ordinals();
  Timer::stop(Timer::analyze);
  mdb_ = std::make_shared<MetaDataBuilder>();
}

template <typename ES>
void QLT<ES>::init_ordinals () {
  gci2lci_ = std::make_shared<Gci2LciMap>();
  for (const auto& idx : ns_->levels[0].nodes) {
    const auto n = ns_->node_h(idx);
    (*gci2lci_)[n->id] = n->offset;
  }
}

template <typename ES>
QLT<ES>::QLT (const Parallel::Ptr& p, const Int& ncells, const tree::Node::Ptr& tree,
              Options options)
  : CDR(options)
{
  init(p, ncells, tree);
  cedr_throw_if(nlclcells() == 0, "QLT does not support 0 cells on a rank.");
}

template <typename ES>
void QLT<ES>::print (std::ostream& os) const {
  ns_->print(os);
}

// Number of cells owned by this rank.
template <typename ES>
Int QLT<ES>::nlclcells () const { return ns_->levels[0].nodes.size(); }

// Cells owned by this rank, in order of local numbering. Thus,
// gci2lci(gcis[i]) == i. Ideally, the caller never actually calls gci2lci(),
// and instead uses the information from get_owned_glblcells to determine
// local cell indices.
template <typename ES>
void QLT<ES>::get_owned_glblcells (std::vector<Long>& gcis) const {
  gcis.resize(ns_->levels[0].nodes.size());
  for (const auto& idx : ns_->levels[0].nodes) {
    const auto n = ns_->node_h(idx);
    gcis[n->offset] = n->id;
  }
}

// For global cell index cellidx, i.e., the globally unique ordinal associated
// with a cell in the caller's tree, return this rank's local index for
// it. This is not an efficient operation.
template <typename ES>
Int QLT<ES>::gci2lci (const Int& gci) const {
  const auto it = gci2lci_->find(gci);
  if (it == gci2lci_->end()) {
    pr(puf(gci));
    std::vector<Long> gcis;
    get_owned_glblcells(gcis);
    mprarr(gcis);
  }
  cedr_throw_if(it == gci2lci_->end(), "gci " << gci << " not in gci2lci map.");
  return it->second;
}

template <typename ES>
void QLT<ES>::declare_tracer (int problem_type, const Int& rhomidx) {
  cedr_throw_if( ! mdb_, "end_tracer_declarations was already called; "
                 "it is an error to call declare_tracer now.");
  cedr_throw_if(rhomidx > 0, "rhomidx > 0 is not supported yet.");
  // For its exception side effect, and to get canonical problem type, since
  // some possible problem types map to the same canonical one:
  problem_type = md_.get_problem_type(md_.get_problem_type_idx(problem_type));
  mdb_->trcr2prob.push_back(problem_type);
}

template <typename ES>
void QLT<ES>::end_tracer_declarations () {
  md_.init(*mdb_);
  mdb_ = nullptr;
}

template <typename ES>
void QLT<ES>::get_buffers_sizes (size_t& buf1, size_t& buf2) {
  const auto nslots = ns_->nslots;
  buf1 = md_.a_h.prob2bl2r[md_.nprobtypes]*nslots;
  buf2 = md_.a_h.prob2br2l[md_.nprobtypes]*nslots;
}

template <typename ES>
void QLT<ES>::set_buffers (Real* buf1, Real* buf2) {
  size_t l2r_sz, r2l_sz;
  get_buffers_sizes(l2r_sz, r2l_sz);
  bd_.init(buf1, l2r_sz, buf2, r2l_sz);
}

template <typename ES>
void QLT<ES>::finish_setup () {
  if (bd_.inited()) return;
  size_t l2r_sz, r2l_sz;
  get_buffers_sizes(l2r_sz, r2l_sz);
  bd_.init(l2r_sz, r2l_sz);
}

template <typename ES>
int QLT<ES>::get_problem_type (const Int& tracer_idx) const {
  cedr_throw_if(tracer_idx < 0 || tracer_idx > md_.a_h.trcr2prob.extent_int(0),
                "tracer_idx is out of bounds: " << tracer_idx);
  return md_.a_h.trcr2prob[tracer_idx];
}

template <typename ES>
Int QLT<ES>::get_num_tracers () const {
  return md_.a_h.trcr2prob.size();
}

template <typename ES> void QLT<ES>
::l2r_recv (const impl::NodeSets::Level& lvl, const Int& l2rndps) const {
  for (size_t i = 0; i < lvl.kids.size(); ++i) {
    const auto& mmd = lvl.kids[i];
    mpi::irecv(*p_, bd_.l2r_data.data() + mmd.offset*l2rndps, mmd.size*l2rndps,
               mmd.rank, impl::NodeSets::mpitag, &lvl.kids_req[i]);
  }
  Timer::start(Timer::waitall);
  mpi::waitall(lvl.kids_req.size(), lvl.kids_req.data());
  Timer::stop(Timer::waitall);
}

template <typename ES> void QLT<ES>
::l2r_combine_kid_data (const Int& lvlidx, const Int& l2rndps) const {
  if (cedr::impl::OnGpu<ES>::value) {
    const auto d = *nsdd_;
    const auto l2r_data = bd_.l2r_data;
    const auto a = md_.a_d;
    const Int ntracer = a.trcr2prob.size();
    const Int nfield = ntracer + 1;
    const Int lvl_os = nshd_->lvlptr(lvlidx);
    const Int N = nfield*(nshd_->lvlptr(lvlidx+1) - lvl_os);
    const auto combine_kid_data = KOKKOS_LAMBDA (const Int& k) {
      const Int il = lvl_os + k / nfield;
      const Int fi = k % nfield;
      const auto node_idx = d.lvl(il);
      const auto& n = d.node(node_idx);
      if ( ! n.nkids) return;
      cedr_kernel_assert(n.nkids == 2);
      if (fi == 0) {
        // Total density.
        l2r_data(n.offset*l2rndps) =
          (l2r_data(d.node(n.kids[0]).offset*l2rndps) +
           l2r_data(d.node(n.kids[1]).offset*l2rndps));
      } else {
        // Tracers. Order by bulk index for efficiency of memory access.
        const Int bi = fi - 1; // bulk index
        const Int ti = a.bidx2trcr(bi); // tracer (user) index
        const Int problem_type = a.trcr2prob(ti);
        const bool nonnegative = problem_type & ProblemType::nonnegative;
        const bool shapepreserve = problem_type & ProblemType::shapepreserve;
        const bool conserve = problem_type & ProblemType::conserve;
        const Int bdi = a.trcr2bl2r(ti);
        Real* const me = &l2r_data(n.offset*l2rndps + bdi);
        const auto& kid0 = d.node(n.kids[0]);
        const auto& kid1 = d.node(n.kids[1]);
        const Real* const k0 = &l2r_data(kid0.offset*l2rndps + bdi);
        const Real* const k1 = &l2r_data(kid1.offset*l2rndps + bdi);
        if (nonnegative) {
          me[0] = k0[0] + k1[0];
          if (conserve) me[1] = k0[1] + k1[1];
        } else {
          me[0] = shapepreserve ? k0[0] + k1[0] : cedr::impl::min(k0[0], k1[0]);
          me[1] = k0[1] + k1[1];
          me[2] = shapepreserve ? k0[2] + k1[2] : cedr::impl::max(k0[2], k1[2]);
          if (conserve) me[3] = k0[3] + k1[3] ;
        }
      }
    };
    Kokkos::parallel_for(Kokkos::RangePolicy<ES>(0, N), combine_kid_data);
    Kokkos::fence();
  } else {
    const auto& lvl = ns_->levels[lvlidx];
    const Int n_lvl_nodes = lvl.nodes.size();
#ifdef KOKKOS_ENABLE_OPENMP
#   pragma omp parallel for
#endif
    for (Int ni = 0; ni < n_lvl_nodes; ++ni) {
      const auto lvlidx = lvl.nodes[ni];
      const auto n = ns_->node_h(lvlidx);
      if ( ! n->nkids) continue;
      cedr_assert(n->nkids == 2);
      // Total density.
      bd_.l2r_data(n->offset*l2rndps) =
        (bd_.l2r_data(ns_->node_h(n->kids[0])->offset*l2rndps) +
         bd_.l2r_data(ns_->node_h(n->kids[1])->offset*l2rndps));
      // Tracers.
      for (Int pti = 0; pti < md_.nprobtypes; ++pti) {
        const Int problem_type = md_.get_problem_type(pti);
        const bool nonnegative = problem_type & ProblemType::nonnegative;
        const bool shapepreserve = problem_type & ProblemType::shapepreserve;
        const bool conserve = problem_type & ProblemType::conserve;
        const Int bis = md_.a_d.prob2trcrptr[pti], bie = md_.a_d.prob2trcrptr[pti+1];
        for (Int bi = bis; bi < bie; ++bi) {
          const Int bdi = md_.a_d.trcr2bl2r(md_.a_d.bidx2trcr(bi));
          Real* const me = &bd_.l2r_data(n->offset*l2rndps + bdi);
          const auto kid0 = ns_->node_h(n->kids[0]);
          const auto kid1 = ns_->node_h(n->kids[1]);
          const Real* const k0 = &bd_.l2r_data(kid0->offset*l2rndps + bdi);
          const Real* const k1 = &bd_.l2r_data(kid1->offset*l2rndps + bdi);
          if (nonnegative) {
            me[0] = k0[0] + k1[0];
            if (conserve) me[1] = k0[1] + k1[1];
          } else {
            me[0] = shapepreserve ? k0[0] + k1[0] : cedr::impl::min(k0[0], k1[0]);
            me[1] = k0[1] + k1[1];
            me[2] = shapepreserve ? k0[2] + k1[2] : cedr::impl::max(k0[2], k1[2]);
            if (conserve) me[3] = k0[3] + k1[3] ;
          }
        }
      }
    }
  }
}

template <typename ES> void QLT<ES>
::l2r_send_to_parents (const impl::NodeSets::Level& lvl, const Int& l2rndps) const {
  for (size_t i = 0; i < lvl.me.size(); ++i) {
    const auto& mmd = lvl.me[i];
    mpi::isend(*p_, bd_.l2r_data.data() + mmd.offset*l2rndps, mmd.size*l2rndps,
               mmd.rank, impl::NodeSets::mpitag);
  }  
}

template <typename ES> void QLT<ES>
::root_compute (const Int& l2rndps, const Int& r2lndps) const {
  if (ns_->levels.empty() || ns_->levels.back().nodes.size() != 1 ||
      ns_->node_h(ns_->levels.back().nodes[0])->parent >= 0)
    return;
  const auto d = *nsdd_;
  const auto l2r_data = bd_.l2r_data;
  const auto r2l_data = bd_.r2l_data;
  const auto a = md_.a_d;
  const Int nlev = nshd_->lvlptr.size() - 1;
  const Int node_idx = nshd_->lvl(nshd_->lvlptr(nlev-1));
  const Int ntracer = a.trcr2prob.size();
  const auto compute = KOKKOS_LAMBDA (const Int& bi) {
    const auto& n = d.node(node_idx);
    const Int ti = a.bidx2trcr(bi);
    const Int problem_type = a.trcr2prob(ti);
    const Int l2rbdi = a.trcr2bl2r(a.bidx2trcr(bi));
    const Int r2lbdi = a.trcr2br2l(a.bidx2trcr(bi));
    // If QLT is enforcing global mass conservation, set root's r2l Qm value to
    // the l2r Qm_prev's sum; otherwise, copy the l2r Qm value to the r2l one.
    const Int os = (problem_type & ProblemType::conserve ?
                    md_.get_problem_type_l2r_bulk_size(problem_type) - 1 :
                    (problem_type & ProblemType::nonnegative ? 0 : 1));
    r2l_data(n.offset*r2lndps + r2lbdi) = l2r_data(n.offset*l2rndps + l2rbdi + os);
    if ((problem_type & ProblemType::consistent) &&
        ! (problem_type & ProblemType::shapepreserve)) {
      // Consistent but not shape preserving, so we're solving a dynamic range
      // preservation problem. We now know the global q_{min,max}. Start
      // propagating it leafward.
      r2l_data(n.offset*r2lndps + r2lbdi + 1) = l2r_data(n.offset*l2rndps + l2rbdi + 0);
      r2l_data(n.offset*r2lndps + r2lbdi + 2) = l2r_data(n.offset*l2rndps + l2rbdi + 2);
    }
  };
  Kokkos::parallel_for(Kokkos::RangePolicy<ES>(0, ntracer), compute);
  Kokkos::fence();
}

template <typename ES> void QLT<ES>
::r2l_recv (const impl::NodeSets::Level& lvl, const Int& r2lndps) const {
  for (size_t i = 0; i < lvl.me.size(); ++i) {
    const auto& mmd = lvl.me[i];
    mpi::irecv(*p_, bd_.r2l_data.data() + mmd.offset*r2lndps, mmd.size*r2lndps,
               mmd.rank, impl::NodeSets::mpitag, &lvl.me_recv_req[i]);
  }
  Timer::start(Timer::waitall);
  mpi::waitall(lvl.me_recv_req.size(), lvl.me_recv_req.data());
  Timer::stop(Timer::waitall);
}

template <typename Data> KOKKOS_INLINE_FUNCTION
void r2l_solve_qp_set_q (
  const Data& l2r_data, const Data& r2l_data,
  const Int& os, const Int& l2rndps, const Int& r2lndps,
  const Int& l2rbdi, const Int& r2lbdi, const Real& q_min, const Real& q_max)
{
  l2r_data(os*l2rndps + l2rbdi + 0) = q_min;
  l2r_data(os*l2rndps + l2rbdi + 2) = q_max;
  r2l_data(os*r2lndps + r2lbdi + 1) = q_min;
  r2l_data(os*r2lndps + r2lbdi + 2) = q_max; 
}

template <typename Data> KOKKOS_INLINE_FUNCTION
void r2l_solve_qp_solve_node_problem (
  const Data& l2r_data, const Data& r2l_data, const Int& problem_type,
  const impl::NodeSets::Node& n,
  const impl::NodeSets::Node& k0, const impl::NodeSets::Node& k1,
  const Int& l2rndps, const Int& r2lndps,
  const Int& l2rbdi, const Int& r2lbdi,
  const bool prefer_mass_con_to_bounds)
{
  impl::solve_node_problem(
    problem_type,
     l2r_data( n.offset*l2rndps),
    &l2r_data( n.offset*l2rndps + l2rbdi),
     r2l_data( n.offset*r2lndps + r2lbdi),
     l2r_data(k0.offset*l2rndps),
    &l2r_data(k0.offset*l2rndps + l2rbdi),
     r2l_data(k0.offset*r2lndps + r2lbdi),
     l2r_data(k1.offset*l2rndps),
    &l2r_data(k1.offset*l2rndps + l2rbdi),
     r2l_data(k1.offset*r2lndps + r2lbdi),
    prefer_mass_con_to_bounds);
}

template <typename ES> void QLT<ES>
::r2l_solve_qp (const Int& lvlidx, const Int& l2rndps, const Int& r2lndps) const {
  Timer::start(Timer::snp);
  const bool prefer_mass_con_to_bounds =
    options_.prefer_numerical_mass_conservation_to_numerical_bounds;
  if (cedr::impl::OnGpu<ES>::value) {
    const auto d = *nsdd_;
    const auto l2r_data = bd_.l2r_data;
    const auto r2l_data = bd_.r2l_data;
    const auto a = md_.a_d;
    const Int ntracer = a.trcr2prob.size();
    const Int lvl_os = nshd_->lvlptr(lvlidx);
    const Int N = ntracer*(nshd_->lvlptr(lvlidx+1) - lvl_os);
    const auto solve_qp = KOKKOS_LAMBDA (const Int& k) {
      const Int il = lvl_os + k / ntracer;
      const Int bi = k % ntracer;
      const auto node_idx = d.lvl(il);
      const auto& n = d.node(node_idx);
      if ( ! n.nkids) return;
      const Int ti = a.bidx2trcr(bi);
      const Int problem_type = a.trcr2prob(ti);
      const Int l2rbdi = a.trcr2bl2r(a.bidx2trcr(bi));
      const Int r2lbdi = a.trcr2br2l(a.bidx2trcr(bi));
      cedr_kernel_assert(n.nkids == 2);
      if ((problem_type & ProblemType::consistent) &&
          ! (problem_type & ProblemType::shapepreserve)) {
        // Pass q_{min,max} info along. l2r data are updated for use in
        // solve_node_problem. r2l data are updated for use in isend.
        const Real q_min = r2l_data(n.offset*r2lndps + r2lbdi + 1);
        const Real q_max = r2l_data(n.offset*r2lndps + r2lbdi + 2);
        l2r_data(n.offset*l2rndps + l2rbdi + 0) = q_min;
        l2r_data(n.offset*l2rndps + l2rbdi + 2) = q_max;
        for (Int k = 0; k < 2; ++k)
          r2l_solve_qp_set_q(l2r_data, r2l_data, d.node(n.kids[k]).offset,
                             l2rndps, r2lndps, l2rbdi, r2lbdi, q_min, q_max);
      }
      r2l_solve_qp_solve_node_problem(
        l2r_data, r2l_data, problem_type, n, d.node(n.kids[0]), d.node(n.kids[1]),
        l2rndps, r2lndps, l2rbdi, r2lbdi, prefer_mass_con_to_bounds);
    };
    Kokkos::parallel_for(Kokkos::RangePolicy<ES>(0, N), solve_qp);
    Kokkos::fence();
  } else {
    const auto& lvl = ns_->levels[lvlidx];
    const Int n_lvl_nodes = lvl.nodes.size();
#ifdef KOKKOS_ENABLE_OPENMP
#   pragma omp parallel for
#endif
    for (Int ni = 0; ni < n_lvl_nodes; ++ni) {
      const auto lvlidx = lvl.nodes[ni];
      const auto n = ns_->node_h(lvlidx);
      if ( ! n->nkids) continue;
      for (Int pti = 0; pti < md_.nprobtypes; ++pti) {
        const Int problem_type = md_.get_problem_type(pti);
        const Int bis = md_.a_d.prob2trcrptr[pti], bie = md_.a_d.prob2trcrptr[pti+1];
        for (Int bi = bis; bi < bie; ++bi) {
          const Int l2rbdi = md_.a_d.trcr2bl2r(md_.a_d.bidx2trcr(bi));
          const Int r2lbdi = md_.a_d.trcr2br2l(md_.a_d.bidx2trcr(bi));
          cedr_assert(n->nkids == 2);
          if ((problem_type & ProblemType::consistent) &&
              ! (problem_type & ProblemType::shapepreserve)) {
            const Real q_min = bd_.r2l_data(n->offset*r2lndps + r2lbdi + 1);
            const Real q_max = bd_.r2l_data(n->offset*r2lndps + r2lbdi + 2);
            bd_.l2r_data(n->offset*l2rndps + l2rbdi + 0) = q_min;
            bd_.l2r_data(n->offset*l2rndps + l2rbdi + 2) = q_max;
            for (Int k = 0; k < 2; ++k)
              r2l_solve_qp_set_q(bd_.l2r_data, bd_.r2l_data,
                                 ns_->node_h(n->kids[k])->offset,
                                 l2rndps, r2lndps, l2rbdi, r2lbdi, q_min, q_max);
          }
          r2l_solve_qp_solve_node_problem(
            bd_.l2r_data, bd_.r2l_data, problem_type, *n, *ns_->node_h(n->kids[0]),
            *ns_->node_h(n->kids[1]), l2rndps, r2lndps, l2rbdi, r2lbdi,
            prefer_mass_con_to_bounds);
        }
      }
    }
  }
  Timer::stop(Timer::snp);
}

template <typename ES> void QLT<ES>
::r2l_send_to_kids (const impl::NodeSets::Level& lvl, const Int& r2lndps) const {
  for (size_t i = 0; i < lvl.kids.size(); ++i) {
    const auto& mmd = lvl.kids[i];
    mpi::isend(*p_, bd_.r2l_data.data() + mmd.offset*r2lndps, mmd.size*r2lndps,
               mmd.rank, impl::NodeSets::mpitag);
  }
}

template <typename ES>
void QLT<ES>::run () {
  cedr_assert(bd_.inited());
  Timer::start(Timer::qltrunl2r);
  // Number of data per slot.
  const Int l2rndps = md_.a_h.prob2bl2r[md_.nprobtypes];
  const Int r2lndps = md_.a_h.prob2br2l[md_.nprobtypes];
  for (size_t il = 0; il < ns_->levels.size(); ++il) {
    auto& lvl = ns_->levels[il];
    if (lvl.kids.size()) l2r_recv(lvl, l2rndps);
    l2r_combine_kid_data(il, l2rndps);    
    if (lvl.me.size()) l2r_send_to_parents(lvl, l2rndps);
  }
  Timer::stop(Timer::qltrunl2r); Timer::start(Timer::qltrunr2l);
  root_compute(l2rndps, r2lndps);
  for (size_t il = ns_->levels.size(); il > 0; --il) {
    auto& lvl = ns_->levels[il-1];
    if (lvl.me.size()) r2l_recv(lvl, r2lndps);
    r2l_solve_qp(il-1, l2rndps, r2lndps);
    if (lvl.kids.size()) r2l_send_to_kids(lvl, r2lndps);
  }
  Timer::stop(Timer::qltrunr2l);
}

namespace test {
using namespace impl;

class TestQLT : public cedr::test::TestRandomized {
public:
  typedef QLT<Kokkos::DefaultExecutionSpace> QLTT;

  TestQLT (const Parallel::Ptr& p, const tree::Node::Ptr& tree,
           const Int& ncells, const bool external_memory, const bool verbose,
           CDR::Options options)
    : TestRandomized("QLT", p, ncells, verbose, options),
      qlt_(p, ncells, tree, options), tree_(tree), external_memory_(external_memory)
  {
    if (verbose) qlt_.print(std::cout);
    init();
  }

private:
  QLTT qlt_;
  tree::Node::Ptr tree_;
  bool external_memory_;
  typename QLTT::RealList buf1_, buf2_;

  CDR& get_cdr () override { return qlt_; }

  void init_numbering () override {
    init_numbering(tree_);
  }

  void init_numbering (const tree::Node::Ptr& node) {
    check(qlt_);
    // TestQLT doesn't actually care about a particular ordering, as there is no
    // geometry to the test problem. However, use *some* ordering to model what
    // a real problem must do.
    if ( ! node->nkids) {
      if (node->rank == p_->rank())
        gcis_.push_back(node->cellidx);
      return;
    }
    for (Int i = 0; i < node->nkids; ++i)
      init_numbering(node->kids[i]);
  }

  static void check (const QLTT& qlt) {
    const Int n = qlt.nlclcells();
    std::vector<Long> gcis;
    qlt.get_owned_glblcells(gcis);
    cedr_assert(static_cast<Int>(gcis.size()) == n);
    for (Int i = 0; i < n; ++i)
      cedr_assert(qlt.gci2lci(gcis[i]) == i);
  }

  void init_tracers () override {
    for (const auto& t : tracers_)
      qlt_.declare_tracer(t.problem_type, 0);
    qlt_.end_tracer_declarations();
    if (external_memory_) {
      size_t l2r_sz, r2l_sz;
      qlt_.get_buffers_sizes(l2r_sz, r2l_sz);
      buf1_ = typename QLTT::RealList("buf1", l2r_sz);
      buf2_ = typename QLTT::RealList("buf2", r2l_sz);
      qlt_.set_buffers(buf1_.data(), buf2_.data());
    }
    qlt_.finish_setup();
    cedr_assert(qlt_.get_num_tracers() == static_cast<Int>(tracers_.size()));
    for (size_t i = 0; i < tracers_.size(); ++i) {
      const auto pt = qlt_.get_problem_type(i);
      cedr_assert((pt == ((tracers_[i].problem_type | ProblemType::consistent) &
                          ~ProblemType::nonnegative)) ||
                  (pt == ((tracers_[i].problem_type | ProblemType::nonnegative) &
                          ~ProblemType::consistent)));
    }
  }
  
  void run_impl (const Int trial) override {
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
};

// Test all QLT variations and situations.
Int test_qlt (const Parallel::Ptr& p, const tree::Node::Ptr& tree,
              const Int& ncells, const Int nrepeat,
              const bool write, const bool external_memory,
              const bool prefer_mass_con_to_bounds, const bool verbose) {
  CDR::Options options;
  options.prefer_numerical_mass_conservation_to_numerical_bounds =
    prefer_mass_con_to_bounds;
  return TestQLT(p, tree, ncells, external_memory, verbose, options)
    .run<TestQLT::QLTT>(nrepeat, write);
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
    cedr_throw_if(nranks_ > nc_, "#GIDs < #ranks is not supported.");
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
} // namespace test
} // namespace oned

tree::Node::Ptr tree::make_tree_over_1d_mesh (const Parallel::Ptr& p, const Int& ncells,
                                              const bool imbalanced) {
  return oned::make_tree(oned::Mesh(ncells, p), imbalanced);
}

namespace test {
Int unittest_NodeSets (const Parallel::Ptr& p) {
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
  for (size_t is = 0, islim = sizeof(szs)/sizeof(*szs); is < islim; ++is) {
    for (size_t id = 0, idlim = sizeof(dists)/sizeof(*dists); id < idlim; ++id) {
      for (bool imbalanced : {false, true}) {
        for (bool prefer_mass_con_to_bounds : {false, true}) {
          const auto external_memory = imbalanced;
          if (p->amroot()) {
            std::cout << " (" << szs[is] << ", " << id << ", " << imbalanced << ", "
                      << prefer_mass_con_to_bounds << ")";
            std::cout.flush();
          }
          Mesh m(szs[is], p, dists[id]);
          tree::Node::Ptr tree = make_tree(m, imbalanced);
          const bool write = (write_requested && m.ncell() < 3000 &&
                              is == islim-1 && id == idlim-1);
          nerr += test::test_qlt(p, tree, m.ncell(), 1, write, external_memory,
                                 prefer_mass_con_to_bounds, false);
        }
      }
    }
  }
  return nerr;
}

Int run_unit_and_randomized_tests (const Parallel::Ptr& p, const Input& in) {
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
  // Performance test.
  if (in.perftest && in.ncells > 0) {
    oned::Mesh m(in.ncells, p,
                 (in.pseudorandom ?
                  oned::Mesh::ParallelDecomp::pseudorandom :
                  oned::Mesh::ParallelDecomp::contiguous));
    Timer::init();
    Timer::start(Timer::total); Timer::start(Timer::tree);
    tree::Node::Ptr tree = make_tree(m, false);
    Timer::stop(Timer::tree);
    test::test_qlt(p, tree, in.ncells, in.nrepeat, false, false, false, in.verbose);
    Timer::stop(Timer::total);
    if (p->amroot()) Timer::print();
  }
  return nerr;
}

} // namespace test
} // namespace qlt
} // namespace cedr

#ifdef KOKKOS_ENABLE_SERIAL
template class cedr::qlt::QLT<Kokkos::Serial>;
#endif
#ifdef KOKKOS_ENABLE_OPENMP
template class cedr::qlt::QLT<Kokkos::OpenMP>;
#endif
#ifdef KOKKOS_ENABLE_CUDA
template class cedr::qlt::QLT<Kokkos::Cuda>;
#endif
#ifdef KOKKOS_ENABLE_THREADS
template class cedr::qlt::QLT<Kokkos::Threads>;
#endif

//>> cedr_caas.cpp
// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

//#include "cedr_caas.hpp"
//#include "cedr_util.hpp"
//#include "cedr_test_randomized.hpp"

namespace Kokkos {
struct Real2 {
  cedr::Real v[2];
  KOKKOS_INLINE_FUNCTION Real2 () { v[0] = v[1] = 0; }

  KOKKOS_INLINE_FUNCTION void operator= (const Real2& s) {
    v[0] = s.v[0];
    v[1] = s.v[1];
  }
  KOKKOS_INLINE_FUNCTION void operator= (const volatile Real2& s) volatile {
    v[0] = s.v[0];
    v[1] = s.v[1];
  }

  KOKKOS_INLINE_FUNCTION Real2& operator+= (const Real2& o) {
    v[0] += o.v[0];
    v[1] += o.v[1];
    return *this;
  }
};

template<> struct reduction_identity<Real2> {
  KOKKOS_INLINE_FUNCTION static Real2 sum() { return Real2(); }
};
} // namespace Kokkos

namespace cedr {
namespace caas {

template <typename ES>
CAAS<ES>::CAAS (const mpi::Parallel::Ptr& p, const Int nlclcells,
                const typename UserAllReducer::Ptr& uar)
  : p_(p), user_reducer_(uar), nlclcells_(nlclcells), nrhomidxs_(0),
    need_conserve_(false), finished_setup_(false)
{
  cedr_throw_if(nlclcells == 0, "CAAS does not support 0 cells on a rank.");
  tracer_decls_ = std::make_shared<std::vector<Decl> >();  
}

template <typename ES>
void CAAS<ES>::declare_tracer(int problem_type, const Int& rhomidx) {
  cedr_throw_if( ! (problem_type & ProblemType::shapepreserve),
                "CAAS does not support ! shapepreserve yet.");
  cedr_throw_if(rhomidx > 0, "rhomidx > 0 is not supported yet.");
  tracer_decls_->push_back(Decl(problem_type, rhomidx));
  if (problem_type & ProblemType::conserve)
    need_conserve_ = true;
  nrhomidxs_ = std::max(nrhomidxs_, rhomidx+1);
}

template <typename ES>
void CAAS<ES>::end_tracer_declarations () {
  cedr_throw_if(tracer_decls_->size() == 0, "#tracers is 0.");
  cedr_throw_if(nrhomidxs_ == 0, "#rhomidxs is 0.");
  probs_ = IntList("CAAS probs", static_cast<Int>(tracer_decls_->size()));
  probs_h_ = Kokkos::create_mirror_view(probs_);
  //t2r_ = IntList("CAAS t2r", static_cast<Int>(tracer_decls_->size()));
  for (Int i = 0; i < probs_.extent_int(0); ++i) {
    probs_h_(i) = (*tracer_decls_)[i].probtype;
    //t2r_(i) = (*tracer_decls_)[i].rhomidx;
  }
  Kokkos::deep_copy(probs_, probs_h_);
  tracer_decls_ = nullptr;
}

template <typename ES>
void CAAS<ES>::get_buffers_sizes (size_t& buf1, size_t& buf2, size_t& buf3) {
  const Int e = need_conserve_ ? 1 : 0;
  const auto nslots = 4*probs_.size();
  buf1 = nlclcells_ * ((3+e)*probs_.size() + 1);
  cedr_assert( ! user_reducer_ || nlclcells_ % user_reducer_->n_accum_in_place() == 0);
  buf2 = nslots*(user_reducer_ ? (nlclcells_ / user_reducer_->n_accum_in_place()) : 1);
  buf3 = nslots;
}

template <typename ES>
void CAAS<ES>::get_buffers_sizes (size_t& buf1, size_t& buf2) {
  size_t buf3;
  get_buffers_sizes(buf1, buf2, buf3);
  buf2 += buf3;
}

template <typename ES>
void CAAS<ES>::set_buffers (Real* buf1, Real* buf2) {
  size_t buf1sz, buf2sz, buf3sz;
  get_buffers_sizes(buf1sz, buf2sz, buf3sz);
  d_ = RealList(buf1, buf1sz);
  send_ = RealList(buf2, buf2sz);
  recv_ = RealList(buf2 + buf2sz, buf3sz);
}

template <typename ES>
void CAAS<ES>::finish_setup () {
  if (recv_.size() > 0) {
    finished_setup_ = true;
    return;
  }
  size_t buf1, buf2, buf3;
  get_buffers_sizes(buf1, buf2, buf3);
  // (rho, Qm, Qm_min, Qm_max, [Qm_prev])
  d_ = RealList("CAAS data", buf1);
  // (e'Qm_clip, e'Qm, e'Qm_min, e'Qm_max, [e'Qm_prev])
  send_ = RealList("CAAS send", buf2);
  recv_ = RealList("CAAS recv", buf3);
  finished_setup_ = true;
}

template <typename ES>
int CAAS<ES>::get_problem_type (const Int& tracer_idx) const {
  cedr_assert(tracer_idx >= 0 && tracer_idx < probs_.extent_int(0));
  return probs_h_[tracer_idx];
}

template <typename ES>
Int CAAS<ES>::get_num_tracers () const {
  return probs_.extent_int(0);
}

template <typename RealList, typename IntList>
KOKKOS_INLINE_FUNCTION static void
calc_Qm_scalars (const RealList& d, const IntList& probs,
                 const Int& nt, const Int& nlclcells,
                 const Int& k, const Int& os, const Int& i,
                 Real& Qm_clip, Real& Qm_term) {
  const Real Qm = d(os+i);
  Qm_term = (probs(k) & ProblemType::conserve ?
             d(os + nlclcells*3*nt + i) /* Qm_prev */ :
             Qm);
  const Real Qm_min = d(os + nlclcells*  nt + i);
  const Real Qm_max = d(os + nlclcells*2*nt + i);
  Qm_clip = cedr::impl::min(Qm_max, cedr::impl::max(Qm_min, Qm));
}

template <typename Func>
void homme_parallel_for (const Int& beg, const Int& end, const Func& f) {
#ifdef HORIZ_OPENMP
# pragma omp for
#endif
  for (Int i = beg; i < end; ++i)
    f(i);
}

template <typename ES>
void CAAS<ES>::reduce_locally () {
  const bool user_reduces = user_reducer_ != nullptr;
  ConstExceptGnu Int nt = probs_.size(), nlclcells = nlclcells_;

  const auto probs = probs_;
  const auto send = send_;
  const auto d = d_;
  if (user_reduces) {
    const Int n_accum_in_place = user_reducer_->n_accum_in_place();
    const Int nlclaccum = nlclcells / n_accum_in_place;
    const auto calc_Qm_clip = KOKKOS_LAMBDA (const Int& j) {
      const auto k = j / nlclaccum;
      const auto bi = j % nlclaccum;
      const auto os = (k+1)*nlclcells;
      Real accum_clip = 0, accum_term = 0;
      for (Int ai = 0; ai < n_accum_in_place; ++ai) {
        const Int i = n_accum_in_place*bi + ai;
        Real Qm_clip, Qm_term;
        calc_Qm_scalars(d, probs, nt, nlclcells, k, os, i, Qm_clip, Qm_term);
        d(os + i) = Qm_clip;
        accum_clip += Qm_clip;
        accum_term += Qm_term;
      }
      send(nlclaccum*      k  + bi) = accum_clip;
      send(nlclaccum*(nt + k) + bi) = accum_term;
    };
    homme_parallel_for(0, nt*nlclaccum, calc_Qm_clip);
    const auto set_Qm_minmax = KOKKOS_LAMBDA (const Int& j) {
      const auto k = 2*nt + j / nlclaccum;
      const auto bi = j % nlclaccum;
      const auto os = (k-nt+1)*nlclcells;
      Real accum_ext = 0;
      for (Int ai = 0; ai < n_accum_in_place; ++ai) {
        const Int i = n_accum_in_place*bi + ai;
        accum_ext += d(os + i);
      }
      send(nlclaccum*k + bi) = accum_ext;
    };
    homme_parallel_for(0, 2*nt*nlclaccum, set_Qm_minmax);
  } else {
    using ESU = cedr::impl::ExeSpaceUtils<ES>;
    const auto calc_Qm_clip = KOKKOS_LAMBDA (const typename ESU::Member& t) {
      const auto k = t.league_rank();
      const auto os = (k+1)*nlclcells;
      const auto reduce = [&] (const Int& i, Kokkos::Real2& accum) {
        Real Qm_clip, Qm_term;
        calc_Qm_scalars(d, probs, nt, nlclcells, k, os, i, Qm_clip, Qm_term);
        d(os+i) = Qm_clip;
        accum.v[0] += Qm_clip;
        accum.v[1] += Qm_term;
      };
      Kokkos::Real2 accum;
      Kokkos::parallel_reduce(Kokkos::TeamThreadRange(t, nlclcells),
                              reduce, Kokkos::Sum<Kokkos::Real2>(accum));
      send(     k) = accum.v[0];
      send(nt + k) = accum.v[1];
    };
    Kokkos::parallel_for(ESU::get_default_team_policy(nt, nlclcells),
                         calc_Qm_clip);
    const auto set_Qm_minmax = KOKKOS_LAMBDA (const typename ESU::Member& t) {
      const auto k = 2*nt + t.league_rank();
      const auto os = (k-nt+1)*nlclcells;
      Real accum = 0;
      Kokkos::parallel_reduce(Kokkos::TeamThreadRange(t, nlclcells),
                              [&] (const Int& i, Real& accum) { accum += d(os+i); },
                              Kokkos::Sum<Real>(accum));
      send(k) = accum;
    };
    Kokkos::parallel_for(ESU::get_default_team_policy(2*nt, nlclcells),
                         set_Qm_minmax);
  }
}

template <typename ES>
void CAAS<ES>::reduce_globally () {
  const int err = mpi::all_reduce(*p_, send_.data(), recv_.data(),
                                  send_.size(), MPI_SUM);
  cedr_throw_if(err != MPI_SUCCESS,
                "CAAS::reduce_globally MPI_Allreduce returned " << err);
}

template <typename ES>
void CAAS<ES>::finish_locally () {
  using ESU = cedr::impl::ExeSpaceUtils<ES>;
  ConstExceptGnu Int nt = probs_.size(), nlclcells = nlclcells_;
  const auto recv = recv_;
  const auto d = d_;
  const auto adjust_Qm = KOKKOS_LAMBDA (const Int& k) {
    const auto os = (k+1)*nlclcells;
    const auto Qm_clip_sum = recv(     k);
    const auto Qm_sum      = recv(nt + k);
    const auto m = Qm_sum - Qm_clip_sum;
    if (m < 0) {
      const auto Qm_min_sum = recv(2*nt + k);
      auto fac = Qm_clip_sum - Qm_min_sum;
      if (fac > 0) {
        fac = m/fac;
        for (Int i = 0; i < nlclcells; ++i) {
          const auto Qm_min = d(os + nlclcells * nt + i);
          auto& Qm = d(os+i);
          Qm += fac*(Qm - Qm_min);
          Qm = impl::max(Qm_min, Qm);
        };
      }
    } else if (m > 0) {
      const auto Qm_max_sum = recv(3*nt + k);
      auto fac = Qm_max_sum - Qm_clip_sum;
      if (fac > 0) {
        fac = m/fac;
        for (Int i = 0; i < nlclcells; ++i) {
          const auto Qm_max = d(os + nlclcells*2*nt + i);
          auto& Qm = d(os+i);
          Qm += fac*(Qm_max - Qm);
          Qm = impl::min(Qm_max, Qm);
        };
      }
    }
  };
  homme_parallel_for(0, nt, adjust_Qm);
}

template <typename ES>
void CAAS<ES>::run () {
  cedr_assert(finished_setup_);
  reduce_locally();
  const bool user_reduces = user_reducer_ != nullptr;
  if (user_reduces)
    (*user_reducer_)(*p_, send_.data(), recv_.data(),
                     nlclcells_ / user_reducer_->n_accum_in_place(),
                     recv_.size(), MPI_SUM);
  else
    reduce_globally();
  finish_locally();
}

namespace test {
struct TestCAAS : public cedr::test::TestRandomized {
  typedef CAAS<Kokkos::DefaultExecutionSpace> CAAST;

  struct TestAllReducer : public CAAST::UserAllReducer {
    int operator() (const mpi::Parallel& p, Real* sendbuf, Real* rcvbuf,
                    int nlcl, int count, MPI_Op op) const override {
      Kokkos::View<Real*> s(sendbuf, nlcl*count), r(rcvbuf, count);
      const auto s_h = Kokkos::create_mirror_view(s);
      Kokkos::deep_copy(s_h, s);
      const auto r_h = Kokkos::create_mirror_view(r);
      for (int i = 1; i < nlcl; ++i)
        s_h(0) += s_h(i);
      for (int k = 1; k < count; ++k) {
        s_h(k) = s_h(nlcl*k);
        for (int i = 1; i < nlcl; ++i)
          s_h(k) += s_h(nlcl*k + i);
      }
      const int err = mpi::all_reduce(p, s_h.data(), r_h.data(), count, op);
      Kokkos::deep_copy(r, r_h);
      return err;
    }
  };

  TestCAAS (const mpi::Parallel::Ptr& p, const Int& ncells,
            const bool use_own_reducer, const bool external_memory,
            const bool verbose)
    : TestRandomized("CAAS", p, ncells, verbose),
      p_(p), external_memory_(external_memory)
  {
    const auto np = p->size(), rank = p->rank();
    nlclcells_ = ncells / np;
    const Int todo = ncells - nlclcells_ * np;
    if (rank < todo) ++nlclcells_;
    caas_ = std::make_shared<CAAST>(
      p, nlclcells_,
      use_own_reducer ? std::make_shared<TestAllReducer>() : nullptr);
    init();
  }

  CDR& get_cdr () override { return *caas_; }

  void init_numbering () override {
    const auto np = p_->size(), rank = p_->rank();
    Int start = 0;
    for (Int lrank = 0; lrank < rank; ++lrank)
      start += get_nllclcells(ncells_, np, lrank);
    gcis_.resize(nlclcells_);
    for (Int i = 0; i < nlclcells_; ++i)
      gcis_[i] = start + i;
  }

  void init_tracers () override {
    // CAAS doesn't yet support everything, so remove a bunch of the tracers.
    std::vector<TestRandomized::Tracer> tracers;
    Int idx = 0;
    for (auto& t : tracers_) {
      if ( ! (t.problem_type & ProblemType::shapepreserve) ||
           ! t.local_should_hold)
        continue;
      t.idx = idx++;
      tracers.push_back(t);
      caas_->declare_tracer(t.problem_type, 0);
    }
    tracers_ = tracers;
    caas_->end_tracer_declarations();
    if (external_memory_) {
      size_t l2r_sz, r2l_sz;
      caas_->get_buffers_sizes(l2r_sz, r2l_sz);
      buf1_ = typename CAAST::RealList("buf1", l2r_sz);
      buf2_ = typename CAAST::RealList("buf2", r2l_sz);
      caas_->set_buffers(buf1_.data(), buf2_.data());
    }
    caas_->finish_setup();
  }

  void run_impl (const Int trial) override {
    caas_->run();
  }

private:
  mpi::Parallel::Ptr p_;
  bool external_memory_;
  Int nlclcells_;
  CAAST::Ptr caas_;
  typename CAAST::RealList buf1_, buf2_;

  static Int get_nllclcells (const Int& ncells, const Int& np, const Int& rank) {
    Int nlclcells = ncells / np;
    const Int todo = ncells - nlclcells * np;
    if (rank < todo) ++nlclcells;
    return nlclcells;
  }
};

Int unittest (const mpi::Parallel::Ptr& p) {
  const auto np = p->size();
  Int nerr = 0;
  for (Int nlclcells : {1, 2, 4, 11}) {
    Long ncells = np*nlclcells;
    if (ncells > np) ncells -= np/2;
    for (const bool own_reducer : {false, true})
      for (const bool external_memory : {false, true})
        nerr += TestCAAS(p, ncells, own_reducer, external_memory, false)
          .run<TestCAAS::CAAST>(1, false);
  }
  return nerr;
}
} // namespace test
} // namespace caas
} // namespace cedr

#ifdef KOKKOS_ENABLE_SERIAL
template class cedr::caas::CAAS<Kokkos::Serial>;
#endif
#ifdef KOKKOS_ENABLE_OPENMP
template class cedr::caas::CAAS<Kokkos::OpenMP>;
#endif
#ifdef KOKKOS_ENABLE_CUDA
template class cedr::caas::CAAS<Kokkos::Cuda>;
#endif
#ifdef KOKKOS_ENABLE_THREADS
template class cedr::caas::CAAS<Kokkos::Threads>;
#endif

//>> cedr_test_randomized.cpp
// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

//#include "cedr_test_randomized.hpp"

namespace cedr {
namespace test {

std::string TestRandomized::Tracer::str () const {
  std::stringstream ss;
  ss << "(ti " << idx;
  if (problem_type & PT::conserve) ss << " c";
  if (problem_type & PT::shapepreserve) ss << " s";
  if (problem_type & PT::consistent) ss << " t";
  if (problem_type & PT::nonnegative) ss << " nn";
  ss << " pt " << perturbation_type << " ssh " << safe_should_hold
     << " lsh " << local_should_hold << ")";
  return ss.str();
}

TestRandomized::Writer::~Writer () {
  if ( ! fh) return;
  fprintf(fh.get(), "  return s\n");
}

void TestRandomized::init_tracers_vector () {
  typedef Tracer::PT PT;
  static const Int pts[] = {
    PT::conserve | PT::shapepreserve | PT::consistent,
    PT::shapepreserve,
    PT::conserve | PT::consistent,
    PT::consistent,
    PT::nonnegative,
    PT::nonnegative | PT::conserve
  };
  Int tracer_idx = 0;
  for (Int perturb = 0; perturb < 6; ++perturb)
    for (size_t ti = 0; ti < sizeof(pts)/sizeof(*pts); ++ti) {
      Tracer t;
      t.problem_type = pts[ti];
      const bool shapepreserve = t.problem_type & PT::shapepreserve;
      const bool nonnegative = t.problem_type & PT::nonnegative;
      t.idx = tracer_idx++;
      t.perturbation_type = perturb;
      t.safe_should_hold = true;
      t.no_change_should_hold = perturb == 0;
      t.local_should_hold = perturb < 4 && (shapepreserve || nonnegative);
      t.write = perturb == 2 && ti == 2;
      tracers_.push_back(t);
    }
}

static Real urand () { return rand() / ((Real) RAND_MAX + 1.0); }

void TestRandomized::generate_rho (Values& v) {
  auto r = v.rhom();
  const Int n = v.ncells();
  for (Int i = 0; i < n; ++i)
    r[i] = 0.5*(1 + urand());
}

void TestRandomized::generate_Q (const Tracer& t, Values& v) {
  Real* rhom = v.rhom(), * Qm_min = v.Qm_min(t.idx), * Qm = v.Qm(t.idx),
    * Qm_max = v.Qm_max(t.idx), * Qm_prev = v.Qm_prev(t.idx);
  const Int n = v.ncells();
  const bool nonneg_only = t.problem_type & ProblemType::nonnegative;
  for (Int i = 0; i < n; ++i) {
    if (nonneg_only) {
      // Make sure the generated Qm is globally positive.
      Qm[i] = (t.no_change_should_hold ?
               urand() :
               ((i % 2 == 0) ? 0.75 : -0.75) + urand());
      // Qm_min,max are unused in QLT, but need them to be set here for
      // bookkeeping in check().
      Qm_min[i] = 0;
      Qm_max[i] = 10;
    } else {
      // The typical use case has 0 <= q_min <= q, but test with general sign.
      const Real
        q_min = -0.75 + urand(),
        q_max = q_min + urand(),
        q = q_min + (q_max - q_min)*urand();
      // Check correctness up to FP.
      cedr_assert(q_min <= q && q <= q_max);
      Qm_min[i] = q_min*rhom[i];
      Qm_max[i] = q_max*rhom[i];
      // Protect against FP error.
      Qm[i] = std::max<Real>(Qm_min[i], std::min(Qm_max[i], q*rhom[i]));
    }
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
void TestRandomized::permute_Q (const Tracer& t, Values& v) {
  Real* const Qm = v.Qm(t.idx);
  const Int N = v.ncells();
  std::vector<Int> p;
  gen_rand_perm(N, p);
  std::vector<Real> Qm_orig(N);
  std::copy(Qm, Qm + N, Qm_orig.begin());
  for (Int i = 0; i < N; ++i)
    Qm[i] = Qm_orig[p[i]];
}

void TestRandomized
::add_const_to_Q (const Tracer& t, Values& v,
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
  if (safety_problem && v.ncells()) {
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

void TestRandomized::perturb_Q (const Tracer& t, Values& v) {
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

std::string TestRandomized::get_tracer_name (const Tracer& t) {
  std::stringstream ss;
  ss << "t" << t.idx;
  return ss.str();
}

void TestRandomized::init_writer () {
  if (p_->amroot()) {
    w_ = std::make_shared<Writer>();
    w_->fh = std::unique_ptr<FILE, cedr::util::FILECloser>(fopen("out_QLT.py", "w"));
    int n = gcis_.size();
    w_->ngcis.resize(p_->size());
    mpi::gather(*p_, &n, 1, w_->ngcis.data(), 1, p_->root());
    w_->displs.resize(p_->size() + 1);
    w_->displs[0] = 0;
    for (size_t i = 0; i < w_->ngcis.size(); ++i)
      w_->displs[i+1] = w_->displs[i] + w_->ngcis[i];
    cedr_assert(w_->displs.back() == ncells_);
    w_->gcis.resize(ncells_);
    mpi::gatherv(*p_, gcis_.data(), gcis_.size(), w_->gcis.data(), w_->ngcis.data(),
                 w_->displs.data(), p_->root());
  } else {
    int n = gcis_.size();
    mpi::gather(*p_, &n, 1, static_cast<int*>(nullptr), 0, p_->root());
    Long* Lnull = nullptr;
    const int* inull = nullptr;
    mpi::gatherv(*p_, gcis_.data(), gcis_.size(), Lnull, inull, inull, p_->root());
  }
  write_inited_ = true;
}

void TestRandomized
::gather_field (const Real* Qm_lcl, std::vector<Real>& Qm_gbl,
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

void TestRandomized
::write_field (const std::string& tracer_name, const std::string& field_name,
               const std::vector<Real>& Qm) {
  if ( ! p_->amroot()) return;
  fprintf(w_->fh.get(), "  s.%s.%s = [", tracer_name.c_str(), field_name.c_str());
  for (const auto& e : Qm)
    fprintf(w_->fh.get(), "%1.15e, ", e);
  fprintf(w_->fh.get(), "]\n");
}

void TestRandomized::write_pre (const Tracer& t, Values& v) {
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

void TestRandomized::write_post (const Tracer& t, Values& v) {
  if ( ! t.write) return;
  const auto name = get_tracer_name(t);
  std::vector<Real> Qm, wrk;
  gather_field(v.Qm(t.idx), Qm, wrk);
  write_field(name, "Qm_qlt", Qm);
}

Int TestRandomized
::check (const std::string& cdr_name, const mpi::Parallel& p,
         const std::vector<Tracer>& ts, const Values& v) {
  static const bool details = true;
  static const Real eps = std::numeric_limits<Real>::epsilon();
  static const Real ulp3 = 3*eps;
  const auto prefer_mass_con_to_bounds =
    options_.prefer_numerical_mass_conservation_to_numerical_bounds;
  const Real lv_tol = prefer_mass_con_to_bounds ? 100*eps : 0;
  const Real safety_tol = prefer_mass_con_to_bounds ? 100*eps : ulp3;
  Int nerr = 0;
  std::vector<Real> lcl_mass(2*ts.size()), q_min_lcl(ts.size()), q_max_lcl(ts.size());
  std::vector<Int> t_ok(ts.size(), 1), local_violated(ts.size(), 0);
  for (size_t ti = 0; ti < ts.size(); ++ti) {
    const auto& t = ts[ti];

    cedr_assert(t.safe_should_hold);
    const bool safe_only = ! t.local_should_hold;
    const bool nonneg_only = t.problem_type & ProblemType::nonnegative;
    const Int n = v.ncells();
    const Real* rhom = v.rhom(), * Qm_min = v.Qm_min(t.idx), * Qm = v.Qm(t.idx),
      * Qm_max = v.Qm_max(t.idx), * Qm_prev = v.Qm_prev(t.idx);

    q_min_lcl[ti] =  1e3;
    q_max_lcl[ti] = -1e3;
    for (Int i = 0; i < n; ++i) {
      const bool lv = (nonneg_only ?
                       Qm[i] < 0 :
                       Qm[i] < Qm_min[i] - lv_tol || Qm[i] > Qm_max[i] + lv_tol);
      if (lv) local_violated[ti] = 1;
      if ( ! safe_only && lv) {
        // If this fails at ~ machine eps, check r2l_nl_adjust_bounds code in
        // solve_node_problem.
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
    const bool nonneg_only = t.problem_type & ProblemType::nonnegative;
    if (safe_only) {
      const Int n = v.ncells();
      const Real* rhom = v.rhom(), * Qm_min = v.Qm_min(t.idx), * Qm = v.Qm(t.idx),
        * Qm_max = v.Qm_max(t.idx);
      const Real q_min = nonneg_only ? 0 : q_min_gbl[ti], q_max = q_max_gbl[ti];
      for (Int i = 0; i < n; ++i) {
        const Real delta = (q_max - q_min)*safety_tol;
        const bool lv = (nonneg_only ?
                         Qm[i] < -ulp3 :
                         (Qm[i] < q_min*rhom[i] - delta ||
                          Qm[i] > q_max*rhom[i] + delta));
        if (lv) {
          if (details)
            pr("check q (safety) " << t.str() << ": " << q_min*rhom[i] << " "
               << Qm_min[i] << " " << Qm[i] << " " << Qm_max[i] << " "
               << q_max*rhom[i] << " | " << (Qm[i] < q_min*rhom[i] ?
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
        rd = cedr::util::reldif(desired_mass, actual_mass);
      const bool mass_failed = rd > tol;
      if (mass_failed) {
        ++nerr;
        t_ok_gbl[ti] = false;
      }
      if ( ! t_ok_gbl[ti]) {
        std::cout << "FAIL " << cdr_name << ": " << ts[ti].str();
        if (mass_failed) std::cout << " mass re " << rd;
        std::cout << "\n";
        //pr(puf(desired_mass) pu(actual_mass));
      }
    }
  }

  return nerr;
}
  
TestRandomized
::TestRandomized (const std::string& name, const mpi::Parallel::Ptr& p,
                  const Int& ncells, const bool verbose,
                  const CDR::Options options)
  : cdr_name_(name), options_(options), p_(p), ncells_(ncells),
    write_inited_(false)
{}

void TestRandomized::init () {
  init_numbering();
  init_tracers_vector();
  init_tracers();
}

} // namespace test
} // namespace cedr

//>> cedr_test_1d_transport.cpp
// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

//#include "cedr_test.hpp"
//#include "cedr_qlt.hpp"
//#include "cedr_caas.hpp"

#include <algorithm>

namespace cedr {
namespace test {
namespace transport1d {

namespace interp {
inline Real to_periodic_core (const Real& xl, const Real& xr, const Real& x) {
  if (x >= xl && x <= xr) return x;
  const Real w = xr - xl, xmxl = x - xl;
  return x - w*std::floor(xmxl / w);
}

inline Real get_slope (const Real x[2], const Real y[2]) {
  return (y[1] - y[0]) / (x[1] - x[0]);
}

inline void
get_cubic (Real dx, Real v1, Real s1, Real v2, Real s2, Real c[4]) {
  Real dx2 = dx*dx;
  Real dx3 = dx2*dx;
  Real den = -dx3;
  Real b1, b2;
  c[2] = s1;
  c[3] = v1;
  b1 = v2 - dx*c[2] - c[3];
  b2 = s2 - c[2];
  c[0] = (2.0*b1 - dx*b2) / den;
  c[1] = (-3.0*dx*b1 + dx2*b2) / den;
}

void cubic_interp_periodic (
  const Real* const x, const Int nx, const Real* const y,
  const Real* const xi, const Int nxi, Real* const yi,
  Int* const dod)
{
  const int nc = nx - 1;
#ifdef _OPENMP
# pragma omp parallel for
#endif
  for (Int j = 0; j < nxi; ++j) {
    const Real xi_per = to_periodic_core(x[0], x[nc], xi[j]);
    Int ip1 = std::upper_bound(x, x + nx, xi_per) - x;
    // Handle numerical issues at boundaries.
    if (ip1 == 0) ++ip1;
    else if (ip1 == nx) --ip1;
    const Int i = ip1 - 1;
    // Domain of dependence.
    Int* dodj = dod + 4*j;
    for (Int k = 0; k < 4; ++k)
      dodj[k] = (i - 1 + k + nc) % nc;
    // Slopes.
    const bool at_start = i == 0, at_end = i == nc - 1;
    const Real smid = get_slope(x+i, y+i);
    Real s1, s2;
    if (at_start) {
      const Real a = (x[nc] - x[nc-1]) / ((x[1] - x[0]) + (x[nc] - x[nc-1]));
      s1 = (1 - a)*get_slope(x+nc-1, y+nc-1) + a*smid;
    } else {
      const Real a = (x[i] - x[i-1]) / (x[ip1] - x[i-1]);
      s1 = (1 - a)*get_slope(x+i-1, y+i-1) + a*smid;
    }
    if (at_end) {
      const Real a = (x[ip1] - x[i]) / ((x[ip1] - x[i]) + (x[1] - x[0]));
      s2 = (1 - a)*smid + a*get_slope(x, y);
    } else {
      const Real a = (x[ip1] - x[i]) / (x[i+2] - x[i]);
      s2 = (1 - a)*smid + a*get_slope(x+ip1, y+ip1);
    }
    // Interp.
    Real c[4];
    get_cubic(x[ip1] - x[i], y[i], s1, y[ip1], s2, c);
    const Real xij = xi_per - x[i];
    yi[j] = (((c[0]*xij + c[1])*xij) + c[2])*xij + c[3];
  }
}
} // namespace interp

class PyWriter {
  typedef std::unique_ptr<FILE, cedr::util::FILECloser> FilePtr;
  FilePtr fh_;
public:
  PyWriter(const std::string& filename);
  void write(const std::string& field_name, const std::vector<Real>& v) const;
};

PyWriter::PyWriter (const std::string& filename) {
  fh_ = FilePtr(fopen((filename + ".py").c_str(), "w"));
  fprintf(fh_.get(), "s = {};\n");
}

void PyWriter::write (const std::string& field_name, const std::vector<Real>& v) const {
  fprintf(fh_.get(), "s['%s'] = [", field_name.c_str());
  for (const auto& e: v)
    fprintf(fh_.get(), " %1.15e,", e);
  fprintf(fh_.get(), "]\n");
}

struct InitialCondition {
  enum Enum { sin, bell, rect, uniform };
  static std::string convert (const Enum& e) {
    switch (e) {
    case Enum::sin: return "sin";
    case Enum::bell: return "bell";
    case Enum::rect: return "rect";
    case Enum::uniform: return "uniform";
    }
    cedr_throw_if(true, "InitialCondition::convert can't convert " << e);
  }
  static Enum convert (const std::string& s) {
    using util::eq;
    if (eq(s, "sin")) return Enum::sin;
    if (eq(s, "bell")) return Enum::bell;
    if (eq(s, "rect")) return Enum::rect;
    if (eq(s, "uniform")) return Enum::uniform;
    cedr_throw_if(true, "InitialCondition::convert can't convert " << s);
  }
  static Real eval (const Enum& ic, const Real x) {
    switch (ic) {
    case Enum::sin: return 0.1 + 0.8*0.5*(1 + std::sin(6*M_PI*x));
    case Enum::bell: return x < 0.5 ? std::sin(2*M_PI*x) : 0;
    case Enum::rect: return x > 0.66 || x < 0.33 ? 0 : 1;
    case Enum::uniform: return 0.42;
    }
    cedr_throw_if(true, "InitialCondition::eval can't convert " << ic);
  }
};

class Problem1D {
  std::vector<Real> xb_, xcp_, rwrk_;
  std::vector<Int> iwrk_;

  void init_mesh (const Int ncells, const bool nonuniform_mesh) {
    xb_.resize(ncells+1);
    xcp_.resize(ncells+1);
    xb_[0] = 0;
    if (nonuniform_mesh) {
      // Large-scale, continuous variation in cell size, plus a huge jump at the
      // periodic boundary.
      for (Int i = 1; i <= ncells; ++i) {
        const Real x = cedr::util::square(Real(i) / ncells);
        xb_[i] = 0.01 + sin(0.5*M_PI*x*x*x*x);
      }
      // Random local cell sizes.
      for (Int i = 1; i <= ncells; ++i)
        xb_[i] *= 0.3 + cedr::util::urand();
      // Cumsum.
      for (Int i = 1; i <= ncells; ++i)
        xb_[i] += xb_[i-1];
      // Normalize.
      for (Int i = 1; i <= ncells; ++i)
        xb_[i] /= xb_[ncells];
    } else {
      xb_.back() = 1;
      for (Int i = 1; i < ncells; ++i)
        xb_[i] = Real(i) / ncells;
    }
    for (Int i = 0; i < ncells; ++i)
      xcp_[i] = 0.5*(xb_[i] + xb_[i+1]);
    xcp_.back() = 1 + xcp_[0];
  }

  static void run_cdr (const Problem1D& p, CDR& cdr,
                       const Real* yp, Real* y, const Int* dods) {
    const Int n = p.ncells();
    for (Int i = 0; i < n; ++i) {
      const Int* dod = dods + 4*i;
      Real min = yp[dod[0]], max = min;
      for (Int j = 1; j < 4; ++j) {
        const Real v = yp[dod[j]];
        min = std::min(min, v);
        max = std::max(max, v);
      }
      const Real area_i = p.area(i);
      cdr.set_Qm(i, 0, y[i]*area_i, min*area_i, max*area_i, yp[i]*area_i);
    }
    cdr.run();
    for (Int i = 0; i < n; ++i)
      y[i] = cdr.get_Qm(i, 0) / p.area(i);
    y[n] = y[0];
  }

  static void run_caas (const Problem1D& p, const Real* yp, Real* y, const Int* dods) {
    const Int n = p.ncells();
    std::vector<Real> lo(n), up(n), w(n);
    Real m = 0;
    for (Int i = 0; i < n; ++i) {
      const Int* dod = dods + 4*i;
      Real min = yp[dod[0]], max = min;
      for (Int j = 1; j < 4; ++j) {
        const Real v = yp[dod[j]];
        min = std::min(min, v);
        max = std::max(max, v);
      }
      const Real area_i = p.area(i);
      lo[i] = min*area_i;
      up[i] = max*area_i;
      y[i] = std::max(min, std::min(max, y[i]));
      m += (yp[i] - y[i])*area_i;
    }
    Real wsum = 0;
    for (Int i = 0; i < n; ++i) {
      w[i] = m >= 0 ? up[i] - y[i]*p.area(i) : y[i]*p.area(i) - lo[i];
      wsum += w[i];
    }
    for (Int i = 0; i < n; ++i)
      y[i] += (m/(wsum*p.area(i)))*w[i];
  }

public:
  Problem1D (const Int ncells, const bool nonuniform_mesh = false) {
    init_mesh(ncells, nonuniform_mesh);
  }

  Int ncells () const { return xb_.size() - 1; }
  Real xb (const Int& i) const { return xb_[i]; }
  Real xcp (const Int& i) const { return xcp_[i]; }
  Real area (const Int& i) const { return xb_[i+1] - xb_[i]; }

  const std::vector<Real> get_xb () const { return xb_; }
  const std::vector<Real> get_xcp () const { return xcp_; }

  void cycle (const Int& nsteps, const Real* y0, Real* yf, CDR* cdr = nullptr) {
    const Int n = xcp_.size();
    rwrk_.resize(2*n);
    iwrk_.resize(4*n);
    Real* xcpi = rwrk_.data();
    Int* dod = iwrk_.data();

    const Real xos = -1.0 / nsteps;
    for (Int i = 0; i < n; ++i)
      xcpi[i] = xcp_[i] + xos;

    Real* ys[] = {xcpi + n, yf};
    std::copy(y0, y0 + n, ys[0]);
    for (Int ti = 0; ti < nsteps; ++ti) {
      interp::cubic_interp_periodic(xcp_.data(), n, ys[0],
                                    xcpi, n, ys[1], dod);
      if (cdr)
        run_cdr(*this, *cdr, ys[0], ys[1], dod);
      else
        run_caas(*this, ys[0], ys[1], dod);
      std::swap(ys[0], ys[1]);
    }
    std::copy(ys[0], ys[0] + n, yf);
  }
};

//todo Clean this up. Right now everything is hardcoded and kludgy.
// - optional write
// - some sort of brief quantitative output
// - better, more canonical IC
// - optional tree imbalance
// - optional mesh nonuniformity
// - parallel?
Int run (const mpi::Parallel::Ptr& parallel, const Input& in) {
  cedr_throw_if(parallel->size() > 1, "run_1d_transport_test runs in serial only.");
  Int nerr = 0;

  Problem1D p(in.ncells, false /* nonuniform_mesh */ );

  auto tree = qlt::tree::make_tree_over_1d_mesh(parallel, in.ncells,
                                                false /* imbalanced */);
  typedef qlt::QLT<Kokkos::DefaultHostExecutionSpace> QLTT;
  QLTT qltnn(parallel, in.ncells, tree), qlt(parallel, in.ncells, tree);

  typedef caas::CAAS<Kokkos::DefaultHostExecutionSpace> CAAST;
  CAAST caas(parallel, in.ncells);

  CDR* cdrs[] = {&qltnn, &qlt, &caas};
  const int ncdrs = sizeof(cdrs)/sizeof(*cdrs);

  bool first = true;
  for (CDR* cdr : cdrs) {
    if (first) {
      QLTT* qlt = dynamic_cast<QLTT*>(cdr);
      cedr_assert(qlt);
      qlt->declare_tracer(cedr::ProblemType::conserve |
                          cedr::ProblemType::nonnegative, 0);
      first = false;
    } else
      cdr->declare_tracer(cedr::ProblemType::conserve |
                          cedr::ProblemType::shapepreserve, 0);
    cdr->end_tracer_declarations();
    cdr->finish_setup();
    for (Int i = 0; i < in.ncells; ++i)
      cdr->set_rhom(i, 0, p.area(i));
    cdr->print(std::cout);
  }

  std::vector<Real> y0(in.ncells+1);
  for (Int i = 0, nc = p.ncells(); i < nc; ++i)
    y0[i] = (p.xcp(i) < 0.4 || p.xcp(i) > 0.9 ?
             InitialCondition::eval(InitialCondition::sin, p.xcp(i)) :
             InitialCondition::eval(InitialCondition::rect, p.xcp(i)));
  y0.back() = y0[0];

  PyWriter w("out_transport1d");
  w.write("xb", p.get_xb());
  w.write("xcp", p.get_xcp());
  w.write("y0", y0);

  std::vector<Real> yf(in.ncells+1);
  const Int nsteps = Int(3.17*in.ncells);
  const Int ncycles = 1;
  
  const char* names[] = {"yqltnn", "yqlt", "ycaas"};
  for (int ic = 0; ic < ncdrs; ++ic) {
    std::copy(y0.begin(), y0.end(), yf.begin());
    for (Int i = 0; i < ncycles; ++i)
      p.cycle(nsteps, yf.data(), yf.data(), cdrs[ic]);
    w.write(names[ic], yf);
  }

  std::copy(y0.begin(), y0.end(), yf.begin());
  for (Int i = 0; i < ncycles; ++i)
    p.cycle(nsteps, yf.data(), yf.data());
  w.write("ylcaas", yf);

  return nerr;
}

} // namespace transport1d
} // namespace test
} // namespace cedr

//>> cedr_test.cpp
// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

//#include "cedr_qlt.hpp"
//#include "cedr_caas.hpp"
//#include "cedr_mpi.hpp"
//#include "cedr_util.hpp"
//#include "cedr_test.hpp"

#include <stdexcept>
#include <sstream>

namespace cedr { 
struct InputParser {
  qlt::test::Input qin;
  test::transport1d::Input tin;

  class ArgAdvancer {
    const int argc_;
    char const* const* argv_;
    int i_;
  public:
    ArgAdvancer (int argc, char** argv) : argc_(argc), argv_(argv), i_(1) {}
    const char* advance () {
      if (i_+1 >= argc_) cedr_throw_if(true, "Command line is missing an argument.");
      return argv_[++i_];
    }
    const char* token () const { return argv_[i_]; }
    void incr () { ++i_; }
    bool more () const { return i_ < argc_; }
  };
  
  InputParser (int argc, char** argv, const qlt::Parallel::Ptr& p) {
    using util::eq;
    qin.unittest = false;
    qin.perftest = false;
    qin.write = false;
    qin.ncells = 0;
    qin.ntracers = 1;
    qin.tracer_type = 0;
    qin.nrepeat = 1;
    qin.pseudorandom = false;
    qin.verbose = false;
    tin.ncells = 0;
    for (ArgAdvancer aa(argc, argv); aa.more(); aa.incr()) {
      const char* token = aa.token();
      if (eq(token, "-t", "--unittest")) qin.unittest = true;
      else if (eq(token, "-pt", "--perftest")) qin.perftest = true;
      else if (eq(token, "-w", "--write")) qin.write = true;
      else if (eq(token, "-nc", "--ncells")) qin.ncells = std::atoi(aa.advance());
      else if (eq(token, "-nt", "--ntracers")) qin.ntracers = std::atoi(aa.advance());
      else if (eq(token, "-tt", "--tracertype")) qin.tracer_type = std::atoi(aa.advance());
      else if (eq(token, "-nr", "--nrepeat")) qin.nrepeat = std::atoi(aa.advance());
      else if (eq(token, "--proc-random")) qin.pseudorandom = true;
      else if (eq(token, "-v", "--verbose")) qin.verbose = true;
      else if (eq(token, "-t1d", "--transport1dtest")) tin.ncells = 1;
      else cedr_throw_if(true, "Invalid token " << token);
    }

    if (tin.ncells) {
      tin.ncells = qin.ncells;
      tin.verbose = qin.verbose;
    }

    cedr_throw_if(qin.tracer_type < 0 || qin.tracer_type >= 4,
                  "Tracer type is out of bounds [0, 3].");
    cedr_throw_if(qin.ntracers < 1, "Number of tracers is < 1.");
  }

  void print (std::ostream& os) const {
    os << "ncells " << qin.ncells
       << " nrepeat " << qin.nrepeat;
    if (qin.pseudorandom) os << " random";
    os << "\n";
  }
};
} // namespace cedr

// -------------------- Homme-specific impl details -------------------- //

#if 0
// Includes if using compose library.
#include "compose/cedr.hpp"
// Use these when rewriting each CDR's run() function to interact nicely with
// Homme's nested OpenMP and top-level horizontal threading scheme.
#include "compose/cedr_qlt.hpp"
#include "compose/cedr_caas.hpp"
#endif

#define THREAD_QLT_RUN
#ifndef QLT_MAIN
# ifdef HAVE_CONFIG_H
#  include "config.h.c"
# endif
#endif

namespace homme {
namespace compose {

template <typename ES>
class QLT : public cedr::qlt::QLT<ES> {
  typedef cedr::Int Int;
  typedef cedr::Real Real;
  typedef cedr::qlt::QLT<ES> Super;
  typedef typename Super::RealList RealList;

  //todo All of this VerticalLevelsData-related code should be impl'ed
  // using a new optional root-node-function QLT registration function.
  struct VerticalLevelsData {
    typedef std::shared_ptr<VerticalLevelsData> Ptr;

    RealList lo, hi, mass, ones, wrk;

    VerticalLevelsData (const cedr::Int n)
      : lo("lo", n), hi("hi", n), mass("mass", n), ones("ones", n), wrk("wrk", n)
    {
      for (cedr::Int k = 0; k < n; ++k) ones(k) = 1;
    }
  };

  typename VerticalLevelsData::Ptr vld_;

  void reconcile_vertical (const Int problem_type, const Int bd_os,
                           const Int bis, const Int bie) {
    using cedr::ProblemType;

    cedr_assert((problem_type & ProblemType::shapepreserve) &&
                (problem_type & ProblemType::conserve));

    auto& md = this->md_;
    auto& bd = this->bd_;
    const auto& vld = *vld_;
    const Int nlev = vld.lo.extent_int(0);
    const Int nprob = (bie - bis)/nlev;

#if defined THREAD_QLT_RUN && defined HORIZ_OPENMP
#   pragma omp master
#endif
    for (Int pi = 0; pi < nprob; ++pi) {
      const Int bd_os_pi = bd_os + md.a_d.trcr2bl2r(md.a_d.bidx2trcr(bis + pi));
#ifdef RV_DIAG
      Real oob = 0;
#endif
      Real tot_mass = 0;
      for (Int k = 0; k < nlev; ++k) {
        const Int bd_os_k = bd_os_pi + nprob*4*k;
        vld.lo  (k) = bd.l2r_data(bd_os_k    );
        vld.hi  (k) = bd.l2r_data(bd_os_k + 2);
        vld.mass(k) = bd.l2r_data(bd_os_k + 3); // previous mass, not current one
        tot_mass += vld.mass(k);
#ifdef RV_DIAG
        if (vld.mass(k) < vld.lo(k)) oob += vld.lo(k) - vld.mass(k);
        if (vld.mass(k) > vld.hi(k)) oob += vld.mass(k) - vld.hi(k);
#endif
      }
      solve(nlev, vld, tot_mass);
      Real tot_mass_slv = 0, oob_slv = 0;
      for (Int k = 0; k < nlev; ++k) {
        const Int bd_os_k = bd_os_pi + nprob*4*k;
        bd.l2r_data(bd_os_k + 3) = vld.mass(k); // previous mass, not current one
#ifdef RV_DIAG
        tot_mass_slv += vld.mass(k);
        if (vld.mass(k) < vld.lo(k)) oob_slv += vld.lo(k) - vld.mass(k);
        if (vld.mass(k) > vld.hi(k)) oob_slv += vld.mass(k) - vld.hi(k);
#endif
      }
#ifdef RV_DIAG
      printf("%2d %9.2e %9.2e %9.2e %9.2e\n", pi,
             oob/tot_mass, oob_slv/tot_mass,
             tot_mass, std::abs(tot_mass_slv - tot_mass)/tot_mass);
#endif
    }
#if defined THREAD_QLT_RUN && defined HORIZ_OPENMP
#   pragma omp barrier
#endif
  }

  static Int solve (const Int n, const Real* a, const Real b,
                    const Real* xlo, const Real* xhi,
                    Real* x, Real* wrk) {
#ifndef NDEBUG
    cedr_assert(b >= 0);
    for (Int i = 0; i < n; ++i) {
      cedr_assert(a[i] > 0);
      cedr_assert(xlo[i] >= 0);
      cedr_assert(xhi[i] >= xlo[i]);
    }
#endif
    Int status = 0;
    Real tot_lo = 0, tot_hi = 0;
    for (Int i = 0; i < n; ++i) tot_lo += a[i]*xlo[i];
    for (Int i = 0; i < n; ++i) tot_hi += a[i]*xhi[i];
    if (b < tot_lo) {
      status = -2;
      for (Int i = 0; i < n; ++i) wrk[i] = 0;
      for (Int i = 0; i < n; ++i) x[i] = xlo[i];
      // Find a new xlo >= 0 minimally far from the current one. This
      // is also the solution x.
      cedr::local::caas(n, a, b, wrk, xlo, x, x, false);
    } else if (b > tot_hi) {
      status = -1;
      const Real f = b/tot_hi;
      // a[i] divides out.
      for (Int i = 0; i < n; ++i) x[i] = f*xhi[i];
    } else {
      cedr::local::caas(n, a, b, xlo, xhi, x, x, false);
    }
    return status;
  }

  static Int solve (const Int nlev, const VerticalLevelsData& vld,
                    const Real& tot_mass) {
    solve(nlev, vld.ones.data(), tot_mass, vld.lo.data(), vld.hi.data(),
          vld.mass.data(), vld.wrk.data());    
  }

  static Int solve_unittest () {
    static const auto eps = std::numeric_limits<Real>::epsilon();
    static const Int n = 7;

    Real a[n], xlo[n], xhi[n], x[n], wrk[n];
    static const Real x0  [n] = { 1.2, 0.5,3  , 2  , 1.5, 1.8,0.2};
    static const Real dxlo[n] = {-0.1,-0.2,0.5,-1.5,-0.1,-1.1,0.1};
    static const Real dxhi[n] = { 0.1,-0.1,1  ,-0.5, 0.1,-0.2,0.5};
    for (Int i = 0; i < n; ++i) a[i] = i+1;
    for (Int i = 0; i < n; ++i) xlo[i] = x0[i] + dxlo[i];
    for (Int i = 0; i < n; ++i) xhi[i] = x0[i] + dxhi[i];
    Real b, b1;
    Int status, nerr = 0;

    const auto check_mass = [&] () {
      b1 = 0;
      for (Int i = 0; i < n; ++i) b1 += a[i]*x[i];
      if (std::abs(b1 - b) >= 10*eps*b) ++nerr;
    };

    for (Int i = 0; i < n; ++i) x[i] = x0[i];
    b = 0;
    for (Int i = 0; i < n; ++i) b += a[i]*xlo[i];
    b *= 0.9;
    status = solve(n, a, b, xlo, xhi, x, wrk);
    if (status != -2) ++nerr;
    check_mass();
    for (Int i = 0; i < n; ++i) if (x[i] > xhi[i]*(1 + 10*eps)) ++nerr;

    for (Int i = 0; i < n; ++i) x[i] = x0[i];
    b = 0;
    for (Int i = 0; i < n; ++i) b += a[i]*xhi[i];
    b *= 1.1;
    status = solve(n, a, b, xlo, xhi, x, wrk);
    if (status != -1) ++nerr;
    check_mass();
    for (Int i = 0; i < n; ++i) if (x[i] < xlo[i]*(1 - 10*eps)) ++nerr;

    for (Int i = 0; i < n; ++i) x[i] = x0[i];
    b = 0;
    for (Int i = 0; i < n; ++i) b += 0.5*a[i]*(xlo[i] + xhi[i]);
    status = solve(n, a, b, xlo, xhi, x, wrk);
    if (status != 0) ++nerr;
    check_mass();
    for (Int i = 0; i < n; ++i) if (x[i] < xlo[i]*(1 - 10*eps)) ++nerr;
    for (Int i = 0; i < n; ++i) if (x[i] > xhi[i]*(1 + 10*eps)) ++nerr;

    return nerr;
  }

public:
  QLT (const cedr::mpi::Parallel::Ptr& p, const cedr::Int& ncells,
       const cedr::qlt::tree::Node::Ptr& tree, const cedr::CDR::Options& options,
       const cedr::Int& vertical_levels)
    : cedr::qlt::QLT<ES>(p, ncells, tree, options)
  {
    if (vertical_levels)
      vld_ = std::make_shared<VerticalLevelsData>(vertical_levels);
  }

  void run () override {
    static const int mpitag = 42;
    using cedr::Int;
    using cedr::Real;
    using cedr::ProblemType;
    using cedr::qlt::impl::NodeSets;
    namespace mpi = cedr::mpi;
    auto& md_ = this->md_;
    auto& bd_ = this->bd_;
    auto& ns_ = this->ns_;
    auto& p_ = this->p_;
#if ! defined THREAD_QLT_RUN && defined HORIZ_OPENMP
#   pragma omp master
    {
#endif
      // Number of data per slot.
      const Int l2rndps = md_.a_d.prob2bl2r[md_.nprobtypes];
      const Int r2lndps = md_.a_d.prob2br2l[md_.nprobtypes];

      // Leaves to root.
      for (size_t il = 0; il < ns_->levels.size(); ++il) {
        auto& lvl = ns_->levels[il];

        // Set up receives.
        if (lvl.kids.size()) {
#if defined THREAD_QLT_RUN && defined HORIZ_OPENMP
#         pragma omp master
#endif
          {
            for (size_t i = 0; i < lvl.kids.size(); ++i) {
              const auto& mmd = lvl.kids[i];
              mpi::irecv(*p_, &bd_.l2r_data(mmd.offset*l2rndps), mmd.size*l2rndps, mmd.rank,
                         mpitag, &lvl.kids_req[i]);
            }
            mpi::waitall(lvl.kids_req.size(), lvl.kids_req.data());
          }
#if defined THREAD_QLT_RUN && defined HORIZ_OPENMP
#         pragma omp barrier
#endif
        }

        // Combine kids' data.
        if (lvl.nodes.size()) {
#if defined THREAD_QLT_RUN && defined HORIZ_OPENMP
#         pragma omp for
#endif
          for (size_t ni = 0; ni < lvl.nodes.size(); ++ni) {
            const auto lvlidx = lvl.nodes[ni];
            const auto n = ns_->node_h(lvlidx);
            if ( ! n->nkids) continue;
            cedr_kernel_assert(n->nkids == 2);
            // Total density.
            bd_.l2r_data(n->offset*l2rndps) =
              (bd_.l2r_data(ns_->node_h(n->kids[0])->offset*l2rndps) +
               bd_.l2r_data(ns_->node_h(n->kids[1])->offset*l2rndps));
            // Tracers.
            for (Int pti = 0; pti < md_.nprobtypes; ++pti) {
              const Int problem_type = md_.get_problem_type(pti);
              const bool nonnegative = problem_type & ProblemType::nonnegative;
              const bool shapepreserve = problem_type & ProblemType::shapepreserve;
              const bool conserve = problem_type & ProblemType::conserve;
              const Int bis = md_.a_d.prob2trcrptr[pti], bie = md_.a_d.prob2trcrptr[pti+1];
#if defined THREAD_QLT_RUN && defined COLUMN_OPENMP
#             pragma omp parallel for
#endif
              for (Int bi = bis; bi < bie; ++bi) {
                const Int bdi = md_.a_d.trcr2bl2r(md_.a_d.bidx2trcr(bi));
                Real* const me = &bd_.l2r_data(n->offset*l2rndps + bdi);
                const auto kid0 = ns_->node_h(n->kids[0]);
                const auto kid1 = ns_->node_h(n->kids[1]);
                const Real* const k0 = &bd_.l2r_data(kid0->offset*l2rndps + bdi);
                const Real* const k1 = &bd_.l2r_data(kid1->offset*l2rndps + bdi);
                if (nonnegative) {
                  me[0] = k0[0] + k1[0];
                  if (conserve) me[1] = k0[1] + k1[1];
                } else {
                  me[0] = shapepreserve ? k0[0] + k1[0] : cedr::impl::min(k0[0], k1[0]);
                  me[1] = k0[1] + k1[1];
                  me[2] = shapepreserve ? k0[2] + k1[2] : cedr::impl::max(k0[2], k1[2]);
                  if (conserve) me[3] = k0[3] + k1[3] ;
                }
              }
            }
          }
        }

        // Send to parents.
        if (lvl.me.size())
#if defined THREAD_QLT_RUN && defined HORIZ_OPENMP
#       pragma omp master
#endif
        {
          for (size_t i = 0; i < lvl.me.size(); ++i) {
            const auto& mmd = lvl.me[i];
            mpi::isend(*p_, &bd_.l2r_data(mmd.offset*l2rndps), mmd.size*l2rndps,
                       mmd.rank, mpitag);
          }
        }

#if defined THREAD_QLT_RUN && defined HORIZ_OPENMP
#       pragma omp barrier
#endif
      }

      // Root.
      if ( ! (ns_->levels.empty() || ns_->levels.back().nodes.size() != 1 ||
              ns_->node_h(ns_->levels.back().nodes[0])->parent >= 0)) {
        const auto n = ns_->node_h(ns_->levels.back().nodes[0]);
        for (Int pti = 0; pti < md_.nprobtypes; ++pti) {
          const Int bis = md_.a_d.prob2trcrptr[pti], bie = md_.a_d.prob2trcrptr[pti+1];
          if (bie == bis) continue;
          const Int problem_type = md_.get_problem_type(pti);
          if (vld_) reconcile_vertical(problem_type, n->offset*l2rndps, bis, bie);
#if defined THREAD_QLT_RUN && defined HORIZ_OPENMP && defined COLUMN_OPENMP
#         pragma omp parallel
#endif
#if defined THREAD_QLT_RUN && (defined HORIZ_OPENMP || defined COLUMN_OPENMP)
#         pragma omp for
#endif
          for (Int bi = bis; bi < bie; ++bi) {
            const Int l2rbdi = md_.a_d.trcr2bl2r(md_.a_d.bidx2trcr(bi));
            const Int r2lbdi = md_.a_d.trcr2br2l(md_.a_d.bidx2trcr(bi));
            // If QLT is enforcing global mass conservation, set the root's r2l Qm
            // value to the l2r Qm_prev's sum; otherwise, copy the l2r Qm value to
            // the r2l one.
            const Int os = (problem_type & ProblemType::conserve ?
                            md_.get_problem_type_l2r_bulk_size(problem_type) - 1 :
                            (problem_type & ProblemType::nonnegative ? 0 : 1));
            bd_.r2l_data(n->offset*r2lndps + r2lbdi) =
              bd_.l2r_data(n->offset*l2rndps + l2rbdi + os);
            if ((problem_type & ProblemType::consistent) &&
                ! (problem_type & ProblemType::shapepreserve)) {
              // Consistent but not shape preserving, so we're solving a dynamic range
              // preservation problem. We now know the global q_{min,max}. Start
              // propagating it leafward.
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

        if (lvl.me.size()) {
#if defined THREAD_QLT_RUN && defined HORIZ_OPENMP
#         pragma omp master
#endif
          {
            for (size_t i = 0; i < lvl.me.size(); ++i) {
              const auto& mmd = lvl.me[i];
              mpi::irecv(*p_, &bd_.r2l_data(mmd.offset*r2lndps), mmd.size*r2lndps, mmd.rank,
                         mpitag, &lvl.me_recv_req[i]);
            }
            mpi::waitall(lvl.me_recv_req.size(), lvl.me_recv_req.data());
          }
#if defined THREAD_QLT_RUN && defined HORIZ_OPENMP
#         pragma omp barrier
#endif
        }

        // Solve QP for kids' values.
#if defined THREAD_QLT_RUN && defined HORIZ_OPENMP
#       pragma omp for
#endif
        for (size_t ni = 0; ni < lvl.nodes.size(); ++ni) {
          const auto lvlidx = lvl.nodes[ni];
          const auto n = ns_->node_h(lvlidx);
          if ( ! n->nkids) continue;
          for (Int pti = 0; pti < md_.nprobtypes; ++pti) {
            const Int problem_type = md_.get_problem_type(pti);
            const Int bis = md_.a_d.prob2trcrptr[pti], bie = md_.a_d.prob2trcrptr[pti+1];
#if defined THREAD_QLT_RUN && defined COLUMN_OPENMP
#           pragma omp parallel for
#endif
            for (Int bi = bis; bi < bie; ++bi) {
              const Int l2rbdi = md_.a_d.trcr2bl2r(md_.a_d.bidx2trcr(bi));
              const Int r2lbdi = md_.a_d.trcr2br2l(md_.a_d.bidx2trcr(bi));
              cedr_assert(n->nkids == 2);
              if ((problem_type & ProblemType::consistent) &&
                  ! (problem_type & ProblemType::shapepreserve)) {
                // Pass q_{min,max} info along. l2r data are updated for use in
                // solve_node_problem. r2l data are updated for use in isend.
                const Real q_min = bd_.r2l_data(n->offset*r2lndps + r2lbdi + 1);
                const Real q_max = bd_.r2l_data(n->offset*r2lndps + r2lbdi + 2);
                bd_.l2r_data(n->offset*l2rndps + l2rbdi + 0) = q_min;
                bd_.l2r_data(n->offset*l2rndps + l2rbdi + 2) = q_max;
                for (Int k = 0; k < 2; ++k) {
                  const auto os = ns_->node_h(n->kids[k])->offset;
                  bd_.l2r_data(os*l2rndps + l2rbdi + 0) = q_min;
                  bd_.l2r_data(os*l2rndps + l2rbdi + 2) = q_max;
                  bd_.r2l_data(os*r2lndps + r2lbdi + 1) = q_min;
                  bd_.r2l_data(os*r2lndps + r2lbdi + 2) = q_max;
                }
              }
              const auto k0 = ns_->node_h(n->kids[0]);
              const auto k1 = ns_->node_h(n->kids[1]);
              cedr::qlt::impl::solve_node_problem(
                problem_type,
                 bd_.l2r_data( n->offset*l2rndps),
                &bd_.l2r_data( n->offset*l2rndps + l2rbdi),
                 bd_.r2l_data( n->offset*r2lndps + r2lbdi),
                 bd_.l2r_data(k0->offset*l2rndps),
                &bd_.l2r_data(k0->offset*l2rndps + l2rbdi),
                 bd_.r2l_data(k0->offset*r2lndps + r2lbdi),
                 bd_.l2r_data(k1->offset*l2rndps),
                &bd_.l2r_data(k1->offset*l2rndps + l2rbdi),
                 bd_.r2l_data(k1->offset*r2lndps + r2lbdi),
                this->options_.prefer_numerical_mass_conservation_to_numerical_bounds);
            }
          }
        }

        // Send.
        if (lvl.kids.size())
#if defined THREAD_QLT_RUN && defined HORIZ_OPENMP
#       pragma omp master
#endif
        {
          for (size_t i = 0; i < lvl.kids.size(); ++i) {
            const auto& mmd = lvl.kids[i];
            mpi::isend(*p_, &bd_.r2l_data(mmd.offset*r2lndps), mmd.size*r2lndps,
                       mmd.rank, mpitag);
          }
        }
      }
#if ! defined THREAD_QLT_RUN && defined HORIZ_OPENMP
    }
#endif
  }

  static Int unittest () {
    return solve_unittest();
  }
};

// We explicitly use Kokkos::Serial here so we can run the Kokkos kernels in the
// super class w/o triggering an expecution-space initialization error in
// Kokkos. This complication results from the interaction of Homme's
// HORIZ_OPENMP threading with Kokkos kernels.
struct CAAS : public cedr::caas::CAAS<Kokkos::Serial> {
  typedef cedr::caas::CAAS<Kokkos::Serial> Super;

  CAAS (const cedr::mpi::Parallel::Ptr& p, const cedr::Int nlclcells,
        const typename Super::UserAllReducer::Ptr& uar)
    : Super(p, nlclcells, uar)
  {}

  void run () override {
    Super::run();
  }
};

} // namespace compose
} // namespace homme

#ifdef QLT_MAIN
int main (int argc, char** argv) {
  int nerr = 0, retval = 0;
  MPI_Init(&argc, &argv);
  auto p = cedr::mpi::make_parallel(MPI_COMM_WORLD);
  srand(p->rank());
  Kokkos::initialize(argc, argv);
#if 0
  try
#endif
  {
    cedr::InputParser inp(argc, argv, p);
    if (p->amroot()) inp.print(std::cout);
    if (inp.qin.unittest) {
      nerr += cedr::local::unittest();
      nerr += cedr::caas::test::unittest(p);
    }
    if (inp.qin.unittest || inp.qin.perftest)
      nerr += cedr::qlt::test::run_unit_and_randomized_tests(p, inp.qin);
    if (inp.tin.ncells > 0)
      nerr += cedr::test::transport1d::run(p, inp.tin);
    {
      int gnerr;
      cedr::mpi::all_reduce(*p, &nerr, &gnerr, 1, MPI_SUM);
      retval = gnerr != 0 ? -1 : 0;
      if (p->amroot())
        std::cout << (gnerr != 0 ? "FAIL" : "PASS") << "\n";
    }
  }
#if 0
  catch (const std::exception& e) {
    if (p->amroot())
      std::cerr << e.what();
    retval = -1;
  }
#endif
  Kokkos::finalize();
  if (nerr) prc(nerr);
  MPI_Finalize();
  return retval;
}
#endif

namespace homme {
namespace qlt = cedr::qlt;
using cedr::Int;
using cedr::Real;

Int rank2sfc_search (const Int* rank2sfc, const Int& nrank, const Int& sfc) {
  Int lo = 0, hi = nrank+1;
  while (hi > lo + 1) {
    const Int mid = (lo + hi)/2;
    if (sfc >= rank2sfc[mid])
      lo = mid;
    else
      hi = mid;
  }
  return lo;
}

// Change leaf node->cellidx from index into space-filling curve to global cell
// index. owned_ids is in SFC index order for this rank.
void renumber (const Int nrank, const Int nelem, const Int my_rank, const Int* owned_ids,
               const Int* rank2sfc, const qlt::tree::Node::Ptr& node) {
  if (node->nkids) {
    for (Int k = 0; k < node->nkids; ++k)
      renumber(nrank, nelem, my_rank, owned_ids, rank2sfc, node->kids[k]);
  } else {
    const Int sfc = node->cellidx;
    node->cellidx = node->rank == my_rank ? owned_ids[sfc - rank2sfc[my_rank]] : -1;
    cedr_assert((node->rank != my_rank && node->cellidx == -1) ||
                (node->rank == my_rank && node->cellidx >= 0 && node->cellidx < nelem));
  }
}

void renumber (const Int* sc2gci, const Int* sc2rank,
               const qlt::tree::Node::Ptr& node) {
  if (node->nkids) {
    for (Int k = 0; k < node->nkids; ++k)
      renumber(sc2gci, sc2rank, node->kids[k]);
  } else {
    const Int ci = node->cellidx;
    node->cellidx = sc2gci[ci];
    node->rank = sc2rank[ci];
  }
}

// Build a subtree over [0, nsublev).
void add_sub_levels (const qlt::tree::Node::Ptr& node, const Int nsublev,
                     const Int gci, const Int my_rank, const Int rank,
                     const bool calc_level, const Int slb, const Int sle) {
  if (slb+1 == sle) {
    node->cellidx = rank == my_rank ? nsublev*gci + slb : -1;
    if (calc_level) node->level = 0;
  } else {
    node->nkids = 2;
    for (Int k = 0; k < 2; ++k) {
      auto kid = std::make_shared<qlt::tree::Node>();
      kid->parent = node.get();
      kid->rank = rank;
      node->kids[k] = kid;
    }
    const Int mid = slb + (sle - slb)/2;
    add_sub_levels(node->kids[0], nsublev, gci, my_rank, rank, calc_level, slb, mid);
    add_sub_levels(node->kids[1], nsublev, gci, my_rank, rank, calc_level, mid, sle);
    if (calc_level)
      node->level = 1 + std::max(node->kids[0]->level, node->kids[1]->level);
  }
}

// Recurse to each leaf and call add_sub_levels above.
void add_sub_levels (const Int my_rank, const qlt::tree::Node::Ptr& node,
                     const Int nsublev, const Int level_offset) {
  if (node->nkids) {
    for (Int k = 0; k < node->nkids; ++k)
      add_sub_levels(my_rank, node->kids[k], nsublev, level_offset);
    node->level += level_offset;
  } else {
    const Int gci = node->cellidx;
    const Int rank = node->rank;
    add_sub_levels(node, nsublev, gci, my_rank, rank, level_offset, 0, nsublev);
    // Level already calculated if requested.
    cedr_assert(level_offset == 0 || node->level == level_offset);
  }
}

// This impl carefully follows the requirements that
// cedr::qlt::impl::init_tree, level_schedule_and_collect
// establish. init_tree has to be modified to have the condition in
// the line
//   if (node->rank < 0) node->rank = node->kids[0]->rank
// since here we're assigning the ranks ourselves. Similarly, it must
// check for node->level >= 0; if the tree is partial, it is unable to
// compute node level.
qlt::tree::Node::Ptr
make_my_tree_part (const qlt::oned::Mesh& m, const Int cs, const Int ce,
                   const qlt::tree::Node* parent,
                   const Int& nrank, const Int* rank2sfc) {
  const auto my_rank = m.parallel()->rank();
  const Int cn = ce - cs, cn0 = cn/2;
  qlt::tree::Node::Ptr n = std::make_shared<qlt::tree::Node>();
  n->parent = parent;
  n->rank = rank2sfc_search(rank2sfc, nrank, cs);
  n->cellidx = n->rank == my_rank ? cs : -1;
  cedr_assert(n->rank >= 0 && n->rank < nrank);
  if (cn == 1) {
    n->nkids = 0;
    n->level = 0;
    return n;
  }
  const auto k1 = make_my_tree_part(m, cs, cs + cn0, n.get(), nrank, rank2sfc);
  const auto k2 = make_my_tree_part(m, cs + cn0, ce, n.get(), nrank, rank2sfc);
  n->level = 1 + std::max(k1->level, k2->level);
  if (n->rank == my_rank) {
    // Need to know both kids for comm.
    n->nkids = 2;
    n->kids[0] = k1;
    n->kids[1] = k2;
  } else {
    // Prune parts of the tree irrelevant to my rank.
    n->nkids = 0;
    if (k1->nkids > 0 || k1->rank == my_rank) n->kids[n->nkids++] = k1;
    if (k2->nkids > 0 || k2->rank == my_rank) n->kids[n->nkids++] = k2;
    if (n->nkids == 0) {
      // Signal a non-leaf node with 0 kids to init_tree.
      n->nkids = -1;
    }
  }
  cedr_assert(n->level > 0 || n->nkids == 0);
  return n;
}

qlt::tree::Node::Ptr
make_my_tree_part (const cedr::mpi::Parallel::Ptr& p, const Int& ncells,
                   const Int& nrank, const Int* rank2sfc) {
  qlt::oned::Mesh m(ncells, p);
  return make_my_tree_part(m, 0, m.ncell(), nullptr, nrank, rank2sfc);
}

static size_t nextpow2 (size_t n) {
  size_t p = 1;
  while (p < n) p <<= 1;
  return p;
}

static size_t get_tree_height (size_t nleaf) {
  size_t height = 0;
  nleaf = nextpow2(nleaf);
  while (nleaf) {
    ++height;
    nleaf >>= 1;
  }
  return height;
}

qlt::tree::Node::Ptr
make_tree_sgi (const cedr::mpi::Parallel::Ptr& p, const Int nelem,
               const Int* owned_ids, const Int* rank2sfc, const Int nsublev) {
  // Partition 0:nelem-1, the space-filling curve space.
  auto tree = make_my_tree_part(p, nelem, p->size(), rank2sfc);
  // Renumber so that node->cellidx records the global element number, and
  // associate the correct rank with the element.
  const auto my_rank = p->rank();
  renumber(p->size(), nelem, my_rank, owned_ids, rank2sfc, tree);
  if (nsublev > 1) {
    const Int level_offset = get_tree_height(nsublev) - 1;
    add_sub_levels(my_rank, tree, nsublev, level_offset);
  }
  return tree;
}

qlt::tree::Node::Ptr
make_tree_non_sgi (const cedr::mpi::Parallel::Ptr& p, const Int nelem,
                   const Int* sc2gci, const Int* sc2rank, const Int nsublev) {
  auto tree = qlt::tree::make_tree_over_1d_mesh(p, nelem);
  renumber(sc2gci, sc2rank, tree);
  const auto my_rank = p->rank();
  if (nsublev > 1) add_sub_levels(my_rank, tree, nsublev, 0);
  return tree;
}

qlt::tree::Node::Ptr
clone (const qlt::tree::Node::Ptr& in, const qlt::tree::Node* parent = nullptr) {
  const auto out = std::make_shared<qlt::tree::Node>(*in);
  cedr_assert(out->rank == in->rank && out->level == in->level &&
              out->nkids == in->nkids && out->cellidx == in->cellidx);
  out->parent = parent;
  for (Int k = 0; k < in->nkids; ++k)
    out->kids[k] = clone(in->kids[k], out.get());
  return out;
}

void renumber_leaves (const qlt::tree::Node::Ptr& node, const Int horiz_nleaf,
                      const Int supidx) {
  if (node->nkids) {
    for (Int k = 0; k < node->nkids; ++k)
      renumber_leaves(node->kids[k], horiz_nleaf, supidx);
  } else {
    if (node->cellidx != -1) {
      cedr_assert(node->cellidx >= 0 && node->cellidx < horiz_nleaf);
      node->cellidx += horiz_nleaf*supidx;
    }
  }
}

void attach_and_renumber_horizontal_trees (const qlt::tree::Node::Ptr& supnode,
                                           const qlt::tree::Node::Ptr& htree,
                                           const Int horiz_nleaf) {
  Int level = -1, rank;
  for (Int k = 0; k < supnode->nkids; ++k) {
    auto& kid = supnode->kids[k];
    if (kid->nkids) {
      attach_and_renumber_horizontal_trees(kid, htree, horiz_nleaf);
    } else {
      const auto supidx = kid->cellidx;
      supnode->kids[k] = clone(htree);
      kid = supnode->kids[k];
      kid->parent = supnode.get();
      kid->cellidx = -1;
      renumber_leaves(kid, horiz_nleaf, supidx);
    }
    rank = kid->rank;
    level = std::max(level, kid->level);
  }
  if (level != -1) ++level;
  supnode->level = level;
  supnode->rank = rank;
}

qlt::tree::Node::Ptr
make_tree_over_index_range (const Int cs, const Int ce,
                            const qlt::tree::Node* parent = nullptr) {
  const Int cn = ce - cs, cn0 = cn/2;
  const auto n = std::make_shared<qlt::tree::Node>();
  n->parent = parent;
  if (cn == 1) {
    n->nkids = 0;
    n->cellidx = cs;
  } else {
    n->nkids = 2;
    n->kids[0] = make_tree_over_index_range(cs, cs + cn0, n.get());
    n->kids[1] = make_tree_over_index_range(cs + cn0, ce, n.get());
  }
  return n;
}

qlt::tree::Node::Ptr
combine_superlevels(const qlt::tree::Node::Ptr& horiz_tree, const Int horiz_nleaf,
                    const Int nsuplev) {
  cedr_assert(horiz_tree->nkids > 0);
  // In this tree, cellidx 0 is the top super level.
  const auto suptree = make_tree_over_index_range(0, nsuplev);
  attach_and_renumber_horizontal_trees(suptree, horiz_tree, horiz_nleaf);
  return suptree;
}

void check_tree (const cedr::mpi::Parallel::Ptr& p, const qlt::tree::Node::Ptr& n,
                 const Int nleaf) {
#ifndef NDEBUG
  cedr_assert(n->nkids >= -1 && n->nkids <= 2);
  cedr_assert(n->rank >= 0);
  cedr_assert(n->reserved == -1);
  if (n->nkids == 2)
    cedr_assert(n->level == 1 + std::max(n->kids[0]->level, n->kids[1]->level));
  if (n->nkids == 1) cedr_assert(n->level >= 1 + n->kids[0]->level);
  if (n->nkids == 0) cedr_assert(n->level == 0);
  if (n->rank != p->rank()) cedr_assert(n->cellidx == -1);
  else cedr_assert(n->cellidx < nleaf);
  for (Int k = 0; k < n->nkids; ++k) {
    cedr_assert(n.get() == n->kids[k]->parent);
    check_tree(p, n->kids[k], nleaf);
  }
#endif
}

qlt::tree::Node::Ptr
make_tree (const cedr::mpi::Parallel::Ptr& p, const Int nelem,
           const Int* gid_data, const Int* rank_data, const Int nsublev,
           const bool use_sgi, const bool cdr_over_super_levels,
           const Int nsuplev) {
  auto tree = use_sgi ?
    make_tree_sgi    (p, nelem, gid_data, rank_data, nsublev) :
    make_tree_non_sgi(p, nelem, gid_data, rank_data, nsublev);
  Int nleaf = nelem*nsublev;
  if (cdr_over_super_levels) {
    tree = combine_superlevels(tree, nleaf, nsuplev);
    nleaf *= nsuplev;
  }
  if (use_sgi) check_tree(p, tree, nleaf);
  return tree;
}

Int test_tree_maker () {
  Int nerr = 0;
  if (nextpow2(3) != 4) ++nerr;
  if (nextpow2(4) != 4) ++nerr;
  if (nextpow2(5) != 8) ++nerr;
  if (get_tree_height(3) != 3) ++nerr;
  if (get_tree_height(4) != 3) ++nerr;
  if (get_tree_height(5) != 4) ++nerr;
  if (get_tree_height(8) != 4) ++nerr;
  return nerr;
}

extern "C"
void compose_repro_sum(const Real* send, Real* recv,
                       Int nlocal, Int nfld, Int fcomm);

struct ReproSumReducer :
    public compose::CAAS::UserAllReducer {
  ReproSumReducer (Int fcomm, Int n_accum_in_place)
    : fcomm_(fcomm), n_accum_in_place_(n_accum_in_place)
  {}

  int n_accum_in_place () const override { return n_accum_in_place_; }

  int operator() (const cedr::mpi::Parallel& p, Real* sendbuf, Real* rcvbuf,
                  int nlocal, int count, MPI_Op op) const override {
    cedr_assert(op == MPI_SUM);
#ifdef HORIZ_OPENMP
#   pragma omp barrier
#   pragma omp master
#endif
    compose_repro_sum(sendbuf, rcvbuf, nlocal, count, fcomm_);
#ifdef HORIZ_OPENMP
#   pragma omp barrier
#endif
    return 0;
  }

private:
  const Int fcomm_, n_accum_in_place_;
};

struct CDR {
  typedef std::shared_ptr<CDR> Ptr;
  typedef compose::QLT<Kokkos::DefaultExecutionSpace> QLTT;
  typedef compose::CAAS CAAST;

  struct Alg {
    enum Enum { qlt, qlt_super_level, qlt_super_level_local_caas, caas, caas_super_level };
    static Enum convert (Int cdr_alg) {
      switch (cdr_alg) {
      case 2:  return qlt;
      case 20: return qlt_super_level;
      case 21: return qlt_super_level_local_caas;
      case 3:  return caas;
      case 30: return caas_super_level;
      case 42: return caas_super_level; // actually none
      default: cedr_throw_if(true,  "cdr_alg " << cdr_alg << " is invalid.");
      }
    }
    static bool is_qlt (Enum e) {
      return (e == qlt || e == qlt_super_level ||
              e == qlt_super_level_local_caas);
    }
    static bool is_caas (Enum e) {
      return e == caas || e == caas_super_level;
    }
    static bool is_suplev (Enum e) {
      return (e == qlt_super_level || e == caas_super_level ||
              e == qlt_super_level_local_caas);
    }
  };

  enum { nsublev_per_suplev = 8 };
  
  const Alg::Enum alg;
  const Int ncell, nlclcell, nlev, nsublev, nsuplev;
  const bool threed, cdr_over_super_levels, caas_in_suplev, hard_zero;
  const cedr::mpi::Parallel::Ptr p;
  qlt::tree::Node::Ptr tree; // Don't need this except for unit testing.
  cedr::CDR::Ptr cdr;
  std::vector<Int> ie2gci; // Map Homme ie to Homme global cell index.
  std::vector<Int> ie2lci; // Map Homme ie to CDR local cell index (lclcellidx).
  std::vector<char> nonneg;

  CDR (Int cdr_alg_, Int ngblcell_, Int nlclcell_, Int nlev_, bool use_sgi,
       bool independent_time_steps, const bool hard_zero_, const Int* gid_data,
       const Int* rank_data, const cedr::mpi::Parallel::Ptr& p_, Int fcomm)
    : alg(Alg::convert(cdr_alg_)),
      ncell(ngblcell_), nlclcell(nlclcell_), nlev(nlev_),
      nsublev(Alg::is_suplev(alg) ? nsublev_per_suplev : 1),
      nsuplev((nlev + nsublev - 1) / nsublev),
      threed(independent_time_steps),
      cdr_over_super_levels(threed && Alg::is_caas(alg)),
      caas_in_suplev(alg == Alg::qlt_super_level_local_caas && nsublev > 1),
      hard_zero(hard_zero_),
      p(p_), inited_tracers_(false)
  {
    const Int n_id_in_suplev = caas_in_suplev ? 1 : nsublev;
    if (Alg::is_qlt(alg)) {
      tree = make_tree(p, ncell, gid_data, rank_data, n_id_in_suplev, use_sgi,
                       cdr_over_super_levels, nsuplev);
      cedr::CDR::Options options;
      options.prefer_numerical_mass_conservation_to_numerical_bounds = true;
      Int nleaf = ncell*n_id_in_suplev;
      if (cdr_over_super_levels) nleaf *= nsuplev;
      cdr = std::make_shared<QLTT>(p, nleaf, tree, options,
                                   threed ? nsuplev : 0);
      tree = nullptr;
    } else if (Alg::is_caas(alg)) {
      const Int n_accum_in_place = n_id_in_suplev*(cdr_over_super_levels ?
                                                   nsuplev : 1);
      const auto caas = std::make_shared<CAAST>(
        p, nlclcell*n_accum_in_place,
        std::make_shared<ReproSumReducer>(fcomm, n_accum_in_place));
      cdr = caas;
    } else {
      cedr_throw_if(true, "Invalid semi_lagrange_cdr_alg " << alg);
    }
    ie2gci.resize(nlclcell);
  }

  void init_tracers (const Int qsize, const bool need_conservation) {
    nonneg.resize(qsize, hard_zero);
    typedef cedr::ProblemType PT;
    const Int nt = cdr_over_super_levels ? qsize : nsuplev*qsize;
    for (Int ti = 0; ti < nt; ++ti)
      cdr->declare_tracer(PT::shapepreserve |
                          (need_conservation ? PT::conserve : 0), 0);
    cdr->end_tracer_declarations();
  }

  void get_buffers_sizes (size_t& s1, size_t &s2) {
    cdr->get_buffers_sizes(s1, s2);
  }

  void set_buffers (Real* b1, Real* b2) {
    cdr->set_buffers(b1, b2);
    cdr->finish_setup();
  }

private:
  bool inited_tracers_;
};

void set_ie2gci (CDR& q, const Int ie, const Int gci) { q.ie2gci[ie] = gci; }

void init_ie2lci (CDR& q) {
  const Int n_id_in_suplev = q.caas_in_suplev ? 1 : q.nsublev;
  const Int nleaf =
    n_id_in_suplev*
    q.ie2gci.size()*
    (q.cdr_over_super_levels ? q.nsuplev : 1);
  q.ie2lci.resize(nleaf);
  if (CDR::Alg::is_qlt(q.alg)) {
    auto qlt = std::static_pointer_cast<CDR::QLTT>(q.cdr);
    if (q.cdr_over_super_levels) {
      const auto nlevwrem = q.nsuplev*n_id_in_suplev;
      for (size_t ie = 0; ie < q.ie2gci.size(); ++ie)
        for (Int spli = 0; spli < q.nsuplev; ++spli)
          for (Int sbli = 0; sbli < n_id_in_suplev; ++sbli)
            //       local indexing is fastest over the whole column
            q.ie2lci[nlevwrem*ie + n_id_in_suplev*spli + sbli] =
              //           but global indexing is organized according to the tree
              qlt->gci2lci(n_id_in_suplev*(q.ncell*spli + q.ie2gci[ie]) + sbli);
    } else {
      for (size_t ie = 0; ie < q.ie2gci.size(); ++ie)
        for (Int sbli = 0; sbli < n_id_in_suplev; ++sbli)
          q.ie2lci[n_id_in_suplev*ie + sbli] =
            qlt->gci2lci(n_id_in_suplev*q.ie2gci[ie] + sbli);
    }
  } else {
    if (q.cdr_over_super_levels) {
      const auto nlevwrem = q.nsuplev*n_id_in_suplev;
      for (size_t ie = 0; ie < q.ie2gci.size(); ++ie)
        for (Int spli = 0; spli < q.nsuplev; ++spli)
          for (Int sbli = 0; sbli < n_id_in_suplev; ++sbli) {
            const Int id = nlevwrem*ie + n_id_in_suplev*spli + sbli;
            q.ie2lci[id] = id;
          }
    } else {
      for (size_t ie = 0; ie < q.ie2gci.size(); ++ie)
        for (Int sbli = 0; sbli < n_id_in_suplev; ++sbli) {
          const Int id = n_id_in_suplev*ie + sbli;
          q.ie2lci[id] = id;
        }
    }
  }
}

void init_tracers (CDR& q, const Int nlev, const Int qsize,
                   const bool need_conservation) {
  q.init_tracers(qsize, need_conservation);
}

namespace sl { // For sl_advection.F90
// Fortran array wrappers.
template <typename T> using FA1 =
  Kokkos::View<T*,     Kokkos::LayoutLeft, Kokkos::HostSpace>;
template <typename T> using FA2 =
  Kokkos::View<T**,    Kokkos::LayoutLeft, Kokkos::HostSpace>;
template <typename T> using FA3 =
  Kokkos::View<T***,   Kokkos::LayoutLeft, Kokkos::HostSpace>;
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
  const Real* dp0;

  struct Check {
    Kokkos::View<Real**, Kokkos::Serial>
      mass_p, mass_c, mass_lo, mass_hi,
      q_lo, q_hi, q_min_l, q_max_l, qd_lo, qd_hi;
    Check (const Int nlev, const Int qsize)
      : mass_p("mass_p", nlev, qsize), mass_c("mass_c", nlev, qsize),
        mass_lo("mass_lo", nlev, qsize), mass_hi("mass_hi", nlev, qsize),
        q_lo("q_lo", nlev, qsize), q_hi("q_hi", nlev, qsize),
        q_min_l("q_min_l", nlev, qsize), q_max_l("q_max_l", nlev, qsize),
        qd_lo("qd_lo", nlev, qsize), qd_hi("qd_hi", nlev, qsize)
    {}
  };
  std::shared_ptr<Check> check;

  Data (Int lcl_ncell, Int np_, Int nlev_, Int qsize_, Int qsize_d_, Int timelevels_)
    : np(np_), nlev(nlev_), qsize(qsize_), qsize_d(qsize_d_), timelevels(timelevels_),
      spheremp(lcl_ncell, nullptr), dp3d_c(lcl_ncell, nullptr), q_c(lcl_ncell, nullptr),
      qdp_pc(lcl_ncell, nullptr)
  {}
};

static void check (const CDR& q, const Data& d) {
  cedr_assert(q.nlclcell == static_cast<Int>(d.spheremp.size()));
}

template <typename T>
void insert (std::vector<T*>& r, const Int i, T* v) {
  cedr_assert(i >= 0 && i < static_cast<int>(r.size()));
  r[i] = v;
}

void insert (const Data::Ptr& d, const Int ie, const Int ptridx, Real* array,
             const Int i0 = 0, const Int i1 = 0) {
  cedr_assert(d);
  switch (ptridx) {
  case 0: insert<const double>(d->spheremp, ie, array); break;
  case 1: insert<      double>(d->qdp_pc,   ie, array); d->n0_qdp = i0; d->n1_qdp = i1; break;
  case 2: insert<const double>(d->dp3d_c,   ie, array); d->tl_np1 = i0; break;
  case 3: insert<      double>(d->q_c,      ie, array); break;
  case 4: d->dp0 = array; break;
  default: cedr_throw_if(true, "Invalid pointer index " << ptridx);
  }
}

static void run_cdr (CDR& q) {
#ifdef HORIZ_OPENMP
# pragma omp barrier
#endif
  q.cdr->run();
#ifdef HORIZ_OPENMP
# pragma omp barrier
#endif
}

template <int np_>
void accum_values (const Int ie, const Int k, const Int q, const Int tl_np1,
                   const Int n0_qdp, const Int np, const bool nonneg,
                   const FA1<const Real>& spheremp, const FA3<const Real>& dp3d_c,
                   const FA4<Real>& q_min, const FA4<const Real>& q_max,
                   const FA4<const Real>& qdp_p, const FA3<const Real>& q_c,
                   Real& volume, Real& rhom, Real& Qm, Real& Qm_prev,
                   Real& Qm_min, Real& Qm_max) {
  cedr_assert(np == np_);
  static const Int np2 = np_*np_;
  for (Int g = 0; g < np2; ++g) {
    volume += spheremp(g); // * dp0[k];
    const Real rhomij = dp3d_c(g,k,tl_np1) * spheremp(g);
    rhom += rhomij;
    Qm += q_c(g,k,q) * rhomij;
    if (nonneg) q_min(g,k,q,ie) = std::max<Real>(q_min(g,k,q,ie), 0);
    Qm_min += q_min(g,k,q,ie) * rhomij;
    Qm_max += q_max(g,k,q,ie) * rhomij;
    Qm_prev += qdp_p(g,k,q,n0_qdp) * spheremp(g);
  }
}

template <int np_>
void run_global (CDR& cdr, const Data& d, Real* q_min_r, const Real* q_max_r,
                 const Int nets, const Int nete) {
  static const Int np2 = np_*np_;
  const Int np = d.np, nlev = d.nlev, qsize = d.qsize,
    nlevwrem = cdr.nsuplev*cdr.nsublev;
  cedr_assert(np == np_);
  
  FA4<      Real> q_min(q_min_r, np2, nlev, qsize, nete+1);
  FA4<const Real> q_max(q_max_r, np2, nlev, qsize, nete+1);

  for (Int ie = nets; ie <= nete; ++ie) {
    FA1<const Real> spheremp(d.spheremp[ie], np2);
    FA4<const Real> qdp_p(d.qdp_pc[ie], np2, nlev, d.qsize_d, 2);
    FA3<const Real> dp3d_c(d.dp3d_c[ie], np2, nlev, d.timelevels);
    FA3<const Real> q_c(d.q_c[ie], np2, nlev, d.qsize_d);
#ifdef COLUMN_OPENMP
#   pragma omp parallel for
#endif
    for (Int q = 0; q < qsize; ++q) {
      for (Int spli = 0; spli < cdr.nsuplev; ++spli) {
        const Int k0 = cdr.nsublev*spli;
        const bool nonneg = cdr.nonneg[q];
        const Int ti = cdr.cdr_over_super_levels ? q : spli*qsize + q;
        Real Qm = 0, Qm_min = 0, Qm_max = 0, Qm_prev = 0, rhom = 0, volume = 0;
        Int ie_idx;
        if (cdr.caas_in_suplev)
          ie_idx = cdr.cdr_over_super_levels ?
            cdr.nsuplev*ie + spli :
            ie;
        for (Int sbli = 0; sbli < cdr.nsublev; ++sbli) {
          const auto k = k0 + sbli;
          if ( ! cdr.caas_in_suplev)
            ie_idx = cdr.cdr_over_super_levels ?
              nlevwrem*ie + k :
              cdr.nsublev*ie + sbli;
          const auto lci = cdr.ie2lci[ie_idx];
          if ( ! cdr.caas_in_suplev) {
            Qm = 0; Qm_min = 0; Qm_max = 0; Qm_prev = 0;
            rhom = 0;
            volume = 0;
          }
          if (k < nlev)
            accum_values<np_>(ie, k, q, d.tl_np1, d.n0_qdp, np, nonneg,
			      spheremp, dp3d_c, q_min, q_max, qdp_p, q_c,
			      volume, rhom, Qm, Qm_prev, Qm_min, Qm_max);
          const bool write = ! cdr.caas_in_suplev || sbli == cdr.nsublev-1;
          if (write) {
            // For now, handle just one rhom. For feasible global problems,
            // it's used only as a weight vector in QLT, so it's fine. In fact,
            // use just the cell geometry, rather than total density, since in QLT
            // this field is used as a weight vector.
            //todo Generalize to one rhom field per level. Until then, we're not
            // getting QLT's safety benefit.
            if (ti == 0) cdr.cdr->set_rhom(lci, 0, volume);
            cdr.cdr->set_Qm(lci, ti, Qm, Qm_min, Qm_max, Qm_prev);
            if (Qm_prev < -0.5) {
              static bool first = true;
              if (first) {
                first = false;
                std::stringstream ss;
                ss << "Qm_prev < -0.5: Qm_prev = " << Qm_prev
                   << " on rank " << cdr.p->rank()
                   << " at (ie,gid,spli,k0,q,ti,sbli,lci,k,n0_qdp,tl_np1) = ("
                   << ie << "," << cdr.ie2gci[ie] << "," << spli << "," << k0 << ","
                   << q << "," << ti << "," << sbli << "," << lci << "," << k << ","
                   << d.n0_qdp << "," << d.tl_np1 << ")\n";
                ss << "Qdp(:,:,k,q,n0_qdp) = [";
                for (Int g = 0; g < np2; ++g)
                  ss << " " << qdp_p(g,k,q,d.n0_qdp);
                ss << "]\n";
                ss << "dp3d(:,:,k,tl_np1) = [";
                for (Int g = 0; g < np2; ++g)
                  ss << " " << dp3d_c(g,k,d.tl_np1);
                ss << "]\n";
                pr(ss.str());
              }
            }
          }
        }
      }
    }
  }

  run_cdr(cdr);
}

template <int np_>
void solve_local (const Int ie, const Int k, const Int q,
                  const Int tl_np1, const Int n1_qdp, const Int np, 
                  const bool scalar_bounds, const Int limiter_option,
                  const FA1<const Real>& spheremp, const FA3<const Real>& dp3d_c,
                  const FA4<const Real>& q_min, const FA4<const Real>& q_max,
                  const Real Qm, FA4<Real>& qdp_c, FA3<Real>& q_c) {
  cedr_assert(np == np_);
  static const Int np2 = np_*np_;

  Real wa[np2], qlo[np2], qhi[np2], y[np2], x[np2];
  Real rhom = 0;
  for (Int g = 0; g < np2; ++g) {
    const Real rhomij = dp3d_c(g,k,tl_np1) * spheremp(g);
    rhom += rhomij;
    wa[g] = rhomij;
    y[g] = q_c(g,k,q);
    x[g] = y[g];
  }

  //todo Replace with ReconstructSafely.
  if (scalar_bounds) {
    qlo[0] = q_min(0,k,q,ie);
    qhi[0] = q_max(0,k,q,ie);
    for (Int i = 1; i < np2; ++i) qlo[i] = qlo[0];
    for (Int i = 1; i < np2; ++i) qhi[i] = qhi[0];
    // We can use either 2-norm minimization or ClipAndAssuredSum for
    // the local filter. CAAS is the faster. It corresponds to limiter
    // = 0. 2-norm minimization is the same in spirit as limiter = 8,
    // but it assuredly achieves the first-order optimality conditions
    // whereas limiter 8 does not.
    if (limiter_option == 8)
      cedr::local::solve_1eq_bc_qp(np2, wa, wa, Qm, qlo, qhi, y, x);
    else {
      // We need to use *some* limiter; if 8 isn't chosen, default to
      // CAAS.
      cedr::local::caas(np2, wa, Qm, qlo, qhi, y, x);
    }
  } else {
    for (Int g = 0; g < np2; ++g) {
      qlo[g] = q_min(g,k,q,ie);
      qhi[g] = q_max(g,k,q,ie);
    }
    for (Int trial = 0; trial < 3; ++trial) {
      int info;
      if (limiter_option == 8) {
        info = cedr::local::solve_1eq_bc_qp(
          np2, wa, wa, Qm, qlo, qhi, y, x);
        if (info == 1) info = 0;
      } else {
        info = 0;
        cedr::local::caas(np2, wa, Qm, qlo, qhi, y, x, false /* clip */);
        // Clip for numerics against the cell extrema.
        Real qlo_s = qlo[0], qhi_s = qhi[0];
        for (Int i = 1; i < np2; ++i) {
          qlo_s = std::min(qlo_s, qlo[i]);
          qhi_s = std::max(qhi_s, qhi[i]);
        }
        for (Int i = 0; i < np2; ++i)
          x[i] = cedr::impl::max(qlo_s, cedr::impl::min(qhi_s, x[i]));
      }
      if (info == 0 || trial == 1) break;
      switch (trial) {
      case 0: {
        Real qlo_s = qlo[0], qhi_s = qhi[0];
        for (Int i = 1; i < np2; ++i) {
          qlo_s = std::min(qlo_s, qlo[i]);
          qhi_s = std::max(qhi_s, qhi[i]);
        }
        for (Int i = 0; i < np2; ++i) qlo[i] = qlo_s;
        for (Int i = 0; i < np2; ++i) qhi[i] = qhi_s;
      } break;
      case 1: {
        const Real q = Qm / rhom;
        for (Int i = 0; i < np2; ++i) qlo[i] = std::min(qlo[i], q);
        for (Int i = 0; i < np2; ++i) qhi[i] = std::max(qhi[i], q);                
      } break;
      }
    }
  }
        
  for (Int g = 0; g < np2; ++g) {
    q_c(g,k,q) = x[g];
    qdp_c(g,k,q,n1_qdp) = q_c(g,k,q) * dp3d_c(g,k,tl_np1);
  }
}

Int vertical_caas_backup (const Int n, Real* rhom,
                          const Real q_min, const Real q_max,
                          Real Qmlo_tot, Real Qmhi_tot, const Real Qm_tot,
                          Real* Qmlo, Real* Qmhi, Real* Qm) {
  Int status = 0;
  if (Qm_tot < Qmlo_tot || Qm_tot > Qmhi_tot) {
    if (Qm_tot < Qmlo_tot) {
      status = -2;
      for (Int i = 0; i < n; ++i) Qmhi[i] = Qmlo[i];
      for (Int i = 0; i < n; ++i) Qmlo[i] = q_min*rhom[i];
      Qmlo_tot = 0;
      for (Int i = 0; i < n; ++i) Qmlo_tot += Qmlo[i];
      if (Qm_tot < Qmlo_tot) status = -4;
    } else {
      status = -1;
      for (Int i = 0; i < n; ++i) Qmlo[i] = Qmhi[i];
      for (Int i = 0; i < n; ++i) Qmhi[i] = q_max*rhom[i];
      Qmhi_tot = 0;
      for (Int i = 0; i < n; ++i) Qmhi_tot += Qmhi[i];
      if (Qm_tot > Qmhi_tot) status = -3;
    }
    if (status < -2) {
      Real rhom_tot = 0;
      for (Int i = 0; i < n; ++i) rhom_tot += rhom[i];
      const Real q = Qm_tot/rhom_tot;
      for (Int i = 0; i < n; ++i) Qm[i] = q*rhom[i];
      return status;
    }
  }
  for (Int i = 0; i < n; ++i) rhom[i] = 1;
  cedr::local::caas(n, rhom, Qm_tot, Qmlo, Qmhi, Qm, Qm, false);
  return status;
}

void accum_values (const Int ie, const Int k, const Int q, const Int tl_np1,
                   const Int n0_qdp, const Int np2,
                   const FA1<const Real>& spheremp, const FA3<const Real>& dp3d_c,
                   const FA4<Real>& q_min, const FA4<const Real>& q_max,
                   const FA4<const Real>& qdp_p, const FA3<const Real>& q_c,
                   Real& rhom, Real& Qm, Real& Qm_min, Real& Qm_max) {
  for (Int g = 0; g < np2; ++g) {
    const Real rhomij = dp3d_c(g,k,tl_np1) * spheremp(g);
    rhom += rhomij;
    Qm += q_c(g,k,q) * rhomij;
    Qm_min += q_min(g,k,q,ie) * rhomij;
    Qm_max += q_max(g,k,q,ie) * rhomij;
  }
}

template <int np_>
void run_local (CDR& cdr, const Data& d, Real* q_min_r, const Real* q_max_r,
                const Int nets, const Int nete, const bool scalar_bounds,
                const Int limiter_option) {
  static const Int np2 = np_*np_;
  const Int np = d.np, nlev = d.nlev, qsize = d.qsize,
    nlevwrem = cdr.nsuplev*cdr.nsublev;
  cedr_assert(np == np_);

  FA4<      Real> q_min(q_min_r, np2, nlev, qsize, nete+1);
  FA4<const Real> q_max(q_max_r, np2, nlev, qsize, nete+1);

  for (Int ie = nets; ie <= nete; ++ie) {
    FA1<const Real> spheremp(d.spheremp[ie], np2);
    FA4<      Real> qdp_c(d.qdp_pc[ie], np2, nlev, d.qsize_d, 2);
    FA3<const Real> dp3d_c(d.dp3d_c[ie], np2, nlev, d.timelevels);
    FA3<      Real> q_c(d.q_c[ie], np2, nlev, d.qsize_d);
#ifdef COLUMN_OPENMP
#   pragma omp parallel for
#endif
    for (Int q = 0; q < qsize; ++q) {
      for (Int spli = 0; spli < cdr.nsuplev; ++spli) {
        const Int k0 = cdr.nsublev*spli;
        const Int ti = cdr.cdr_over_super_levels ? q : spli*qsize + q;
        if (cdr.caas_in_suplev) {
          const auto ie_idx = cdr.cdr_over_super_levels ?
            cdr.nsuplev*ie + spli :
            ie;
          const auto lci = cdr.ie2lci[ie_idx];
          const Real Qm_tot = cdr.cdr->get_Qm(lci, ti);
          Real Qm_min_tot = 0, Qm_max_tot = 0;
          Real rhom[CDR::nsublev_per_suplev], Qm[CDR::nsublev_per_suplev],
            Qm_min[CDR::nsublev_per_suplev], Qm_max[CDR::nsublev_per_suplev];
          // Redistribute mass in the vertical direction of the super level.
          Int n = cdr.nsublev;
          for (Int sbli = 0; sbli < cdr.nsublev; ++sbli) {
            const Int k = k0 + sbli;
            if (k >= nlev) {
              n = sbli;
              break;
            }
            rhom[sbli] = 0; Qm[sbli] = 0; Qm_min[sbli] = 0; Qm_max[sbli] = 0;
            accum_values(ie, k, q, d.tl_np1, d.n0_qdp, np2,
                         spheremp, dp3d_c, q_min, q_max, qdp_c, q_c,
                         rhom[sbli], Qm[sbli], Qm_min[sbli], Qm_max[sbli]);
            Qm_min_tot += Qm_min[sbli];
            Qm_max_tot += Qm_max[sbli];
          }
          if (Qm_tot >= Qm_min_tot && Qm_tot <= Qm_max_tot) {
            for (Int i = 0; i < n; ++i) rhom[i] = 1;
            cedr::local::caas(n, rhom, Qm_tot, Qm_min, Qm_max, Qm, Qm, false);
          } else {
            Real q_min_s, q_max_s;
            bool first = true;
            for (Int sbli = 0; sbli < n; ++sbli) {
              const Int k = k0 + sbli;
              for (Int g = 0; g < np2; ++g) {
                if (first) {
                  q_min_s = q_min(g,k,q,ie);
                  q_max_s = q_max(g,k,q,ie);
                  first = false;
                } else {
                  q_min_s = std::min(q_min_s, q_min(g,k,q,ie));
                  q_max_s = std::max(q_max_s, q_max(g,k,q,ie));
                }
              }
            }
            vertical_caas_backup(n, rhom, q_min_s, q_max_s,
                                 Qm_min_tot, Qm_max_tot, Qm_tot,
                                 Qm_min, Qm_max, Qm);
          }
          // Redistribute mass in the horizontal direction of each level.
          for (Int i = 0; i < n; ++i) {
            const Int k = k0 + i;
            solve_local<np_>(ie, k, q, d.tl_np1, d.n1_qdp, np,
			     scalar_bounds, limiter_option,
			     spheremp, dp3d_c, q_min, q_max, Qm[i], qdp_c, q_c);
          }
        } else {
          for (Int sbli = 0; sbli < cdr.nsublev; ++sbli) {
            const Int k = k0 + sbli;
            if (k >= nlev) break;
            const auto ie_idx = cdr.cdr_over_super_levels ?
              nlevwrem*ie + k :
              cdr.nsublev*ie + sbli;
            const auto lci = cdr.ie2lci[ie_idx];
            const Real Qm = cdr.cdr->get_Qm(lci, ti);
            solve_local<np_>(ie, k, q, d.tl_np1, d.n1_qdp, np,
			     scalar_bounds, limiter_option,
			     spheremp, dp3d_c, q_min, q_max, Qm, qdp_c, q_c);
          }
        }
      }
    }
  }
}

void check (CDR& cdr, Data& d, const Real* q_min_r, const Real* q_max_r,
            const Int nets, const Int nete) {
  using cedr::mpi::reduce;

  const Int np = d.np, nlev = d.nlev, nsuplev = cdr.nsuplev, qsize = d.qsize,
    nprob = cdr.threed ? 1 : nsuplev;

  Kokkos::View<Real**, Kokkos::Serial>
    mass_p("mass_p", nprob, qsize), mass_c("mass_c", nprob, qsize),
    mass_lo("mass_lo", nprob, qsize), mass_hi("mass_hi", nprob, qsize),
    q_lo("q_lo", nprob, qsize), q_hi("q_hi", nprob, qsize),
    q_min_l("q_min_l", nprob, qsize), q_max_l("q_max_l", nprob, qsize),
    qd_lo("qd_lo", nprob, qsize), qd_hi("qd_hi", nprob, qsize);
  FA5<const Real>
    q_min(q_min_r, np, np, nlev, qsize, nete+1),
    q_max(q_max_r, np, np, nlev, qsize, nete+1);
  Kokkos::deep_copy(q_lo,  1e200);
  Kokkos::deep_copy(q_hi, -1e200);
  Kokkos::deep_copy(q_min_l,  1e200);
  Kokkos::deep_copy(q_max_l, -1e200);
  Kokkos::deep_copy(qd_lo, 0);
  Kokkos::deep_copy(qd_hi, 0);

  Int iprob = 0;

  bool fp_issue = false; // Limit output once the first issue is seen.
  for (Int ie = nets; ie <= nete; ++ie) {
    FA2<const Real> spheremp(d.spheremp[ie], np, np);
    FA5<const Real> qdp_pc(d.qdp_pc[ie], np, np, nlev, d.qsize_d, 2);
    FA4<const Real> dp3d_c(d.dp3d_c[ie], np, np, nlev, d.timelevels);
    FA4<const Real> q_c(d.q_c[ie], np, np, nlev, d.qsize_d);
    for (Int spli = 0; spli < nsuplev; ++spli) {
      if (nprob > 1) iprob = spli;
      for (Int k = spli*cdr.nsublev; k < (spli+1)*cdr.nsublev; ++k) {
        if (k >= nlev) continue;
        if ( ! fp_issue) {
          for (Int j = 0; j < np; ++j)
            for (Int i = 0; i < np; ++i) {
              // FP issues.
              if (std::isnan(dp3d_c(i,j,k,d.tl_np1)))
              { pr("dp3d NaN:" pu(k) pu(i) pu(j)); fp_issue = true; }
              if (std::isinf(dp3d_c(i,j,k,d.tl_np1)))
              { pr("dp3d Inf:" pu(k) pu(i) pu(j)); fp_issue = true; }
            }
        }
        for (Int q = 0; q < qsize; ++q) {
          Real qlo_s = q_min(0,0,k,q,ie), qhi_s = q_max(0,0,k,q,ie);
          for (Int j = 0; j < np; ++j)
            for (Int i = 0; i < np; ++i) {
              qlo_s = std::min(qlo_s, q_min(i,j,k,q,ie));
              qhi_s = std::max(qhi_s, q_max(i,j,k,q,ie));
            }
          for (Int j = 0; j < np; ++j)
            for (Int i = 0; i < np; ++i) {
              // FP issues.
              if ( ! fp_issue) {
                for (Int i_qdp : {0, 1}) {
                  const Int n_qdp = i_qdp == 0 ? d.n0_qdp : d.n1_qdp;
                  if (std::isnan(qdp_pc(i,j,k,q,n_qdp)))
                  { pr("qdp NaN:" puf(i_qdp) pu(q) pu(k) pu(i) pu(j)); fp_issue = true; }
                  if (std::isinf(qdp_pc(i,j,k,q,n_qdp)))
                  { pr("qdp Inf:" puf(i_qdp) pu(q) pu(k) pu(i) pu(j)); fp_issue = true; }
                }
                if (std::isnan(q_c(i,j,k,q)))
                { pr("q NaN:" pu(q) pu(k) pu(i) pu(j)); fp_issue = true; }
                if (std::isinf(q_c(i,j,k,q)))
                { pr("q Inf:" pu(q) pu(k) pu(i) pu(j)); fp_issue = true; }
              }
              // Mass conservation.
              mass_p(iprob,q) += qdp_pc(i,j,k,q,d.n0_qdp) * spheremp(i,j);
              mass_c(iprob,q) += qdp_pc(i,j,k,q,d.n1_qdp) * spheremp(i,j);
              // Local bound constraints w.r.t. cell-local extrema.
              if (q_c(i,j,k,q) < qlo_s)
                qd_lo(iprob,q) = std::max(qd_lo(iprob,q), qlo_s - q_c(i,j,k,q));
              if (q_c(i,j,k,q) > qhi_s)
                qd_hi(iprob,q) = std::max(qd_hi(iprob,q), q_c(i,j,k,q) - qhi_s);
              // Safety problem bound constraints.
              mass_lo(iprob,q) += (q_min(i,j,k,q,ie) * dp3d_c(i,j,k,d.tl_np1) *
                                  spheremp(i,j));
              mass_hi(iprob,q) += (q_max(i,j,k,q,ie) * dp3d_c(i,j,k,d.tl_np1) *
                                  spheremp(i,j));
              q_lo(iprob,q) = std::min(q_lo(iprob,q), q_min(i,j,k,q,ie));
              q_hi(iprob,q) = std::max(q_hi(iprob,q), q_max(i,j,k,q,ie));
              q_min_l(iprob,q) = std::min(q_min_l(iprob,q), q_min(i,j,k,q,ie));
              q_max_l(iprob,q) = std::max(q_max_l(iprob,q), q_max(i,j,k,q,ie));
            }
        }
      }
    }
  }

#ifdef HORIZ_OPENMP
# pragma omp barrier
# pragma omp master
#endif
  {
    if ( ! d.check)
      d.check = std::make_shared<Data::Check>(nprob, qsize);
    auto& c = *d.check;
    Kokkos::deep_copy(c.mass_p, 0);
    Kokkos::deep_copy(c.mass_c, 0);
    Kokkos::deep_copy(c.mass_lo, 0);
    Kokkos::deep_copy(c.mass_hi, 0);
    Kokkos::deep_copy(c.q_lo,  1e200);
    Kokkos::deep_copy(c.q_hi, -1e200);
    Kokkos::deep_copy(c.q_min_l,  1e200);
    Kokkos::deep_copy(c.q_max_l, -1e200);
    Kokkos::deep_copy(c.qd_lo, 0);
    Kokkos::deep_copy(c.qd_hi, 0);
  }

#ifdef HORIZ_OPENMP
# pragma omp barrier
# pragma omp critical
#endif
  {
    auto& c = *d.check;
    for (Int spli = 0; spli < nprob; ++spli) {
      if (nprob > 1) iprob = spli;
      for (Int q = 0; q < qsize; ++q) {
        c.mass_p(iprob,q) += mass_p(iprob,q);
        c.mass_c(iprob,q) += mass_c(iprob,q);
        c.qd_lo(iprob,q) = std::max(c.qd_lo(iprob,q), qd_lo(iprob,q));
        c.qd_hi(iprob,q) = std::max(c.qd_hi(iprob,q), qd_hi(iprob,q));
        c.mass_lo(iprob,q) += mass_lo(iprob,q);
        c.mass_hi(iprob,q) += mass_hi(iprob,q);
        c.q_lo(iprob,q) = std::min(c.q_lo(iprob,q), q_lo(iprob,q));
        c.q_hi(iprob,q) = std::max(c.q_hi(iprob,q), q_hi(iprob,q));
        c.q_min_l(iprob,q) = std::min(c.q_min_l(iprob,q), q_min_l(iprob,q));
        c.q_max_l(iprob,q) = std::max(c.q_max_l(iprob,q), q_max_l(iprob,q));
      }
    }
  }

#ifdef HORIZ_OPENMP
# pragma omp barrier
# pragma omp master
#endif
  {
    Kokkos::View<Real**, Kokkos::Serial>
      mass_p_g("mass_p_g", nprob, qsize), mass_c_g("mass_c_g", nprob, qsize),
      mass_lo_g("mass_lo_g", nprob, qsize), mass_hi_g("mass_hi_g", nprob, qsize),
      q_lo_g("q_lo_g", nprob, qsize), q_hi_g("q_hi_g", nprob, qsize),
      q_min_g("q_min_g", nprob, qsize), q_max_g("q_max_g", nprob, qsize),
      qd_lo_g("qd_lo_g", nprob, qsize), qd_hi_g("qd_hi_g", nprob, qsize);

    const auto& p = *cdr.p;
    const auto& c = *d.check;
    const auto root = cdr.p->root();
    const auto N = nprob*qsize;

    reduce(p, c.mass_p.data(), mass_p_g.data(), N, MPI_SUM, root);
    reduce(p, c.mass_c.data(), mass_c_g.data(), N, MPI_SUM, root);
    reduce(p, c.qd_lo.data(), qd_lo_g.data(), N, MPI_MAX, root);
    reduce(p, c.qd_hi.data(), qd_hi_g.data(), N, MPI_MAX, root);
    // Safety problem.
    reduce(p, c.mass_lo.data(), mass_lo_g.data(), N, MPI_SUM, root);
    reduce(p, c.mass_hi.data(), mass_hi_g.data(), N, MPI_SUM, root);
    reduce(p, c.q_lo.data(), q_lo_g.data(), N, MPI_MIN, root);
    reduce(p, c.q_hi.data(), q_hi_g.data(), N, MPI_MAX, root);
    reduce(p, c.q_min_l.data(), q_min_g.data(), N, MPI_MIN, root);
    reduce(p, c.q_max_l.data(), q_max_g.data(), N, MPI_MAX, root);

    if (cdr.p->amroot()) {
      const Real tol = 1e4*std::numeric_limits<Real>::epsilon();
      for (Int k = 0; k < nprob; ++k)
        for (Int q = 0; q < qsize; ++q) {
          const Real rd = cedr::util::reldif(mass_p_g(k,q), mass_c_g(k,q));
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
}

} // namespace sl
} // namespace homme

// Interface for Homme, through compose_mod.F90.
extern "C" void kokkos_init () {
  Kokkos::InitArguments args;
  args.disable_warnings = true;
  Kokkos::initialize(args);
}

extern "C" void kokkos_finalize () {
  Kokkos::finalize();
}

static homme::CDR::Ptr g_cdr;

extern "C" void
cedr_init_impl (const homme::Int fcomm, const homme::Int cdr_alg, const bool use_sgi,
                const homme::Int* gid_data, const homme::Int* rank_data,
                const homme::Int gbl_ncell, const homme::Int lcl_ncell,
                const homme::Int nlev, const bool independent_time_steps, const bool hard_zero,
                const homme::Int, const homme::Int) {
  const auto p = cedr::mpi::make_parallel(MPI_Comm_f2c(fcomm));
  g_cdr = std::make_shared<homme::CDR>(
    cdr_alg, gbl_ncell, lcl_ncell, nlev, use_sgi, independent_time_steps, hard_zero,
    gid_data, rank_data, p, fcomm);
}

extern "C" void cedr_query_bufsz (homme::Int* sendsz, homme::Int* recvsz) {
  cedr_assert(g_cdr);
  size_t s1, s2;
  g_cdr->get_buffers_sizes(s1, s2);
  *sendsz = static_cast<homme::Int>(s1);
  *recvsz = static_cast<homme::Int>(s2);
}

extern "C" void cedr_set_bufs (homme::Real* sendbuf, homme::Real* recvbuf,
                               homme::Int, homme::Int) {
  g_cdr->set_buffers(sendbuf, recvbuf);
}

extern "C" void cedr_unittest (const homme::Int fcomm, homme::Int* nerrp) {
#if 0
  cedr_assert(g_cdr);
  cedr_assert(g_cdr->tree);
  auto p = cedr::mpi::make_parallel(MPI_Comm_f2c(fcomm));
  if (homme::CDR::Alg::is_qlt(g_cdr->alg))
    *nerrp = cedr::qlt::test::test_qlt(p, g_cdr->tree, g_cdr->nsublev*g_cdr->ncell,
                                       1, false, false, true, false);
  else
    *nerrp = cedr::caas::test::unittest(p);
#endif
  *nerrp += homme::test_tree_maker();
  *nerrp += homme::CDR::QLTT::unittest();
}

extern "C" void cedr_set_ie2gci (const homme::Int ie, const homme::Int gci) {
  cedr_assert(g_cdr);
  // Now is a good time to drop the tree, whose persistence was used for unit
  // testing if at all.
  g_cdr->tree = nullptr;
  homme::set_ie2gci(*g_cdr, ie - 1, gci - 1);
}

static homme::sl::Data::Ptr g_sl;

extern "C" homme::Int cedr_sl_init (
  const homme::Int np, const homme::Int nlev, const homme::Int qsize,
  const homme::Int qsized, const homme::Int timelevels,
  const homme::Int need_conservation)
{
  cedr_assert(g_cdr);
  g_sl = std::make_shared<homme::sl::Data>(g_cdr->nlclcell, np, nlev, qsize,
                                           qsized, timelevels);
  homme::init_ie2lci(*g_cdr);
  homme::init_tracers(*g_cdr, nlev, qsize, need_conservation);
  homme::sl::check(*g_cdr, *g_sl);
  return 1;
}

extern "C" void cedr_sl_set_pointers_begin (homme::Int nets, homme::Int nete) {}
extern "C" void cedr_sl_set_spheremp (homme::Int ie, homme::Real* v)
{ homme::sl::insert(g_sl, ie - 1, 0, v); }
extern "C" void cedr_sl_set_qdp (homme::Int ie, homme::Real* v, homme::Int n0_qdp,
                                 homme::Int n1_qdp)
{ homme::sl::insert(g_sl, ie - 1, 1, v, n0_qdp - 1, n1_qdp - 1); }
extern "C" void cedr_sl_set_dp3d (homme::Int ie, homme::Real* v, homme::Int tl_np1)
{ homme::sl::insert(g_sl, ie - 1, 2, v, tl_np1 - 1); }
extern "C" void cedr_sl_set_dp (homme::Int ie, homme::Real* v)
{ homme::sl::insert(g_sl, ie - 1, 2, v, 0); }
extern "C" void cedr_sl_set_q (homme::Int ie, homme::Real* v)
{ homme::sl::insert(g_sl, ie - 1, 3, v); }
extern "C" void cedr_sl_set_dp0 (homme::Real* v)
{ homme::sl::insert(g_sl, 0, 4, v); }
extern "C" void cedr_sl_set_pointers_end () {}

// Run QLT.
extern "C" void cedr_sl_run (homme::Real* minq, const homme::Real* maxq,
                             homme::Int nets, homme::Int nete) {
  cedr_assert(minq != maxq);
  cedr_assert(g_cdr);
  cedr_assert(g_sl);
  homme::sl::run_global<4>(*g_cdr, *g_sl, minq, maxq, nets-1, nete-1);
}

// Run the cell-local limiter problem.
extern "C" void cedr_sl_run_local (homme::Real* minq, const homme::Real* maxq,
                                   homme::Int nets, homme::Int nete, homme::Int use_ir,
                                   homme::Int limiter_option) {
  cedr_assert(minq != maxq);
  cedr_assert(g_cdr);
  cedr_assert(g_sl);
  homme::sl::run_local<4>(*g_cdr, *g_sl, minq, maxq, nets-1, nete-1, use_ir,
                          limiter_option);
}

// Check properties for this transport step.
extern "C" void cedr_sl_check (const homme::Real* minq, const homme::Real* maxq,
                               homme::Int nets, homme::Int nete) {
  cedr_assert(g_cdr);
  cedr_assert(g_sl);
  homme::sl::check(*g_cdr, *g_sl, minq, maxq, nets-1, nete-1);
}

extern "C" void cedr_finalize () {
  g_sl = nullptr;
  g_cdr = nullptr;
}
