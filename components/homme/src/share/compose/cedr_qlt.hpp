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

#include "cedr_cdr.hpp"
#include "cedr_util.hpp"
#include "cedr_tree.hpp"

namespace cedr {
// QLT: Quasi-local tree-based non-iterative tracer density reconstructor for
//      mass conservation, shape preservation, and tracer consistency.
namespace qlt {
using cedr::mpi::Parallel;

template <typename ExeSpace = Kokkos::DefaultExecutionSpace>
class QLT : public cedr::CDR {
public:
  typedef typename cedr::impl::DeviceType<ExeSpace>::type Device;
  typedef QLT<ExeSpace> Me;
  typedef std::shared_ptr<Me> Ptr;
  typedef Kokkos::View<Real*, Device> RealList;
  typedef Kokkos::View<Int*, Device> IntList;
  typedef cedr::impl::Const<IntList> ConstIntList;
  typedef cedr::impl::ConstUnmanaged<IntList> ConstUnmanagedIntList;

protected:
  struct MetaDataBuilder {
    typedef std::shared_ptr<MetaDataBuilder> Ptr;
    std::vector<int> trcr2prob;
  };

public:
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
  
  // Set up QLT topology and communication data structures based on a tree. Both
  // ncells and tree refer to the global mesh, not just this processor's part.
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

  struct DeviceOp : public CDR::DeviceOp {
    // lclcellidx is gci2lci(cellidx).
    KOKKOS_INLINE_FUNCTION
    void set_rhom(const Int& lclcellidx, const Int& rhomidx, const Real& rhom) const override;

    // lclcellidx is gci2lci(cellidx).
    KOKKOS_INLINE_FUNCTION
    void set_Qm(const Int& lclcellidx, const Int& tracer_idx,
                const Real& Qm, const Real& Qm_min, const Real& Qm_max,
                const Real Qm_prev = cedr::impl::TypeTraits<Real>::infinity) const override;

    KOKKOS_INLINE_FUNCTION
    Real get_Qm(const Int& lclcellidx, const Int& tracer_idx) const override;

    /// View data for host and device computation.
    // Constructed in end_tracer_declarations().
    MetaData md_;
    BulkData bd_;
  };

  const DeviceOp& get_device_op() override;

  void run() override;

protected:
  static void init(const std::string& name, IntList& d,
                   typename IntList::HostMirror& h, size_t n);

  void init(const Parallel::Ptr& p, const Int& ncells, const tree::Node::Ptr& tree);

  void init_ordinals();

  /// Pointer data for initialization and host computation.
  Parallel::Ptr p_;
  // Tree and communication topology.
  std::shared_ptr<const tree::NodeSets> ns_;
  // Data extracted from ns_ for use in run() on device.
  std::shared_ptr<tree::NodeSetsDeviceData<ExeSpace> > nsdd_;
  std::shared_ptr<tree::NodeSetsHostData> nshd_;
  // Globally unique cellidx -> rank-local index.
  typedef std::map<Int,Int> Gci2LciMap;
  std::shared_ptr<Gci2LciMap> gci2lci_;
  // Temporary to collect caller's tracer information prior to calling
  // end_tracer_declarations().
  typename MetaDataBuilder::Ptr mdb_;
  DeviceOp o;

PRIVATE_CUDA:
  void l2r_recv(const tree::NodeSets::Level& lvl, const Int& l2rndps) const;
  void l2r_combine_kid_data(const Int& lvlidx, const Int& l2rndps) const;
  void l2r_send_to_parents(const tree::NodeSets::Level& lvl, const Int& l2rndps) const;
  void root_compute(const Int& l2rndps, const Int& r2lndps) const;
  void r2l_recv(const tree::NodeSets::Level& lvl, const Int& r2lndps) const;
  void r2l_solve_qp(const Int& lvlidx, const Int& l2rndps, const Int& r2lndps) const;
  void r2l_send_to_kids(const tree::NodeSets::Level& lvl, const Int& r2lndps) const;
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
#include "cedr_qlt_inl.hpp"

#endif
