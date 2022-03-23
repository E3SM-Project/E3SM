// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

#ifndef INCLUDE_CEDR_TEST_RANDOMIZED_HPP
#define INCLUDE_CEDR_TEST_RANDOMIZED_HPP

#include "cedr_cdr.hpp"
#include "cedr_mpi.hpp"
#include "cedr_util.hpp"

#include <vector>

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

#include "cedr_test_randomized_inl.hpp"

#endif
