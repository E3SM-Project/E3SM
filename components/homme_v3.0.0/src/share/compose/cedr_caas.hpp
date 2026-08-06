// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

#ifndef INCLUDE_CEDR_CAAS_HPP
#define INCLUDE_CEDR_CAAS_HPP

#include "cedr_cdr.hpp"

#include <vector>

namespace cedr {
// ClipAndAssuredSum.
namespace caas {

template <typename ExeSpace = Kokkos::DefaultExecutionSpace>
class CAAS : public CDR {
public:
  typedef typename cedr::impl::DeviceType<ExeSpace>::type Device;
  typedef CAAS<ExeSpace> Me;
  typedef std::shared_ptr<Me> Ptr;
  typedef Kokkos::View<Real*, Kokkos::LayoutRight, Device> RealList;
  typedef Kokkos::View<Int*, Kokkos::LayoutRight, Device> IntList;

public:

  // The caller may optionally provide its own all-reduce implementation.
  struct UserAllReducer {
    typedef std::shared_ptr<const UserAllReducer> Ptr;

    // MPI_Allreduce-like interface.
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

    // Reductions sometimes must be done BFB-invariant to rank decomposition;
    // that is the likeliest motivation for creating a UserAllReducer object. By
    // default, CAAS provides each DOF to the UserAllReducer. However, DOF order
    // may allow n DOFs in a row to be accumulated in a BFB-invariant way, e.g.,
    // if those DOFs are guaranteed always to be on the same processor. If so,
    // expose that value n here.
    virtual int n_accum_in_place () const { return 1; }
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

  struct DeviceOp : public CDR::DeviceOp {
    // lclcellidx is trivial; it is the user's index for the cell.
    KOKKOS_INLINE_FUNCTION
    void set_rhom(const Int& lclcellidx, const Int& rhomidx, const Real& rhom) const override;

    KOKKOS_INLINE_FUNCTION
    void set_Qm(const Int& lclcellidx, const Int& tracer_idx,
                const Real& Qm, const Real& Qm_min, const Real& Qm_max,
                const Real Qm_prev = cedr::impl::TypeTraits<Real>::infinity) const override;

    KOKKOS_INLINE_FUNCTION
    Real get_Qm(const Int& lclcellidx, const Int& tracer_idx) const override;

    Int nlclcells_, nrhomidxs_;
    bool need_conserve_;
    IntList probs_;
    RealList d_;
  };

  const DeviceOp& get_device_op() override;

  void run() override;

protected:
  typedef cedr::impl::Unmanaged<RealList> UnmanagedRealList;

  struct Decl {
    int probtype;
    Int rhomidx;
    Decl (const int probtype_, const Int rhomidx_)
      : probtype(probtype_), rhomidx(rhomidx_) {}
  };

  mpi::Parallel::Ptr p_;
  typename UserAllReducer::Ptr user_reducer_;
  std::shared_ptr<std::vector<Decl> > tracer_decls_;
  typename IntList::HostMirror probs_h_;
  IntList t2r_;
  RealList send_, recv_;
  bool finished_setup_;
  DeviceOp o;

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

#include "cedr_caas_inl.hpp"

#endif
