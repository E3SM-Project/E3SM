// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

#include "cedr_caas.hpp"
#include "cedr_util.hpp"
#include "cedr_test_randomized.hpp"

namespace Kokkos {
struct ComposeReal2 {
  cedr::Real v[2];
  KOKKOS_INLINE_FUNCTION ComposeReal2 () { v[0] = v[1] = 0; }

  KOKKOS_INLINE_FUNCTION void operator= (const ComposeReal2& s) {
    v[0] = s.v[0];
    v[1] = s.v[1];
  }
  KOKKOS_INLINE_FUNCTION void operator= (const volatile ComposeReal2& s) volatile {
    v[0] = s.v[0];
    v[1] = s.v[1];
  }

  KOKKOS_INLINE_FUNCTION ComposeReal2& operator+= (const ComposeReal2& o) {
    v[0] += o.v[0];
    v[1] += o.v[1];
    return *this;
  }
};

template<> struct reduction_identity<ComposeReal2> {
  KOKKOS_INLINE_FUNCTION static ComposeReal2 sum() { return ComposeReal2(); }
};
} // namespace Kokkos

namespace cedr {
namespace caas {

template <typename ES>
CAAS<ES>::CAAS (const mpi::Parallel::Ptr& p, const Int nlclcells,
                const typename UserAllReducer::Ptr& uar) {
  p_ = p;
  user_reducer_ = uar;
  o.nlclcells_ = nlclcells;
  o.nrhomidxs_ = 0;
  o.need_conserve_ = false;
  finished_setup_ = false;
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
    o.need_conserve_ = true;
  o.nrhomidxs_ = std::max(o.nrhomidxs_, rhomidx+1);
}

template <typename ES>
void CAAS<ES>::end_tracer_declarations () {
  cedr_throw_if(tracer_decls_->size() == 0, "#tracers is 0.");
  cedr_throw_if(o.nrhomidxs_ == 0, "#rhomidxs is 0.");
  o.probs_ = IntList("CAAS probs", static_cast<Int>(tracer_decls_->size()));
  probs_h_ = Kokkos::create_mirror_view(o.probs_);
  //t2r_ = IntList("CAAS t2r", static_cast<Int>(tracer_decls_->size()));
  for (Int i = 0; i < o.probs_.extent_int(0); ++i) {
    probs_h_(i) = (*tracer_decls_)[i].probtype;
    //t2r_(i) = (*tracer_decls_)[i].rhomidx;
  }
  Kokkos::deep_copy(o.probs_, probs_h_);
  tracer_decls_ = nullptr;
}

template <typename ES>
void CAAS<ES>::get_buffers_sizes (size_t& buf1, size_t& buf2, size_t& buf3) {
  const Int e = o.need_conserve_ ? 1 : 0;
  const auto nslots = 4*o.probs_.size();
  buf1 = o.nlclcells_ * ((3+e)*o.probs_.size() + 1);
  cedr_assert( ! user_reducer_ || o.nlclcells_ % user_reducer_->n_accum_in_place() == 0);
  buf2 = nslots*(user_reducer_ ? (o.nlclcells_ / user_reducer_->n_accum_in_place()) : 1);
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
  o.d_ = RealList(buf1, buf1sz);
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
  o.d_ = RealList("CAAS data", buf1);
  // (e'Qm_clip, e'Qm, e'Qm_min, e'Qm_max, [e'Qm_prev])
  send_ = RealList("CAAS send", buf2);
  recv_ = RealList("CAAS recv", buf3);
  finished_setup_ = true;
}

template <typename ES>
int CAAS<ES>::get_problem_type (const Int& tracer_idx) const {
  cedr_assert(tracer_idx >= 0 && tracer_idx < o.probs_.extent_int(0));
  return probs_h_[tracer_idx];
}

template <typename ES>
Int CAAS<ES>::get_num_tracers () const {
  return o.probs_.extent_int(0);
}

template <typename ES>
void CAAS<ES>::reduce_locally () {
  const bool user_reduces = user_reducer_ != nullptr;
  ConstExceptGnu Int nt = o.probs_.size(), nlclcells = o.nlclcells_;

  const auto probs = o.probs_;
  const auto send = send_;
  const auto d = o.d_;
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
    Kokkos::parallel_for(Kokkos::RangePolicy<ES>(0, nt*nlclaccum), calc_Qm_clip);
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
    Kokkos::parallel_for(Kokkos::RangePolicy<ES>(0, 2*nt*nlclaccum), set_Qm_minmax);
  } else {
    using ESU = cedr::impl::ExeSpaceUtils<ES>;
    const auto calc_Qm_clip = KOKKOS_LAMBDA (const typename ESU::Member& t) {
      const auto k = t.league_rank();
      const auto os = (k+1)*nlclcells;
      const auto reduce = [&] (const Int& i, Kokkos::ComposeReal2& accum) {
        Real Qm_clip, Qm_term;
        calc_Qm_scalars(d, probs, nt, nlclcells, k, os, i, Qm_clip, Qm_term);
        d(os+i) = Qm_clip;
        accum.v[0] += Qm_clip;
        accum.v[1] += Qm_term;
      };
      Kokkos::ComposeReal2 accum;
      Kokkos::parallel_reduce(Kokkos::TeamThreadRange(t, nlclcells),
                              reduce, Kokkos::Sum<Kokkos::ComposeReal2>(accum));
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
  ConstExceptGnu Int nt = o.probs_.size(), nlclcells = o.nlclcells_;
  const auto recv = recv_;
  const auto d = o.d_;
  const auto adjust_Qm = KOKKOS_LAMBDA (const typename ESU::Member& t) {
    const auto k = t.league_rank();
    const auto os = (k+1)*nlclcells;
    const auto Qm_clip_sum = recv(     k);
    const auto Qm_sum      = recv(nt + k);
    const auto m = Qm_sum - Qm_clip_sum;
    if (m < 0) {
      const auto Qm_min_sum = recv(2*nt + k);
      auto fac = Qm_clip_sum - Qm_min_sum;
      if (fac > 0) {
        fac = m/fac;
        const auto adjust = [&] (const Int& i) {
          const auto Qm_min = d(os + nlclcells * nt + i);
          auto& Qm = d(os+i);
          Qm += fac*(Qm - Qm_min);
          Qm = impl::max(Qm_min, Qm);
        };
        Kokkos::parallel_for(Kokkos::TeamThreadRange(t, nlclcells), adjust);
      }
    } else if (m > 0) {
      const auto Qm_max_sum = recv(3*nt + k);
      auto fac = Qm_max_sum - Qm_clip_sum;
      if (fac > 0) {
        fac = m/fac;
        const auto adjust = [&] (const Int& i) {
          const auto Qm_max = d(os + nlclcells*2*nt + i);
          auto& Qm = d(os+i);
          Qm += fac*(Qm_max - Qm);
          Qm = impl::min(Qm_max, Qm);
        };
        Kokkos::parallel_for(Kokkos::TeamThreadRange(t, nlclcells), adjust);
      }
    }
  };
  Kokkos::parallel_for(ESU::get_default_team_policy(nt, nlclcells),
                       adjust_Qm);
}

template <typename ES>
const typename CAAS<ES>::DeviceOp& CAAS<ES>::get_device_op() { return o; }

template <typename ES>
void CAAS<ES>::run () {
  cedr_assert(finished_setup_);
  reduce_locally();
  const bool user_reduces = user_reducer_ != nullptr;
  if (user_reduces)
    (*user_reducer_)(*p_, send_.data(), recv_.data(),
                     o.nlclcells_ / user_reducer_->n_accum_in_place(),
                     recv_.size(), MPI_SUM);
  else
    reduce_globally();
  finish_locally();
}

namespace test {
struct TestCAAS : public cedr::test::TestRandomized {
  typedef CAAS<Kokkos::DefaultExecutionSpace> CAAST;

  struct TestAllReducer : public CAAST::UserAllReducer {
    TestAllReducer (const Int n_accum) : n_(n_accum) {}

    Int n_accum_in_place () const override { return n_; }

    int operator() (const mpi::Parallel& p, Real* sendbuf, Real* rcvbuf,
                    int nlcl, int count, MPI_Op op) const override {
      Kokkos::View<Real*> s(sendbuf, nlcl*count), r(rcvbuf, count);
      const auto s_h = Kokkos::create_mirror_view(s);
      Kokkos::deep_copy(s_h, s);
      const auto r_h = Kokkos::create_mirror_view(r);
      for (int k = 0; k < count; ++k) {
        // When k == 0, s_h(0:nlcl-1) is summed. Then s_h(1:nlcl-1) is
        // free to be overwritten.
        s_h(k) = s_h(nlcl*k);
        for (int i = 1; i < nlcl; ++i)
          s_h(k) += s_h(nlcl*k + i);
      }
      const int err = mpi::all_reduce(p, s_h.data(), r_h.data(), count, op);
      Kokkos::deep_copy(r, r_h);
      return err;
    }

  private:
    Int n_;
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
    CAAST::UserAllReducer::Ptr reducer;
    if (use_own_reducer) {
      // Although it doesn't make sense to do this in practice, nothing prevents
      // us from using different n_accum_in_place values in different ranks. So
      // just compute a local value.
      const Int n_accum = (nlclcells_ % 3 == 0 ? 3 :
                           nlclcells_ % 2 == 0 ? 2 : 1);
      reducer = std::make_shared<TestAllReducer>(n_accum);
    }
    caas_ = std::make_shared<CAAST>( p, nlclcells_, reducer);
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
#ifdef CEDR_ENABLE_GPU
template class cedr::caas::CAAS<CedrGpuExeSpace>;
#endif
#ifdef KOKKOS_ENABLE_THREADS
template class cedr::caas::CAAS<Kokkos::Threads>;
#endif
