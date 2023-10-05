// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

#include "cedr_qlt.hpp"
#include "cedr_test_randomized.hpp"

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
void init_device_data (const tree::NodeSets& ns, tree::NodeSetsHostData& h,
                       tree::NodeSetsDeviceData<ES>& d) {
  typedef tree::NodeSetsDeviceData<ES> NSDD;
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
  ns_ = tree::analyze(p, ncells, tree);
  nshd_ = std::make_shared<tree::NodeSetsHostData>();
  nsdd_ = std::make_shared<tree::NodeSetsDeviceData<ES> >();
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
  problem_type = o.md_.get_problem_type(o.md_.get_problem_type_idx(problem_type));
  mdb_->trcr2prob.push_back(problem_type);
}

template <typename ES>
void QLT<ES>::end_tracer_declarations () {
  o.md_.init(*mdb_);
  mdb_ = nullptr;
}

template <typename ES>
void QLT<ES>::get_buffers_sizes (size_t& buf1, size_t& buf2) {
  const auto nslots = ns_->nslots;
  buf1 = o.md_.a_h.prob2bl2r[o.md_.nprobtypes]*nslots;
  buf2 = o.md_.a_h.prob2br2l[o.md_.nprobtypes]*nslots;
}

template <typename ES>
void QLT<ES>::set_buffers (Real* buf1, Real* buf2) {
  size_t l2r_sz, r2l_sz;
  get_buffers_sizes(l2r_sz, r2l_sz);
  o.bd_.init(buf1, l2r_sz, buf2, r2l_sz);
}

template <typename ES>
void QLT<ES>::finish_setup () {
  if (o.bd_.inited()) return;
  size_t l2r_sz, r2l_sz;
  get_buffers_sizes(l2r_sz, r2l_sz);
  o.bd_.init(l2r_sz, r2l_sz);
}

template <typename ES>
int QLT<ES>::get_problem_type (const Int& tracer_idx) const {
  cedr_throw_if(tracer_idx < 0 || tracer_idx > o.md_.a_h.trcr2prob.extent_int(0),
                "tracer_idx is out of bounds: " << tracer_idx);
  return o.md_.a_h.trcr2prob[tracer_idx];
}

template <typename ES>
Int QLT<ES>::get_num_tracers () const {
  return o.md_.a_h.trcr2prob.size();
}

template <typename ES> void QLT<ES>
::l2r_recv (const tree::NodeSets::Level& lvl, const Int& l2rndps) const {
  for (size_t i = 0; i < lvl.kids.size(); ++i) {
    const auto& mmd = lvl.kids[i];
    mpi::irecv(*p_, o.bd_.l2r_data.data() + mmd.offset*l2rndps, mmd.size*l2rndps,
               mmd.rank, tree::NodeSets::mpitag, &lvl.kids_req[i]);
  }
  Timer::start(Timer::waitall);
  mpi::waitall(lvl.kids_req.size(), lvl.kids_req.data());
  Timer::stop(Timer::waitall);
}

template <typename ES> void QLT<ES>
::l2r_combine_kid_data (const Int& lvlidx, const Int& l2rndps) const {
  if (cedr::impl::OnGpu<ES>::value) {
    const auto d = *nsdd_;
    const auto l2r_data = o.bd_.l2r_data;
    const auto a = o.md_.a_d;
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
      o.bd_.l2r_data(n->offset*l2rndps) =
        (o.bd_.l2r_data(ns_->node_h(n->kids[0])->offset*l2rndps) +
         o.bd_.l2r_data(ns_->node_h(n->kids[1])->offset*l2rndps));
      // Tracers.
      for (Int pti = 0; pti < o.md_.nprobtypes; ++pti) {
        const Int problem_type = o.md_.get_problem_type(pti);
        const bool nonnegative = problem_type & ProblemType::nonnegative;
        const bool shapepreserve = problem_type & ProblemType::shapepreserve;
        const bool conserve = problem_type & ProblemType::conserve;
        const Int bis = o.md_.a_d.prob2trcrptr[pti], bie = o.md_.a_d.prob2trcrptr[pti+1];
        for (Int bi = bis; bi < bie; ++bi) {
          const Int bdi = o.md_.a_d.trcr2bl2r(o.md_.a_d.bidx2trcr(bi));
          Real* const me = &o.bd_.l2r_data(n->offset*l2rndps + bdi);
          const auto kid0 = ns_->node_h(n->kids[0]);
          const auto kid1 = ns_->node_h(n->kids[1]);
          const Real* const k0 = &o.bd_.l2r_data(kid0->offset*l2rndps + bdi);
          const Real* const k1 = &o.bd_.l2r_data(kid1->offset*l2rndps + bdi);
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
::l2r_send_to_parents (const tree::NodeSets::Level& lvl, const Int& l2rndps) const {
  for (size_t i = 0; i < lvl.me.size(); ++i) {
    const auto& mmd = lvl.me[i];
    mpi::isend(*p_, o.bd_.l2r_data.data() + mmd.offset*l2rndps, mmd.size*l2rndps,
               mmd.rank, tree::NodeSets::mpitag);
  }  
}

template <typename ES> void QLT<ES>
::root_compute (const Int& l2rndps, const Int& r2lndps) const {
  if (ns_->levels.empty() || ns_->levels.back().nodes.size() != 1 ||
      ns_->node_h(ns_->levels.back().nodes[0])->parent >= 0)
    return;
  const auto d = *nsdd_;
  const auto l2r_data = o.bd_.l2r_data;
  const auto r2l_data = o.bd_.r2l_data;
  const auto a = o.md_.a_d;
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
                    MetaData::get_problem_type_l2r_bulk_size(problem_type) - 1 :
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
::r2l_recv (const tree::NodeSets::Level& lvl, const Int& r2lndps) const {
  for (size_t i = 0; i < lvl.me.size(); ++i) {
    const auto& mmd = lvl.me[i];
    mpi::irecv(*p_, o.bd_.r2l_data.data() + mmd.offset*r2lndps, mmd.size*r2lndps,
               mmd.rank, tree::NodeSets::mpitag, &lvl.me_recv_req[i]);
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
  const tree::NodeSets::Node& n,
  const tree::NodeSets::Node& k0, const tree::NodeSets::Node& k1,
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
    const auto l2r_data = o.bd_.l2r_data;
    const auto r2l_data = o.bd_.r2l_data;
    const auto a = o.md_.a_d;
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
      for (Int pti = 0; pti < o.md_.nprobtypes; ++pti) {
        const Int problem_type = o.md_.get_problem_type(pti);
        const Int bis = o.md_.a_d.prob2trcrptr[pti], bie = o.md_.a_d.prob2trcrptr[pti+1];
        for (Int bi = bis; bi < bie; ++bi) {
          const Int l2rbdi = o.md_.a_d.trcr2bl2r(o.md_.a_d.bidx2trcr(bi));
          const Int r2lbdi = o.md_.a_d.trcr2br2l(o.md_.a_d.bidx2trcr(bi));
          cedr_assert(n->nkids == 2);
          if ((problem_type & ProblemType::consistent) &&
              ! (problem_type & ProblemType::shapepreserve)) {
            const Real q_min = o.bd_.r2l_data(n->offset*r2lndps + r2lbdi + 1);
            const Real q_max = o.bd_.r2l_data(n->offset*r2lndps + r2lbdi + 2);
            o.bd_.l2r_data(n->offset*l2rndps + l2rbdi + 0) = q_min;
            o.bd_.l2r_data(n->offset*l2rndps + l2rbdi + 2) = q_max;
            for (Int k = 0; k < 2; ++k)
              r2l_solve_qp_set_q(o.bd_.l2r_data, o.bd_.r2l_data,
                                 ns_->node_h(n->kids[k])->offset,
                                 l2rndps, r2lndps, l2rbdi, r2lbdi, q_min, q_max);
          }
          r2l_solve_qp_solve_node_problem(
            o.bd_.l2r_data, o.bd_.r2l_data, problem_type, *n, *ns_->node_h(n->kids[0]),
            *ns_->node_h(n->kids[1]), l2rndps, r2lndps, l2rbdi, r2lbdi,
            prefer_mass_con_to_bounds);
        }
      }
    }
  }
  Timer::stop(Timer::snp);
}

template <typename ES> void QLT<ES>
::r2l_send_to_kids (const tree::NodeSets::Level& lvl, const Int& r2lndps) const {
  for (size_t i = 0; i < lvl.kids.size(); ++i) {
    const auto& mmd = lvl.kids[i];
    mpi::isend(*p_, o.bd_.r2l_data.data() + mmd.offset*r2lndps, mmd.size*r2lndps,
               mmd.rank, tree::NodeSets::mpitag);
  }
}

template <typename ES>
const typename QLT<ES>::DeviceOp& QLT<ES>::get_device_op() { return o; }

template <typename ES>
void QLT<ES>::run () {
  cedr_assert(o.bd_.inited());
  Timer::start(Timer::qltrunl2r);
  // Number of data per slot.
  const Int l2rndps = o.md_.a_h.prob2bl2r[o.md_.nprobtypes];
  const Int r2lndps = o.md_.a_h.prob2br2l[o.md_.nprobtypes];
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

Int unittest_QLT (const Parallel::Ptr& p, const bool write_requested=false) {
  using Mesh = tree::oned::Mesh;
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

namespace test {
Int run_unit_and_randomized_tests (const Parallel::Ptr& p, const Input& in) {
  Int nerr = 0;
  if (in.unittest) {
    Int ne = unittest_QLT(p, in.write);
    if (ne && p->amroot()) std::cerr << "FAIL: tree::oned::unittest_QLT()\n";
    nerr += ne;
    if (p->amroot()) std::cout << "\n";
  }
  // Performance test.
  if (in.perftest && in.ncells > 0) {
    tree::oned::Mesh m(in.ncells, p,
                       (in.pseudorandom ?
                        tree::oned::Mesh::ParallelDecomp::pseudorandom :
                        tree::oned::Mesh::ParallelDecomp::contiguous));
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
#ifdef CEDR_ENABLE_GPU
template class cedr::qlt::QLT<CedrGpuExeSpace>;
#endif
#ifdef KOKKOS_ENABLE_THREADS
template class cedr::qlt::QLT<Kokkos::Threads>;
#endif
