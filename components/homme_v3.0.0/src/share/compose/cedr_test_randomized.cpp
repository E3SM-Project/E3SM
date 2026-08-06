// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

#include "cedr_test_randomized.hpp"

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
  std::vector<Real> lcl_mass(3*ts.size()), q_min_lcl(ts.size()), q_max_lcl(ts.size());
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
      lcl_mass[3*ti    ] += Qm_prev[i];
      lcl_mass[3*ti + 1] += Qm[i];
      // Because the randomized test permits q < 0, in the global mass relative
      // error we need to be careful to choose a denominator that doesn't have
      // cancellation of masses from the signs. So accumulate abs as the
      // denominator.
      lcl_mass[3*ti + 2] += std::abs(Qm_prev[i]);
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

  std::vector<Real> glbl_mass(3*ts.size(), 0);
  mpi::reduce(p, lcl_mass.data(), glbl_mass.data(), lcl_mass.size(), MPI_SUM,
              p.root());
  std::vector<Int> t_ok_gbl(ts.size(), 0);
  mpi::reduce(p, t_ok.data(), t_ok_gbl.data(), t_ok.size(), MPI_MIN, p.root());
  // Right now we're not using these:
  std::vector<Int> local_violated_gbl(ts.size(), 0);
  mpi::reduce(p, local_violated.data(), local_violated_gbl.data(),
              local_violated.size(), MPI_MAX, p.root());

  if (p.amroot()) {
    const Real tol = 5e2*std::numeric_limits<Real>::epsilon();
    for (size_t ti = 0; ti < ts.size(); ++ti) {
      // Check mass conservation.
      const Real desired_mass = glbl_mass[3*ti], actual_mass = glbl_mass[3*ti+1],
        den_mass = glbl_mass[3*ti+2],
        rd = std::abs(actual_mass - desired_mass)/std::abs(den_mass);
      const bool mass_failed = rd > tol;
      if (mass_failed) {
        ++nerr;
        t_ok_gbl[ti] = false;
      }
      if ( ! t_ok_gbl[ti]) {
        std::cout << "FAIL " << cdr_name << ": " << ts[ti].str();
        if (mass_failed) std::cout << " mass re " << rd;
        std::cout << "\n";
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
