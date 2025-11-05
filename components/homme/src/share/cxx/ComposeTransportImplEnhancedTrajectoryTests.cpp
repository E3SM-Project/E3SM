/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "Config.hpp"
#ifdef HOMME_ENABLE_COMPOSE

#include "ComposeTransportImplEnhancedTrajectoryImpl.hpp"

#include <random>

namespace Homme {
namespace { // anon

Kokkos::TeamPolicy<ExecSpace>
get_test_team_policy (const int nelem, const int nlev, const int ncol=NP*NP) {
  ThreadPreferences tp;
  tp.max_threads_usable = ncol;
  tp.max_vectors_usable = nlev;
  tp.prefer_threads = true;
  tp.prefer_larger_team = true;
  return Homme::get_default_team_policy<ExecSpace>(nelem, tp);
}

struct TestData {
  const ComposeTransportImpl& cti;
  std::mt19937_64 engine;
  static const Real eps;

  TestData (const CTI& cti_, const int seed = 0)
    : cti(cti_), engine(seed == 0 ? std::random_device()() : seed)
  {}

  Real urand (const Real lo = 0, const Real hi = 1) {
    std::uniform_real_distribution<Real> urb(lo, hi);
    return urb(engine);
  }
};

// Data to deal with views of packs easily in tests.
struct ColData {
  int npack;
  ExecView<Scalar*> d;
  ExecView<Scalar*>::HostMirror h;
  ExecView<Real*>::HostMirror r;

  ColData (const std::string& name, const int nlev) {
    npack = calc_npack(nlev);
    d = decltype(d)(name, npack);
    h = Kokkos::create_mirror_view(d);
    r = decltype(r)(cti::pack2real(h), calc_nscal(npack));
  }

  void h2d () { Kokkos::deep_copy(d, h); }
};

struct ElData {
  int npack;
  ExecView<Scalar***> d;
  ExecView<Scalar***>::HostMirror h;
  ExecView<Real***>::HostMirror r;

  ElData (const std::string& name, const int nlev) {
    npack = calc_npack(nlev);
    d = decltype(d)(name, NP, NP, npack);
    h = Kokkos::create_mirror_view(d);
    r = decltype(r)(cti::pack2real(h), NP, NP, calc_nscal(npack));
  }

  void d2h () { Kokkos::deep_copy(h, d); }
  void h2d () { Kokkos::deep_copy(d, h); }
};

const Real TestData::eps = std::numeric_limits<Real>::epsilon();

int test_deriv (const TestData& td) {
  int ne = 0;
  { // The two impls are equivalent in infinite precision.
    const int n = 3;
    const Real x[n] = {-0.3, 0.1, 1.1};
    Real y[n] = {-0.1, 0.1, 0};
    for (int t = 0; t < 10; ++t) {
      y[2] = 0.2*(t - 9.5);
      Real g1, g2;
      g1 = cti::approx_derivative1(x[0], x[1], x[2], y[0], y[1], y[2]);
      g2 = cti::approx_derivative (x[0], x[1], x[2], y[0], y[1], y[2]);
      if (std::abs(g1 - g2) > 100*td.eps) ++ne;
    }
  }
  { // Exactly recovers quadratic.
    const int n = 3;
    const Real x[n] = {-0.3, 0.1, 1.1};
    Real y[n];
    for (int i = 0; i < n; ++i)
      y[i] = 0.7*x[i]*x[i] - 1.2*x[i] + 0.7;
    const Real g_true = 1.4*x[1] - 1.2;
    const Real g_est = cti::approx_derivative(x[0], x[1], x[2], y[0], y[1], y[2]);
    if (std::abs(g_est - g_true) > 100*td.eps) ++ne;
  }
  return ne;
}

int test_find_support (TestData&) {
  int ne = 0;
  const int n = 97;
  std::vector<Real> x(n);
  for (int i = 0; i < n; ++i) x[i] = -11.7 + (i*i)/n;
  const int ntest = 10000;
  for (int i = 0; i < ntest; ++i) {
    const Real xi = x[0] + (Real(i)/ntest)*(x[n-1] - x[0]);
    for (int x_idx : {0, 1, n/3, n/2, n-2, n-1}) {
      const int sup = find_support(n, x.data(), x_idx, xi);
      if (sup > n-2) ++ne;
      else if (xi < x[sup] or xi > x[sup+1]) ++ne;
    }
  }
  return ne;
}

void todev (const std::vector<Real>& h, const RnV& d) {
  assert(h.size() <= d.size());
  const auto m = Kokkos::create_mirror_view(d);
  for (size_t i = 0; i < h.size(); ++i) m(i) = h[i];
  Kokkos::deep_copy(d, m);
}

void fillcols (const int n, const Real* const h, const RelnV::HostMirror& a) {
  assert(n <= a.extent_int(2));
  for (int i = 0; i < a.extent_int(0); ++i)
    for (int j = 0; j < a.extent_int(1); ++j)
      for (int k = 0; k < n; ++k)
        a(i,j,k) = h[k];  
}

void todev (const int n, const Real* const h, const RelnV& d) {
  const auto m = Kokkos::create_mirror_view(d);
  fillcols(n, h, m)  ;
  Kokkos::deep_copy(d, m);
}

void todev (const std::vector<Real>& h, const RelnV& d) {
  todev(h.size(), h.data(), d);
}

void tohost (const ExecView<const Real*>& d, std::vector<Real>& h) {
  assert(h.size() <= d.size());
  const auto m = Kokkos::create_mirror_view(d);
  Kokkos::deep_copy(m, d);
  for (size_t i = 0; i < h.size(); ++i) h[i] = m(i);
}

void run_linterp (const std::vector<Real>& x, const std::vector<Real>& y,
                  std::vector<Real>& xi, std::vector<Real>& yi) {
  const auto n = x.size(), ni = xi.size();
  assert(y.size() == n); assert(yi.size() == ni);
  // input -> device (test different sizes >= n)
  ExecView<Real*> xv("xv", n), yv("yv", n+1), xiv("xiv", ni+2), yiv("yiv", ni+3);
  todev(x, xv);
  todev(y, yv);
  todev(xi, xiv);
  // call linterp
  const auto f = KOKKOS_LAMBDA(const cti::MT& team) {
    const auto range = Kokkos::TeamVectorRange(team, ni);
    linterp(range, n, xv, yv, ni, xiv, yiv, 0, "unittest");
  };
  Homme::ThreadPreferences tp;
  tp.max_threads_usable = 1;
  tp.max_vectors_usable = ni;
  tp.prefer_threads = false;
  tp.prefer_larger_team = true;
  const auto policy = get_test_team_policy(1, n);
  Kokkos::parallel_for(policy, f);
  Kokkos::fence();
  // output -> host
  tohost(yiv, yi);
}

void make_random_sorted (TestData& td, const int n, const Real xlo, const Real xhi,
                         std::vector<Real>& x) {
  assert(n >= 2);
  x.resize(n);
  x[0] = xlo;
  for (int i = 1; i < n-1; ++i) x[i] = td.urand(xlo, xhi);
  x[n-1] = xhi;
  std::sort(x.begin(), x.end());
}

int test_linterp (TestData& td) {
  int nerr = 0;
  { // xi == x => yi == y.
    int ne = 0;
    const int n = 30;
    std::vector<Real> x(n), y(n), xi(n), yi(n);
    make_random_sorted(td, n, -0.1, 1.2, x);
    make_random_sorted(td, n, -3, -1, y);
    for (int i = 0; i < n; ++i) xi[i] = x[i];
    run_linterp(x, y, xi, yi);
    for (int i = 0; i < n; ++i)
      if (yi[i] != y[i])
        ++ne;
    nerr += ne;
  }
  { // Reconstruct a linear function exactly.
    int ne = 0;
    const int n = 56, ni = n-3;
    const Real xlo = -1.2, xhi = 3.1;
    const auto f = [&] (const Real x) { return -0.7 + 1.3*x; };
    std::vector<Real> x(n), y(n), xi(ni), yi(ni);
    for (int trial = 0; trial < 4; ++trial) {
      make_random_sorted(td, n, xlo, xhi, x);
      make_random_sorted(td, ni,
                         xlo + (trial == 1 or trial == 3 ?  0.1 : 0),
                         xhi + (trial == 2 or trial == 3 ? -0.5 : 0),
                         xi);
      for (int i = 0; i < n; ++i) y[i] = f(x[i]);
      run_linterp(x, y, xi, yi);
      for (int i = 0; i < ni; ++i)
        if (std::abs(yi[i] - f(xi[i])) > 100*td.eps)
          ++ne;
    }
    nerr += ne;
  }
  return nerr;
}

int make_random_deta (TestData& td, const Real deta_tol, const int nlev,
                      Real* const deta) {
  int nerr = 0;
  Real sum = 0;
  for (int k = 0; k < nlev; ++k) {
    deta[k] = td.urand(0, 1) + 0.1;
    sum += deta[k];
  }
  for (int k = 0; k < nlev; ++k) {
    deta[k] /= sum;
    if (deta[k] < deta_tol) ++nerr;
  }
  return nerr;
}

int make_random_deta (TestData& td, const Real deta_tol, const RnV& deta) {
  int nerr = 0;
  const int nlev = deta.extent_int(0);
  const auto m = Kokkos::create_mirror_view(deta);
  nerr = make_random_deta(td, deta_tol, nlev, &m(0));
  Kokkos::deep_copy(deta, m);
  return nerr;  
}

int make_random_deta (TestData& td, const Real deta_tol, const RelnV& deta) {
  int nerr = 0;
  const int nlev = deta.extent_int(2);
  const auto m = Kokkos::create_mirror_view(deta);
  for (int i = 0; i < NP; ++i)
    for (int j = 0; j < NP; ++j)
      nerr += make_random_deta(td, deta_tol, nlev, &m(i,j,0));
  Kokkos::deep_copy(deta, m);
  return nerr;
}

// Wrapper to main deta_caas routine; this wrapper is used in the unit test of
// the main routine.
KOKKOS_FUNCTION void
deta_caas (const KernelVariables& kv, const int nlevp, const CRnV& deta_ref,
           const Real low, const RelnV& wrk, const RelnV& deta) {
  assert(deta_ref.extent_int(0) >= nlevp);
  assert_eln(wrk, nlevp);
  assert_eln(deta, nlevp);
  const auto ttr = Kokkos::TeamThreadRange(kv.team, NP*NP);
  const auto tvr = Kokkos::ThreadVectorRange(kv.team, nlevp);
  const auto f = [&] (const int idx) {
    const int i = idx / NP, j = idx % NP;
    deta_caas(kv, tvr, deta_ref, low, getcol(wrk,i,j), getcol(deta,i,j));
  };
  Kokkos::parallel_for(ttr, f);
}

int test_deta_caas (TestData& td) {
  int nerr = 0;
  const Real tol = 100*td.eps;

  for (const int nlev : {15, 128, 161}) {
    const Real deta_tol = 10*td.eps/nlev;
    const auto err = [&] (const char* lbl) {
      ++nerr;
      printf("test_deta_caa nlev %d: %s\n", nlev, lbl);
    };

    // nlev+1 deltas: deta = diff([0, etam, 1])
    ExecView<Real*> deta_ref("deta_ref", nlev+1);
    ExecView<Real***> deta("deta",NP,NP,nlev+1), wrk("wrk",NP,NP,nlev+1);
    nerr += make_random_deta(td, deta_tol, deta_ref);

    const auto policy = get_test_team_policy(1, nlev);
    const auto run = [&] (const RelnV& deta) {
      const auto f = KOKKOS_LAMBDA(const cti::MT& team) {
        KernelVariables kv(team);
        deta_caas(kv, nlev+1, deta_ref, deta_tol, wrk, deta);
      };
      Kokkos::parallel_for(policy, f);
      Kokkos::fence();
    };

    { // Test that if all is OK, the input is not altered.
      nerr += make_random_deta(td, deta_tol, deta);
      ExecView<Real***>::HostMirror copy("copy",NP,NP,nlev+1);
      Kokkos::deep_copy(copy, deta);
      run(deta);
      const auto m = cti::cmvdc(deta);
      bool diff = false;
      for (int i = 0; i < NP; ++i)
        for (int j = 0; j < NP; ++j)
          for (int k = 0; k <= nlev; ++k)
            if (m(i,j,k) != copy(i,j,k))
              diff = true;
      if (diff) err("input not altered");
    }

    { // Modify one etam and test that only adjacent intervals change beyond eps.
      // nlev midpoints
      ExecView<Real*> etam_ref("etam_ref",nlev);
      const auto her = Kokkos::create_mirror_view(etam_ref);
      const auto hder = cti::cmvdc(deta_ref);
      {
        her(0) = hder(0);
        for (int k = 1; k < nlev; ++k)
          her(k) = her(k-1) + hder(k);
        Kokkos::deep_copy(etam_ref, her);
      }
      std::vector<Real> etam(nlev);
      const auto hde = Kokkos::create_mirror_view(deta);
      const auto get_idx = [&] (const int i, const int j) {
        const int idx = static_cast<int>(0.15*nlev);
        return std::max(1, std::min(nlev-2, idx+NP*i+j));
      };
      for (int trial = 0; trial < 2; ++trial) {
        for (int i = 0; i < NP; ++i)
          for (int j = 0; j < NP; ++j) {
            for (int k = 0; k < nlev; ++k) etam[k] = her(k);
            // Perturb one level.
            const int idx = get_idx(i,j);
            etam[idx] += trial == 0 ? 1.1 : -13.1;
            hde(i,j,0) = etam[0];
            for (int k = 1; k < nlev; ++k) hde(i,j,k) = etam[k] - etam[k-1];
            hde(i,j,nlev) = 1 - etam[nlev-1];
            // Make sure we have a meaningful test.
            Real minval = 1;
            for (int k = 0; k <= nlev; ++k) minval = std::min(minval, hde(i,j,k));
            if (minval >= deta_tol) err("meaningful test");
          }
        Kokkos::deep_copy(deta, hde);
        run(deta);
        Kokkos::deep_copy(hde, deta);
        for (int i = 0; i < NP; ++i)
          for (int j = 0; j < NP; ++j) {
            const int idx = get_idx(i,j);
            // Min val should be deta_tol.
            Real minval = 1;
            for (int k = 0; k <= nlev; ++k) minval = std::min(minval, hde(i,j,k));
            if (minval != deta_tol) err("min val");
            // Sum of levels should be 1.
            Real sum = 0;
            for (int k = 0; k <= nlev; ++k) sum += hde(i,j,k);
            if (std::abs(sum - 1) > tol) err("sum 1");
            // Only two deltas should be affected.
            Real maxdiff = 0;
            for (int k = 0; k <= nlev; ++k) {
              const auto diff = std::abs(hde(i,j,k) - hder(k));
              if (k == idx or k == idx+1) {
                if (diff <= deta_tol) err("2 deltas a");
              } else {
                maxdiff = std::max(maxdiff, diff);
              }
            }
            if (maxdiff > tol) err("2 deltas b");
          }
      }
    }

    { // Test generally (and highly) perturbed levels.
      const auto hde = Kokkos::create_mirror_view(deta);
      for (int i = 0; i < NP; ++i)
        for (int j = 0; j < NP; ++j) {
          Real sum = 0;
          for (int k = 0; k <= nlev; ++k) {
            hde(i,j,k) = td.urand(-0.5, 0.5);
            sum += hde(i,j,k);
          }
          // Make the column sum to 0.2 for safety in the next step.
          const Real colsum = 0.2;
          for (int k = 0; k <= nlev; ++k) hde(i,j,k) += (colsum - sum)/(nlev+1);
          for (int k = 0; k <= nlev; ++k) hde(i,j,k) /= colsum;
          sum = 0;
          for (int k = 0; k <= nlev; ++k) sum += hde(i,j,k);
          if (std::abs(sum - 1) > 10*tol) err("general sum 1");
        }
      Kokkos::deep_copy(deta, hde);
      run(deta);
      Kokkos::deep_copy(hde, deta);
      for (int i = 0; i < NP; ++i)
        for (int j = 0; j < NP; ++j) {
          Real sum = 0, minval = 1;
          for (int k = 0; k <= nlev; ++k) sum += hde(i,j,k);
          for (int k = 0; k <= nlev; ++k) minval = std::min(minval, hde(i,j,k));
          if (std::abs(sum - 1) > 1e3*td.eps) ++nerr;
          if (minval != deta_tol) err("general minval");
        }
    }
  }
  
  return nerr;
}

struct HybridLevels {
  Real ps0, a_eta, b_eta;
  std::vector<Real> ai, dai, bi, dbi, am, bm, etai, etam, detai, detam;
};

// Follow DCMIP2012 3D tracer transport specification for a, b, eta.
void fill (HybridLevels& h, const int n) {
  h.ai.resize(n+1); h.bi.resize(n+1);
  h.am.resize(n  ); h.bm.resize(n  );
  h.etai.resize(n+1); h.etam.resize(n);

  const auto Rd = PhysicalConstants::Rgas;
  const auto T0 = 300; // K
  const auto p0 = PhysicalConstants::p0;
  const auto g = PhysicalConstants::g;
  const Real ztop = 12e3; // m

  h.ps0 = p0;

  const auto calc_pressure = [&] (const Real z) {
    return p0*std::exp(-g*z/(Rd*T0));
  };

  const Real eta_top = calc_pressure(ztop)/p0;
  assert(eta_top > 0);
  for (int i = 0; i <= n; ++i) {
    const auto z = (Real(n - i)/n)*ztop;
    h.etai[i] = calc_pressure(z)/p0;
    h.bi[i] = i == 0 ? 0 : (h.etai[i] - eta_top)/(1 - eta_top);
    h.ai[i] = h.etai[i] - h.bi[i];
    assert(i == 0 or h.etai[i] > h.etai[i-1]);
  }
  assert(h.bi  [0] == 0); // Real(n - i)/n is exactly 1, so exact = holds
  assert(h.bi  [n] == 1); // exp(0) is exactly 0, so exact = holds
  assert(h.etai[n] == 1); // same
  // b = (eta - eta_top)/(1 - eta_top) => b_eta = 1/(1 - eta_top)
  // a = eta - b => a_eta = 1 - b_eta = -eta_top/(1 - eta_top)
  // p_eta = a_eta p0 + b_eta ps
  h.b_eta = 1/(1 - eta_top);
  h.a_eta = 1 - h.b_eta;

  const auto tomid = [&] (const std::vector<Real>& in, std::vector<Real>& mi) {
    for (int i = 0; i < n; ++i) mi[i] = (in[i] + in[i+1])/2;
  };
  tomid(h.ai, h.am);
  tomid(h.bi, h.bm);
  tomid(h.etai, h.etam);

  const auto diff = [&] (const std::vector<Real>& ai, std::vector<Real>& dai) {
    dai.resize(n);
    for (int i = 0; i < n; ++i) dai[i] = ai[i+1] - ai[i];
  };
  diff(h.ai, h.dai);
  diff(h.bi, h.dbi);
  diff(h.etai, h.detai);

  h.detam.resize(n+1);
  h.detam[0] = h.etam[0] - h.etai[0];
  for (int i = 1; i < n; ++i) h.detam[i] = h.etam[i] - h.etam[i-1];
  h.detam[n] = h.etai[n] - h.etam[n-1];
}

int test_limit_etai (TestData& td) {
  int nerr = 0;
  const Real tol = 100*td.eps;

  for (const int nlev : {143, 128, 81}) {
    const Real deta_tol = 1e5*td.eps/nlev;

    ExecView<Real*> hy_etai("hy_etai",nlev+1), detai("detai",nlev);
    ExecView<Real***> wrk1("wrk1",NP,NP,nlev), wrk2("wrk2",NP,NP,nlev);
    ExecView<Real***> etai("etai",NP,NP,nlev);

    HybridLevels h;
    fill(h, nlev);
    todev(h.etai, hy_etai);
    todev(h.detai, detai);

    const auto he = Kokkos::create_mirror_view(etai);

    const auto policy = get_test_team_policy(1, nlev);
    const auto run = [&] () {
      Kokkos::deep_copy(etai, he);
      const auto f = KOKKOS_LAMBDA(const cti::MT& team) {
        KernelVariables kv(team);
        limit_etai(kv, nlev, hy_etai, detai, deta_tol, wrk1, wrk2, etai);
      };
      Kokkos::parallel_for(policy, f);
      Kokkos::fence();
      Kokkos::deep_copy(he, etai);
    };

    fillcols(h.etai.size()-1, h.etai.data(), he);
    // Col 0 should be untouched. Cols 1 and 2 should have very specific
    // changes. etai[0] is intended to be constant, so don't perturb that entry.
    const int col1_idx = std::max(1, static_cast<int>(0.25*nlev));
    he(0,1,col1_idx) += 0.3;
    const int col2_idx = std::max(1, static_cast<int>(0.8*nlev));
    he(0,2,col2_idx) -= 5.3;
    // The rest of the columns get wild changes.
    for (int idx = 3; idx < NP*NP; ++idx) {
      const int i = idx / NP, j = idx % NP;
      for (int k = 1; k < nlev; ++k)
        he(i,j,k) += td.urand(-1, 1)*(h.etam[k] - h.etam[k-1]);
    }
    run();
    bool ok = true;
    for (int k = 0; k < nlev; ++k)
      if (he(0,0,k) != h.etai[k]) ok = false;
    if (he(0,1,0) != h.etai[0]) ok = false;
    for (int k = 0; k < nlev; ++k) {
      if (k == col1_idx) continue;
      if (std::abs(he(0,1,k) - h.etai[k]) > tol) ok = false;
    }
    if (he(0,2,0) != h.etai[0]) ok = false;
    for (int k = 0; k < nlev; ++k) {
      if (k == col2_idx) continue;
      if (std::abs(he(0,2,k) - h.etai[k]) > tol) ok = false;
    }
    Real mingap = 1;
    for (int i = 0; i < NP; ++i)
      for (int j = 0; j < NP; ++j) {
        if (he(i,j,0) != h.etai[0]) ok = false;
        for (int k = 1; k < nlev; ++k)
          mingap = std::min(mingap, he(i,j,k) - he(i,j,k-1));
        mingap = std::min(mingap, h.etai[nlev] - he(i,j,nlev-1));
      }
    // Test minimum level delta, with room for numerical error.
    if (mingap < 0.8*deta_tol) ok = false;
    if (not ok) ++nerr;
  }
  
  return nerr;
}

int test_eta_interp (TestData& td) {
  int nerr = 0;
  const Real tol = 100*td.eps;

  for (const int nlev : {15, 128, 161}) {
    HybridLevels h;
    fill(h, nlev);

    ExecView<Real*> hy_etai("hy_etai",nlev+1);
    ExecView<Real***> x("x",NP,NP,nlev), y("y",NP,NP,nlev);
    ExecView<Real***> xi("xi",NP,NP,nlev+1), yi("yi",NP,NP,nlev+1);
    ExecView<Real***> xwrk("xwrk",NP,NP,nlev+2), ywrk("ywrk",NP,NP,nlev+2);

    todev(h.etai, hy_etai);

    const auto xh  = Kokkos::create_mirror_view(x );
    const auto yh  = Kokkos::create_mirror_view(y );
    const auto xih = Kokkos::create_mirror_view(xi);
    const auto yih = Kokkos::create_mirror_view(yi);

    const auto policy = get_test_team_policy(1, nlev);
    const auto run_eta = [&] (const int ni) {
      Kokkos::deep_copy(x, xh); Kokkos::deep_copy(y, yh);
      Kokkos::deep_copy(xi, xih);
      const auto f = KOKKOS_LAMBDA(const cti::MT& team) {
        KernelVariables kv(team);
        eta_interp_eta(kv, nlev, hy_etai,
                       nlev, 0, x, getcolc(y,0,0),
                       xwrk, getcol(ywrk,0,0),
                       ni, 0, getcolc(xi,0,0), yi);
      };
      Kokkos::parallel_for(policy, f);
      Kokkos::fence();
      Kokkos::deep_copy(yih, yi);
    };
    const auto run_horiz = [&] () {
      Kokkos::deep_copy(x, xh); Kokkos::deep_copy(y, yh);
      Kokkos::deep_copy(xi, xih);
      const auto f = KOKKOS_LAMBDA(const cti::MT& team) {
        KernelVariables kv(team);
        eta_interp_horiz(kv, nlev, hy_etai,
                         getcolc(x,0,0), y,
                         getcol(xwrk,0,0), ywrk,
                         xi, yi);
      };
      Kokkos::parallel_for(policy, f);
      Kokkos::fence();
      Kokkos::deep_copy(yih, yi);
    };

    std::vector<Real> v;
    const Real d = 1e-6, vlo = h.etai[0]+d, vhi = h.etai[nlev]-d;

    for (const int ni : {int(0.7*nlev), nlev-1, nlev, nlev+1}) {
      make_random_sorted(td, nlev, vlo, vhi, v);
      fillcols(nlev, v.data(), xh);
      fillcols(nlev, v.data(), yh);
      make_random_sorted(td, ni, vlo, vhi, v);
      fillcols(ni, v.data(), xih);
      run_eta(ni);
      bool ok = true;
      for (int i = 0; i < NP; ++i)
        for (int j = 0; j < NP; ++j)
          for (int k = 0; k < ni; ++k)
            if (std::abs(yih(i,j,k) - xih(i,j,k)) > tol)
              ok = false;
      if (not ok) ++nerr;
    }

    { // Test exact interp of line in the interior, const interp near the bdys.
      make_random_sorted(td, nlev, vlo+0.05, vhi-0.1, v);
      fillcols(nlev, v.data(), xh);
      for (int i = 0; i < NP; ++i)
        for (int j = 0; j < NP; ++j) {
          for (int k = 0; k < nlev; ++k)
            yh(i,j,k) = i*xh(0,0,k) - j;
          make_random_sorted(td, nlev, vlo, vhi, v);
          for (int k = 0; k < nlev; ++k)
            xih(i,j,k) = v[k];
        }
      run_horiz();
      bool ok = true;
      for (int i = 0; i < NP; ++i)
        for (int j = 0; j < NP; ++j)
          for (int k = 0; k < nlev; ++k) {
            if (xih(i,j,k) < xh(0,0,0)) {
              if (std::abs(yih(i,j,k) - yih(i,j,0)) > tol)
                ok = false;
            } else if (xih(i,j,k) > xh(0,0,nlev-1)) {
              if (std::abs(yih(i,j,k) - yih(i,j,nlev-1)) > tol)
                ok = false;            
            } else {
              if (std::abs(yih(i,j,k) - (i*xih(i,j,k) - j)) > tol)
                ok = false;
            }
          }
      if (not ok) ++nerr;
    }
  }
  
  return nerr;
}

int test_eta_to_dp (TestData& td) {
  int nerr = 0;
  const Real tol = 100*td.eps;

  for (const int nlev : {143, 128, 81}) {
    const auto err = [&] (const char* lbl) {
      ++nerr;
      printf("test_eta_to_dp nlev %d: %s\n", nlev, lbl);
    };

    HybridLevels h;
    fill(h, nlev);

    ExecView<Real*> hy_bi("hy_bi",nlev+1), hy_etai("hy_etai",nlev+1);
    ExecView<Real***> etai("etai",NP,NP,nlev+1), wrk("wrk",NP,NP,nlev+1);
    ExecView<Real***> dp("dp",NP,NP,nlev);
    ExecView<Real[NP][NP]> ps("ps");
    const Real hy_ps0 = h.ps0;

    todev(h.bi, hy_bi);
    todev(h.etai, hy_etai);

    const auto psm = Kokkos::create_mirror_view(ps);
    HostView<Real***> dp1("dp1",NP,NP,nlev);
    for (int i = 0; i < NP; ++i)
      for (int j = 0; j < NP; ++j)
        psm(i,j) = (1 + 0.1*td.urand(-1, 1))*h.ps0;
    Kokkos::deep_copy(ps, psm);

    const auto policy = get_test_team_policy(1, nlev);
    const auto run = [&] () {
      const auto f = KOKKOS_LAMBDA(const cti::MT& team) {
        KernelVariables kv(team);
        eta_to_dp(kv, nlev, hy_ps0, hy_bi, hy_etai, ps, etai, wrk, dp);
      };
      Kokkos::parallel_for(policy, f);
      Kokkos::fence();
    };

    { // Test that for etai_ref we get the same as the usual formula.
      todev(h.etai, etai);
      HostView<Real***> dp1("dp1",NP,NP,nlev);
      Real dp1_max = 0;
      for (int i = 0; i < NP; ++i)
        for (int j = 0; j < NP; ++j)
          for (int k = 0; k < nlev; ++k) {
            dp1(i,j,k) = ((h.ai[k+1] - h.ai[k])*h.ps0 +
                          (h.bi[k+1] - h.bi[k])*psm(i,j));
            dp1_max = std::max(dp1_max, std::abs(dp1(i,j,k)));
          }
      run();
      const auto dph = cti::cmvdc(dp);
      Real err_max = 0;
      for (int i = 0; i < NP; ++i)
        for (int j = 0; j < NP; ++j)
          for (int k = 0; k < nlev; ++k)
            err_max = std::max(err_max, std::abs(dph(i,j,k) - dp1(i,j,k)));
      if (err_max > tol*dp1_max) err("t1");
    }

    { // Test that sum(dp) = ps for random input etai.
      std::vector<Real> etai_r;
      make_random_sorted(td, nlev+1, h.etai[0], h.etai[nlev], etai_r);
      todev(etai_r, etai);
      run();
      const auto dph1 = cti::cmvdc(dp);
      for (int i = 0; i < NP; ++i)
        for (int j = 0; j < NP; ++j) {
          Real ps = h.ai[0]*h.ps0;
          for (int k = 0; k < nlev; ++k)
            ps += dph1(i,j,k);
          if (std::abs(ps - psm(i,j)) > tol*psm(i,j)) err("t2");
        }    
      // Test that values on input don't affect solution.
      Kokkos::deep_copy(wrk, 0);
      Kokkos::deep_copy(dp, 0);
      run();
      const auto dph2 = cti::cmvdc(dp);
      bool alleq = true;
      for (int i = 0; i < NP; ++i)
        for (int j = 0; j < NP; ++j)
          for (int k = 0; k < nlev; ++k)
            if (dph2(i,j,k) != dph1(i,j,k))
              alleq = false;
      if (not alleq) err("t3");
    }
  }

  return nerr;
}

struct Snapshots {
  const Real alpha[2];
  const ExecViewUnmanaged<Scalar***> dps[2];

  struct Element {
    const ExecViewUnmanaged<Scalar***> dps[2];

    KOKKOS_INLINE_FUNCTION
    Element (const Snapshots& s, const int ie)
      : dps{s.dps[0], s.dps[1]}
    {}

    KOKKOS_INLINE_FUNCTION
    Real get_dp_real(const int t, const int i, const int j, const int k) const {
      return dps[t](i,j, k / VECTOR_SIZE)[k % VECTOR_SIZE];
    }
  };

  Snapshots (const Real as[2], const ExecViewUnmanaged<Scalar***>& dp1,
             const ExecViewUnmanaged<Scalar***>& dp2)
    : alpha{as[0], as[1]}, dps{dp1, dp2}
  {}

  KOKKOS_INLINE_FUNCTION Real get_alpha (const int t) const { return alpha[t]; }

  KOKKOS_INLINE_FUNCTION Element get_element (const int ie) const {
    return Element(*this, ie);
  }
};

int test_calc_ps (TestData& td) {
  int nerr = 0;
  const Real tol = 100*td.eps;

  for (const int nlev : {15, 128, 161}) {
    HybridLevels h;
    fill(h, nlev);
    const auto ps0 = h.ps0, hyai0 = h.ai[0];

    ElData dp1("dp1", nlev), dp2("dp2", nlev);
    for (int i = 0; i < NP; ++i)
      for (int j = 0; j < NP; ++j)
        for (int k = 0; k < nlev; ++k) {
          dp1.r(i,j,k) = td.urand(0, 1000);
          dp2.r(i,j,k) = td.urand(0, 1000);
        }
    dp1.h2d();
    dp2.h2d();

    const Real alpha[] = {td.urand(0,1), td.urand(0,1)};

    ExecView<Real[NP][NP]> ps("ps");
    ExecView<Real[2][NP][NP]> ps2("ps2");
    const auto policy = get_test_team_policy(1, nlev);
    Snapshots snaps(alpha, dp1.d, dp2.d);
    const auto f = KOKKOS_LAMBDA(const cti::MT& team) {
      KernelVariables kv(team);
      calc_ps(kv, nlev, ps0, hyai0, snaps, ps2);
      calc_ps(kv, nlev, ps0, hyai0, dp1.d, ps);
    };
    Kokkos::parallel_for(policy, f);
    Kokkos::fence();

    const auto ps_h = cti::cmvdc(ps);
    const auto ps2_h = cti::cmvdc(ps2);
    for (int i = 0; i < NP; ++i)
      for (int j = 0; j < NP; ++j) {
        {
          Real ps = h.ai[0]*h.ps0;
          for (int k = 0; k < nlev; ++k)
            ps += dp1.r(i,j,k);
          if (std::abs(ps_h(i,j) - ps) > tol*ps) ++nerr;
        }
        for (int t = 0; t < 2; ++t) {
          Real ps = h.ai[0]*h.ps0;
          for (int k = 0; k < nlev; ++k)
            ps += (1 - alpha[t])*dp1.r(i,j,k) + alpha[t]*dp2.r(i,j,k);
          if (std::abs(ps2_h(t,i,j) - ps) > tol*ps) ++nerr;
        }
      }
  }

  return nerr;
}

// Transform eta_dot_dpdn at interfaces to eta_dot at midpoints using the
// formula
//     eta_dot = eta_dot_dpdn/(A_eta p0 + B_eta ps)
//            a= eta_dot_dpdn diff(eta)/(diff(A) p0 + diff(B) ps).
// I'm keeping this because it's unit tested and might be used in the future.
// But it's not currently used, so it's in this unit-test file.
KOKKOS_FUNCTION void calc_etadotmid_from_etadotdpdnint (
  const KernelVariables& kv, const int nlev,
  const Real& ps0, const CSnV& hydai, const CSnV& hydbi,
  const CSnV& hydetai, const CRelV& ps, const SelnV& wrk,
  //  in: eta_dot_dpdn at interfaces
  // out: eta_dot at midpoints, final slot unused
  const SelnV& ed)
{
  assert(calc_nscal(hydai.extent_int(0)) >= nlev);
  assert(calc_nscal(hydbi.extent_int(0)) >= nlev);
  assert(calc_nscal(hydetai.extent_int(0)) >= nlev);
  assert_eln(wrk, nlev+1);
  assert_eln(ed, nlev+1);
  const auto& edd_mid = wrk;
  {
    const CRelnV edd(elp2r(ed));
    const RelnV tmp(elp2r(wrk));
    const auto f = [&] (const int i, const int j, const int k) {
      tmp(i,j,k) = (edd(i,j,k) + edd(i,j,k+1))/2;
    };
    cti::loop_ijk(nlev, kv, f);
  }
  kv.team_barrier();
  {
    const auto f = [&] (const int i, const int j, const int kp) {
      ed(i,j,kp) = (edd_mid(i,j,kp)
                    * hydetai(kp)
                    / (hydai(kp)*ps0 + hydbi(kp)*ps(i,j)));
    };
    cti::loop_ijk(calc_npack(nlev), kv, f);
  }
}

int test_calc_etadot_from_etadotdpdnint (TestData& td) {
  int nerr = 0;
  const Real tol = 100*td.eps;

  for (const int nlev : {143, 128, 81}) {
    HybridLevels h;
    fill(h, nlev);

    // Test function:
    //     eta_dot_dpdn(eta) = c eta + d.
    // Then
    //     eta_dot = eta_dot_dpdn(eta)/dpdn(eta)
    //             = (c eta + d)/(a_eta p0 + b_eta ps).
    // Since a_eta, b_eta are constants independent of eta in this test, eta_dot
    // is then also a linear function of eta. Thus, we can test for exact
    // agreement with the true solution.

    ColData hydai("hydai",nlev), hydbi("hydbi",nlev), hydetai("hydetai",nlev);
    ElData wrk("wrk",nlev+1), ed("ed",nlev+1);
    ExecView<Real[NP][NP]> ps("ps");
    const Real ps0 = h.ps0;
    ExecView<Scalar*> db_deta_i("db_deta_i", calc_npack(nlev+1));
    Kokkos::deep_copy(db_deta_i, h.b_eta);

    for (int trial = 0; trial < 2; ++trial) {
      const auto ps_m = Kokkos::create_mirror_view(ps);
      for (int i = 0; i < NP; ++i)
        for (int j = 0; j < NP; ++j) {
          ps_m(i,j) = td.urand(0.5, 1.2)*ps0;
          for (int k = 0; k < nlev; ++k) {
            hydai.r[k] = h.dai[k];
            hydbi.r[k] = h.dbi[k];
            hydetai.r[k] = h.detai[k];
          }
          for (int k = 0; k <= nlev; ++k)
            ed.r(i,j,k) = (i-j)*h.etai[k] + 0.3;
          if (trial == 1) {
            ed.r(i,j,0) = 0;
            ed.r(i,j,nlev) = 0;
          }
        }
      Kokkos::deep_copy(ps, ps_m);
      hydai.h2d(); hydbi.h2d(); hydetai.h2d();
      ed.h2d();

      const auto policy = get_test_team_policy(1, nlev);
      const auto fmid = KOKKOS_LAMBDA(const cti::MT& team) {
        KernelVariables kv(team);
        calc_etadotmid_from_etadotdpdnint(
          kv, nlev, ps0, hydai.d, hydbi.d, hydetai.d, ps, wrk.d, ed.d);
      };
      const auto fint = KOKKOS_LAMBDA(const cti::MT& team) {
        KernelVariables kv(team);
        calc_etadotint_from_etadotdpdnint(
          kv, nlev, ps0, db_deta_i, ps, ed.d);
      };
      if (trial == 0)
        Kokkos::parallel_for(policy, fmid);
      else
        Kokkos::parallel_for(policy, fint);
      Kokkos::fence();
      ed.d2h();

      for (int i = 0; i < NP; ++i)
        for (int j = 0; j < NP; ++j) {
          const auto den = h.a_eta*h.ps0 + h.b_eta*ps_m(i,j);
          for (int k = 0; k <= nlev; ++k) {
            if (trial == 0 and k == nlev) continue;
            if (trial == 1 and (k == 0 or k == nlev)) {
              if (ed.r(i,j,k) != 0) ++nerr;
              continue;
            }
            const auto eta = trial == 0 ? h.etam[k] : h.etai[k];
            const auto ed_true = ((i-j)*eta + 0.3)/den;
            if (std::abs(ed.r(i,j,k) - ed_true) > tol*(10/den)) ++nerr;
          }
        }
    }
  }

  return nerr;
}

int test1_init_velocity_record (
  const int dtf, const int drf, const int nsub, const int nvel)
{
  const auto eps = std::numeric_limits<Real>::epsilon();
  int e = 0;

  if (dtf % drf != 0) {
    printf("Testing erro: dtf %% drf == 0 is required: %d %d\n", dtf, drf);
    ++e; 
  }

  const CTI::VelocityRecord v(dtf, drf, nsub, nvel);
  if (v.dtf() != dtf) ++e;
  if (v.nsub() != nsub) ++e;

  // Check that t_vel is monotonically increasing.
  for (int n = 1; n < v.nvel(); ++n)
    if (v.t_vel(n) <= v.t_vel(n-1))
      ++e;

  // Check that obs_slots does not reference end points. This should not happen
  // b/c nvel <= navail and observations are uniformly spaced.
  for (int n = 0; n < dtf; ++n)
    for (int i = 0; i < 2; ++i)
      if (v.obs_slots(n,i) == 0 or v.obs_slots(n,i) == dtf)
        ++e;

  // Check that weights sum to 1.
  std::vector<Real> ys(dtf);
  for (int n = 0; n < dtf; ++n)
    for (int i = 0; i < 2; ++i)
      if (v.obs_slots(n,i) >= 0)
        ys[v.obs_slots(n,i)] += v.obs_wts(n,i);
  for (int i = 1; i < v.nvel()-1; ++i)
    if (std::abs(ys[i] - 1) > 1e3*eps)
      ++e;

  // Test for exact interp of an affine function.
  const auto tfn = [] (const Real x) { return 7.1*x - 11.5; };
  //   Observe data forward in time.
  Real endslots[2];
  endslots[0] = tfn(0);
  endslots[1] = tfn(dtf);
  ys[0] = -1000; // unused
  for (int i = 1; i < dtf; ++i) ys[i] = 0;
  for (int n = 0; n < dtf; ++n) {
    if ((n+1) % drf != 0) continue;
    const Real y = tfn(n+1);
    for (int i = 0; i < 2; ++i) {
      if (v.obs_slots(n,i) < 0) continue;
      ys[v.obs_slots(n,i)] += v.obs_wts(n,i)*y;
    }
  }
  //   Use the data backward in time.
  for (int n = 0; n < nsub; ++n) {
    // Each segment orders the data forward in time. Thus, data are always
    // ordered forward in time but used backward.
    Real xsup[2], ysup[2];
    for (int i = 0; i < 2; ++i) {
      int k = nsub - (n+1) + i;
      xsup[i] = (k*v.t_vel(v.nvel()-1))/nsub;
      k = v.run_step(n);
      const Real
        y0 = k == 1 ? endslots[0] : ys[k-1],
        y1 = k == v.nvel()-1 ? endslots[1] : ys[k];
      ysup[i] = (((v.t_vel(k) - xsup[i])*y0 + (xsup[i] - v.t_vel(k-1))*y1) /
                 (v.t_vel(k) - v.t_vel(k-1)));
    }
    for (int i = 0; i <= 10; ++i) {
      const Real
        a = Real(i)/10,
        x = (1-a)*xsup[0] + a*xsup[1],
        y = (1-a)*ysup[0] + a*ysup[1];
      if (std::abs(y - tfn(x)) > 1e3*eps) {
        printf("n %d i %2d x %7.3f y %7.3f t %7.3f\n", n, i, x, y, tfn(x));
        ++e;
      }
    }
  }
  
  if (e) {
    printf("ERROR e %d\n", e);
    printf("dtf %d drf %d nsub %d nvel %d v.nvel %d\n",
           dtf, drf, nsub, nvel, v.nvel());
    printf("  t_vel:");
    for (int i = 0; i < v.nvel(); ++i) printf(" %1.3f", v.t_vel(i));
    printf("\n  obs:\n");
    for (int n = 0; n <= dtf; ++n)
      printf("    %2d %2d %2d %1.3f %1.3f\n",
             n, v.obs_slots(n,0), v.obs_slots(n,1), v.obs_wts(n,0),
             v.obs_wts(n,1));
    printf("  run_step:\n");
    for (int n = 0; n <= nsub; ++n) printf("    %2d %2d\n", n, v.run_step(n));
  }

  return e;
}

int test_init_velocity_record (TestData& td) {
  int dtf, drf, nsub, nvel, error;

  error = 0;

  const auto f = [&] () {
    const int e = test1_init_velocity_record(dtf, drf, nsub, nvel);
    if (e > 0) error = 1;
  };

  error = 0;
  dtf = 6;
  drf = 2;
  nsub = 3;
  nvel = 4;
  f();
  nvel = 3;
  f();
  drf = 3;
  nvel = 6;
  f();
  drf = 1;
  nsub = 5;
  f();
  dtf = 12;
  drf = 2;
  nsub = 3;
  nvel = -1;
  f();
  nsub = 5;
  nvel = 5;
  f();
  dtf = 27;
  drf = 3;
  nsub = 51;
  nvel = 99;
  f();

  return error;
}

} // namespace anon

#define comunittest(f) do {                     \
    ne = f(td);                                 \
    if (ne) printf(#f " ne %d\n", ne);          \
    nerr += ne;                                 \
  } while (0)

int ComposeTransportImpl::run_enhanced_trajectory_unit_tests () {
  int nerr = 0, ne;
  TestData td(*this);
  comunittest(test_deriv);
  comunittest(test_find_support);
  comunittest(test_linterp);
  comunittest(test_eta_interp);
  comunittest(test_eta_to_dp);
  comunittest(test_deta_caas);
  comunittest(test_limit_etai);
  comunittest(test_calc_ps);
  comunittest(test_calc_etadot_from_etadotdpdnint);
  comunittest(test_init_velocity_record);
  return nerr;
}

#undef comunittest

} // namespace Homme

#endif
