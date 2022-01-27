#include <catch2/catch.hpp>

#include "DirkFunctorImpl.hpp"

#include <random>

#include "Types.hpp"
#include "Context.hpp"
#include "mpi/Comm.hpp"
#include "mpi/Connectivity.hpp"
#include "SimulationParams.hpp"
#include "Elements.hpp"
#include "PhysicalConstants.hpp"

#include "utilities/TestUtils.hpp"
#include "utilities/SyncUtils.hpp"
#include "utilities/ViewUtils.hpp"

using namespace Homme;

extern "C" {
  void init_dirk_f90(int ne, const Real* hyai, const Real* hybi, const Real* hyam,
                     const Real* hybm, Real ps0);
  void cleanup_f90();
  void pnh_and_exner_from_eos_f90(const Real* vtheta_dp, const Real* dp3d, const Real* dphi,
                                  Real* pnh, Real* exner, Real* dpnh_dp_i);
  void compute_gwphis_f90(Real* gwh_i, const Real* dp3d, const Real* v, const Real* gradphis);
  void tridiag_diagdom_bfb_a1x1(int n, void* dl, void* d, void* du, void* x, int real_size);
  void c2f_f90(int nelem, int nlev, int nlevp, const Real* dp3d, const Real* w_i, const Real* v,
               const Real* vtheta_dp, const Real* phinh_i, const Real* gradphis, const Real* phis);
  void f2c_f90(int nelem, int nlev, int nlevp, Real* dp3d, Real* w_i, Real* v,
               Real* vtheta_dp, Real* phinh_i, Real* gradphis);
  void get_dirk_jacobian_f90(Real* dl, Real* d, Real* du, Real dt, const Real* dp3d,
                             const Real* dphi, const Real* pnh);
  void compute_stage_value_dirk_f90(int nm1, Real alphadt_nm1, int n0, Real alphadt_n0,
                                    int np1, Real dt2);
  void phi_from_eos_f90(const Real* phis, const Real* vtheta_dp, const Real* dp, Real* phi_i);
} // extern "C"

using FA3 = Kokkos::View<Real*[NP][NP], Kokkos::LayoutRight, Kokkos::HostSpace>;
using FA3d = Kokkos::View<Real***, Kokkos::LayoutRight, Kokkos::HostSpace>;
using FA4 = Kokkos::View<Real*[2][NP][NP], Kokkos::LayoutRight, Kokkos::HostSpace>;
using dfi = DirkFunctorImpl;

template <typename V>
decltype(Kokkos::create_mirror_view(V())) cmvdc (const V& v) {
  const auto h = Kokkos::create_mirror_view(v);
  deep_copy(h, v);
  return h;
}

class Random {
  using rngalg = std::mt19937_64;
  using rpdf = std::uniform_real_distribution<Real>;
  using ipdf = std::uniform_int_distribution<int>;
  std::random_device rd;
  unsigned int seed;
  rngalg engine;
public:
  Random (unsigned int seed_ = Catch::rngSeed()) : seed(seed_ == 0 ? rd() : seed_), engine(seed) {}
  unsigned int gen_seed () { return seed; }
  Real urrng (const Real lo = 0, const Real hi = 1) { return rpdf(lo, hi)(engine); }
  int  uirng (const int lo, const int hi) { return ipdf(lo, hi)(engine); }
};

struct Session {
  HybridVCoord h;
  Random r;
  Elements e;
  const int ne = 2;
  int nelemd;

  //Session () : r(269041989) {}

  void init () {
    printf("seed %u\n", r.gen_seed());
    auto& c = Context::singleton();
    c.create<HybridVCoord>().random_init(r.gen_seed());
    h = c.get<HybridVCoord>();

    const auto hyai = cmvdc(h.hybrid_ai);
    const auto hybi = cmvdc(h.hybrid_bi);
    const auto hyam = cmvdc(h.hybrid_am);
    const auto hybm = cmvdc(h.hybrid_bm);

    init_dirk_f90(ne, hyai.data(), hybi.data(), &hyam(0)[0], &hybm(0)[0], h.ps0);

    nelemd = c.get<Connectivity>().get_num_local_elements();
    e = c.create<Elements>();
  }

  void cleanup () {
    cleanup_f90();
    auto& c = Context::singleton();
    c.finalize_singleton();
  }

  static Session& singleton () {
    if ( ! s_session) {
      s_session = std::make_shared<Session>();
      s_session->init();
    }
    return *s_session;
  }

  // Call only in last line of last TEST_CASE.
  static void delete_singleton () {
    if (s_session) s_session->cleanup();
    s_session = nullptr;
  }

private:
  static std::shared_ptr<Session> s_session;
};

std::shared_ptr<Session> Session::s_session;

static bool almost_equal (const Real& a, const Real& b,
                          const Real tol = 0) {
  const auto re = std::abs(a-b)/(1 + std::abs(a));
  const bool good = re <= tol;
  if ( ! good)
    printf("equal: a,b = %23.16e %23.16e re = %23.16e tol %9.2e\n",
           a, b, re, tol);
  return good;
}

static bool equal (const Real& a, const Real& b,
                   // Used only if not defined HOMMEXX_BFB_TESTING.
                   const Real tol = 0) {
#ifdef HOMMEXX_BFB_TESTING
  if (a != b)
    printf("equal: a,b = %23.16e %23.16e re = %23.16e\n",
           a, b, std::abs((a-b)/a));
  return a == b;
#else
  return almost_equal(a, b, tol);
#endif
}

template <typename V>
void fill (Random& r, const V& a, const Real scale = 1,
           typename std::enable_if<V::rank == 3>::type* = 0) {
  const auto am = cmvdc(a);
  for (int i = 0; i < a.extent_int(0); ++i)
    for (int j = 0; j < a.extent_int(1); ++j)
      for (int k = 0; k < a.extent_int(2); ++k)
        for (int s = 0; s < dfi::packn; ++s)
          am(i,j,k)[s] = scale*r.urrng(-1,1); 
  deep_copy(a, am);
}

template <typename V>
void fill (Random& r, const V& a,
           typename std::enable_if<V::rank == 4>::type* = 0) {
  const auto am = cmvdc(a);
  for (int i = 0; i < a.extent_int(0); ++i)
    for (int j = 0; j < a.extent_int(1); ++j)
      for (int k = 0; k < a.extent_int(2); ++k)
        for (int l = 0; l < a.extent_int(3); ++l)
          for (int s = 0; s < dfi::packn; ++s)
            am(i,j,k,l)[s] = r.urrng(-1,1); 
  deep_copy(a, am);
}

template <typename V>
void fill_inc (Random& r, const int nlev, const Real top, const Real bottom,
               const V& a, typename std::enable_if<V::rank == 3>::type* = 0) {
  const auto am = cmvdc(a);
  for (int i = 0; i < a.extent_int(0); ++i)
    for (int j = 0; j < a.extent_int(1); ++j) {
      Real* const p = &am(i,j,0)[0];
      Real sum = 0;
      for (int k = 0; k < nlev; ++k) {
        p[k] = r.urrng(); 
        sum += p[k];
      }
      const Real f = (bottom - top)/sum;
      for (int k = 0; k < nlev; ++k)
        p[k] *= f;
    }
  deep_copy(a, am);
}

template <typename V>
void fill_mid (Random& r, const int nlev, const Real top, const Real bottom,
               const V& a, typename std::enable_if<V::rank == 3>::type* = 0) {
  const auto am = cmvdc(a);
  fill_inc(r, nlev, top, bottom, am);
  for (int i = 0; i < a.extent_int(0); ++i)
    for (int j = 0; j < a.extent_int(1); ++j) {
      Real* const p = &am(i,j,0)[0];
      auto delta = p[0];
      p[0] = top + 0.5*delta;
      for (int k = 1; k < nlev; ++k) {
        const auto tmp = p[k];
        p[k] = p[k-1] + 0.5*(delta + p[k]);
        delta = tmp;
      }
    }
  deep_copy(a, am);
}

template <typename V>
void require_all_equal (const int nlev, const V& a, const V& b,
                        typename std::enable_if<V::rank == 3>::type* = 0) {
  const auto am = cmvdc(a);
  const auto bm = cmvdc(b);
  for (int i = 0; i < a.extent_int(0); ++i)
    for (int j = 0; j < a.extent_int(1); ++j)
      for (int k = 0; k < a.extent_int(2); ++k)
        for (int s = 0; s < dfi::packn; ++s)
          if (dfi::packn*k + s < nlev)
            REQUIRE(am(i,j,k)[s] == bm(i,j,k)[s]);
}

template <typename V>
void c2f (const V& c, const FA3& f) {
  const auto cm = cmvdc(c);
  for (int i = 0; i < c.extent_int(0); ++i)
    for (int j = 0; j < c.extent_int(1); ++j) {
      Real* const p = &cm(i,j,0)[0];
      for (int k = 0; k < f.extent_int(0); ++k)
        f(k,i,j) = p[k];
    }
}

template <typename V>
void f2c (const FA3& f, const V& c) {
  const auto cm = Kokkos::create_mirror_view(c);
  for (int i = 0; i < c.extent_int(0); ++i)
    for (int j = 0; j < c.extent_int(1); ++j) {
      Real* const p = &cm(i,j,0)[0];
      for (int k = 0; k < f.extent_int(0); ++k)
        p[k] = f(k,i,j);
    }
  deep_copy(c, cm);
}

template <typename V>
void c2f (const V& c, const FA4& f) {
  const auto cm = cmvdc(c);
  for (int d = 0; d < c.extent_int(0); ++d)
    for (int i = 0; i < c.extent_int(1); ++i)
      for (int j = 0; j < c.extent_int(2); ++j) {
        Real* const p = &cm(d,i,j,0)[0];
        for (int k = 0; k < f.extent_int(0); ++k)
          f(k,d,i,j) = p[k];
      }
}

static void init (dfi& d, FunctorsBuffersManager& fbm) {
  fbm.request_size(d.requested_buffer_size());
  fbm.allocate();
  d.init_buffers(fbm);
}

TEST_CASE ("dirk_pieces_testing") {
  using Kokkos::parallel_for;
  using Kokkos::fence;

  auto& s = Session::singleton();
  const auto& hvcoord = s.h;
  auto& r = s.r;
  const auto eps = std::numeric_limits<Real>::epsilon();

  const int nelem = 1;
  dfi d1(nelem);
  FunctorsBuffersManager fbm;
  init(d1, fbm);

  SECTION ("transpose") {
    ExecView<Scalar[NP][NP][NUM_LEV  ]> ham0("ham0"), ham1("ham1");
    ExecView<Scalar[NP][NP][NUM_LEV_P]> hai0("hai0"), hai1("hai1");
    const int nlev = NUM_LEV*VECTOR_SIZE - 2; // test with some remainder
    fill(r, ham0);
    fill(r, hai0);
    const auto w = d1.m_work;
    const auto a = dfi::get_work_slot(w, 0, 0);
    const auto f1 = KOKKOS_LAMBDA(const dfi::MT& t) {
      KernelVariables kv(t);
      dfi::transpose(kv, nlev, ham0, a);
      kv.team_barrier();
      dfi::transpose(kv, nlev, a, ham1);
    };
    parallel_for(d1.m_policy, f1); fence();
    require_all_equal(nlev, ham0, ham1);
    const auto f2 = KOKKOS_LAMBDA(const dfi::MT& t) {
      KernelVariables kv(t);
      dfi::transpose(kv, nlev+1, hai0, a);
      kv.team_barrier();
      dfi::transpose(kv, nlev+1, a, hai1);
    };
    parallel_for(d1.m_policy, f2); fence();
    require_all_equal(nlev+1, hai0, hai1);
  }

  SECTION ("jacobian") {
    const int nlev = NUM_PHYSICAL_LEV, np = NP;
    const Real dt = r.urrng(32, 64);

    ExecView<Scalar[NP][NP][NUM_LEV]> dp3d("dp3d"), dphi("dphi"), pnh("pnh");
    fill_inc(r, nlev, 5, 10000, dp3d);
    fill_mid(r, nlev, 5, 10000, pnh);
    fill_inc(r, nlev, 500000, 100, dphi);

    const auto w = d1.m_work;
    const auto
      dp3dw = dfi::get_work_slot(w, 0, 0),
      dphiw = dfi::get_work_slot(w, 0, 1),
      pnhw  = dfi::get_work_slot(w, 0, 2);
    const auto ls = d1.m_ls;
    const auto
      dl = dfi::get_ls_slot(ls, 0, 0),
      d  = dfi::get_ls_slot(ls, 0, 1),
      du = dfi::get_ls_slot(ls, 0, 2);

    const auto f1 = KOKKOS_LAMBDA(const dfi::MT& t) {
      KernelVariables kv(t);
      dfi::transpose(kv, nlev, dp3d, dp3dw);
      dfi::transpose(kv, nlev, dphi, dphiw);
      dfi::transpose(kv, nlev, pnh, pnhw);
      kv.team_barrier();
      dfi::calc_jacobian(kv, dt, dp3dw, dphiw, pnhw, dl, d, du);
    };
    parallel_for(d1.m_policy, f1); fence();

    const auto dlm = cmvdc(dl), dm = cmvdc(d), dum = cmvdc(du);

    FA3 dp3df("dp3df", nlev), dphif("dphif", nlev), pnhf("pnhf", nlev);
    FA3 dlf("dlf", nlev-1), df("df", nlev), duf("duf", nlev-1);
    c2f(dp3d, dp3df);
    c2f(dphi, dphif);
    c2f(pnh, pnhf);
    get_dirk_jacobian_f90(dlf.data(), df.data(), duf.data(), dt,
                          dp3df.data(), dphif.data(), pnhf.data());

    const auto eq = [&] (const decltype(dm)& dc, const int os, const FA3& df,
                         const int n) {
      for (int i = 0; i < np; ++i)
        for (int j = 0; j < np; ++j) {
          const auto
          idx = np*i + j,
          pi = idx / dfi::packn,
          si = idx % dfi::packn;
          for (int k = 0; k < n; ++k)
            REQUIRE(equal(dc(k+os,pi)[si], df(k,i,j), 1e3*eps));
        }
    };
    eq(dlm, 1, dlf, nlev-1); eq(dm, 0, df, nlev); eq(dum, 0, duf, nlev-1);

    // Verify diagonal dominance.
    for (int i = 0; i < np; ++i)
      for (int j = 0; j < np; ++j) {
        const auto
          idx = np*i + j,
          pi = idx / dfi::packn,
          si = idx % dfi::packn;
        for (int k = 0; k < nlev; ++k) {
          const auto v = -((k > 0      ? dlm(k,pi)[si] : 0) +
                           (k < nlev-1 ? dum(k,pi)[si] : 0));
          REQUIRE(v >= 0);
          REQUIRE(dm(k,pi)[si] > v);
        }
      }

    // Test solvers.
    dfi::LinearSystem ls1("w1", 1);
    const auto
      x1  = dfi::get_ls_slot(ls , 0, 3),
      x2  = dfi::get_ls_slot(ls1, 0, 0),
      dl2 = dfi::get_ls_slot(ls1, 0, 1),
      d2  = dfi::get_ls_slot(ls1, 0, 2),
      du2 = dfi::get_ls_slot(ls1, 0, 3);
    FA3d x3("x3", np, np, nlev);
    const auto x1m = create_mirror_view(x1);
    // Fill RHS with random numbers.
    for (int k = 0; k < nlev; ++k)
      for (int i = 0; i < np; ++i)
        for (int j = 0; j < np; ++j) {
          const auto
            idx = i*np + j,
            pi = idx / dfi::packn,
            si = idx % dfi::packn;
          const auto v = r.urrng(-1,1);
          x1m(k,pi)[si] = v;
          x3(i,j,k) = v;
        }
    deep_copy(x1, x1m);
    deep_copy(x2, x1);
    deep_copy(dl2, dl); deep_copy(d2, d); deep_copy(du2, du);
    const auto f2 = KOKKOS_LAMBDA(const dfi::MT& t) {
      KernelVariables kv(t);
      dfi::solve   (kv, dl , d , du , x1);
      dfi::solvebfb(kv, dl2, d2, du2, x2);
    };
    parallel_for(d1.m_policy, f2); fence();
    // Test that BFB and non-BFB solvers give nearly the same answers.
    deep_copy(x1m, x1);
    const auto x2m = cmvdc(x2);
    for (int i = 0; i < np; ++i)
      for (int j = 0; j < np; ++j) {
        const auto
          idx = np*i + j,
          pi = idx / dfi::packn,
          si = idx % dfi::packn;
        for (int k = 0; k < nlev; ++k) {
          REQUIRE(almost_equal(x1m(k,pi)[si], x2m(k,pi)[si], 1e4*eps));
        }
      }
    // Test BFB F90 and C++.
    for (int i = 0; i < np; ++i)
      for (int j = 0; j < np; ++j) {
        static const int n = NUM_PHYSICAL_LEV;
        Real dl[n], d[n], du[n];
        for (int k = 0; k < nlev-1; ++k) dl[k] = dlf(k,i,j);
        for (int k = 0; k < nlev  ; ++k) d [k] = df (k,i,j);
        for (int k = 0; k < nlev-1; ++k) du[k] = duf(k,i,j);
        tridiag_diagdom_bfb_a1x1(nlev, dl, d, du, &x3(i,j,0), sizeof(Real));
        const auto
          idx = np*i + j,
          pi = idx / dfi::packn,
          si = idx % dfi::packn;
        for (int k = 0; k < nlev; ++k)
          REQUIRE(equal(x2m(k,pi)[si], x3(i,j,k), 1e3*eps));
      }
  }

  SECTION ("pnh_and_exner_from_eos") {
    const int nlev = NUM_PHYSICAL_LEV, np = NP;

    ExecView<Scalar[NP][NP][NUM_LEV]> dp3d("dp3d"), dphi("dphi"), vtheta_dp("vtheta_dp");
    fill_inc(r, nlev, 5, 10000, dp3d);
    fill_inc(r, nlev, 500000, 100, dphi);
    fill_inc(r, nlev, 5*300, 10000*300, vtheta_dp);

    const auto w = d1.m_work;
    const auto
      dp3dw = dfi::get_work_slot(w, 0, 0),
      dphiw = dfi::get_work_slot(w, 0, 1),
      vtheta_dpw = dfi::get_work_slot(w, 0, 2),
      pnhw = dfi::get_work_slot(w, 0, 3),
      exnerw = dfi::get_work_slot(w, 0, 4),
      dpnh_dp_iw = dfi::get_work_slot(w, 0, 5);

    const auto f = KOKKOS_LAMBDA(const dfi::MT& t) {
      KernelVariables kv(t);
      dfi::transpose(kv, nlev, vtheta_dp, vtheta_dpw);
      dfi::transpose(kv, nlev, dp3d, dp3dw);
      dfi::transpose(kv, nlev, dphi, dphiw);
      kv.team_barrier();
      dfi::pnh_and_exner_from_eos(kv, hvcoord, vtheta_dpw, dp3dw, dphiw,
                                  pnhw, exnerw, dpnh_dp_iw);
    };
    parallel_for(d1.m_policy, f); fence();
    
    const auto dpnh_dp_im = cmvdc(dpnh_dp_iw);

    FA3 dp3df("dp3df", nlev), dphif("dphif", nlev), vtheta_dpf("vtheta_dpf", nlev),
      pnhf("pnhf", nlev), exnerf("exner", nlev), dpnh_dp_if("dpnh_dp_if", nlev+1);
    c2f(dp3d, dp3df);
    c2f(dphi, dphif);
    c2f(vtheta_dp, vtheta_dpf);
    pnh_and_exner_from_eos_f90(vtheta_dpf.data(), dp3df.data(), dphif.data(),
                               pnhf.data(), exnerf.data(), dpnh_dp_if.data());

    for (int i = 0; i < np; ++i)
      for (int j = 0; j < np; ++j) {
        const auto
          idx = np*i + j,
          pi = idx / dfi::packn,
          si = idx % dfi::packn;
        for (int k = 0; k < nlev; ++k)
          REQUIRE(equal(dpnh_dp_im(k,pi)[si], dpnh_dp_if(k,i,j), 1e2*eps));
      }
  }

  SECTION ("calc_gwphis") {
    const int nlev = NUM_PHYSICAL_LEV, np = NP;

    ExecView<Scalar[NP][NP][NUM_LEV]> dp3d("dp3d");
    fill_inc(r, nlev, 5, 10000, dp3d);
    ExecView<Scalar[2][NP][NP][NUM_LEV]> v("v");
    fill(r, v);
    ExecView<Real[2][NP][NP]> gradphis("gradphis");
    const auto gpm = create_mirror_view(gradphis);
    for (int d = 0; d < 2; ++d)
      for (int i = 0; i < np; ++i)
        for (int j = 0; j < np; ++j)
          gpm(d,i,j) = r.urrng(-1,1);
    deep_copy(gradphis, gpm);
    const auto hybi = hvcoord.hybrid_bi;
    
    const auto w = d1.m_work;
    const auto gwh_i = dfi::get_work_slot(w, 0, 0);

    const auto f = KOKKOS_LAMBDA(const dfi::MT& t) {
      KernelVariables kv(t);
      dfi::calc_gwphis(kv, dp3d, v, gradphis, hybi, gwh_i);
    };
    parallel_for(d1.m_policy, f); fence();
    const auto gwh_im = cmvdc(gwh_i);

    FA3 dp3df("dp3df", nlev), gwh_if("gwh_if", nlev+1);
    c2f(dp3d, dp3df);
    FA4 vf("vf", nlev);
    c2f(v, vf);
    compute_gwphis_f90(gwh_if.data(), dp3df.data(), vf.data(), gpm.data());

    for (int i = 0; i < np; ++i)
      for (int j = 0; j < np; ++j) {
        const auto
          idx = np*i + j,
          pi = idx / dfi::packn,
          si = idx % dfi::packn;
        for (int k = 0; k < nlev; ++k)
          REQUIRE(equal(gwh_im(k,pi)[si], gwh_if(k,i,j), 1e3*eps));
      }
  }

  SECTION ("threshold") {
    const int nlev = dfi::num_phys_lev, nvec = dfi::npack, ilim = nvec*dfi::packn;
    const auto w = d1.m_work;
    const auto
      a = dfi::get_work_slot(w, 0, 0),
      b = dfi::get_work_slot(w, 0, 1),
      wa = dfi::get_work_slot(w, 0, 2),
      wb = dfi::get_work_slot(w, 0, 3);
    const Real threshold = 0.42;

    const auto am = create_mirror_view(a);
    const auto bm = create_mirror_view(b);
    for (int k = 0; k < nlev; ++k) {
      Real* const ap = &am(k,0)[0];
      Real* const bp = &bm(k,0)[0];
      for (int i = 0; i < ilim; ++i)
        ap[i] = bp[i] = r.urrng(-2,0);
      if (k == 1) ap[2] = bp[2] = threshold;
      if (k == 2) ap[3] = bp[3] = threshold + 1e-10;
    }
    deep_copy(a, am);
    deep_copy(b, bm);

    const auto f1 = KOKKOS_LAMBDA(const dfi::MT& t) {
      KernelVariables kv(t);
      dfi::calc_whether_gt_and_set(kv, nlev, nvec, threshold, a, wa);
    };
    parallel_for(d1.m_policy, f1); fence();

    const auto f2 = KOKKOS_LAMBDA(const dfi::MT& t) {
      KernelVariables kv(t);
      dfi::calc_whether_ge(kv, nlev, nvec, threshold, b, wb);
    };
    parallel_for(d1.m_policy, f2); fence();

    deep_copy(am, a);
    deep_copy(bm, b);
    const auto wam = cmvdc(wa);
    const auto wbm = cmvdc(wb);
    Real* wp = &wam(0,0)[0];
    REQUIRE(wp[3] == 1);
    for (int i = 0; i < ilim; ++i) REQUIRE((i == 3 || wp[i] == 0));
    wp = &wbm(0,0)[0];
    REQUIRE(wp[2] == 1);
    REQUIRE(wp[3] == 1);
    for (int i = 0; i < ilim; ++i) REQUIRE((i == 2 || i == 3 || wp[i] == 0));
    int thr_cnt = 0;
    for (int k = 0; k < nlev; ++k) {
      Real* const ap = &am(k,0)[0];
      Real* const bp = &bm(k,0)[0];
      if (k == 2) REQUIRE(ap[3] == threshold);
      for (int i = 0; i < ilim; ++i) {
        if (ap[i] >= threshold) ++thr_cnt;
        if (bp[i] >= threshold) ++thr_cnt;
      }
    }
    REQUIRE(thr_cnt == 4);
  }
}

static void c2f (const Elements& e) {
  const auto dp3d = cmvdc(e.m_state.m_dp3d);
  const auto w_i = cmvdc(e.m_state.m_w_i);
  const auto v = cmvdc(e.m_state.m_v);
  const auto vtheta_dp = cmvdc(e.m_state.m_vtheta_dp);
  const auto phinh_i = cmvdc(e.m_state.m_phinh_i);
  const auto gradphis = cmvdc(e.m_geometry.m_gradphis);
  const auto phis = cmvdc(e.m_geometry.m_phis);
  c2f_f90(e.num_elems(), VECTOR_SIZE*NUM_LEV, VECTOR_SIZE*NUM_LEV_P,
    reinterpret_cast<Real*>(dp3d.data()),
    reinterpret_cast<Real*>(w_i.data()),
    reinterpret_cast<Real*>(v.data()),
    reinterpret_cast<Real*>(vtheta_dp.data()),
    reinterpret_cast<Real*>(phinh_i.data()),
    reinterpret_cast<Real*>(gradphis.data()),
    reinterpret_cast<Real*>(phis.data()));
}

static void f2c (Elements& e) {
  using Kokkos::create_mirror_view;
  using Kokkos::deep_copy;
  const auto dp3d = create_mirror_view(e.m_state.m_dp3d);
  const auto w_i = create_mirror_view(e.m_state.m_w_i);
  const auto v = create_mirror_view(e.m_state.m_v);
  const auto vtheta_dp = create_mirror_view(e.m_state.m_vtheta_dp);
  const auto phinh_i = create_mirror_view(e.m_state.m_phinh_i);
  const auto gradphis = create_mirror_view(e.m_geometry.m_gradphis);
  f2c_f90(e.num_elems(), VECTOR_SIZE*NUM_LEV, VECTOR_SIZE*NUM_LEV_P,
    reinterpret_cast<Real*>(dp3d.data()),
    reinterpret_cast<Real*>(w_i.data()),
    reinterpret_cast<Real*>(v.data()),
    reinterpret_cast<Real*>(vtheta_dp.data()),
    reinterpret_cast<Real*>(phinh_i.data()),
    reinterpret_cast<Real*>(gradphis.data()));
  deep_copy(e.m_state.m_dp3d, dp3d);
  deep_copy(e.m_state.m_w_i, w_i);
  deep_copy(e.m_state.m_v, v); 
  deep_copy(e.m_state.m_vtheta_dp, vtheta_dp );
  deep_copy(e.m_state.m_phinh_i, phinh_i);
  deep_copy(e.m_geometry.m_gradphis, gradphis);
}

static void init_elems (int, int nelemd, Random& r, const HybridVCoord& hvcoord,
                        Elements& e) {
  using Kokkos::create_mirror_view;
  using Kokkos::deep_copy;
  using Kokkos::subview;

  const int nlev = NUM_PHYSICAL_LEV, np = NP;
  const auto all = Kokkos::ALL();

  e.init(nelemd, false, true, PhysicalConstants::rearth0);
  const auto max_pressure = 1000 + hvcoord.ps0;
  auto& geo = e.m_geometry;
  e.m_geometry.randomize(r.gen_seed());
  e.m_state.randomize(r.gen_seed(), max_pressure, hvcoord.ps0, hvcoord.hybrid_ai0,
                      geo.m_phis);
  e.m_derived.randomize(r.gen_seed(), 10);

  // We want mixed signs in w_i, v, and gradphis.
  for (int ie = 0; ie < e.num_elems(); ++ie)
    for (int i = 0; i < 3; ++i) {
      fill(r, subview(e.m_state.m_w_i,ie,i,all,all,all), 5 ); // Pa/s
      fill(r, subview(e.m_state.m_v,ie,i,0,all,all,all), 10); // m/s
      fill(r, subview(e.m_state.m_v,ie,i,1,all,all,all), 10);
    }
  const auto gpm = create_mirror_view(e.m_geometry.m_gradphis);
  for (int ie = 0; ie < e.num_elems(); ++ie)
    for (int d = 0; d < 2; ++d)
      for (int i = 0; i < np; ++i)
        for (int j = 0; j < np; ++j)
          gpm(ie,d,i,j) = r.urrng(-1,1);
  deep_copy(e.m_geometry.m_gradphis, gpm);

  // Make sure dphi <= -g.
  const auto phis = cmvdc(geo.m_phis);
  const auto phinh_i = cmvdc(e.m_state.m_phinh_i);
  for (int ie = 0; ie < e.num_elems(); ++ie)
    for (int t = 0; t < NUM_TIME_LEVELS; ++t)
      for (int i = 0; i < np; ++i)
        for (int j = 0; j < np; ++j) {
          Real* const phi = &phinh_i(ie,t,i,j,0)[0];
          phi[nlev] = phis(ie,i,j);
          for (int k = nlev-1; k >= 0; --k)
            if (phi[k] - phi[k+1] < PhysicalConstants::g)
              for (int k1 = k; k1 >= 0; --k1)
                phi[k1] += PhysicalConstants::g;
        }
  deep_copy(e.m_state.m_phinh_i, phinh_i);
}

TEST_CASE ("dirk_toplevel_testing") {
  using Kokkos::create_mirror_view;
  using Kokkos::parallel_for;
  using Kokkos::fence;
  using Kokkos::subview;
  using Kokkos::deep_copy;

  const int np = NP, n0 = 1, np1 = 2, ne = 2;
  const int nlev = dfi::num_phys_lev, nvec = dfi::npack;
  const auto eps = std::numeric_limits<Real>::epsilon();
  Real dt2 = 0.15;

  auto& s = Session::singleton();
  const auto& hvcoord = s.h;
  auto& r = s.r;
  auto& e = s.e;
  const auto nelemd = s.nelemd;

  DirkFunctorImpl d(nelemd);
  FunctorsBuffersManager fbm;
  init(d, fbm);

  { // Test initial guess function.
    init_elems(ne, nelemd, r, hvcoord, e);
    { // C++ version with DIRK-newton-loop policy.
      const auto e_phis = e.m_geometry.m_phis;
      const auto e_vtheta_dp = e.m_state.m_vtheta_dp;
      const auto e_dp3d = e.m_state.m_dp3d;
      const auto e_phinh_i = e.m_state.m_phinh_i;
      const auto w = d.m_work;
      const auto f = KOKKOS_LAMBDA(const dfi::MT& t) {
        KernelVariables kv(t);
        const auto ie = kv.ie;
        const auto a = Kokkos::ALL();
        const auto
        vtheta_dp = dfi::get_work_slot(w, kv.team_idx, 0),
        dp3d = dfi::get_work_slot(w, kv.team_idx, 1),
        phi_np1 = dfi::get_work_slot(w, kv.team_idx, 2);
        dfi::transpose(kv, nlev, subview(e_vtheta_dp,ie,np1,a,a,a), vtheta_dp);
        dfi::transpose(kv, nlev, subview(e_dp3d,ie,np1,a,a,a), dp3d);
        kv.team_barrier();
        dfi::phi_from_eos(kv, nlev, nvec, hvcoord, subview(e_phis,ie,a,a),
                          vtheta_dp, dp3d, phi_np1);
        kv.team_barrier();
        dfi::transpose(kv, nlev, phi_np1, subview(e_phinh_i,ie,np1,a,a,a));
      };
      parallel_for(d.m_policy, f); fence();
    }
    const auto phic1_m = cmvdc(e.m_state.m_phinh_i);
    decltype(phic1_m) phic1("phic1", phic1_m.extent_int(0));
    deep_copy(phic1, phic1_m);
    { // C++ version with separate dispatch and Hommexx patterns.
      d.run_initial_guess(np1, e, hvcoord);
    }
    const auto phic2_m = cmvdc(e.m_derived.m_divdp_proj);
    decltype(phic2_m) phic2("phic2", phic2_m.extent_int(0));
    deep_copy(phic2, phic2_m);
    { // F90 version.
      const auto e_phis = cmvdc(e.m_geometry.m_phis);
      const auto e_vtheta_dp = cmvdc(e.m_state.m_vtheta_dp);
      const auto e_dp3d = cmvdc(e.m_state.m_dp3d);
      const auto e_phinh_i = cmvdc(e.m_state.m_phinh_i);
      const auto a = Kokkos::ALL();
      for (int ie = 0; ie < nelemd; ++ie) {
        const auto phis = subview(e_phis,ie,a,a);
        const auto vtheta_dp = subview(e_vtheta_dp,ie,np1,a,a,a);
        const auto dp3d = subview(e_dp3d,ie,np1,a,a,a);
        const auto phi_i = subview(e_phinh_i,ie,np1,a,a,a);
        FA3 vtheta_dpf("vtheta_dpf", nlev), dp3df("dp3df", nlev), phif("phif", nlev+1);
        c2f(vtheta_dp, vtheta_dpf); c2f(dp3d, dp3df); c2f(phi_i, phif);
        phi_from_eos_f90(phis.data(), vtheta_dpf.data(), dp3df.data(), phif.data());
        f2c(phif, phi_i);
      }
      deep_copy(e.m_state.m_phinh_i, e_phinh_i);
    }
    const auto phif = cmvdc(e.m_state.m_phinh_i);
    for (int ie = 0; ie < nelemd; ++ie)
      for (int i = 0; i < np; ++i)
        for (int j = 0; j < np; ++j) {
          Real* pf = &phif(ie,np1,i,j,0)[0];
          Real* pc1 = &phic1(ie,np1,i,j,0)[0];
          Real* pc2 = &phic2(ie,i,j,0)[0];
          for (int k = 0; k < nlev; ++k) REQUIRE(equal(pf[k], pc1[k], 1e6*eps));
          for (int k = 0; k < nlev; ++k) REQUIRE(equal(pf[k], pc2[k], 1e6*eps));
        }
  }

  for (Real alphadtwt_nm1 : {0.0, 0.3}) {
    const int nm1 = alphadtwt_nm1 == 0.0 ? -1 : 0;
    for (Real alphadtwt_n0 : {0.0, 0.7}) {
      decltype(ElementsState::m_w_i) w_i("w_i", nelemd),
        w_i1("w_i1", nelemd), w_i2("w_i2", nelemd);
      decltype(ElementsState::m_phinh_i) phinh_i("phinh_i", nelemd),
        phinh_i1("phinh_i1", nelemd), phinh_i2("phinh_i2", nelemd);

      bool good = false;
      for (int trial = 0; trial < 100 /* don't enter an inf loop */; ++trial) {
        init_elems(ne, nelemd, r, hvcoord, e);

        deep_copy(w_i, e.m_state.m_w_i);
        deep_copy(phinh_i, e.m_state.m_phinh_i);

        // Run C++ with BFB solver.
        d.run(nm1, alphadtwt_nm1*dt2, n0, alphadtwt_n0*dt2, np1, dt2,
              e, hvcoord, true /* BFB solver */);
        fence();
        deep_copy(w_i2, e.m_state.m_w_i);
        deep_copy(phinh_i2, e.m_state.m_phinh_i);
        // Restore state.
        deep_copy(e.m_state.m_w_i, w_i);
        deep_copy(e.m_state.m_phinh_i, phinh_i);

        // The F90 code will hard-abort if it can't solve the equations, which
        // interferes with unit testing. To work around this, detect NaNs in the
        // C++ code and reset if any exists. This of course could permit the C++
        // code to NaN when the F90 wouldn't, but BFB testing strongly mitigates
        // that risk. Moreover, being able to generate hard problems but back
        // off it can't be solved (by the specific F90 alg) makes for better
        // testing.
        const auto w2m = cmvdc(w_i2);
        const auto phinh2m = cmvdc(phinh_i2);
        bool ok = true;
        for (int ie = 0; ie < nelemd; ++ie)
          for (int i = 0; i < np; ++i)
            for (int j = 0; j < np; ++j) {
              for (int f = 0; f < 2; ++f) {
                Real* p = f == 0 ? &phinh2m(ie,np1,i,j,0)[0] : &w2m(ie,np1,i,j,0)[0];
                for (int k = 0; k < nlev+1; ++k)
                  if (std::isnan(p[k]) || std::isinf(p[k]))
                    ok = false;
                if (f == 0)
                  for (int k = 0; k < nlev; ++k)
                    if (p[k] <= p[k+1])
                      ok = false;
              }
            }
        if ( ! ok) {
          // Make the problems a little easier.
          const Real f = 0.99;
          dt2 *= f;
          alphadtwt_nm1 *= f;
          alphadtwt_n0 *= f;
          continue;
        }
        good = true;

        // Run C++ with non-BFB solver.
        d.run(nm1, alphadtwt_nm1*dt2, n0, alphadtwt_n0*dt2, np1, dt2,
              e, hvcoord, false /* non-BFB solver */);
        fence();
        deep_copy(w_i1, e.m_state.m_w_i);
        deep_copy(phinh_i1, e.m_state.m_phinh_i);
        // Restore state.
        deep_copy(e.m_state.m_w_i, w_i);
        deep_copy(e.m_state.m_phinh_i, phinh_i);

        break;
      }

      if ( ! good) std::cout << "After 100 attempts, did not generate a good problem.\n";

      const auto w1m = cmvdc(w_i1);
      const auto w2m = cmvdc(w_i2);
      const auto phinh1m = cmvdc(phinh_i1);
      const auto phinh2m = cmvdc(phinh_i2);

      // Test that running with BFB and non-BFB solvers produces similar answers.
      for (int ie = 0; ie < nelemd; ++ie)
        for (int i = 0; i < np; ++i)
          for (int j = 0; j < np; ++j)
            for (int f = 0; f < 2; ++f) {
              Real* p1 = f == 0 ? &w1m(ie,np1,i,j,0)[0] : &phinh1m(ie,np1,i,j,0)[0];
              Real* p2 = f == 0 ? &w2m(ie,np1,i,j,0)[0] : &phinh2m(ie,np1,i,j,0)[0];
              for (int k = 0; k < nlev+1; ++k)
                REQUIRE(almost_equal(p1[k], p2[k], 1e6*eps));
            }

      // Run F90 with BFB solver.
      c2f(e);
      compute_stage_value_dirk_f90(nm1+1, alphadtwt_nm1*dt2, n0+1, alphadtwt_n0*dt2, np1+1, dt2);
      Elements ef90;
      ef90.init(nelemd, false, true, PhysicalConstants::rearth0);
      f2c(ef90);

      const auto phif = cmvdc(ef90.m_state.m_phinh_i);
      const auto wif  = cmvdc(ef90.m_state.m_w_i);

      // Test that C++ and F90 with BFB solvers produce the same answer.
      for (int ie = 0; ie < 1/*nelemd*/; ++ie)
        for (int i = 0; i < np; ++i)
          for (int j = 0; j < np; ++j) {
            for (int f = 0; f < 2; ++f) {
              Real* pf = f == 0 ? &phif   (ie,np1,i,j,0)[0] : &wif(ie,np1,i,j,0)[0];
              Real* pc = f == 0 ? &phinh2m(ie,np1,i,j,0)[0] : &w2m(ie,np1,i,j,0)[0];
              for (int k = 0; k < nlev; ++k)
                REQUIRE(equal(pf[k], pc[k], 1e8*eps));
            }
          }
    }
  }

  Session::delete_singleton();
}
