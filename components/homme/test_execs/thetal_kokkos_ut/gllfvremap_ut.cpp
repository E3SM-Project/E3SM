#include "GllFvRemapImpl.hpp"

#include "Types.hpp"
#include "Context.hpp"
#include "mpi/Comm.hpp"
#include "mpi/Connectivity.hpp"
#include "mpi/MpiBuffersManager.hpp"
#include "FunctorsBuffersManager.hpp"
#include "SimulationParams.hpp"
#include "Elements.hpp"
#include "Tracers.hpp"
#include "TimeLevel.hpp"
#include "HybridVCoord.hpp"
#include "KernelVariables.hpp"
#include "ElementOps.hpp"
#include "EquationOfState.hpp"
#include "profiling.hpp"
#include "ErrorDefs.hpp"

#include "utilities/TestUtils.hpp"
#include "utilities/SyncUtils.hpp"
#include "utilities/ViewUtils.hpp"

#include <catch2/catch.hpp>
#include <random>

using namespace Homme;

extern int hommexx_catch2_argc;
extern char** hommexx_catch2_argv;

extern "C" {
  void limiter1_clip_and_sum_f90(int n, Real* spheremp, Real* qmin, Real* qmax, Real* dp, Real* q);
  void calc_dp_fv_f90(int nf, Real* ps, Real* dp_fv);
  void get_temperature_f90(int ie, int nf, bool theta_hydrostatic_mode, Real* T);

  void init_gllfvremap_f90(int ne, const Real* hyai, const Real* hybi, const Real* hyam,
                           const Real* hybm, Real ps0, Real* dvv, Real* mp, int qsize,
                           bool is_sphere);
  void init_geometry_f90();

  void run_gfr_test(int* nerr);
  void run_gfr_check_api(bool theta_hydrostatic_mode, int* nerr);

  void gfr_init_f90(int nf, bool theta_hydrostatic_mode);
  void gfr_init_hxx();
  void gfr_finish_f90();

  void init_dyn_data_f90(int nlev_align, int nlevp_align, int nq, Real* ps, Real* phis, Real* dp3d,
                         Real* vthdp, Real* uv, Real* omega, Real* q, Real* phinh_i);
  void gfr_dyn_to_fv_phys_f90(int nf, int nt, Real* ps, Real* phis, Real* T, Real* uv,
                              Real* omega_p, Real* q);

  void gfr_fv_phys_to_dyn_f90(int nf, int nt, Real* T, Real* uv, Real* q);
  void cmp_dyn_data_f90(int nlev_align, int nq, Real* ft, Real* fm, Real* q, Real* fq, int* nerr);
} // extern "C"

using CA1d = Kokkos::View<Real*,     Kokkos::LayoutRight, Kokkos::HostSpace>;
using CA2d = Kokkos::View<Real**,    Kokkos::LayoutRight, Kokkos::HostSpace>;
using CA3d = Kokkos::View<Real***,   Kokkos::LayoutRight, Kokkos::HostSpace>;
using CA4d = Kokkos::View<Real**** , Kokkos::LayoutRight, Kokkos::HostSpace>;
using CA5d = Kokkos::View<Real*****, Kokkos::LayoutRight, Kokkos::HostSpace>;

template <typename V>
decltype(Kokkos::create_mirror_view(V())) cmv (const V& v) {
  return Kokkos::create_mirror_view(v);
}

template <typename V>
decltype(Kokkos::create_mirror_view(V())) cmvdc (const V& v) {
  const auto h = Kokkos::create_mirror_view(v);
  deep_copy(h, v);
  return h;
}

// Make a,b a little better behaved so levels don't get too thin.
static void clean (HybridVCoord& h) {
  static const int n = NUM_INTERFACE_LEV, nh = (n+1)/2, nh0 = n-nh;
  const auto ai = cmvdc(h.hybrid_ai);
  const auto bi = cmvdc(h.hybrid_bi);
  const auto amp = cmvdc(h.hybrid_am);
  const auto bmp = cmvdc(h.hybrid_bm);
  const CA1d am(reinterpret_cast<Real*>(amp.data()), n-1);
  const CA1d bm(reinterpret_cast<Real*>(bmp.data()), n-1);

  for (int i = 0; i < n; ++i) bi(i) = 0;
  for (int i = 0; i < nh; ++i) {
    assert(nh0+i < n);
    const Real a = Real(i)/(nh-1);
    bi(nh0+i) = (1-a)*0.02 + a*1;
  }
  assert(bi(n-1) == 1);

  Real etai[n];
  for (int i = 0; i < n; ++i) {
    const Real a = Real(i)/(n-1);
    etai[i] = (1-a)*0.0001 + a*1;
  }

  for (int i = 0; i < n; ++i) {
    ai(i) = etai[i] - bi(i);
    assert(ai(i) >= 0);
  }

  for (int i = 0; i < n-1; ++i) am(i) = (ai(i) + ai(i+1))/2;
  for (int i = 0; i < n-1; ++i) bm(i) = (bi(i) + bi(i+1))/2;

  deep_copy(h.hybrid_ai, ai);
  deep_copy(h.hybrid_bi, bi);
  deep_copy(h.hybrid_am, amp);
  deep_copy(h.hybrid_bm, bmp);

  h.hybrid_ai0 = ai(0);
  h.compute_deltas();
  h.compute_eta();
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
  //Random () : seed(346068100), engine(seed) {}
  unsigned int gen_seed () { return seed; }
  Real urrng (const Real lo = 0, const Real hi = 1) { return rpdf(lo, hi)(engine); }
  int  uirng (const int lo, const int hi) { return ipdf(lo, hi)(engine); }
};

template <typename V>
void fill (Random& r, const V& a, const Real scale = 1,
           typename std::enable_if<V::rank == 3>::type* = 0) {
  const auto am = cmvdc(a);
  for (int i = 0; i < a.extent_int(0); ++i)
    for (int j = 0; j < a.extent_int(1); ++j)
      for (int k = 0; k < a.extent_int(2); ++k)
        for (int s = 0; s < VECTOR_SIZE; ++s)
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
          for (int s = 0; s < VECTOR_SIZE; ++s)
            am(i,j,k,l)[s] = r.urrng(-1,1); 
  deep_copy(a, am);
}

struct Session {
  int ne;
  bool is_sphere;
  HybridVCoord h;
  Random r;
  std::shared_ptr<Elements> e;
  int nelemd, qsize, nlev, np;
  FunctorsBuffersManager fbm;

  //Session () : r(269041989) {}

  void init () {
    const auto seed = r.gen_seed();
    printf("seed %u\n", seed);

    parse_command_line();
    assert(is_sphere); // planar isn't available in Hxx yet

    auto& c = Context::singleton();

    c.create<HybridVCoord>().random_init(seed);
    auto& h_ref = c.get<HybridVCoord>();
    clean(h_ref);
    h = h_ref;

    auto& p = c.create<SimulationParams>();
    p.qsize = qsize;
    p.hypervis_scaling = 0;
    p.transport_alg = 0;
    p.moisture = MoistDry::MOIST;
    p.theta_hydrostatic_mode = false;
    p.scale_factor = is_sphere ? PhysicalConstants::rearth0 : 1;
    p.laplacian_rigid_factor = is_sphere ? 1/p.scale_factor : 0;
    p.params_set = true;

    const auto hyai = cmvdc(h.hybrid_ai);
    const auto hybi = cmvdc(h.hybrid_bi);
    const auto hyam = cmvdc(h.hybrid_am);
    const auto hybm = cmvdc(h.hybrid_bm);
    auto& ref_FE = c.create<ReferenceElement>();
    std::vector<Real> dvv(NP*NP), mp(NP*NP);
    init_gllfvremap_f90(ne, hyai.data(), hybi.data(), &hyam(0)[0], &hybm(0)[0], h.ps0,
                        dvv.data(), mp.data(), qsize, is_sphere);
    ref_FE.init_mass(mp.data());
    ref_FE.init_deriv(dvv.data());

    nelemd = c.get<Connectivity>().get_num_local_elements();
    auto& bmm = c.create<MpiBuffersManagerMap>();
    bmm.set_connectivity(c.get_ptr<Connectivity>());
    e = c.get_ptr<Elements>();
    c.create<TimeLevel>();

    init_geometry_f90();    
    auto& geo = c.get<ElementsGeometry>();

    auto& sphop = c.create<SphereOperators>();
    sphop.setup(geo, ref_FE);

    auto& gfr = c.create<GllFvRemap>();
    gfr.reset(p);
    fbm.request_size(gfr.requested_buffer_size());
    fbm.allocate();
    gfr.init_buffers(fbm);
    gfr.init_boundary_exchanges();

    nlev = NUM_PHYSICAL_LEV;
    assert(nlev > 0);
    np = NP;
    assert(np == 4);
  }

  void cleanup () {
    const auto& c = Context::singleton();
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

  Comm& get_comm () const { return Context::singleton().get<Comm>(); }

private:
  static std::shared_ptr<Session> s_session;

  // gllfvremap_ut hommexx -ne NE -qsize QSIZE
  void parse_command_line () {
    const bool am_root = get_comm().root();
    ne = 2;
    qsize = QSIZE_D;
    is_sphere = true;
    bool ok = true;
    int i;
    for (i = 0; i < hommexx_catch2_argc; ++i) {
      const std::string tok(hommexx_catch2_argv[i]);
      if (tok == "-ne") {
        if (i+1 == hommexx_catch2_argc) { ok = false; break; }
        ne = std::atoi(hommexx_catch2_argv[++i]);
      } else if (tok == "-qsize") {
        if (i+1 == hommexx_catch2_argc) { ok = false; break; }
        qsize = std::atoi(hommexx_catch2_argv[++i]);
      } else if (tok == "-planar") {
        is_sphere = false;
      }
    }
    ne = std::max(2, std::min(128, ne));
    qsize = std::max(1, std::min(QSIZE_D, qsize));
    if ( ! ok && am_root)
      printf("gllfvremap_ut> Failed to parse command line, starting with: %s\n",
             hommexx_catch2_argv[i]);
    if (am_root) {
      const int bfb =
#ifdef HOMMEXX_BFB_TESTING
        1;
#else
        0;
#endif
      printf("gllfvremap_ut> bfb %d ne %d qsize %d\n", bfb, ne, qsize);
    }
  }
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
                   const Real tol = 1e4*std::numeric_limits<Real>::epsilon()) {
#ifdef HOMMEXX_BFB_TESTING
  if (a != b)
    printf("equal: a,b = %23.16e %23.16e re = %23.16e\n",
           a, b, std::abs((a-b)/a));
  return a == b;
#else
  return almost_equal(a, b, tol);
#endif
}

static void test_calc_dp_fv (Random& r, const HybridVCoord& hvcoord) {
  using Kokkos::deep_copy;
  using g = GllFvRemapImpl;
  
  const int nf = 3, ncol = nf*nf;
  
  CA1d ps_f90("ps", ncol);
  CA2d dp_fv_f90("dp_fv", g::num_phys_lev, ncol);
  for (int i = 0; i < ncol; ++i) ps_f90(i) = r.urrng(0.9e5, 1.05e5);
  calc_dp_fv_f90(nf, ps_f90.data(), dp_fv_f90.data());

  const ExecView<Real*> ps_d("ps", ncol);
  const ExecView<Scalar**> dp_fv_p("dp_fv", ncol, g::num_lev_pack);
  const auto ps_h = cmv(ps_d);
  for (int i = 0; i < ncol; ++i) ps_h(i) = ps_f90(i);
  deep_copy(ps_d, ps_h);
  Kokkos::parallel_for(
    Homme::get_default_team_policy<ExecSpace>(1),
    KOKKOS_LAMBDA (const g::MT& team) {
      g::calc_dp_fv(team, hvcoord, ncol, g::num_lev_pack, ps_d, dp_fv_p);
    });
  const ExecViewUnmanaged<Real**> dp_fv_d(g::pack2real(dp_fv_p), ncol, g::num_lev_aligned);
  const auto dp_fv_h = cmvdc(dp_fv_d);
  
  for (int i = 0; i < ncol; ++i)
    for (int k = 0; k < g::num_phys_lev; ++k)
      REQUIRE(equal(dp_fv_f90(k,i), dp_fv_h(i,k)));
}

static void sfwd_remapd (const int m, const int n, const Real* A,
                         const Real* d1, const Real s2, const Real* d2,
                         const Real* x, Real* y) {
  for (int i = 0; i < m; ++i) {
    y[i] = 0;
    for (int j = 0; j < n; ++j)
      y[i] += A[n*i + j] * (x[j] * d1[j]);
    y[i] /= s2 * d2[i];
  }
}

static void sfwd_remapd (const int m, const int n, const Real* A,
                         const Real* Dinv, const Real s2, const Real* D,
                         const Real* x, Real* wx, Real* wy, Real* y) {
  static const int nones = 128;
  Real ones[nones]; assert(std::max(m,n) < nones);
  for (int i = 0; i < nones; ++i) ones[i] = 1;
  for (int d = 0; d < 2; ++d) {
    const int os = d*n;
    for (int i = 0; i < n; ++i)
      wx[os+i] = Dinv[4*i+2*d]*x[i] + Dinv[4*i+2*d+1]*x[n+i];
  }
  for (int d = 0; d < 2; ++d)
    sfwd_remapd(m, n, A, ones, s2, ones, wx + d*n, wy + d*m);
  for (int d = 0; d < 2; ++d) {
    const int os = d*m;
    for (int i = 0; i < m; ++i)
      y[os+i] = D[4*i+2*d]*wy[i] + D[4*i+2*d+1]*wy[m+i];
  }  
}

// Comparison of straightforwardly computed remapd vs GllFvRemapImpl version.
static void test_remapds (Random& r, const int m, const int n, const int nlev) {
  using Kokkos::deep_copy;
  using g = GllFvRemapImpl;

  // Random data for both scalar and vector remapd.

  std::vector<Real> A(m*n), d1(n), d2(m), x(2*n), y(2*m), Dinv(4*n), D(4*m);
  for (int i = 0; i < m*n; ++i) A [i] = r.urrng(0.1, 1.3);
  for (int i = 0; i <   n; ++i) d1[i] = r.urrng(0.2, 1.2);
  for (int i = 0; i < m  ; ++i) d2[i] = r.urrng(0.3, 1.1);
  for (int i = 0; i < 2*n; ++i) x [i] = r.urrng(0.4, 0.9);
  for (int i = 0; i < 4*n; ++i) Dinv[i] = r.urrng(0.5, 1.1);
  for (int i = 0; i < 4*m; ++i) D   [i] = r.urrng(0.6, 1.2);

  const int nlevpk = (nlev + g::packn - 1)/g::packn;
  const int nlevsk = nlevpk*g::packn;
  // Size arrays larger than needed because we want to support that case.
  const ExecView<Real*[2][2]> Dinv_d("Dinv", n), D_d("D", m);
  const ExecView<Real**> A_d("A", m+1, n+2);
  const ExecView<Real*> d1_d("d1", n+1), d2_d("d2", m+3);
  const ExecView<Scalar**> x_p("x", n+1, nlevpk), y_p("y", m+2, nlevpk);
  const ExecView<Real**> x_d(g::pack2real(x_p), n+1, nlevsk),
    y_d(g::pack2real(y_p), m+2, nlevsk);
  const auto Dinv_h = cmv(Dinv_d), D_h = cmv(D_d);
  const auto A_h = cmv(A_d);
  const auto d1_h = cmv(d1_d), d2_h = cmv(d2_d);
  const auto x_h = cmv(x_d), y_h = cmv(y_d);
    for (int d1 = 0; d1 < 2; ++d1)
      for (int d2 = 0; d2 < 2; ++d2) {
        for (int i = 0; i < n; ++i) Dinv_h(i,d1,d2) = Dinv[4*i + 2*d1 + d2];
        for (int i = 0; i < m; ++i) D_h   (i,d1,d2) = D   [4*i + 2*d1 + d2];
      }
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j)
      A_h(i,j) = A[n*i + j];
  for (int i = 0; i < n; ++i) d1_h(i) = d1[i];
  for (int i = 0; i < m; ++i) d2_h(i) = d2[i];
  for (int k = 0; k < nlev; ++k)
    for (int i = 0; i < n; ++i)
      x_h(i,k) = x[i];
  deep_copy(Dinv_d, Dinv_h); deep_copy(D_d, D_h);
  deep_copy(A_d, A_h); deep_copy(d1_d, d1_h); deep_copy(d2_d, d2_h);
  deep_copy(x_d, x_h);

  // Scalar remapd.
  sfwd_remapd(m, n, A.data(), d1.data(), 5, d2.data(), x.data(), y.data());
  Kokkos::parallel_for(
    Homme::get_default_team_policy<ExecSpace>(1),
    KOKKOS_LAMBDA (const g::MT& team) {
      g::remapd(team, m, n, nlevpk, A_d, d1_d, 5, d2_d, x_p, x_p, y_p); });
  deep_copy(x_h, x_d); deep_copy(y_h, y_d);
  for (int k = 0; k < nlev; ++k)
    for (int i = 0; i < m; ++i)
      REQUIRE(equal(y_h(i,k), y[i]));

  // Vector remapd.
  std::vector<Real> wx(2*n), wy(2*m);
  sfwd_remapd(m, n, A.data(), Dinv.data(), 5, D.data(), x.data(), wx.data(), wy.data(), y.data());
  { // x(i,d), y(d,i)
    const ExecView<Scalar***> x2_p("x2", n+1, 2, nlevpk), y2_p("y", 2, m+2, nlevpk);
    const ExecView<Real***> x2_d(g::pack2real(x2_p), n+1, 2, nlevsk),
      y2_d(g::pack2real(y2_p), 2, m+2, nlevsk);
    const auto x2_h = cmv(x2_d), y2_h = cmv(y2_d);
    for (int k = 0; k < nlev; ++k)
      for (int i = 0; i < n; ++i)
        for (int d = 0; d < 2; ++d)
          x2_h(i,d,k) = x[d*n+i];
    deep_copy(x2_d, x2_h);
    Kokkos::parallel_for(
      Homme::get_default_team_policy<ExecSpace>(1),
      KOKKOS_LAMBDA (const g::MT& team) {
        g::remapd<true>(team, m, n, nlevpk, A_d, Dinv_d, 5, D_d, x2_p, x2_p, y2_p); });
    deep_copy(y2_h, y2_d);
    for (int k = 0; k < nlev; ++k)
      for (int i = 0; i < m; ++i)
        for (int d = 0; d < 2; ++d)
          REQUIRE(equal(y2_h(d,i,k), y[d*m+i]));
  }
  { // x(d,i), y(i,d)
    const ExecView<Scalar***> x2_p("x2", 2, n+1, nlevpk), y2_p("y", m+2, 2, nlevpk);
    const ExecView<Real***> x2_d(g::pack2real(x2_p), 2, n+1, nlevsk),
      y2_d(g::pack2real(y2_p), m+2, 2, nlevsk);
    const auto x2_h = cmv(x2_d), y2_h = cmv(y2_d);
    for (int k = 0; k < nlev; ++k)
      for (int i = 0; i < n; ++i)
        for (int d = 0; d < 2; ++d)
          x2_h(d,i,k) = x[d*n+i];
    deep_copy(x2_d, x2_h);
    Kokkos::parallel_for(
      Homme::get_default_team_policy<ExecSpace>(1),
      KOKKOS_LAMBDA (const g::MT& team) {
        g::remapd<false>(team, m, n, nlevpk, A_d, Dinv_d, 5, D_d, x2_p, x2_p, y2_p); });
    deep_copy(y2_h, y2_d);
    for (int k = 0; k < nlev; ++k)
      for (int i = 0; i < m; ++i)
        for (int d = 0; d < 2; ++d)
          REQUIRE(equal(y2_h(i,d,k), y[d*m+i]));
  }
}

template <typename V1O, typename V2O, typename V1s, typename V1, typename V2, typename V2q>
static void
assert_limiter_properties (const int n, const int nlev, const V1s& spheremp,
                           const V1O& qmin_orig, const V1O& qmax_orig,
                           const V1& qmin, const V1& qmax,
                           const V2& dp, const V2O& qorig, const V2q& q,
                           const bool too_tight) {
  static const auto eps = std::numeric_limits<Real>::epsilon();
  const int n2 = n*n;
  int noteq = 0;
  for (int k = 0; k < nlev; ++k) {
    Real qm0 = 0, qm = 0;
    for (int i = 0; i < n2; ++i) {
      REQUIRE(q(i,k) >= (1 - 5e1*eps)*qmin(k));
      REQUIRE(q(i,k) <= (1 + 5e1*eps)*qmax(k));
      qm0 += spheremp(i)*dp(i,k)*qorig(i,k);
      qm  += spheremp(i)*dp(i,k)*q(i,k);
      if (q(i,k) != qorig(i,k)) ++noteq;
    }
    REQUIRE(almost_equal(qm0, qm, 1e2*eps));
    if (too_tight && k % 2 == 1) {
      for (int i = 1; i < n2; ++i)
        REQUIRE(almost_equal(q(i,k), q(0,k), 1e2*eps));
    } else {
      REQUIRE(equal(qmin(k), qmin_orig(k)));
      REQUIRE(equal(qmax(k), qmax_orig(k)));
    }
  }
  REQUIRE(noteq > 0);
}

static void test_limiter (const int nlev, const int n, Random& r, const bool too_tight) {
  using Kokkos::deep_copy;
  using g = GllFvRemapImpl;

  const int n2 = n*n;
  const int nlevpk = (nlev + g::packn - 1)/g::packn;
  const int nlevsk = nlevpk*g::packn;

  const ExecView<Real*> spheremp_d("spheremp", n2);
  const ExecView<Scalar*> qmin_p("qmin", nlevpk), qmax_p("qmax", nlevpk);
  const ExecView<Scalar**> dp_p("dp", n2, nlevpk), q_p("q", n2, nlevpk), wrk_p("wrk", n2, nlevpk);
  const ExecView<Real*> qmin_d(g::pack2real(qmin_p), nlevsk), qmax_d(g::pack2real(qmax_p), nlevsk);
  const ExecView<Real**> dp_d(g::pack2real(dp_p), n2, nlevsk), q_d(g::pack2real(q_p), n2, nlevsk);

  const HostViewManaged<Real*> qmin_orig("qmin_orig", nlevsk), qmax_orig("qmax_orig", nlevsk);
  const HostViewManaged<Real**> qorig("qorig", n2, nlevsk);

  const auto spheremp = cmv(spheremp_d);
  const auto qmin = cmv(qmin_d), qmax = cmv(qmax_d);
  const auto dp = cmv(dp_d), q = cmv(q_d);

  for (int k = 0; k < nlev; ++k) {
    Real mass = 0, qmass = 0;
    qmin(k) = 1; qmax(k) = 0;
    for (int i = 0; i < n2; ++i) {
      if (k == 0) spheremp(i) = r.urrng(0.5, 1.5);
      dp(i,k) = r.urrng(0.5, 1.5);
      qorig(i,k) = q(i,k) = r.urrng(1e-2, 3e-2);
      mass += spheremp(i)*dp(i,k);
      qmass += spheremp(i)*dp(i,k)*q(i,k);
      qmin(k) = std::min(qmin(k), q(i,k));
      qmax(k) = std::max(qmax(k), q(i,k));
    }
    const auto q0 = qmass/mass;
    if (too_tight && k % 2 == 1) {
      qmin(k) = 0.5*(q0 + qmin(k));
      qmax(k) = 0.1*qmin(k) + 0.9*q0;
    } else {
      qmin(k) = 0.5*(q0 + qmin(k));
      qmax(k) = 0.5*(q0 + qmax(k));
    }
    qmin_orig(k) = qmin(k);
    qmax_orig(k) = qmax(k);
  }

  deep_copy(spheremp_d, spheremp);
  deep_copy(qmin_d, qmin); deep_copy(qmax_d, qmax);
  deep_copy(dp_d, dp); deep_copy(q_d, q);

  // F90 limiter properties.
  CA1d qminf90("qminf90", nlevsk), qmaxf90("qmaxf90", nlevsk);
  CA2d qf90("q", n2, nlevsk);
  deep_copy(qf90, qorig); deep_copy(qminf90, qmin_orig); deep_copy(qmaxf90, qmax_orig);
  CA1d dpk("dpk", n2), qk("qk", n2);
  for (int k = 0; k < nlev; ++k) {
    for (int i = 0; i < n2; ++i) dpk(i) = dp(i,k);
    for (int i = 0; i < n2; ++i) qk(i) = qf90(i,k);
    limiter1_clip_and_sum_f90(n, spheremp.data(), &qminf90(k), &qmaxf90(k), dpk.data(),
                              qk.data());
    for (int i = 0; i < n2; ++i) qf90(i,k) = qk(i);
  }
  assert_limiter_properties(n, nlev, spheremp,
                            qmin_orig, qmax_orig, qminf90, qmaxf90,
                            dp, qorig, qf90, too_tight);

  // C++ limiter properties.
  Kokkos::parallel_for(
    Homme::get_default_team_policy<ExecSpace>(1),
    KOKKOS_LAMBDA (const g::MT& team) {
      g::limiter_clip_and_sum(team, n2, nlevpk, 1, spheremp_d, qmin_p, qmax_p, dp_p,
                              wrk_p, q_p); });
  deep_copy(q, q_d); deep_copy(qmin, qmin_d); deep_copy(qmax, qmax_d);
  assert_limiter_properties(n, nlev, spheremp,
                            qmin_orig, qmax_orig, qmin, qmax,
                            dp, qorig, q, too_tight);

  // BFB C++ vs F90.
  for (int k = 0; k < nlev; ++k)
    for (int i = 0; i < n2; ++i)
      REQUIRE(equal(qf90(i,k), q(i,k)));

  { // BFB C++ real1 vs pack

    // Run the pack version with dp_p = ones b/c that is what the real1 supports.
    const ExecView<Scalar**> ones("ones", n2, nlevpk);
    deep_copy(ones, 1);
    deep_copy(q_d, qorig); deep_copy(qmin_d, qmin_orig); deep_copy(qmax_d, qmax_orig);
    Kokkos::parallel_for(
    Homme::get_default_team_policy<ExecSpace>(1),
    KOKKOS_LAMBDA (const g::MT& team) {
      g::limiter_clip_and_sum(team, n2, nlevpk, 1, spheremp_d, qmin_p, qmax_p, ones,
                              wrk_p, q_p); });
    deep_copy(q, q_d); deep_copy(qmin, qmin_d); deep_copy(qmax, qmax_d);

    // Run the real1 version and BFB cmp with the pack version.
    const ExecView<Real*> q1("q1", n2), qmin1("qmin1", nlevsk), qmax1("qmax1", nlevsk),
      wrk("wrk", n2);
    deep_copy(qmin1, qmin_orig); deep_copy(qmax1, qmax_orig);
    const auto q1h = cmv(q1);
    const auto qmin1h = cmv(qmin1), qmax1h = cmv(qmax1);
    for (int k = 0; k < nlev; ++k) {
      for (int i = 0; i < n2; ++i) q1h(i) = qorig(i,k);
      deep_copy(q1, q1h);
      Kokkos::parallel_for(
        Homme::get_default_team_policy<ExecSpace>(1),
        KOKKOS_LAMBDA (const g::MT& team) {
          g::limiter_clip_and_sum_real1(team, n2, 1, spheremp_d, qmin1(k), qmax1(k),
                                        wrk, q1); });
      deep_copy(q1h, q1); deep_copy(qmin1h, qmin1); deep_copy(qmax1h, qmax1);
      for (int i = 0; i < n2; ++i) REQUIRE(equal(q1h(i), q(i,k)));
      REQUIRE(equal(qmin1h(k), qmin(k)));
      REQUIRE(equal(qmax1h(k), qmax(k)));
    }
  }
}

static void init_dyn_data (Session& s) {
  using Kokkos::parallel_for;
  using Kokkos::deep_copy;
  using g = GllFvRemapImpl;

  // randomize C++ data
  auto& c = Context::singleton();
  const auto hai = cmvdc(s.h.hybrid_ai);
  const auto bai = cmvdc(s.h.hybrid_bi);
  const auto state = c.get<ElementsState>();
  const auto derived = c.get<ElementsDerivedState>();
  const auto tracers = c.get<Tracers>();
  const auto geometry = c.get<ElementsGeometry>();
  const int nt = NUM_TIME_LEVELS;
  const auto ps = cmv(state.m_ps_v);

  const g::EVU<Real****>
    omega_s(g::pack2real(derived.m_omega_p), s.nelemd, NP, NP, g::num_lev_aligned);
  const g::EVU<Real*****>
    dp3d_s(g::pack2real(state.m_dp3d), s.nelemd, nt, NP, NP, g::num_lev_aligned),
    vthdp_s(g::pack2real(state.m_vtheta_dp), s.nelemd, nt, NP, NP, g::num_lev_aligned),
    q_s(g::pack2real(tracers.Q), s.nelemd, tracers.Q.extent_int(1), NP, NP, g::num_lev_aligned);
  const g::EVU<Real******>
    uv_s(g::pack2real(state.m_v), s.nelemd, nt, 2, NP, NP, g::num_lev_aligned);

  const auto phis = cmv(geometry.m_phis);
  const auto omega = cmv(omega_s);
  const auto dp3d = cmv(dp3d_s);
  const auto vthdp = cmv(vthdp_s);
  const auto q = cmv(q_s);
  const auto uv = cmv(uv_s);
  for (int ie = 0; ie < s.nelemd; ++ie)
    for (int i = 0; i < s.np; ++i)
      for (int j = 0; j < s.np; ++j)
        for (int k = 0; k < s.nlev; ++k) {
          for (int t = 0; t < nt; ++t) {
            if (k == 0) ps(ie,t,i,j) = s.h.ps0*(s.r.urrng(0.9, 1.1));
            dp3d(ie,t,i,j,k) = ((hai(k+1) - hai(k))*s.h.ps0 +
                                (bai(k+1) - bai(k))*ps(ie,t,i,j));
            vthdp(ie,t,i,j,k) = s.r.urrng(150, 320)*dp3d(ie,t,i,j,k);
            for (int d = 0; d < 2; ++d) uv(ie,t,d,i,j,k) = s.r.urrng(-30, 30);
          }
          if (k == 0) phis(ie,i,j) = s.r.urrng(-1000, 10000);
          omega(ie,i,j,k) = s.r.urrng(-10, 10);
          for (int iq = 0; iq < s.qsize; ++iq)
            q(ie,iq,i,j,k) = s.r.urrng(0, 0.1);
        }
  deep_copy(geometry.m_phis, phis);
  deep_copy(omega, omega_s);
  deep_copy(state.m_ps_v, ps);
  deep_copy(dp3d_s, dp3d);
  deep_copy(vthdp_s, vthdp);
  deep_copy(q_s, q);
  deep_copy(uv_s, uv);

  {
    EquationOfState eos; eos.init(c.get<SimulationParams>().theta_hydrostatic_mode, s.h);
    ElementOps ops; ops.init(s.h);
    const auto phis = geometry.m_phis;
    const auto dp3d = state.m_dp3d;
    const auto vthdp = state.m_vtheta_dp;
    const auto phi_i = state.m_phinh_i;
    const ExecViewManaged<Scalar*[NP][NP][NUM_LEV_P]> p_i("p_i", s.nelemd);
    const ExecViewManaged<Scalar*[NP][NP][NUM_LEV]> p_m("p_m", s.nelemd);
    for (int t = 0; t < nt; ++t) {
      Kokkos::parallel_for(
        Homme::get_default_team_policy<ExecSpace>(s.nelemd),
        KOKKOS_LAMBDA (const g::MT& team) {
          KernelVariables kv(team);
          const auto ie = kv.ie;
          parallel_for(
            Kokkos::TeamThreadRange(team, NP*NP),
            [&] (const int k) {
              const auto i = k / NP, j = k % NP;
              ops.compute_hydrostatic_p(kv, Homme::subview(dp3d,ie,t,i,j),
                                        Homme::subview(p_i,ie,i,j), Homme::subview(p_m,ie,i,j));
            });
          kv.team_barrier();
          eos.compute_phi_i(kv, Homme::subview(phis,ie), Homme::subview(vthdp,ie,t),
                            Homme::subview(p_m,ie), Homme::subview(phi_i,ie,t));
        });
    }
  }
  const g::EVU<Real*****>
    phinh_i_s(g::pack2real(state.m_phinh_i), s.nelemd, nt, NP, NP, g::num_levp_aligned);
  const auto phinh_i = cmvdc(phinh_i_s);
  
  init_dyn_data_f90(g::num_lev_aligned, g::num_levp_aligned, q.extent_int(1),
                    ps.data(), phis.data(), dp3d.data(), vthdp.data(), uv.data(),
                    omega.data(), q.data(), phinh_i.data());
}

static void test_get_temperature (Session& s) {
  using g = GllFvRemapImpl;
  const int nt = NUM_TIME_LEVELS;

  for (const bool theta_hydrostatic_mode : {false, true}) {
    init_dyn_data(s);
  
    const ExecViewManaged<Scalar*[NUM_TIME_LEVELS][NP][NP][NUM_LEV]> Td("T", s.nelemd); {
      auto& c = Context::singleton();
      const auto& sp = c.get<SimulationParams>();
      EquationOfState eos; eos.init(theta_hydrostatic_mode, s.h);
      ElementOps ops; ops.init(s.h);
      const bool use_moisture = sp.moisture == MoistDry::MOIST;
      const auto state = c.get<ElementsState>();
      const auto tracers = c.get<Tracers>();
      const auto dp3d = state.m_dp3d;
      const auto vthdp = state.m_vtheta_dp;
      const auto phi_i = state.m_phinh_i;
      const auto q = tracers.Q;
      const ExecView<Scalar*[NP][NP][NUM_LEV]> wrk("wrk", s.nelemd), exner("exner", s.nelemd);
      const ExecView<Scalar*[NP][NP][NUM_LEV_P]> pi("pi", s.nelemd);
      for (int t = 0; t < nt; ++t) {
        Kokkos::parallel_for(
          Homme::get_default_team_policy<ExecSpace>(s.nelemd),
          KOKKOS_LAMBDA (const g::MT& team) {
            KernelVariables kv(team);
            const auto ie = kv.ie;
            Kokkos::parallel_for(
              Kokkos::TeamThreadRange(team, NP*NP),
              [&] (const int ij) {
                const auto i = ij / NP, j = ij % NP;
                const auto dp3d_ij = Homme::subview(dp3d,ie,t,i,j);
                const auto vthdp_ij = Homme::subview(vthdp,ie,t,i,j);
                const auto wrk_ij = Homme::subview(wrk,ie,i,j);
                const auto exner_ij = Homme::subview(exner,ie,i,j);
                if (theta_hydrostatic_mode) {
                  ops.compute_hydrostatic_p(kv, dp3d_ij, Homme::subview(pi,ie,i,j), exner_ij);
                  eos.compute_exner(kv, exner_ij, exner_ij);
                } else {
                  eos.compute_pnh_and_exner(kv, vthdp_ij, Homme::subview(phi_i,ie,t,i,j), wrk_ij,
                                            exner_ij);
                }
                ops.get_temperature(kv, eos, use_moisture, dp3d_ij, exner_ij,
                                    vthdp_ij, Homme::subview(q,ie,0,i,j), wrk_ij,
                                    Homme::subview(Td,ie,t,i,j));
              });
          });
      }
    }

    const CA5d Tf90("Tf90", s.nelemd, nt, s.nlev, s.np, s.np);
    for (int ie = 0; ie < s.nelemd; ++ie)
      for (int t = 0; t < nt; ++t)
        get_temperature_f90(ie+1, t+1, theta_hydrostatic_mode, &Tf90(ie,t,0,0,0));

    const g::EVU<Real*****> Ts(g::pack2real(Td), s.nelemd, nt, NP, NP, g::num_lev_aligned);
    const auto T = cmvdc(Ts);
    for (int ie = 0; ie < s.nelemd; ++ie)
      for (int t = 0; t < nt; ++t)
        for (int i = 0; i < s.np; ++i)
          for (int j = 0; j < s.np; ++j)
            for (int k = 0; k < s.nlev; ++k)
              REQUIRE(equal(T(ie,t,i,j,k), Tf90(ie,t,k,i,j)));
  }
}

static void
test_dyn_to_fv_phys (Session& s, const int nf, const bool theta_hydrostatic_mode) {
  using g = GllFvRemapImpl;

  const int nf2 = nf*nf;

  gfr_init_f90(nf, theta_hydrostatic_mode);
  gfr_init_hxx();

  init_dyn_data(s);

  CA2d fps("fps", s.nelemd, nf2), fphis("fphis", s.nelemd, nf2);
  CA3d fT("fT", s.nelemd, s.nlev, nf2), fomega("fomega", s.nelemd, s.nlev, nf2);
  CA4d fuv("fuv", s.nelemd, s.nlev, 2, nf2), fq("fq", s.nelemd, s.qsize, s.nlev, nf2);

  const int nq = s.qsize > 1 ? s.qsize-1 : s.qsize;

  const ExecView<Real**> dps("dps", s.nelemd, nf2), dphis("dphis", s.nelemd, nf2);
  const ExecView<Real***> dT("dT", s.nelemd, nf2, g::num_lev_aligned),
    domega("domega", s.nelemd, nf2, g::num_lev_aligned);
  const ExecView<Real****> duv("duv", s.nelemd, nf2, 2, g::num_lev_aligned),
    dq("dq", s.nelemd, nf2, s.qsize, g::num_lev_aligned),
    dq1("dq", s.nelemd, nf2, nq, g::num_lev_aligned);

  const auto& c = Context::singleton();
  auto& gfr = c.get<GllFvRemap>();
  const auto tracers = c.get<Tracers>();
  const auto& q = tracers.Q;

  const g::EVU<const Real****>
    dq1_dyn(g::cpack2real(q), q.extent_int(0), q.extent_int(1),
            q.extent_int(2)*q.extent_int(3), g::num_lev_aligned);

  for (int nt = 0; nt < NUM_TIME_LEVELS; ++nt) {
    gfr_dyn_to_fv_phys_f90(nf, nt+1, fps.data(), fphis.data(), fT.data(), fuv.data(),
                           fomega.data(), fq.data());

    gfr.run_dyn_to_fv_phys(nt, dps, dphis, dT, domega, duv, dq);

    gfr.remap_tracer_dyn_to_fv_phys(nt, nq, dq1_dyn, dq1);

    const auto ps = cmvdc(dps);
    const auto phis = cmvdc(dphis);
    const auto T = cmvdc(dT);
    const auto omega = cmvdc(domega);
    const auto uv = cmvdc(duv);
    const auto q = cmvdc(dq);
    const auto q1 = cmvdc(dq1);

    for (int ie = 0; ie < s.nelemd; ++ie)
      for (int i = 0; i < nf2; ++i) {
        REQUIRE(equal(ps(ie,i), fps(ie,i)));
        REQUIRE(equal(phis(ie,i), fphis(ie,i)));
        for (int k = 0; k < s.nlev; ++k) {
          REQUIRE(equal(omega(ie,i,k), fomega(ie,k,i)));
          REQUIRE(equal(T(ie,i,k), fT(ie,k,i)));
          for (int iq = 0; iq < s.qsize; ++iq)
            REQUIRE(equal(q (ie,i,iq,k), fq(ie,iq,k,i)));
          for (int iq = 0; iq < nq; ++iq)
            REQUIRE(equal(q1(ie,i,iq,k), fq(ie,iq,k,i)));
          for (int d = 0; d < 2; ++d)
            REQUIRE(equal(uv(ie,i,d,k), fuv(ie,k,d,i)));
        }
      }
  }

  gfr_finish_f90();
}

static void
test_fv_phys_to_dyn (Session& s, const int nf, const bool theta_hydrostatic_mode) {
  using Kokkos::deep_copy;
  using g = GllFvRemapImpl;
  
  const int nf2 = nf*nf;

  gfr_init_f90(nf, theta_hydrostatic_mode);
  gfr_init_hxx();

  {
    CA3d fT("fT", s.nelemd, s.nlev, nf2);
    CA4d fuv("fuv", s.nelemd, s.nlev, 2, nf2), ffq("fq", s.nelemd, s.qsize, s.nlev, nf2);

    const ExecView<Real***> dT("dT", s.nelemd, nf2, g::num_lev_aligned);
    const ExecView<Real****> duv("duv", s.nelemd, nf2, 2, g::num_lev_aligned),
      dq("dq", s.nelemd, nf2, s.qsize, g::num_lev_aligned),
      dfq("dq", s.nelemd, nf2, s.qsize, g::num_lev_aligned);

    const auto T = cmv(dT);
    const auto uv = cmv(duv);
    const auto fq = cmv(dfq);

    for (int ie = 0; ie < s.nelemd; ++ie)
      for (int k = 0; k < s.nlev; ++k)
        for (int i = 0; i < nf2; ++i) {
          fT(ie,k,i) = T(ie,i,k) = s.r.urrng(100, 300);
          for (int d = 0; d < 2; ++d)
            fuv(ie,k,d,i) = uv(ie,i,d,k) = s.r.urrng(-30, 30);
          for (int iq = 0; iq < s.qsize; ++iq)
            ffq(ie,iq,k,i) = fq(ie,i,iq,k) = s.r.urrng(0, 0.1);
        }
  
    deep_copy(dT, T);
    deep_copy(duv, uv);
    deep_copy(dfq, fq);

    const auto& c = Context::singleton();
    auto& gfr = c.get<GllFvRemap>();

    const int nt = 1;
    gfr_fv_phys_to_dyn_f90(nf, nt+1, fT.data(), fuv.data(), ffq.data());
    gfr.run_fv_phys_to_dyn(nt, dT, duv, dfq);
    gfr.run_fv_phys_to_dyn_dss();
  }

  {
    auto& c = Context::singleton();
    const auto forcing = c.get<ElementsForcing>();
    const auto tracers = c.get<Tracers>();

    const g::EVU<Real****>
      ft_s(g::pack2real(forcing.m_ft), s.nelemd, NP, NP, g::num_lev_aligned);
    const g::EVU<Real*****>
      q_s( g::pack2real(tracers.Q), s.nelemd, tracers.Q.extent_int(1), NP, NP, g::num_lev_aligned),
      fq_s(g::pack2real(tracers.fq), s.nelemd, tracers.fq.extent_int(1), NP, NP, g::num_lev_aligned),
      fm_s(g::pack2real(forcing.m_fm), s.nelemd, forcing.m_fm.extent_int(1), NP, NP,
           g::num_lev_aligned);

    const auto q = cmvdc(q_s);
    const auto ft = cmvdc(ft_s);
    const auto fq = cmvdc(fq_s);
    const auto fm = cmvdc(fm_s);

    int nerr = 0;
    cmp_dyn_data_f90(g::num_lev_aligned, q.extent_int(1), ft.data(), fm.data(), q.data(),
                     fq.data(), &nerr);
    REQUIRE(nerr == 0);
  }

  gfr_finish_f90();
}

TEST_CASE ("gllfvremap_testing") {
  auto& s = Session::singleton(); try {
    test_get_temperature(s);
    test_calc_dp_fv(s.r, s.h);
    
    // Core scalar and vector remapd routines.
    test_remapds(s.r,  7, 11, 13);
    test_remapds(s.r, 11,  7, 13);
    test_remapds(s.r, 16,  4,  8);
    test_remapds(s.r,  4, 16,  8);

    // Limiter.
    for (const auto too_tight : {false, true}) {
      test_limiter(16, 7, s.r, too_tight);
      test_limiter(16, 4, s.r, too_tight);
    }

    // Existing F90 gllfvremap unit tests.
    {
      int nerr;
      run_gfr_test(&nerr);
      REQUIRE(nerr == 0);
    }
    {
      int nerr;
      const bool thm = true;
      Context::singleton().get<SimulationParams>().theta_hydrostatic_mode = thm;
      init_dyn_data(s);
      run_gfr_check_api(thm, &nerr);
      REQUIRE(nerr == 0);
    }

    // Main remap routines.
    for (const bool theta_hydrostatic_mode : {false, true}) {
      for (const int nf : {2,3,4}) {
        printf("ut> g2f nf %d thm %d\n", nf, (int) theta_hydrostatic_mode);
        test_dyn_to_fv_phys(s, nf, theta_hydrostatic_mode);
      }

      for (const int nf : {2,3,4}) {
        printf("ut> f2g nf %d thm %d\n", nf, (int) theta_hydrostatic_mode);
        test_fv_phys_to_dyn(s, nf, theta_hydrostatic_mode);
      }
    }
  } catch (...) {}
  Session::delete_singleton();
}
