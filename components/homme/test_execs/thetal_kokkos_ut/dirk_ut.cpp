#include <catch2/catch.hpp>

//#include "/home/ambradl/climate/sik/hommexx/dbg.hpp"
// make -j8 dirk_ut;if [ $? == 0 ]; then OMP_NUM_THREADS=1 ./test_execs/thetal_kokkos_ut/dirk_ut; fi
// OMP_NUM_THREADS=2 ctest -R dirk_ut -VV

#include "DirkFunctorImpl.hpp"

#include <random>

#include "Types.hpp"
#include "Context.hpp"
#include "SimulationParams.hpp"
#include "Tracers.hpp"

#include "utilities/TestUtils.hpp"
#include "utilities/SyncUtils.hpp"
#include "utilities/ViewUtils.hpp"

using namespace Homme;
using Kokkos::create_mirror_view;
using Kokkos::deep_copy;
using FA3 = Kokkos::View<Real*[NP][NP], Kokkos::LayoutRight, Kokkos::HostSpace>;
using FA3d = Kokkos::View<Real***, Kokkos::LayoutRight, Kokkos::HostSpace>;
using dfi = DirkFunctorImpl;

extern "C" {
  void init_dirk_f90(const Real* hyai, const Real* hybi, const Real* hyam, const Real* hybm,
                     Real ps0);
  void pnh_and_exner_from_eos_f90(const Real* vtheta_dp, const Real* dp3d, const Real* dphi,
                                  Real* pnh, Real* exner, Real* dpnh_dp_i);
  void get_dirk_jacobian_f90(Real* dl, Real* d, Real* du, Real dt, const Real* dp3d,
                             const Real* dphi, const Real* pnh);
} // extern "C"

class Random {
  using rngalg = std::mt19937_64;
  using rpdf = std::uniform_real_distribution<Real>;
  using ipdf = std::uniform_int_distribution<int>;
  std::random_device rd;
  rngalg engine;
public:
  Random (int seed = -1) : engine(seed == -1 ? rd() : seed) {}
  std::random_device& random_device () { return rd; }
  Real urrng (const Real lo = 0, const Real hi = 1) { return rpdf(lo, hi)(engine); }
  int  uirng (const int lo, const int hi) { return ipdf(lo, hi)(engine); }
};

static void init (Random& r, HybridVCoord& h) {
  h.random_init(r.random_device()());
  const auto hyai = create_mirror_view(h.hybrid_ai); deep_copy(hyai, h.hybrid_ai);
  const auto hybi = create_mirror_view(h.hybrid_bi); deep_copy(hybi, h.hybrid_bi);
  const auto hyam = create_mirror_view(h.hybrid_am); deep_copy(hyam, h.hybrid_am);
  const auto hybm = create_mirror_view(h.hybrid_bm); deep_copy(hybm, h.hybrid_bm);
  init_dirk_f90(hyai.data(), hybi.data(), &hyam(0)[0], &hybm(0)[0], h.ps0);
}

static bool equal (const Real& a, const Real& b) {
  if (a != b)
    printf("equal: a,b = %23.16e %23.16e re = %23.16e\n",
           a, b, std::abs((a-b)/a));
  return a == b;
}

template <typename V>
void fill (Random& r, const V& a,
           typename std::enable_if<V::rank == 3>::type* = 0) {
  const auto am = create_mirror_view(a);
  deep_copy(am, a);
  for (int i = 0; i < a.extent_int(0); ++i)
    for (int j = 0; j < a.extent_int(1); ++j)
      for (int k = 0; k < a.extent_int(2); ++k)
        for (int s = 0; s < dfi::packn; ++s)
          am(i,j,k)[s] = r.urrng(); 
  deep_copy(a, am);
}

template <typename V>
void fill_inc (Random& r, const int nlev, const Real top, const Real bottom,
               const V& a, typename std::enable_if<V::rank == 3>::type* = 0) {
  const auto am = create_mirror_view(a);
  deep_copy(am, a);
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
  const auto am = create_mirror_view(a);
  deep_copy(am, a);
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
  const auto am = create_mirror_view(a); deep_copy(am, a);
  const auto bm = create_mirror_view(b); deep_copy(bm, b);
  for (int i = 0; i < a.extent_int(0); ++i)
    for (int j = 0; j < a.extent_int(1); ++j)
      for (int k = 0; k < a.extent_int(2); ++k)
        for (int s = 0; s < dfi::packn; ++s)
          if (dfi::packn*k + s < nlev)
            REQUIRE(am(i,j,k)[s] == bm(i,j,k)[s]);
}

template <typename V>
void c2f(const V& c, const FA3& f) {
  const auto cm = create_mirror_view(c);
  deep_copy(cm, c);
  for (int i = 0; i < c.extent_int(0); ++i)
    for (int j = 0; j < c.extent_int(1); ++j) {
      Real* const p = &cm(i,j,0)[0];
      for (int k = 0; k < f.extent(0); ++k)
        f(k,i,j) = p[k];
    }
}

TEST_CASE("dirk", "dirk_testing") {
  using Kokkos::parallel_for;
  using Kokkos::fence;

  Random r;

  HybridVCoord hvcoord;
  init(r, hvcoord);

  SECTION ("transpose") {
    const int nelem = 1;
    ExecView<Scalar[NP][NP][NUM_LEV  ]> ham0("ham0"), ham1("ham1");
    ExecView<Scalar[NP][NP][NUM_LEV_P]> hai0("hai0"), hai1("hai1");
    const int nlev = NUM_LEV*VECTOR_SIZE - 2; // test with some remainder
    fill(r, ham0);
    fill(r, hai0);
    dfi d(nelem);
    const auto w = d.m_work;
    const auto a = dfi::get_slot(w, 0, 0);
    const auto f1 = KOKKOS_LAMBDA(const dfi::MT& t) {
      KernelVariables kv(t);
      dfi::transpose(kv, nlev, ham0, a);
      kv.team_barrier();
      dfi::transpose(kv, nlev, a, ham1);
    };
    parallel_for(d.m_policy, f1); fence();
    require_all_equal(nlev, ham0, ham1);
    const auto f2 = KOKKOS_LAMBDA(const dfi::MT& t) {
      KernelVariables kv(t);
      dfi::transpose(kv, nlev+1, hai0, a);
      kv.team_barrier();
      dfi::transpose(kv, nlev+1, a, hai1);
    };
    parallel_for(d.m_policy, f2); fence();
    require_all_equal(nlev+1, hai0, hai1);
  }

  SECTION ("jacobian") {
    const int nlev = NUM_PHYSICAL_LEV, np = NP;
    const Real dt = r.urrng(32, 64);

    ExecView<Scalar[NP][NP][NUM_LEV]> dp3d("dp3d"), dphi("dphi"), pnh("pnh");
    fill_inc(r, nlev, 5, 10000, dp3d);
    fill_mid(r, nlev, 5, 10000, pnh);
    fill_inc(r, nlev, 500000, 100, dphi);

    dfi d1(1);
    const auto w = d1.m_work;
    const auto
      dp3dw = dfi::get_slot(w, 0, 0),
      dphiw = dfi::get_slot(w, 0, 1),
      pnhw  = dfi::get_slot(w, 0, 2),
      dl = dfi::get_slot(w, 0, 3),
      d  = dfi::get_slot(w, 0, 4),
      du = dfi::get_slot(w, 0, 5);

    const auto f = KOKKOS_LAMBDA(const dfi::MT& t) {
      KernelVariables kv(t);
      dfi::transpose(kv, nlev, dp3d, dp3dw);
      dfi::transpose(kv, nlev, dphi, dphiw);
      dfi::transpose(kv, nlev, pnh, pnhw);
      kv.team_barrier();
      dfi::calc_jacobian(kv, dt, dp3dw, dphiw, pnhw, dl, d, du);
    };
    parallel_for(d1.m_policy, f); fence();

    const auto dlm = create_mirror_view(dl), dm = create_mirror_view(d),
      dum = create_mirror_view(du);
    deep_copy(dlm, dl); deep_copy(dm, d); deep_copy(dum, du);

    FA3 dp3df("dp3df", nlev), dphif("dphif", nlev), pnhf("pnhf", nlev);
    FA3d dlf("dlf", np, np, nlev-1), df("df", np, np, nlev), duf("duf", np, np, nlev-1);
    c2f(dp3d, dp3df);
    c2f(dphi, dphif);
    c2f(pnh, pnhf);
    get_dirk_jacobian_f90(dlf.data(), df.data(), duf.data(), dt,
                          dp3df.data(), dphif.data(), pnhf.data());

    const auto eq = [&] (const decltype(dm)& dc, const int os, const FA3d& df,
                         const int n) {
      for (int i = 0; i < np; ++i)
        for (int j = 0; j < np; ++j) {
          const auto
          idx = np*i + j,
          pi = idx / dfi::packn,
          si = idx % dfi::packn;
          for (int k = 0; k < n; ++k)
            REQUIRE(equal(dc(k+os,pi)[si], df(i,j,k)));
        }
    };
    eq(dlm, 1, dlf, nlev-1); eq(dm, 0, df, nlev); eq(dum, 0, duf, nlev-1);
  }

  SECTION ("pnh_and_exner_from_eos") {
    const int nlev = NUM_PHYSICAL_LEV, np = NP;

    ExecView<Scalar[NP][NP][NUM_LEV]> dp3d("dp3d"), dphi("dphi"), vtheta_dp("vtheta_dp");
    fill_inc(r, nlev, 5, 10000, dp3d);
    fill_inc(r, nlev, 500000, 100, dphi);
    fill_inc(r, nlev, 5*300, 10000*300, vtheta_dp);

    dfi d1(1);
    const auto w = d1.m_work;
    const auto
      dp3dw = dfi::get_slot(w, 0, 0),
      dphiw = dfi::get_slot(w, 0, 1),
      vtheta_dpw = dfi::get_slot(w, 0, 2),
      pnhw = dfi::get_slot(w, 0, 3),
      exnerw = dfi::get_slot(w, 0, 4),
      dpnh_dp_iw = dfi::get_slot(w, 0, 5);

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
    
    const auto dpnh_dp_im = create_mirror_view(dpnh_dp_iw);
    deep_copy(dpnh_dp_im, dpnh_dp_iw);

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
          REQUIRE(equal(dpnh_dp_im(k,pi)[si], dpnh_dp_if(k,i,j)));
      }
  }
}
