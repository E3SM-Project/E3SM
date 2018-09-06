#include <catch2/catch.hpp>

#include "utilities/TestUtils.hpp"
#include "EulerStepFunctorImpl.hpp"

#include <limits>
#include <random>

using namespace Homme;
using rngAlg = std::mt19937_64;

extern "C" void limiter_optim_iter_full_c_callable(
  Real* ptens, const Real* sphweights, Real* minp, Real* maxp,
  const Real* dpmass);
extern "C" void limiter_clip_and_sum_c_callable(
  Real* ptens, const Real* sphweights, Real* minp, Real* maxp,
  const Real* dpmass);

struct LimiterTester {
  static constexpr Real eps { std::numeric_limits<Real>::epsilon() };

  struct FortranData {
    typedef Kokkos::LayoutLeft Layout;
    using ArrayPlvl = Kokkos::View<Real[NUM_PHYSICAL_LEV],
                                   Layout, HostMemSpace>;
    using ArrayGll = Kokkos::View<Real[NP][NP],
                                  Layout, HostMemSpace>;
    using ArrayGllPlvl = Kokkos::View<Real[NP][NP][NUM_PHYSICAL_LEV],
                                      Layout, HostMemSpace>;
    ArrayPlvl minp, maxp;
    ArrayGll sphweights;
    ArrayGllPlvl ptens, dpmass;

    FortranData ()
      : minp("minp"), maxp("maxp"), sphweights("sphweights"), ptens("ptens"),
        dpmass("dpmass")
    {}
  };

  using HostGll = HostViewManaged<Real[NP][NP]>;
  using HostGllLvl = HostViewManaged<Scalar[NP][NP][NUM_LEV]>;
  using Host2Lvl = HostViewManaged<Scalar[2][NUM_LEV]>;
  using HostPlvl = HostViewManaged<Real[NUM_PHYSICAL_LEV]>;

  HostGll sphweights;
  HostGllLvl dpmass, ptens;
  Host2Lvl qlim;

  using DevGll = ExecViewManaged<Real[NP][NP]>;
  using DevGllLvl = ExecViewManaged<Scalar[NP][NP][NUM_LEV]>;
  using Dev2Lvl = ExecViewManaged<Scalar[2][NUM_LEV]>;
  using Dev2GllLvl = ExecViewManaged<Scalar[2][NP][NP][NUM_LEV]>;

  DevGll sphweights_d;
  DevGllLvl dpmass_d, ptens_d;
  Dev2Lvl qlim_d;
  Dev2GllLvl rwrk_d;

  // For correctness checking.
  HostGllLvl ptens_orig;
  HostPlvl Qmass;

  void init () {
    sphweights = Kokkos::create_mirror_view(sphweights_d);
    dpmass = Kokkos::create_mirror_view(dpmass_d);
    ptens = Kokkos::create_mirror_view(ptens_d);
    qlim = Kokkos::create_mirror_view(qlim_d);
  }

  LimiterTester ()
    : sphweights_d("sphweights"), dpmass_d("dpmass"),
      ptens_d("ptens"), qlim_d("qlim"), rwrk_d("rwrk"),
      ptens_orig("ptens_orig"), Qmass("Qmass")
  { init(); }

  void deep_copy (const LimiterTester& lv) {
    Kokkos::deep_copy(sphweights, lv.sphweights);
    Kokkos::deep_copy(dpmass, lv.dpmass);
    Kokkos::deep_copy(ptens, lv.ptens);
    Kokkos::deep_copy(qlim, lv.qlim);
    Kokkos::deep_copy(Qmass, lv.Qmass);
    Kokkos::deep_copy(ptens_orig, lv.ptens_orig);
    Kokkos::deep_copy(sphweights_d, lv.sphweights_d);
    Kokkos::deep_copy(dpmass_d, lv.dpmass_d);
    Kokkos::deep_copy(ptens_d, lv.ptens_d);
    Kokkos::deep_copy(qlim_d, lv.qlim_d);
  }

  void todevice () {
    Kokkos::deep_copy(sphweights_d, sphweights);
    Kokkos::deep_copy(dpmass_d, dpmass);
    Kokkos::deep_copy(ptens_d, ptens);
    Kokkos::deep_copy(qlim_d, qlim);
  }

  void fromdevice () {
    Kokkos::deep_copy(sphweights, sphweights_d);
    Kokkos::deep_copy(dpmass, dpmass_d);
    Kokkos::deep_copy(ptens, ptens_d);
    Kokkos::deep_copy(qlim, qlim_d);
  }

  template <typename Array>
  static void urand (Array& a, Real lo, Real hi) {
    std::random_device rd;
    rngAlg engine(rd());
    genRandArray(a, engine, std::uniform_real_distribution<Real>(lo, hi));
  }

  void init_feasible () {
    urand(sphweights, 1.0/16, 2.0/16);
    urand(dpmass, 0.5, 1);
    urand(ptens, 0, 1);
    urand(qlim, 0, 1);

    for (int k = 0; k < NUM_PHYSICAL_LEV; ++k) {
      const int vi = k / VECTOR_SIZE, si = k % VECTOR_SIZE;

      // Order q limits.
      auto q0 = qlim(0,vi)[si], q1 = qlim(1,vi)[si];
      if (q1 < q0) {
        std::swap(q0, q1);
        qlim(0,vi)[si] = q0;
        qlim(1,vi)[si] = q1;
      }

      // Turn ptens into density.
      for (int i = 0; i < NP; ++i)
        for (int j = 0; j < NP; ++j)
          ptens(i,j,vi)[si] *= dpmass(i,j,vi)[si];

      // Make problem feasible by adjusting cell mass.
      Real m = 0, lo = 0, hi = 0;
      for (int i = 0; i < NP; ++i)
        for (int j = 0; j < NP; ++j) {
          m += sphweights(i,j) * ptens(i,j,vi)[si];
          lo += sphweights(i,j) * qlim(0,vi)[si] * dpmass(i,j,vi)[si];
          hi += sphweights(i,j) * qlim(1,vi)[si] * dpmass(i,j,vi)[si];
        }
      const Real dm = m < lo ? lo - m : m > hi ? hi - m : 0;
      if (dm)
        for (int i = 0; i < NP; ++i)
          for (int j = 0; j < NP; ++j)
            ptens(i,j,vi)[si] += (1 + 1e2*eps)*dm/(NP*NP*sphweights(i,j));
      m = 0;
      for (int i = 0; i < NP; ++i)
        for (int j = 0; j < NP; ++j)
          m += sphweights(i,j) * ptens(i,j,vi)[si];
      REQUIRE(m >= (1 - 1e1*eps)*lo);
      REQUIRE(m <= (1 + 1e1*eps)*hi);

      // For checking, record original values and true mass.
      Kokkos::deep_copy(ptens_orig, ptens);
      Qmass(k) = m;
    }

    todevice();
  }

  // Trigger the .not.modified cycle line in limiter_clip_and_sum.
  void init_to_trigger_cycle () {
    init_feasible();

    const int k = 0;
    const int vi = k / VECTOR_SIZE, si = k % VECTOR_SIZE;

    Real minq = 2, maxq = 0;
    for (int i = 0; i < NP; ++i)
      for (int j = 0; j < NP; ++j) {
        const Real q = ptens(i,j,vi)[si] / dpmass(i,j,vi)[si];
        minq = std::min(minq, q);
        maxq = std::max(maxq, q);
      }
    qlim(0,vi)[si] = std::max(0.0, minq - 1e-2);
    qlim(1,vi)[si] = std::min(1.0, maxq + 1e-2);

    Kokkos::deep_copy(qlim_d, qlim);
  }

  void fill_fortran (FortranData& d) {
    for (int k = 0; k < NUM_PHYSICAL_LEV; ++k) {
      const int vi = k / VECTOR_SIZE, si = k % VECTOR_SIZE;
      d.minp(k) = qlim(0,vi)[si];
      d.maxp(k) = qlim(1,vi)[si];
      for (int i = 0; i < NP; ++i)
        for (int j = 0; j < NP; ++j) {
          d.ptens(j,i,k) = ptens(i,j,vi)[si];
          d.dpmass(j,i,k) = dpmass(i,j,vi)[si];
        }
    }
    for (int i = 0; i < NP; ++i)
      for (int j = 0; j < NP; ++j)
        d.sphweights(j,i) = sphweights(i,j);
  }

  void compare_with_fortran (const int limiter_option, FortranData& d) {
    if (limiter_option == 9)
      limiter_clip_and_sum_c_callable(d.ptens.data(), d.sphweights.data(),
                                      d.minp.data(), d.maxp.data(),
                                      d.dpmass.data());
    else
      limiter_optim_iter_full_c_callable(d.ptens.data(), d.sphweights.data(),
                                         d.minp.data(), d.maxp.data(),
                                         d.dpmass.data());      
    for (int k = 0; k < NUM_PHYSICAL_LEV; ++k) {
      const int vi = k / VECTOR_SIZE, si = k % VECTOR_SIZE;
      for (int i = 0; i < NP; ++i)
        for (int j = 0; j < NP; ++j)
          REQUIRE(d.ptens(j,i,k) == ptens(i,j,vi)[si]);
    }
  }

  size_t team_shmem_size (const int team_size) const {
    return Homme::EulerStepFunctorImpl::limiter_team_shmem_size(team_size);
  }

  struct Lim8 {};
  KOKKOS_INLINE_FUNCTION void operator() (const Lim8&, const Homme::TeamMember& team) const {
    Homme::EulerStepFunctorImpl
      ::limiter_optim_iter_full(team, sphweights_d, dpmass_d, qlim_d, ptens_d);
  }

  struct CAAS {};
  KOKKOS_INLINE_FUNCTION void operator() (const CAAS&, const Homme::TeamMember& team) const {
    Homme::EulerStepFunctorImpl
      ::limiter_clip_and_sum(team, sphweights_d, dpmass_d, qlim_d, ptens_d);
  }

  struct SerLim8 {};
  KOKKOS_INLINE_FUNCTION void operator() (const SerLim8&, const Homme::TeamMember& team) const {
    Homme::SerialLimiter<ExecSpace>
      ::run<8>(sphweights_d, dpmass_d, qlim_d, ptens_d, rwrk_d);
  }

  struct SerCAAS {};
  KOKKOS_INLINE_FUNCTION void operator() (const SerCAAS&, const Homme::TeamMember& team) const {
    Homme::SerialLimiter<ExecSpace>
      ::run<9>(sphweights_d, dpmass_d, qlim_d, ptens_d, rwrk_d);
  }

  void check () {
    for (int k = 0; k < NUM_PHYSICAL_LEV; ++k) {
      const int vi = k / VECTOR_SIZE, si = k % VECTOR_SIZE;
      Real m = 0;
      for (int i = 0; i < NP; ++i)
        for (int j = 0; j < NP; ++j) {
          // Check that the mixing ratio is limited.
          REQUIRE(ptens(i,j,vi)[si] >=
                  (1 - 1e1*eps)*qlim(0,vi)[si]*dpmass(i,j,vi)[si]);
          REQUIRE(ptens(i,j,vi)[si] <=
                  (1 + 1e1*eps)*qlim(1,vi)[si]*dpmass(i,j,vi)[si]);
          m += sphweights(i,j) * ptens(i,j,vi)[si];
        }
      // Check mass conservation.
      REQUIRE(std::abs(m - Qmass(k)) <= 1e2*eps*Qmass(k));
    }
  }

  void check_same (const LimiterTester& ref) {
    for (int k = 0; k < NUM_PHYSICAL_LEV; ++k) {
      const int vi = k / VECTOR_SIZE, si = k % VECTOR_SIZE;
      for (int i = 0; i < NP; ++i)
        for (int j = 0; j < NP; ++j) {
          // Check BFB.
          REQUIRE(ptens(i,j,vi)[si] == ref.ptens(i,j,vi)[si]);
        }
    }
  }
};

void test_limiter (const int limiter_option, const int impl, const int init) {
  if (impl == 1 && OnGpu<ExecSpace>::value) return;
  std::cout << "test limiter " << limiter_option << ","
            << (impl == 0 ? " Kokkos impl," : " serial impl,")
            << (init == 0 ? " with feasible\n" : " with early exit\n");
  LimiterTester lv;
  if (init == 0)
    lv.init_feasible();
  else
    lv.init_to_trigger_cycle();
  LimiterTester::FortranData fd;
  lv.fill_fortran(fd);
  if (limiter_option == 8) {
    if (impl == 0)
      Kokkos::parallel_for(Homme::get_default_team_policy<ExecSpace, LimiterTester::Lim8>(1), lv);
    else
      Kokkos::parallel_for(Kokkos::TeamPolicy<ExecSpace, LimiterTester::SerLim8>(1, 1, 1), lv);    
  } else {
    if (impl == 0)
      Kokkos::parallel_for(Homme::get_default_team_policy<ExecSpace, LimiterTester::CAAS>(1), lv);
    else
      Kokkos::parallel_for(Kokkos::TeamPolicy<ExecSpace, LimiterTester::SerCAAS>(1, 1, 1), lv);
  }
  lv.fromdevice();
  lv.check();
  lv.compare_with_fortran(limiter_option, fd);
}

TEST_CASE("lim=8 math correctness", "limiter") {
  for (int init = 0; init < 2; ++init)
    test_limiter(8, 0, init);
}

TEST_CASE("CAAS math correctness", "limiter") {
  for (int init = 0; init < 2; ++init)
    test_limiter(9, 0, init);
}

TEST_CASE("serial lim=8 math correctness", "limiter") {
  for (int init = 0; init < 2; ++init)
    test_limiter(8, 1, init);
}

TEST_CASE("serial CAAS math correctness", "limiter") {
  for (int init = 0; init < 2; ++init)
    test_limiter(9, 1, init);
}

// The best thing to do here would be to check the KKT conditions for the
// 1-norm-minimization problem. That's a fair bit of programming. Instead, check
// that two very different methods that ought to provide 1-norm-minimal
// corrections have corrections that have the same 1-norm.
TEST_CASE("1-norm minimal", "limiters") {
  LimiterTester lv;
  lv.init_feasible();
  LimiterTester lv_deepcopy;
  lv_deepcopy.deep_copy(lv);
  std::vector<Real> lim8norm1(NUM_PHYSICAL_LEV), othernorm1(NUM_PHYSICAL_LEV);

  auto get_norm1 = [&] (const LimiterTester& lv, std::vector<Real>& n1s) {
    for (int k = 0; k < NUM_PHYSICAL_LEV; ++k) {
      const int vi = k / VECTOR_SIZE, si = k % VECTOR_SIZE;
      Real n1 = 0;
      for (int i = 0; i < NP; ++i)
        for (int j = 0; j < NP; ++j)
          n1 += std::abs(lv.ptens(i,j,vi)[si] - lv.ptens_orig(i,j,vi)[si])*lv.sphweights(i,j);
      n1s[k] = n1;
    }
  };

  Kokkos::parallel_for(Homme::get_default_team_policy<ExecSpace, LimiterTester::Lim8>(1), lv);
  lv.fromdevice();
  get_norm1(lv, lim8norm1);

  std::vector<LimiterTester> lts;
  lts.push_back(lv);
  for (int other = 0; other < (OnGpu<ExecSpace>::value ? 1 : 3); ++other) {
    lts.push_back(LimiterTester());
    auto& lv_other = lts.back();
    lv_other.deep_copy(lv_deepcopy);
    switch (other) {
    case 0:
      Kokkos::parallel_for(Homme::get_default_team_policy<ExecSpace, LimiterTester::CAAS>(1), lv_other);
      break;
    case 1:
      Kokkos::parallel_for(Kokkos::TeamPolicy<ExecSpace, LimiterTester::SerLim8>(1, 1, 1), lv_other);
      break;
    case 2:
      Kokkos::parallel_for(Kokkos::TeamPolicy<ExecSpace, LimiterTester::SerCAAS>(1, 1, 1), lv_other);
      break;
    }
    lv_other.fromdevice();
    get_norm1(lv_other, othernorm1);
    for (int k = 0; k < NUM_PHYSICAL_LEV; ++k)
      REQUIRE(std::abs(lim8norm1[k] - othernorm1[k]) <= 1e3*LimiterTester::eps*othernorm1[k]);
  }

  if ( ! OnGpu<ExecSpace>::value) {
    lts[2].check_same(lts[0]);
    lts[3].check_same(lts[1]);
  }
}
