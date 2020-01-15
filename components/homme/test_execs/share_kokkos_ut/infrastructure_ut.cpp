#include <catch2/catch.hpp>

#include <iostream>

#include "ExecSpaceDefs.hpp"
#include "Types.hpp"

#include "utilities/SubviewUtils.hpp"
#include "utilities/TestUtils.hpp"
#include "utilities/BfbUtils.hpp"
#include "utilities/VectorUtils.hpp"

#include <random>

using namespace Homme;

extern "C" {
void run_bfb_pow_f90(Real* x, Real* e, Real* y, int& n);
}

// ====================== EXECUTION SPACE SETUP ====================== //

TEST_CASE("ExecSpaceDefs",
          "Test parallel machine parameterization.") {
  const int plev = 72;

  const auto test_basics = [=] (const ThreadPreferences& tp,
                                const std::pair<int, int>& tv) {
    REQUIRE(tv.first >= 1);
    REQUIRE(tv.second >= 1);
    if (tp.prefer_threads)
      REQUIRE((tv.first == tp.max_threads_usable || tv.second == 1));
    else
      REQUIRE((tv.first == 1 || tv.second > 1));
    if (tv.first * tv.second <= tp.max_threads_usable * tp.max_vectors_usable) {
      REQUIRE(tv.first <= tp.max_threads_usable);
      REQUIRE(tv.second <= tp.max_vectors_usable);
    }
  };

  SECTION("CPU/KNL") {
    for (int num_elem = 0; num_elem <= 3; ++num_elem) {
      for (int qsize = 1; qsize <= 30; ++qsize) {
        for (int pool = 1; pool <= 64; ++pool) {
          for (bool prefer_threads : {true, false}) {
            Homme::ThreadPreferences tp;
            tp.max_vectors_usable = plev;
            tp.prefer_threads = prefer_threads;
            const int npi = num_elem*qsize;
            const auto tv = Homme::Parallel::team_num_threads_vectors_from_pool(
              pool, npi, tp);
            // Correctness tests.
            test_basics(tp, tv);
            REQUIRE(tv.first * tv.second <= pool);
            // Tests for good behavior.
            if (npi >= pool)
              REQUIRE(tv.first * tv.second == 1);
          }
        }
      }
    }
    // Spot check some cases. Numbers are
    //     {#elem, qsize, pool, prefer_threads, #thread, #vector}.
    static const int cases[][6] = {{1, 1, 8, 1, 8, 1},
                                   {1, 30, 8, 1, 1, 1},
                                   {1, 1, 32, 1, 16, 2},
                                   {1, 1, 32, 0, 2, 16}};
    for (unsigned int i = 0; i < sizeof(cases)/sizeof(*cases); ++i) {
      const auto& c = cases[i];
      Homme::ThreadPreferences tp;
      tp.max_vectors_usable = plev/2;
      tp.prefer_threads = c[3];
      const auto tv = Homme::Parallel::team_num_threads_vectors_from_pool(
        c[2], c[0]*c[1], tp);
      REQUIRE(tv.first == c[4]);
      REQUIRE(tv.second == c[5]);
    }
  }

  SECTION("GPU") {
    static const int num_device_warps = 1792;
    static const int min_warps_per_team = 4, max_warps_per_team = 16;
    static const int num_threads_per_warp = 32;
    for (int num_elem = 0; num_elem <= 10000; num_elem += 10) {
      for (int qsize = 1; qsize <= 30; ++qsize) {
        for (bool prefer_threads : {true, false}) {
          Homme::ThreadPreferences tp;
          tp.prefer_threads = prefer_threads;
          tp.max_vectors_usable = plev;
          const int npi = num_elem*qsize;
          const auto tv = Homme::Parallel::team_num_threads_vectors_for_gpu(
            num_device_warps, num_threads_per_warp,
            min_warps_per_team, max_warps_per_team,
            npi, tp);
          // Correctness tests.
          test_basics(tp, tv);
          REQUIRE(Homme::Parallel::prevpow2(tv.second) == tv.second);
          REQUIRE(tv.first * tv.second >= num_threads_per_warp*min_warps_per_team);
          REQUIRE(tv.first * tv.second <= num_threads_per_warp*max_warps_per_team);
          // Tests for good behavior.
          REQUIRE(tv.first * tv.second >= min_warps_per_team * num_threads_per_warp);
          if (npi >= num_device_warps*num_threads_per_warp)
            REQUIRE(tv.first * tv.second == min_warps_per_team * num_threads_per_warp);
        }
      }
    }
    // Numbers are
    //     {#elem, qsize, prefer_threads, #thread, #vector}.
    static const int cases[][5] = {{1, 1, 1, 16, 32},
                                   {96, 30, 1, 16, 8},
                                   {96, 30, 0, 4, 32}};
    for (unsigned int i = 0; i < sizeof(cases)/sizeof(*cases); ++i) {
      const auto& c = cases[i];
      Homme::ThreadPreferences tp;
      tp.max_vectors_usable = plev/2;
      tp.prefer_threads = c[2];
      const auto tv = Homme::Parallel::team_num_threads_vectors_for_gpu(
        num_device_warps, num_threads_per_warp,
        min_warps_per_team, max_warps_per_team,
        c[0]*c[1], tp);
      REQUIRE(tv.first == c[3]);
      REQUIRE(tv.second == c[4]);
    }
  }
}

template <typename Dispatcher, int num_points, int scan_length>
void test_parallel_scan(
    Kokkos::TeamPolicy<ExecSpace,void> policy,
    ExecViewManaged<const Real * [num_points][scan_length]> input,
    ExecViewManaged<Real * [num_points][scan_length]> output) {
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const TeamMember & team) {
    assert(team.league_rank() < input.extent_int(0));
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, input.extent_int(1)),
                         [&](const int ip) {
      assert(ip < input.extent_int(1));
      auto input_vector = Homme::subview(input, team.league_rank(), ip);
      auto output_vector = Homme::subview(output, team.league_rank(), ip);

      auto inclusive_sum = [&](const int i, Real &accum, bool last) {
        accum += input_vector(i);
        if (last) {
          output_vector(i) = accum;
        }
      };
      Dispatcher::parallel_scan(team, input.extent_int(2), inclusive_sum);
    });
  });
}

template <typename ExeSpace> struct KokkosDispatcher {
  template <class Lambda>
  static KOKKOS_FORCEINLINE_FUNCTION void
  parallel_scan(const typename Kokkos::TeamPolicy<ExeSpace,void>::member_type &team,
                const int num_iters, const Lambda &lambda) {
    Kokkos::parallel_scan(Kokkos::ThreadVectorRange(team, num_iters), lambda);
  }
};

TEST_CASE("Parallel_scan",
          "Test parallel vector scan at ThreadVectorRange level.") {
  constexpr int num_elems = 10;
  constexpr int num_points = 16;
  constexpr int scan_length = 97;
  ExecViewManaged<Real * [num_points][scan_length]> input("", num_elems);
  ExecViewManaged<Real * [num_points][scan_length]> output_dispatch("",
                                                                    num_elems);
  ExecViewManaged<Real * [num_points][scan_length]> output_kokkos("",
                                                                  num_elems);

  // Policy used in all parallel for's below
  auto policy = Homme::get_default_team_policy<ExecSpace>(num_elems);
  policy.set_chunk_size(1);

  // Fill the input view.
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const TeamMember& team){
    const int ie = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, num_points),
                         [&](const int ip) {
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, scan_length),
                           [&](const int iv) {
        // We fill the vector entry iv with 1./(iv+1) - 1./(iv+2).
        // This way, the sum of entries [0,N] should return 1 - 1./(N+2).
        input(ie, ip, iv) = 1.0/(iv+1) - 1.0/(iv+2);
      });
    });
  });

  // Note: we want to check whether Dispatch<ExecSpace>::parallel_scan
  //       computes a correct parallel scan. We do not know if the macro
  //       HOMMEXX_BFB_TESTING is defined or not, so we do not know
  //       if we rely on the BFB version or the Kokkos one (which computes
  //       the sum with a parallel-friendly reorganization of sums).
  //       However, we know exactly the mathematical scan sum (see above),
  //       so we can check if, regardless of the HOMMEXX_BFB_TESTING
  //       macro, Dispatch<ExecSpace>::parallel_scan computes a scan sum
  //       within a given tolerance of the correct value

  // Computing a scan sum with Dispatch<ExecSpace> version
  // Note: we do not know if the BFB version is used or not. It depends
  //       on whether HOMMEXX_BFB_TESTING was defined BEFORE including
  //       the ExecSpaceDefs.hpp header.
  test_parallel_scan<Dispatch<ExecSpace>, num_points, scan_length>(
      policy, input, output_dispatch);

  test_parallel_scan<KokkosDispatcher<ExecSpace>, num_points, scan_length>(
      policy, input, output_kokkos);

  // The two versions *should* give the exact answer on CPU,
  // but differ on GPU. Either way, we should be within a good tolerance
  auto output_dispatch_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), output_dispatch);
  auto output_kokkos_h   = Kokkos::create_mirror_view_and_copy(HostMemSpace(), output_kokkos);

  const Real rel_threshold = 4.0 * std::numeric_limits<Real>::epsilon();
  for (int ie=0; ie<num_elems; ++ie) {
    for (int ip=0; ip<num_points; ++ip) {
      for (int iv=0; iv<scan_length; ++iv) {
        const Real computed_kokkos   = output_kokkos_h(ie, ip, iv);
        const Real computed_dispatch = output_dispatch_h(ie, ip, iv);
        const Real exact = 1.0 - 1.0/(iv+2);
        Real rel_error;

        rel_error = compare_answers(computed_kokkos, computed_dispatch);
        REQUIRE(rel_error<=rel_threshold);
        rel_error = compare_answers(exact, computed_kokkos);
        REQUIRE(rel_error<=rel_threshold);
        rel_error = compare_answers(exact, computed_dispatch);
        REQUIRE(rel_error<=rel_threshold);
      }
    }
  }
}

TEST_CASE("bfb_pow", "test taylor appx of pow function") {
  using Kokkos::create_mirror_view;

  using rngAlg = std::mt19937_64;
  std::random_device rd;
  const unsigned int catchRngSeed = Catch::rngSeed();
  const unsigned int seed = catchRngSeed==0 ? rd() : catchRngSeed;
  std::cout << "seed: " << seed << (catchRngSeed==0 ? " (catch rng seed was 0)\n" : "\n");
  rngAlg engine(seed);

  constexpr int num_tests = 10000;

  std::uniform_real_distribution<Real> pdf_exp(0.2, 2.0);
  std::uniform_real_distribution<Real> pdf_base(0.2, 2e6);

  auto cmvdc = [](const ExecViewManaged<Scalar> xd)->ExecViewManaged<Scalar>::HostMirror {
    auto xh = Kokkos::create_mirror_view(xd);
    Kokkos::deep_copy(xh,xd);
    return xh;
  };

  auto achoosek = [](const Real a, const int k) -> Real {
    // a choose k = [ a(a-1)(a-2)...(a-k+1) ] / k!
    //            =  a/k * [(a-1)/(k-1)] * ... * (a-k+2)/2 * (a-k+1)/1

    Real res = 1;
    for (int n=1; n<=k; ++n) {
      res *= (a-n+1)/(k-n+1);
    }
    return res;
  };

  ExecViewManaged<Scalar> b("b"), ycxx("ycxx");
  Real e;
  HostViewManaged<Scalar> yf90("yf90");
  for (int itest=0; itest<num_tests; ++itest) {
    genRandArray(b,engine,pdf_base);
    genRandArray(&e,1,engine,pdf_exp);

    // Run cxx version
    const auto f = KOKKOS_LAMBDA (const int /* idx */) {
      ycxx() = bfb_pow(b(),e);
    };

    Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0,1), f);
    Kokkos::fence();

    // Run f90 version
    auto bh = cmvdc(b);
    auto ycxxh = cmvdc(ycxx);

    Real* bhptr = reinterpret_cast<Real*>(bh.data());
    Real* yf90ptr = reinterpret_cast<Real*>(yf90.data());
    int vs = VECTOR_SIZE;
    run_bfb_pow_f90(bhptr,&e,yf90ptr,vs);

    for (int iv=0; iv<VECTOR_SIZE; ++iv) {
      if(yf90ptr[iv]!=ycxxh()[iv]) {
        printf("x: %3.17f\n",bh()[iv]);
        printf("a: %3.17f\n",e);
        printf("[cxx] x^a = %3.17f\n",ycxx()[iv]);
        printf("[f90] x^a = %3.17f\n",yf90ptr[iv]);
      }

      // Check f90 and cxx impl agree
      REQUIRE(yf90ptr[iv]==ycxxh()[iv]);

      Real acn = std::fabs(achoosek(e,11));
      // Check error bound on computed value via Lagrange remainder formula.
      // If yn ~ (1+x)^a, then x^a-R <= y^n <= x^a+R, where R is the bound on
      // the remainder, given by R = a choose N+1, with N being the number of
      // terms in the Taylor series that were kept (currently 10).
      if (bh()[iv]<2.0) {
        Real tmp = std::pow(bh()[iv],e);
        REQUIRE (ycxxh()[iv]>= tmp - acn);
        REQUIRE (ycxxh()[iv]>= tmp + acn);
      } else {
        // We computed yn = x^a = 1 / (1/x)^a, so the error bound is
        //          1/x^a-R <= yinv <= 1/x^a+R
        // To account for roundings, we throw in a factor of 2.0 on
        // the ubound and 0.5 on the lbound.

        constexpr auto tol = 2.0;
        Real tmp = std::pow(1.0/bh()[iv],e);
        if (1.0/ycxxh()[iv] > tol*(tmp + acn) ||
            1.0/ycxxh()[iv] < (tmp - acn)/tol) {
          printf("x: %3.17f\n",bh()[iv]);
          printf("a: %3.17f\n",e);
          printf("yn: %3.17f\n",ycxx()[iv]);
          printf("x^a: %3.17f\n",std::pow(bh()[iv],e));
          printf("tmp: %3.17f\n",tmp);
          printf("acn: %3.17f\n",acn);
        }
        REQUIRE (1.0/ycxxh()[iv] <= tol*(tmp + acn));
        REQUIRE (1.0/ycxxh()[iv] >= (tmp - acn)/tol);
      }
    }
  }
}
