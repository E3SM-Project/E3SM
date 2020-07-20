#include <catch2/catch.hpp>

#include "tridiag_tests.hpp"

extern "C" {
  void tridiag_diagdom_bfb_a1x1(int n, scream::Real* dl, scream::Real* d,
                                scream::Real* du, scream::Real* x);
  void tridiag_diagdom_bfb_a1xm(int n, int nrhs, scream::Real* dl, scream::Real* d,
                                scream::Real* du, scream::Real* x);
}

namespace scream {
namespace tridiag {
namespace test {
namespace correct {

struct Solver {
  enum Enum { thomas_team_scalar, thomas_team_pack,
              thomas_scalar, thomas_pack,
              cr_scalar, bfb, bfbf90,
              error };

  static std::string convert (Enum e) {
    switch (e) {
    case thomas_team_scalar: return "thomas_team_scalar";
    case thomas_team_pack: return "thomas_team_pack";
    case thomas_scalar: return "thomas_scalar";
    case thomas_pack: return "thomas_pack";
    case cr_scalar: return "cr_scalar";
    case bfb: return "bfb";
    case bfbf90: return "bfbf90";
    default: scream_require_msg(false, "Not a valid solver: " << e);
    }
  }

  static Enum convert (const std::string& s) {
    if (s == "thomas_team_scalar") return thomas_team_scalar;
    if (s == "thomas_team_pack") return thomas_team_pack;
    if (s == "thomas_scalar") return thomas_scalar;
    if (s == "thomas_pack") return thomas_pack;
    if (s == "cr_scalar") return cr_scalar;
    if (s == "bfb") return bfb;
    if (s == "bfbf90") return bfbf90;
    return error;
  }

  static Enum all[];
};

Solver::Enum Solver::all[] = { thomas_team_scalar, thomas_team_pack,
                               thomas_scalar, thomas_pack,
                               cr_scalar, bfb, bfbf90 };

struct TestConfig {
  using TeamLayout = Kokkos::LayoutRight;

  Solver::Enum solver;
  int n_kokkos_thread, n_kokkos_vec;
};

template <typename TridiagArray>
KOKKOS_INLINE_FUNCTION
Kokkos::View<typename TridiagArray::value_type*>
get_diag (const TridiagArray& A, const int diag_idx) {
  assert(A.extent_int(2) == 1);
  return Kokkos::View<typename TridiagArray::value_type*>(
    &A.impl_map().reference(diag_idx, 0, 0),
    A.extent_int(1));
}

template <typename TridiagArray>
Kokkos::View<typename TridiagArray::value_type**, TestConfig::TeamLayout>
KOKKOS_INLINE_FUNCTION
get_diags (const TridiagArray& A, const int diag_idx) {
  return Kokkos::View<typename TridiagArray::value_type**, TestConfig::TeamLayout>(
    &A.impl_map().reference(diag_idx, 0, 0),
    A.extent_int(1), A.extent_int(2));
}

template <typename DataArray>
KOKKOS_INLINE_FUNCTION
Kokkos::View<typename DataArray::value_type*>
get_x (const DataArray& X) {
  assert(X.extent_int(1) == 1);
  return Kokkos::View<typename DataArray::value_type*>(
    &X.impl_map().reference(0, 0), X.extent_int(0));
}

template <typename Scalar>
using TridiagArray = Kokkos::View<Scalar***, TestConfig::TeamLayout>;
template <typename Scalar>
using DataArray = Kokkos::View<Scalar**, TestConfig::TeamLayout>;

template <bool same_pack_size, typename APack, typename DataPack>
struct Solve;

template <typename APack, typename DataPack>
struct Solve<true, APack, DataPack> {
  static void run (const TestConfig& tc, TridiagArray<APack>& A, DataArray<DataPack>& X,
                   const int nprob, const int nrhs) {
    using Kokkos::subview;
    using Kokkos::ALL;
    using scream::pack::scalarize;

    assert(nrhs > 1 || DataPack::n == 1);
    assert(nrhs > 1 || X.extent_int(2) == 1);
    assert(nprob > 1 || APack::n == 1);
    assert(nprob > 1 || A.extent_int(2) == 1);

    using TeamPolicy = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>;
    using MT = typename TeamPolicy::member_type;
    TeamPolicy policy(1, tc.n_kokkos_thread, tc.n_kokkos_vec);

    switch (tc.solver) {
    case Solver::thomas_team_scalar:
    case Solver::thomas_team_pack: {
      const auto As = scalarize(A);
      const auto Xs = scalarize(X);
      const auto f = KOKKOS_LAMBDA (const MT& team) {
        const auto dl = get_diag(As, 0);
        const auto d  = get_diag(As, 1);
        const auto du = get_diag(As, 2);
        if (tc.solver == Solver::thomas_team_scalar)
          scream::tridiag::thomas(team, dl, d, du, Xs);
        else
          scream::tridiag::thomas(team, dl, d, du, X);
      };
      Kokkos::parallel_for(policy, f);
    } break;
    case Solver::thomas_scalar: {
      if (nprob == 1) {
        if (nrhs == 1) {
          const auto As = scalarize(A);
          const auto Xs = scalarize(X);
          const auto f = KOKKOS_LAMBDA (const MT& team) {
            const auto single = [&] () {
              const auto dl = get_diag(As, 0);
              const auto d  = get_diag(As, 1);
              const auto du = get_diag(As, 2);
              const auto x  = get_x(Xs);
              scream::tridiag::thomas(dl, d, du, x);
            };
            Kokkos::single(Kokkos::PerTeam(team), single);
          };
          Kokkos::parallel_for(policy, f);
        } else {
          const auto As = scalarize(A);
          const auto Xs = scalarize(X);
          const auto f = KOKKOS_LAMBDA (const MT& team) {
            const auto single = [&] () {
              const auto dl = get_diag(As, 0);
              const auto d  = get_diag(As, 1);
              const auto du = get_diag(As, 2);
              scream::tridiag::thomas(dl, d, du, Xs);
            };
            Kokkos::single(Kokkos::PerTeam(team), single);
          };
          Kokkos::parallel_for(policy, f);
        }
      } else {
        const auto As = scalarize(A);
        const auto Xs = scalarize(X);
        const auto f = KOKKOS_LAMBDA (const MT& team) {
          const auto single = [&] () {
            const auto dl = get_diags(As, 0);
            const auto d  = get_diags(As, 1);
            const auto du = get_diags(As, 2);
            scream::tridiag::thomas(dl, d, du, Xs);
          };
          Kokkos::single(Kokkos::PerTeam(team), single);
        };
        Kokkos::parallel_for(policy, f);
      }
    } break;
    case Solver::thomas_pack: {
      if (nprob == 1) {
        const auto As = scalarize(A);
        const auto f = KOKKOS_LAMBDA (const MT& team) {
          const auto single = [&] () {
            const auto dl = get_diag(As, 0);
            const auto d  = get_diag(As, 1);
            const auto du = get_diag(As, 2);
            scream::tridiag::thomas(dl, d, du, X);
          };
          Kokkos::single(Kokkos::PerTeam(team), single);
        };
        Kokkos::parallel_for(policy, f);
      } else {
        const auto f = KOKKOS_LAMBDA (const MT& team) {
          const auto single = [&] () {
            const auto dl = get_diags(A, 0);
            const auto d  = get_diags(A, 1);
            const auto du = get_diags(A, 2);
            scream::tridiag::thomas(dl, d, du, X);
          };
          Kokkos::single(Kokkos::PerTeam(team), single);
        };
        Kokkos::parallel_for(policy, f);
      }
    } break;
    case Solver::cr_scalar: {
      if (nprob == 1) {
        if (nrhs == 1) {
          const auto As = scalarize(A);
          const auto Xs = scalarize(X);
          const auto f = KOKKOS_LAMBDA (const MT& team) {
            const auto dl = get_diag(As, 0);
            const auto d  = get_diag(As, 1);
            const auto du = get_diag(As, 2);
            const auto x  = get_x(Xs);
            scream::tridiag::cr(team, dl, d, du, x);
          };
          Kokkos::parallel_for(policy, f);
        } else {
          const auto As = scalarize(A);
          const auto Xs = scalarize(X);
          const auto f = KOKKOS_LAMBDA (const MT& team) {
            const auto dl = get_diag(As, 0);
            const auto d  = get_diag(As, 1);
            const auto du = get_diag(As, 2);
            scream::tridiag::cr(team, dl, d, du, Xs);
          }; 
         Kokkos::parallel_for(policy, f);
        }
      } else {
        const auto As = scalarize(A);
        const auto Xs = scalarize(X);
        const auto f = KOKKOS_LAMBDA (const MT& team) {
          const auto dl = get_diags(As, 0);
          const auto d  = get_diags(As, 1);
          const auto du = get_diags(As, 2);
          scream::tridiag::cr(team, dl, d, du, Xs);
        };
        Kokkos::parallel_for(policy, f);
      }
    } break;
    case Solver::bfb: {
      const auto f = KOKKOS_LAMBDA (const MT& team) {
        const auto dl = get_diags(A, 0);
        const auto d  = get_diags(A, 1);
        const auto du = get_diags(A, 2);
        scream::tridiag::bfb(team, dl, d, du, X);
      };
      Kokkos::parallel_for(policy, f);
    } break;
    case Solver::bfbf90: {
      const auto Am = create_mirror_view(A);
      const auto Xm = create_mirror_view(X);
      deep_copy(Am, A);
      deep_copy(Xm, X);
      const auto As = scalarize(Am);
      const auto Xs = scalarize(Xm);
      const auto dl = get_diag(As, 0);
      const auto d  = get_diag(As, 1);
      const auto du = get_diag(As, 2);
      if (nprob == 1) {
        if (nrhs == 1) {
          const auto x  = get_x(Xs);
          tridiag_diagdom_bfb_a1x1(d.extent_int(0), dl.data(), d.data(),
                                   du.data(), x.data());
        } else {
          tridiag_diagdom_bfb_a1xm(d.extent_int(0), Xs.extent_int(1),
                                   dl.data(), d.data(), du.data(), Xs.data());
        }
      } else {
        scream_require_msg(false, "bfbf90 does not support nprob > 1");
      }
      deep_copy(A, Am);
      deep_copy(X, Xm);
    } break;
    default:
      scream_require_msg(false, "Same pack size: " << Solver::convert(tc.solver));
    }
  }
};

template <typename APack, typename DataPack>
struct Solve<false, APack, DataPack> {
  static void run (const TestConfig& tc, TridiagArray<APack>& A, DataArray<DataPack>& X,
                   const int nprob, const int nrhs) {
    using Kokkos::subview;
    using Kokkos::ALL;
    using scream::pack::scalarize;

    using TeamPolicy = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>;
    using MT = typename TeamPolicy::member_type;
    TeamPolicy policy(1, tc.n_kokkos_thread, tc.n_kokkos_vec);

    assert(nrhs > 1 || DataPack::n == 1);
    assert(nrhs > 1 || X.extent_int(2) == 1);
    assert(nprob == 1);

    switch (tc.solver) {
    case Solver::thomas_team_pack: {
      const auto As = scalarize(A);
      const auto Xs = scalarize(X);
      const auto f = KOKKOS_LAMBDA (const MT& team) {
        const auto dl = get_diag(As, 0);
        const auto d  = get_diag(As, 1);
        const auto du = get_diag(As, 2);
        if (tc.solver == Solver::thomas_team_scalar)
          scream::tridiag::thomas(team, dl, d, du, Xs);
        else
          scream::tridiag::thomas(team, dl, d, du, X);
      };
      Kokkos::parallel_for(policy, f);
    } break;
    case Solver::thomas_pack: {
      const auto As = scalarize(A);
      const auto f = KOKKOS_LAMBDA (const MT& team) {
        const auto single = [&] () {
          const auto dl = get_diag(As, 0);
          const auto d  = get_diag(As, 1);
          const auto du = get_diag(As, 2);
          scream::tridiag::thomas(dl, d, du, X);
        };
        Kokkos::single(Kokkos::PerTeam(team), single);
      };
      Kokkos::parallel_for(policy, f);
    } break;
    case Solver::bfb: {
      const auto As = scalarize(A);
      const auto f = KOKKOS_LAMBDA (const MT& team) {
        const auto dl = get_diags(As, 0);
        const auto d  = get_diags(As, 1);
        const auto du = get_diags(As, 2);
        scream::tridiag::bfb(team, dl, d, du, X);
      };
      Kokkos::parallel_for(policy, f);
    } break;
    default:
      scream_require_msg(false, "Different pack size: " << Solver::convert(tc.solver));
    }
  }
};

template <typename APack, typename DataPack>
struct Data {
  const int nprob, nrhs;
  TridiagArray<APack> A, Acopy;
  DataArray<DataPack> B, X, Y;

  Data (const int nrow, const int nprob_, const int nrhs_)
    : nprob(nprob_), nrhs(nrhs_),
      A("A", 3, nrow, scream::pack::npack<APack>(nprob)),
      Acopy("A", A.extent(0), A.extent(1), A.extent(2)),
      B("B", nrow, scream::pack::npack<DataPack>(nrhs)),
      X("X", B.extent(0), B.extent(1)),
      Y("Y", X.extent(0), X.extent(1))
  {}
};

template <typename APack, typename DataPack>
void fill (Data<APack, DataPack>& dt) {
  using Kokkos::create_mirror_view;
  using Kokkos::deep_copy;
  using Kokkos::subview;
  using Kokkos::ALL;

  const auto Am = create_mirror_view(dt.A);
  const auto Bm = create_mirror_view(dt.B);
  {
    const auto As = scalarize(Am);
    const auto dl = subview(As, 0, ALL(), ALL());
    const auto d  = subview(As, 1, ALL(), ALL());
    const auto du = subview(As, 2, ALL(), ALL());
    fill_tridiag_matrix(dl, d, du, dt.nprob, dt.nrhs /* seed */);
  }
  fill_data_matrix(scalarize(Bm), dt.nrhs);
  deep_copy(dt.A, Am);
  deep_copy(dt.B, Bm);
  deep_copy(dt.Acopy, dt.A);
  deep_copy(dt.X, dt.B);
}

template <typename APack, typename DataPack>
Real relerr (Data<APack, DataPack>& dt) {
  using Kokkos::create_mirror_view;
  using Kokkos::deep_copy;
  using Kokkos::subview;
  using Kokkos::ALL;

  const auto Acopym = create_mirror_view(dt.Acopy);
  const auto As = scalarize(Acopym);
  const auto dl = subview(As, 0, ALL(), ALL());
  const auto d  = subview(As, 1, ALL(), ALL());
  const auto du = subview(As, 2, ALL(), ALL());
  const auto Xm = create_mirror_view(dt.X);
  const auto Ym = create_mirror_view(dt.Y);
  deep_copy(Acopym, dt.Acopy);
  deep_copy(Xm, dt.X);
  matvec(dl, d, du, scalarize(Xm), scalarize(Ym), dt.nprob, dt.nrhs);
  const auto Bm = create_mirror_view(dt.B);
  deep_copy(Bm, dt.B);
  const auto re = reldif(scalarize(Bm), scalarize(Ym), dt.nrhs);
  return re;
}

template <typename Fn>
void run_test_configs (Fn& fn) {
  for (const auto solver : Solver::all) {
    TestConfig tc;
    tc.solver = solver;
    if (scream::util::OnGpu<Kokkos::DefaultExecutionSpace>::value) {
      tc.n_kokkos_vec = 1;
      for (const int n_kokkos_thread : {128, 256, 512}) {
        tc.n_kokkos_thread = n_kokkos_thread;
        fn(tc);
      }
      tc.n_kokkos_vec = 32;
      for (const int n_kokkos_thread : {4}) {
        tc.n_kokkos_thread = n_kokkos_thread;
        fn(tc);
      }
      tc.n_kokkos_vec = 8;
      for (const int n_kokkos_thread : {16}) {
        tc.n_kokkos_thread = n_kokkos_thread;
        fn(tc);
      }
    } else {
      const int concurrency = Kokkos::DefaultExecutionSpace::concurrency();
      const int n_kokkos_thread = concurrency;
      for (const int n_kokkos_vec : {1, 2}) {
        tc.n_kokkos_thread = n_kokkos_thread;
        tc.n_kokkos_vec = n_kokkos_vec;
        fn(tc);
      }
    }
  }
}

template <int A_pack_size, int data_pack_size>
void run_property_test_on_config (const TestConfig& tc) {
  using namespace scream::tridiag::test;

  using scream::Real;
  using APack = scream::pack::Pack<Real, A_pack_size>;
  using DataPack = scream::pack::Pack<Real, data_pack_size>;

#if 1
  const int nrows[] = {1,2,3,4,5, 8,10,16, 32,43, 63,64,65, 111,128,129, 2048};
  const int nrhs_max = 60;
  const int nrhs_inc = 11;
#else
  const int nrows[] = {10};
  const int nrhs_max = 6;
  const int nrhs_inc = 5;  
#endif

  for (const int nrow : nrows) {
    for (int nrhs = 1; nrhs <= nrhs_max; nrhs += nrhs_inc) {
      for (const bool A_many : {false, true}) {
        if (nrhs == 1 && A_many) continue;
        const int nprob = A_many ? nrhs : 1;

        // Skip unsupported solver-problem format combinations.
        if ((tc.solver == Solver::thomas_team_scalar ||
             tc.solver == Solver::thomas_team_pack)
            && nprob > 1)
          continue;

        // Skip combinations generated at this and higher levels that Solve::run
        // doesn't support to reduce redundancies.
        if ((nrhs  == 1 && data_pack_size > 1) ||
            (nprob == 1 && A_pack_size    > 1))
          continue;
        if ((tc.solver == Solver::thomas_team_scalar ||
             tc.solver == Solver::thomas_scalar ||
             tc.solver == Solver::cr_scalar ||
             tc.solver == Solver::bfbf90) &&
            data_pack_size > 1)
          continue;
        if (tc.solver == Solver::bfbf90 && nprob > 1)
          continue;
        if (static_cast<int>(APack::n) != static_cast<int>(DataPack::n) && nprob > 1)
          continue;

        Data<APack, DataPack> dt(nrow, nprob, nrhs);
        fill(dt);

        Solve<A_pack_size == data_pack_size, APack, DataPack>
          ::run(tc, dt.A, dt.X, nprob, nrhs);

        const auto re = relerr(dt);
        const bool pass = re <= 50*std::numeric_limits<Real>::epsilon();
        if ( ! pass) {
          std::stringstream ss;
          ss << Solver::convert(tc.solver) << " " << tc.n_kokkos_thread
             << " " << tc.n_kokkos_vec << " | " << nrow << " " << nrhs << " "
             << A_many << " | log10 reldif " << std::log10(re);
          std::cout << "FAIL: " << ss.str() << "\n";
        }
        REQUIRE(pass);
      }
    }
  }
}

template <int A_pack_size, int data_pack_size>
void run_property_test () {
  run_test_configs(run_property_test_on_config<A_pack_size, data_pack_size>);
}

template <int A_pack_size, int data_pack_size>
void run_bfb_test_on_config (TestConfig& tc) {
  using namespace scream::tridiag::test;

  using scream::Real;
  using APack = scream::pack::Pack<Real, A_pack_size>;
  using DataPack = scream::pack::Pack<Real, data_pack_size>;

  if (tc.solver != Solver::bfb) return;

  const int nrows[] = {1,2,3,4,5, 8,10,16, 32,43, 63,64,65, 111,128,129, 2048};
  const int nrhs_max = 60;
  const int nrhs_inc = 11;

  for (const int nrow : nrows) {
    for (int nrhs = 1; nrhs <= nrhs_max; nrhs += nrhs_inc) {
      for (const bool A_many : {false}) {
        if (nrhs == 1 && A_many) continue;
        const int nprob = A_many ? nrhs : 1;

        if ((nrhs  == 1 && data_pack_size > 1) ||
            (nprob == 1 && A_pack_size    > 1))
          continue;

        Data<APack, DataPack> dt1(nrow, nprob, nrhs);
        fill(dt1);
        tc.solver = Solver::bfb;
        Solve<A_pack_size == data_pack_size, APack, DataPack>
          ::run(tc, dt1.A, dt1.X, nprob, nrhs);

        Data<APack, DataPack> dt2(nrow, nprob, nrhs);
        fill(dt2);
        tc.solver = Solver::bfbf90;
        Solve<A_pack_size == data_pack_size, APack, DataPack>
          ::run(tc, dt2.A, dt2.X, nprob, nrhs);

        const auto X1s = scalarize(dt1.X);
        const auto X2s = scalarize(dt2.X);
        const auto X1 = create_mirror_view(X1s);
        const auto X2 = create_mirror_view(X2s);
        deep_copy(X1, X1s);
        deep_copy(X2, X2s);

        int nerr = 0; // don't blow up the assertion count
        for (int i = 0; i < X1.extent_int(0); ++i)
          for (int j = 0; j < X1.extent_int(1); ++j)
            if (X1(i,j) != X2(i,j)) ++nerr;
        REQUIRE(nerr == 0);
      }
    }
  }
}

template <int A_pack_size, int data_pack_size>
void run_bfb_test () {
  run_test_configs(run_bfb_test_on_config<A_pack_size, data_pack_size>);
}

} // namespace correct
} // namespace test
} // namespace tridiag
} // namespace scream

TEST_CASE("property", "tridiag") {
  scream::tridiag::test::correct::run_property_test<1,1>();
  if (EKAT_PACK_SIZE > 1) {
    scream::tridiag::test::correct::run_property_test<1, EKAT_PACK_SIZE>();
    scream::tridiag::test::correct::run_property_test<EKAT_PACK_SIZE, EKAT_PACK_SIZE>();
  }
}

TEST_CASE("bfb", "tridiag") {
  scream::tridiag::test::correct::run_bfb_test<1,1>();
  if (EKAT_PACK_SIZE > 1)
    scream::tridiag::test::correct::run_bfb_test<EKAT_PACK_SIZE, EKAT_PACK_SIZE>();
}
