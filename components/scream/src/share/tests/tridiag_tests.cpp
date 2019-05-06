#define CATCH_CONFIG_RUNNER
#include "catch2/catch.hpp"

#include "share/util/scream_tridiag.hpp"

#include "share/scream_session.hpp"
#include "share/scream_pack.hpp"
#include "share/scream_pack_kokkos.hpp"
#include "share/util/scream_arch.hpp"

template <typename TridiagDiag>
KOKKOS_INLINE_FUNCTION
void fill_tridiag_matrix (TridiagDiag dl, TridiagDiag d, TridiagDiag du,
                          const int& nprob, const int& seed) {
  const int nrow = d.extent_int(0);

  assert(dl.extent_int(0) == nrow);
  assert(du.extent_int(0) == nrow);

  for (int p = 0; p < nprob; ++p)
    for (int i = 0; i < nrow; ++i) {
      const int k = seed + p + i;
      dl(i,p) = (k % 5 == 0 ? -1 : 1) * 1.3 * (0.1 + ((k*k) % 11));
      du(i,p) = (k % 7 == 0 ? -1 : 1) * 1.7 * (0.2 + ((k*k) % 13));
      d (i,p) = ((k % 3 == 0 ? -1 : 1) *
                 (0.7 + std::abs(dl(i,p)) + std::abs(du(i,p)) + (k % 17)));
    }
}

template <typename DataArray>
KOKKOS_INLINE_FUNCTION
void fill_data_matrix (DataArray X, const int& seed) {
  const int nrow = X.extent_int(0);
  const int nrhs = X.extent_int(1);

  for (int i = 0; i < nrow; ++i)
    for (int j = 0; j < nrhs; ++j)
      X(i,j) = (((7*i + 11*j + 3*i*j) % 3 == 0 ? -1 : 1) *
                1.7 * ((17*(i - 19) + 13*(j - 11) + 5*(i - 5)*(j - 7) + seed) % 47));
}

template <typename TridiagDiag, typename XArray, typename YArray>
KOKKOS_INLINE_FUNCTION
int matvec (TridiagDiag dl, TridiagDiag d, TridiagDiag du, XArray X, YArray Y,
            const int nprob, const int nrhs) {
  const int nrow = d.extent_int(0);

  assert(dl.extent_int(0) == nrow);
  assert(du.extent_int(0) == nrow);
  assert(X .extent_int(0) == nrow);
  assert(X .extent_int(1) >= nrhs);
  assert(Y .extent_int(0) == nrow);
  assert(Y .extent_int(1) == X.extent_int(1));
  assert(dl.extent_int(1) >= nprob);
  assert(nprob == 1 || nprob == nrhs);
  assert(d .extent_int(1) == dl.extent_int(1));
  assert(du.extent_int(1) == dl.extent_int(1));

  const auto dcol = [&] (const int& j) -> int { return nprob > 1 ? j : 0; };

  if (nrow == 1) {
    for (int j = 0; j < nrhs; ++j) {
      const int aj = dcol(j);
      Y(0,j) = d(0,aj) * X(0,j);
    }
    return 0;
  }

  for (int j = 0; j < nrhs; ++j) {
    const int aj = dcol(j);
    Y(0,j) = d(0,aj) * X(0,j) + du(0,aj) * X(1,j);
  }
  for (int i = 1; i < nrow-1; ++i)
    for (int j = 0; j < nrhs; ++j) {
      const int aj = dcol(j);
      Y(i,j) = (dl(i,aj) * X(i-1,j) +
                d (i,aj) * X(i  ,j) +
                du(i,aj) * X(i+1,j));
    }
  const int i = nrow-1;
  for (int j = 0; j < nrhs; ++j) {
    const int aj = dcol(j);
    Y(i,j) = dl(i,aj) * X(i-1,j) + d(i,aj) * X(i,j);
  }

  return 0;  
}

template <typename Array>
scream::Real reldif (const Array& a, const Array& b, const int nrhs) {
  assert(a.extent_int(0) == b.extent_int(0));
  assert(a.extent_int(1) == b.extent_int(1));
  assert(a.rank == 2);
  assert(b.rank == 2);
  scream::Real num = 0, den = 0;
  for (int i = 0; i < a.extent_int(0); ++i)
    for (int j = 0; j < nrhs; ++j) {
      if (std::isnan(a(i,j)) || std::isnan(b(i,j)) ||
          std::isinf(a(i,j)) || std::isinf(b(i,j))) {
        return std::numeric_limits<scream::Real>::infinity();
      }
      num = std::max(num, std::abs(a(i,j) - b(i,j)));
      den = std::max(den, std::abs(a(i,j)));
    }
  return num/den;
}

struct Solver {
  enum Enum { thomas_team_scalar, thomas_team_pack,
              thomas_scalar, thomas_pack,
              cr_scalar,
              error };
  static std::string convert (Enum e) {
    switch (e) {
    case thomas_team_scalar: return "thomas_team_scalar";
    case thomas_team_pack: return "thomas_team_pack";
    case thomas_scalar: return "thomas_scalar";
    case thomas_pack: return "thomas_pack";
    case cr_scalar: return "cr_scalar";
    default: scream_require_msg(false, "Not a valid solver: " << e);
    }
  }
  static Enum convert (const std::string& s) {
    if (s == "thomas_team_scalar") return thomas_team_scalar;
    if (s == "thomas_team_pack") return thomas_team_pack;
    if (s == "thomas_scalar") return thomas_scalar;
    if (s == "thomas_pack") return thomas_pack;
    if (s == "cr_scalar") return cr_scalar;
    return error;
  }

  static Enum all[];
};

Solver::Enum Solver::all[] = { thomas_team_scalar, thomas_team_pack,
                               thomas_scalar, thomas_pack,
                               cr_scalar };

namespace test_correct {
struct TestConfig {
  using BulkLayout = Kokkos::LayoutRight;
  using TeamLayout = Kokkos::LayoutRight;

  Solver::Enum solver;
  int n_kokkos_thread, n_kokkos_vec;
};

template <typename Scalar>
using TridiagArray = Kokkos::View<Scalar***, TestConfig::TeamLayout>;
template <typename Scalar>
using DataArray = Kokkos::View<Scalar**, TestConfig::TeamLayout>;

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
    default:
      scream_require_msg(false, "Different pack size: " << Solver::convert(tc.solver));
    }
  }
};

template <int A_pack_size, int data_pack_size>
void run_test (const TestConfig& tc) {
  using Kokkos::create_mirror_view;
  using Kokkos::deep_copy;
  using Kokkos::subview;
  using Kokkos::ALL;
  using scream::Real;
  using scream::pack::npack;
  using scream::pack::scalarize;
  using APack = scream::pack::Pack<Real, A_pack_size>;
  using DataPack = scream::pack::Pack<Real, data_pack_size>;

#if 1
  const int nrows[] = {1,2,3,4, 8,10,16, 32,43, 63,64,65, 111,128,129, 8192};
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
             tc.solver == Solver::cr_scalar) &&
            data_pack_size > 1)
          continue;
        if (static_cast<int>(APack::n) != static_cast<int>(DataPack::n) && nprob > 1)
          continue;

        const int prob_npack = npack<APack>(nprob);
        const int rhs_npack = npack<DataPack>(nrhs);

        TridiagArray<APack>
          A("A", 3, nrow, prob_npack),
          Acopy("A", A.extent(0), A.extent(1), A.extent(2));
        DataArray<DataPack>
          B("B", nrow, rhs_npack), X("X", B.extent(0), B.extent(1)),
          Y("Y", X.extent(0), X.extent(1));
        const auto Am = create_mirror_view(A);
        const auto Bm = create_mirror_view(B);
        {
          const auto As = scalarize(Am);
          const auto dl = subview(As, 0, ALL(), ALL());
          const auto d  = subview(As, 1, ALL(), ALL());
          const auto du = subview(As, 2, ALL(), ALL());
          fill_tridiag_matrix(dl, d, du, nprob, nrhs /* seed */);
        }
        fill_data_matrix(scalarize(Bm), nrhs);
        deep_copy(A, Am);
        deep_copy(B, Bm);
        deep_copy(Acopy, A);
        deep_copy(X, B);

        Solve<A_pack_size == data_pack_size, APack, DataPack>
          ::run(tc, A, X, nprob, nrhs);

        Real re; {
          const auto Acopym = create_mirror_view(Acopy);
          const auto As = scalarize(Acopym);
          const auto dl = subview(As, 0, ALL(), ALL());
          const auto d  = subview(As, 1, ALL(), ALL());
          const auto du = subview(As, 2, ALL(), ALL());
          const auto Xm = create_mirror_view(X);
          const auto Ym = create_mirror_view(Y);
          deep_copy(Acopym, Acopy);
          deep_copy(Xm, X);
          matvec(dl, d, du, scalarize(Xm), scalarize(Ym), nprob, nrhs);
          re = reldif(scalarize(Bm), scalarize(Ym), nrhs);
        }
        const bool pass = re <= 50*std::numeric_limits<Real>::epsilon();
        std::stringstream ss;
        ss << Solver::convert(tc.solver) << " " << tc.n_kokkos_thread
           << " " << tc.n_kokkos_vec << " | " << nrow << " " << nrhs << " "
           << A_many << " | log10 reldif " << std::log10(re);
        if ( ! pass) std::cout << "FAIL: " << ss.str() << "\n";
        REQUIRE(pass);
        //std::cout << "PASS: " << ss.str() << "\n";
      }
    }
  }
}

template <int A_pack_size, int data_pack_size>
void run_test () {
  for (const auto solver : Solver::all) {
    TestConfig tc;
    tc.solver = solver;
    if (scream::util::OnGpu<Kokkos::DefaultExecutionSpace>::value) {
      tc.n_kokkos_vec = 1;
      for (const int n_kokkos_thread : {128, 256, 512}) {
        tc.n_kokkos_thread = n_kokkos_thread;
        run_test<A_pack_size, data_pack_size>(tc);
      }
      tc.n_kokkos_vec = 32;
      for (const int n_kokkos_thread : {4}) {
        tc.n_kokkos_thread = n_kokkos_thread;
        run_test<A_pack_size, data_pack_size>(tc);
      }
    } else {
      const int concurrency = Kokkos::DefaultExecutionSpace::concurrency();
      const int n_kokkos_thread = concurrency;
      for (const int n_kokkos_vec : {1, 2}) {
        tc.n_kokkos_thread = n_kokkos_thread;
        tc.n_kokkos_vec = n_kokkos_vec;
        run_test<A_pack_size, data_pack_size>(tc);
      }
    }
  }
}
} // namespace test_correct

TEST_CASE("tridiag", "correctness") {
  test_correct::run_test<1,1>();
  if (SCREAM_PACK_SIZE > 1) {
    test_correct::run_test<1, SCREAM_PACK_SIZE>();
    test_correct::run_test<SCREAM_PACK_SIZE, SCREAM_PACK_SIZE>();
  }
}

int main (int argc, char **argv) {
  int num_failed = 0;
  scream::initialize_scream_session(argc, argv); {
    if (argc > 1) {
      // Performance test.
      // TODO
    } else {
      // Correctness tests.
      num_failed = Catch::Session().run(argc, argv);
    }
  } scream::finalize_scream_session();
  return num_failed != 0 ? 1 : 0;
}
