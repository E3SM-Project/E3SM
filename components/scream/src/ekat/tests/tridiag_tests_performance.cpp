#include <chrono>

#include "tridiag_tests.hpp"

namespace scream {
namespace tridiag {
namespace test {
namespace perf {

void expect_another_arg (int i, int argc) {
  if (i == argc-1)
    throw std::runtime_error("Expected another cmd-line arg.");
}

std::string Solver::convert (Enum e) {
  switch (e) {
  case thomas: return "thomas";
  case cr: return "cr";
  default: scream_require_msg(false, "Not a valid solver: " << e);
  }
}

Solver::Enum Solver::convert (const std::string& s) {
  if (s == "thomas") return thomas;
  if (s == "cr") return cr;
  return error;
}

Input::Input ()
  : method(Solver::cr), nprob(2048), nrow(128), nrhs(43), nwarp(-1),
    pack( ! scream::util::OnGpu<Kokkos::DefaultExecutionSpace>::value),
    oneA(false)
{}

bool Input::parse (int argc, char** argv) {
  using scream::util::eq;
  for (int i = 1; i < argc; ++i) {
    if (eq(argv[i], "-m", "--method")) {
      expect_another_arg(i, argc);
      method = Solver::convert(argv[++i]);
      if (method == Solver::error) {
        std::cout << "Not a solver: " << argv[i] << "\n";
        return false;
      }
    } else if (eq(argv[i], "-np", "--nprob")) {
      expect_another_arg(i, argc);
      nprob = std::atoi(argv[++i]);
    } else if (eq(argv[i], "-nr", "--nrow")) {
      expect_another_arg(i, argc);
      nrow = std::atoi(argv[++i]);
    } else if (eq(argv[i], "-nc", "--nrhs")) {
      expect_another_arg(i, argc);
      nrhs = std::atoi(argv[++i]);
    } else if (eq(argv[i], "-1a", "--oneA")) {
      oneA = true;
    } else if (eq(argv[i], "-nw", "--nwarp")) {
      expect_another_arg(i, argc);
      nwarp = std::atoi(argv[++i]);
    } else if (eq(argv[i], "-nop", "--nopack")) {
      pack = false;
    } else {
      std::cout << "Unexpected arg: " << argv[i] << "\n";
      return false;
    }
  }
  if (nrhs == 1) oneA = true;
  if (method == Solver::cr) pack = false;
  return true;
}

std::string string (const Input& in, const int& nwarp) {
  std::stringstream ss;
  ss << "run: solver " << Solver::convert(in.method)
     << " pack " << in.pack
     << " nprob " << in.nprob
     << " nrow " << in.nrow
     << " nA " << (in.oneA ? 1 : in.nrhs)
     << " nrhs " << in.nrhs
     << " nwarp " << nwarp << "\n";
  return ss.str();
}

using BulkLayout = Kokkos::LayoutRight;
using TeamLayout = Kokkos::LayoutRight;

template <typename TridiagArray>
KOKKOS_INLINE_FUNCTION
Kokkos::View<typename TridiagArray::value_type*, TeamLayout, Kokkos::MemoryUnmanaged>
get_diag (const TridiagArray& A, const int& ip, const int& diag_idx) {
  assert(A.extent_int(3) == 1);
  return Kokkos::View<typename TridiagArray::value_type*, TeamLayout, Kokkos::MemoryUnmanaged>(
    &A.impl_map().reference(ip, diag_idx, 0, 0),
    A.extent_int(2));
}

template <typename TridiagArray>
KOKKOS_INLINE_FUNCTION
Kokkos::View<typename TridiagArray::value_type**, TeamLayout, Kokkos::MemoryUnmanaged>
get_diags (const TridiagArray& A, const int& ip, const int& diag_idx) {
  return Kokkos::View<typename TridiagArray::value_type**, TeamLayout, Kokkos::MemoryUnmanaged>(
    &A.impl_map().reference(ip, diag_idx, 0, 0),
    A.extent_int(2), A.extent_int(3));
}

template <typename DataArray>
KOKKOS_INLINE_FUNCTION
Kokkos::View<typename DataArray::value_type*, TeamLayout, Kokkos::MemoryUnmanaged>
get_x (const DataArray& X, const int& ip) {
  assert(X.extent_int(2) == 1);
  return Kokkos::View<typename DataArray::value_type*, TeamLayout, Kokkos::MemoryUnmanaged>(
    &X.impl_map().reference(ip, 0, 0), X.extent_int(1));
}

template <typename DataArray>
KOKKOS_INLINE_FUNCTION
Kokkos::View<typename DataArray::value_type**, TeamLayout, Kokkos::MemoryUnmanaged>
get_xs (const DataArray& X, const int& ip) {
  return Kokkos::View<typename DataArray::value_type**, TeamLayout, Kokkos::MemoryUnmanaged>(
    &X.impl_map().reference(ip, 0, 0), X.extent_int(1), X.extent_int(2));
}

template <typename Scalar>
using TridiagArrays = Kokkos::View<Scalar****, BulkLayout>;
template <typename Scalar>
using DataArrays = Kokkos::View<Scalar***, BulkLayout>;

template <typename Real>
void run (const Input& in) {
  using Kokkos::create_mirror_view;
  using Kokkos::deep_copy;
  using Kokkos::subview;
  using Kokkos::ALL;
  using scream::pack::scalarize;
  using scream::pack::npack;
  using TeamPolicy = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>;
  using MT = typename TeamPolicy::member_type;
  using APack = scream::pack::Pack<Real, SCREAM_SMALL_PACK_SIZE>;
  using DataPack = scream::pack::Pack<Real, SCREAM_SMALL_PACK_SIZE>;

  const auto gettime = [&] () {
    return std::chrono::steady_clock::now();
  };
  using TimePoint = decltype(gettime());
  const auto duration = [&] (const TimePoint& t0, const TimePoint& tf) -> double {
    return 1e-6*std::chrono::duration_cast<std::chrono::microseconds>(tf - t0).count();
  };

  const bool on_gpu = scream::util::OnGpu<Kokkos::DefaultExecutionSpace>::value;
  const int nA = in.oneA ? 1 : in.nrhs;

  scream_require_msg( ! in.pack || in.method != Solver::cr, "CR has no pack version.");

  TridiagArrays<Real> A, Acopy;
  DataArrays<Real> B, X, Y;
  TridiagArrays<APack> Ap, Apcopy;
  DataArrays<DataPack> Bp, Xp, Yp;
  if (in.pack) {
    if (in.oneA) {
      A = TridiagArrays<Real>("A", in.nprob, 3, in.nrow, nA);
      Acopy = TridiagArrays<Real>("Acopy", in.nprob, 3, in.nrow, nA);
    } else {
      Ap = TridiagArrays<APack>("A", in.nprob, 3, in.nrow, npack<APack>(nA));
      Apcopy = TridiagArrays<APack>("Acopy", in.nprob, 3, in.nrow, npack<APack>(nA));
      A = scalarize(Ap);
      Acopy = scalarize(Apcopy);
    }
    const int nrhs = npack<DataPack>(in.nrhs);
    Bp = DataArrays<DataPack>("B", in.nprob, in.nrow, nrhs);
    Xp = DataArrays<DataPack>("X", in.nprob, in.nrow, nrhs);
    Yp = DataArrays<DataPack>("Y", in.nprob, in.nrow, nrhs);
    B = scalarize(Bp);
    X = scalarize(Xp);
    Y = scalarize(Yp);
  } else {
    A = TridiagArrays<Real>("A", in.nprob, 3, in.nrow, nA);
    Acopy = TridiagArrays<Real>("Acopy", in.nprob, 3, in.nrow, nA);
    B = DataArrays<Real>("B", in.nprob, in.nrow, in.nrhs);
    X = DataArrays<Real>("X", in.nprob, in.nrow, in.nrhs);
    Y = DataArrays<Real>("Y", in.nprob, in.nrow, in.nrhs);
  }

  auto Am = create_mirror_view(A);
  auto Bm = create_mirror_view(B);
  const auto fill = [&] (const int i) {
    const auto dl = subview(Am, i, 0, ALL(), ALL());
    const auto d  = subview(Am, i, 1, ALL(), ALL());
    const auto du = subview(Am, i, 2, ALL(), ALL());
    fill_tridiag_matrix(dl, d, du, nA, i);
    fill_data_matrix(subview(Bm, i, ALL(), ALL()), in.nrhs);
  };
  Kokkos::parallel_for(
    Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, in.nprob), fill);
  deep_copy(A, Am);
  deep_copy(B, Bm);
  deep_copy(Acopy, A);
  deep_copy(X, B);

  TeamPolicy policy(in.nprob,
                    on_gpu ? (in.nwarp < 0 ? 128 : 32*in.nwarp) : 1,
                    1);
  assert(in.nwarp < 0 || ! on_gpu || policy.team_size() == 32*in.nwarp);
  std::cout << string(in, policy.team_size()/32);

  Kokkos::fence();
  TimePoint t0, t1;
  switch (in.method) {
  case Solver::thomas: {
    if (on_gpu) {
      scream_require_msg(
        in.oneA, "On GPU, only 1 A/team is supported in the Thomas algorithm.");
      t0 = gettime();
      const auto f = KOKKOS_LAMBDA (const MT& team) {
        const int ip = team.league_rank();
        const auto dl = get_diag(A, ip, 0);
        const auto d  = get_diag(A, ip, 1);
        const auto du = get_diag(A, ip, 2);
        const auto x = get_xs(X, ip);
        scream::tridiag::thomas(team, dl, d, du, x);
      };
      Kokkos::parallel_for(policy, f);
      Kokkos::fence();
      t1 = gettime();
    } else {
      if (in.pack) {
        for (int trial = 0; trial < 2; ++trial) {
          // Manual deep copy b/c Kokkos::deep_copy has some side effect that
          // affects performance on SKX and especially KNL.
          const auto Adc = KOKKOS_LAMBDA (const MT& team) {
            const auto s = [&] () {
              const int i = team.league_rank();
              for (int k = 0; k < 3; ++k)
                for (int r = 0; r < in.nrow; ++r)
                  for (int c = 0; c < nA; ++c)
                    A(i,k,r,c) = Acopy(i,k,r,c);
            };
            Kokkos::single(Kokkos::PerTeam(team), s);
          };
          Kokkos::parallel_for(policy, Adc);
          const int nrhs = npack<DataPack>(in.nrhs);
          const auto Xdc = KOKKOS_LAMBDA (const MT& team) {
            const auto s = [&] () {
              const int i = team.league_rank();
              for (int r = 0; r < in.nrow; ++r)
                for (int c = 0; c < nrhs; ++c)
                  Xp(i,r,c) = Bp(i,r,c);
            };
            Kokkos::single(Kokkos::PerTeam(team), s);
          };
          Kokkos::parallel_for(policy, Xdc);
          Kokkos::fence();
          t0 = gettime();
          if (in.nrhs == 1) {
            assert(in.oneA);
            const auto f = KOKKOS_LAMBDA (const MT& team) {
              const auto single = [&] () {
                const int ip = team.league_rank();
                const auto dl = get_diag(A, ip, 0);
                const auto d  = get_diag(A, ip, 1);
                const auto du = get_diag(A, ip, 2);
                const auto x  = get_x(Xp, ip);
                scream::tridiag::thomas(dl, d, du, x);
              };
              Kokkos::single(Kokkos::PerTeam(team), single);
            };
            Kokkos::parallel_for(policy, f);
          } else {
            if (in.oneA) {
              const auto f = KOKKOS_LAMBDA (const MT& team) {
                const auto single = [&] () {
                  const int ip = team.league_rank();
                  const auto dl = get_diag(A, ip, 0);
                  const auto d  = get_diag(A, ip, 1);
                  const auto du = get_diag(A, ip, 2);
                  const auto x  = get_xs(Xp, ip);
                  scream::tridiag::thomas(dl, d, du, x);
                };
                Kokkos::single(Kokkos::PerTeam(team), single);
              };
              Kokkos::parallel_for(policy, f);
            } else {
              const auto f = KOKKOS_LAMBDA (const MT& team) {
                const auto single = [&] () {
                  const int ip = team.league_rank();
                  const auto dl = get_diags(Ap, ip, 0);
                  const auto d  = get_diags(Ap, ip, 1);
                  const auto du = get_diags(Ap, ip, 2);
                  const auto x  = get_xs(Xp, ip);
                  assert(x.extent_int(1) == nrhs);
                  assert(d.extent_int(1) == nrhs);
                  scream::tridiag::thomas(dl, d, du, x);
                };
                Kokkos::single(Kokkos::PerTeam(team), single);
              };
              Kokkos::parallel_for(policy, f);
            }
          }
          Kokkos::fence();
          t1 = gettime();
        }
      } else {
        for (int trial = 0; trial < 2; ++trial) {
          const auto Adc = KOKKOS_LAMBDA (const MT& team) {
            const auto s = [&] () {
              const int i = team.league_rank();
              for (int k = 0; k < 3; ++k)
                for (int r = 0; r < in.nrow; ++r)
                  for (int c = 0; c < nA; ++c)
                    A(i,k,r,c) = Acopy(i,k,r,c);
            };
            Kokkos::single(Kokkos::PerTeam(team), s);
          };
          Kokkos::parallel_for(policy, Adc);
          const auto Xdc = KOKKOS_LAMBDA (const MT& team) {
            const auto s = [&] () {
              const int i = team.league_rank();
              for (int r = 0; r < in.nrow; ++r)
                for (int c = 0; c < in.nrhs; ++c)
                  X(i,r,c) = B(i,r,c);
            };
            Kokkos::single(Kokkos::PerTeam(team), s);
          };
          Kokkos::parallel_for(policy, Xdc);
          Kokkos::fence();
          t0 = gettime();
          if (in.nrhs == 1) {
            assert(in.oneA);
            const auto f = KOKKOS_LAMBDA (const MT& team) {
              const auto single = [&] () {
                const int ip = team.league_rank();
                const auto dl = get_diag(A, ip, 0);
                const auto d  = get_diag(A, ip, 1);
                const auto du = get_diag(A, ip, 2);
                const auto x  = get_x(X, ip);
                scream::tridiag::thomas(dl, d, du, x);
              };
              Kokkos::single(Kokkos::PerTeam(team), single);
            };
            Kokkos::parallel_for(policy, f);
          } else {
            if (in.oneA) {
              const auto f = KOKKOS_LAMBDA (const MT& team) {
                const auto single = [&] () {
                  const int ip = team.league_rank();
                  const auto dl = get_diag(A, ip, 0);
                  const auto d  = get_diag(A, ip, 1);
                  const auto du = get_diag(A, ip, 2);
                  const auto x  = get_xs(X, ip);
                  scream::tridiag::thomas(dl, d, du, x);
                };
                Kokkos::single(Kokkos::PerTeam(team), single);
              };
              Kokkos::parallel_for(policy, f);
            } else {
              const auto f = KOKKOS_LAMBDA (const MT& team) {
                const auto single = [&] () {
                  const int ip = team.league_rank();
                  const auto dl = get_diags(A, ip, 0);
                  const auto d  = get_diags(A, ip, 1);
                  const auto du = get_diags(A, ip, 2);
                  const auto x  = get_xs(X, ip);
                  assert(x.extent_int(1) == in.nrhs);
                  assert(d.extent_int(1) == in.nrhs);
                  scream::tridiag::thomas(dl, d, du, x);
                };
                Kokkos::single(Kokkos::PerTeam(team), single);
              };
              Kokkos::parallel_for(policy, f);
            }
          }
          Kokkos::fence();
          t1 = gettime();
        }
      }
    }
  } break;
  case Solver::cr: {
    assert( ! in.pack);
    t0 = gettime();
    if (in.nrhs == 1) {
      assert(in.oneA);
      const auto f = KOKKOS_LAMBDA (const MT& team) {
        const int ip = team.league_rank();
        const auto dl = get_diag(A, ip, 0);
        const auto d  = get_diag(A, ip, 1);
        const auto du = get_diag(A, ip, 2);
        const auto x  = get_x(X, ip);
        scream::tridiag::cr(team, dl, d, du, x);
      };
      Kokkos::parallel_for(policy, f);
    } else {
      if (in.oneA) {
        const auto f = KOKKOS_LAMBDA (const MT& team) {
          const int ip = team.league_rank();
          const auto dl = get_diag(A, ip, 0);
          const auto d  = get_diag(A, ip, 1);
          const auto du = get_diag(A, ip, 2);
          const auto x  = get_xs(X, ip);
          scream::tridiag::cr(team, dl, d, du, x);
        };
        Kokkos::parallel_for(policy, f);
      } else {
        const auto f = KOKKOS_LAMBDA (const MT& team) {
          const int ip = team.league_rank();
          const auto dl = get_diags(A, ip, 0);
          const auto d  = get_diags(A, ip, 1);
          const auto du = get_diags(A, ip, 2);
          const auto x  = get_xs(X, ip);
          assert(x.extent_int(1) == in.nrhs);
          assert(d.extent_int(1) == in.nrhs);
          scream::tridiag::cr(team, dl, d, du, x);
        };
        Kokkos::parallel_for(policy, f);
      }
    }
    Kokkos::fence();
    t1 = gettime();    
  } break;
  default:
    std::cout << "run does not support "
              << Solver::convert(in.method) << "\n";
  }

  const auto et = duration(t0, t1);
  printf("run: et %1.3e et/datum %1.3e\n", et, et/(in.nprob*in.nrow*in.nrhs));

  Real re; {
    auto Acopym = create_mirror_view(Acopy);
    auto Xm = create_mirror_view(X);
    auto Ym = create_mirror_view(Y);
    deep_copy(Acopym, Acopy);
    deep_copy(Xm, X);
    const auto ip = std::max(0, in.nprob-1);
    const auto dl = subview(Acopym, ip, 0, ALL(), ALL());
    const auto d  = subview(Acopym, ip, 1, ALL(), ALL());
    const auto du = subview(Acopym, ip, 2, ALL(), ALL());
    matvec(dl, d, du,
           subview(Xm, in.nprob-1, ALL(), ALL()),
           subview(Ym, in.nprob-1, ALL(), ALL()),
           nA, in.nrhs);
    re = reldif(subview(Bm, in.nprob-1, ALL(), ALL()),
                subview(Ym, in.nprob-1, ALL(), ALL()),
                in.nrhs);
  }
  if (re > 50*std::numeric_limits<Real>::epsilon())
    std::cout << "run: " << " re " << re << "\n";
}

template void run<scream::Real>(const Input&);

} // namespace perf
} // namespace test
} // namespace tridiag
} // namespace scream
