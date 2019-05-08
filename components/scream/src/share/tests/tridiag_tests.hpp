#ifndef INCLUDE_SCREAM_TRIDIAG_TESTS
#define INCLUDE_SCREAM_TRIDIAG_TESTS

#include "share/util/scream_tridiag.hpp"
#include "share/util/scream_utils.hpp"
#include "share/util/scream_arch.hpp"
#include "share/scream_pack.hpp"
#include "share/scream_pack_kokkos.hpp"

namespace scream {
namespace tridiag {
namespace test {

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

namespace perf {
struct Solver {
  enum Enum { thomas, cr, error };

  static std::string convert(Enum e);
  static Enum convert(const std::string& s);
};

struct Input {
  Solver::Enum method;
  int nprob, nrow, nrhs, nwarp;
  bool pack, oneA;

  Input();
  bool parse(int argc, char** argv);
};

template <typename Real>
void run(const Input& in);
}

} // namespace test
} // namespace tridiag
} // namespace scream

#endif // INCLUDE_SCREAM_TRIDIAG_TESTS
