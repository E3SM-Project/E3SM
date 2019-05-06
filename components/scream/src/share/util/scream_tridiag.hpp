#ifndef INCLUDE_SCREAM_TRIDIAG
#define INCLUDE_SCREAM_TRIDIAG

#include <cassert>

#include "share/scream_kokkos.hpp"
#include "share/util/scream_utils.hpp"

namespace scream {
namespace tridiag {

/* Team-level solvers for diagonally dominant, scalar tridiagonal systems.

   This file is a header-only library to solve the equation
       A x = b,
   where A is a scalar, tridiagonal, diagonally dominant matrix, on GPU and
   non-GPU architectures, within a Kokkos team. The library supports three
   problem formats:
       1. A x = b: 1 matrix A, one L,RHS x, b;
       2. A X = B: 1 matrix A, multiple L,RHS X, B;
       3. A_i x_i = b_i, i = 1..n: Multiple matrices A, each associated with 1
          L,RHS x, b.

   The nxn matrix A is formatted as (dl, d, du). Using 0-based indexing,
       * the lower diagonal of A is in dl(1:n-1);
       * the diagonal is in d(0:n-1);
       * the upper diagonal is in du(0:n-2).
   In problem format 3, matrix i is stored in (dl(:,i), d(:,i), du(:,i)).

   In problem formats 2, 3, the i'th L,RHS in X, B is stored as X(:,i),
   B(:,i). Otherwise, X(:) = x, B(:) = b.

   In all cases, arrays should have layout LayoutRight. LayoutStride is not
   supported because it degrades performance. This library is intended to help
   you get the best performance, so it disallows nonperformant types.

   This library is intended to be efficient when there are many problems of the
   forms (1,2,3) to solve simultaneously, with one per Kokkos team. For example,
   in a physics parameterization, a device may have 2000 physics columns, and
   each has a problem to solve.

   There are three interface-level functions, and each supports the three
   problem formats. In all functions,
       * X = B on input and X = A \ B on output;
       * (dl, d, du) are overwritten;
       * the value type of each of (dl, d, du) must be the same;
       * the value type of X can differ from that of (dl, d, du).
   The functions do not use any temporary workspace other than registers. The
   interface functions are as follows:

   a. Use the Thomas algorithm to solve a problem within a Kokkos team:

        template <typename TeamMember, typename TridiagDiag, typename DataArray>
        void thomas(const TeamMember& team,
                    TridiagDiag dl, TridiagDiag d, TridiagDiag du, DataArray X);

   b. Use the Thomas algorithm to solve the problem in serial. Here no reference
      is made to Kokkos parallel constructs. The call must be protected by
      Kokkos::single(Kokkos::PerTeam).

        template <typename TridiagDiag, typename DataArray>
        void thomas(TridiagDiag dl, TridiagDiag d, TridiagDiag du, DataArray X);
   
   c. Cyclic reduction at the Kokkos team level. Any Kokkos (thread, vector)
      parameterization works.

        template <typename TeamMember, typename TridiagDiag, typename DataArray>
        void cr(const TeamMember& team,
                TridiagDiag dl, TridiagDiag d, TridiagDiag du, DataArray X);

   In practice, (a, b) are used on a non-GPU computer, and (c) is used on the
   GPU. On a non-GPU computer, the typical use case is that a team has just one
   thread. On a GPU, the typical use case is that a team has 128 to 1024 threads
   (4 to 32 warps).

   On a non-GPU computer, in the case of multiple A or L,RHS per team,
   scream::pack::Pack may be used as the value type. On GPU, as usual, only
   scream::pack::Pack<scalar_type, 1> makes sense, so it also likely makes sense
   that the value type is just the POD (plain-old data) scalar_type.

   The rest of this file contains implementation details. Each of (a, b, c) is
   specialized to the various problem formats. This header documentation is the
   interface, and nothing further needs to be read.
*/

namespace impl {

template <typename Array>
using EnableIfCanUsePointer =
  typename std::enable_if<std::is_same<typename Array::array_layout,
                                       Kokkos::LayoutRight>::value ||
                          (std::is_same<typename Array::array_layout,
                                        Kokkos::LayoutLeft>::value
                           && Array::rank == 1)>::type;

template <typename TeamMember>
KOKKOS_INLINE_FUNCTION
int get_thread_id_within_team (const TeamMember& team) {
  return team.team_rank();
}

template <typename TeamMember>
KOKKOS_INLINE_FUNCTION
int get_team_nthr (const TeamMember& team) {
  return team.team_size();
}

#ifdef KOKKOS_ENABLE_CUDA
KOKKOS_INLINE_FUNCTION
int get_thread_id_within_team (const Kokkos::Impl::CudaTeamMember& team) {
#ifdef __CUDA_ARCH__
  // Can't use team.team_rank() here because vector direction also uses physical
  // threads but TeamMember types don't expose that information.
  return blockDim.x * threadIdx.y + threadIdx.x;
#else
  assert(0);
  return -1;
#endif
}

KOKKOS_INLINE_FUNCTION
int get_team_nthr (const Kokkos::Impl::CudaTeamMember& team) {
#ifdef __CUDA_ARCH__
  return blockDim.x * blockDim.y;
#else
  assert(0);
  return -1;
#endif
}
#endif

// The caller must provide the team_barrier after this function returns before A
// is accessed.
template <typename TeamMember, typename TridiagDiag>
KOKKOS_INLINE_FUNCTION
void thomas_factorize (const TeamMember& team,
                       TridiagDiag dl, TridiagDiag d, TridiagDiag du) {
  const int nrow = d.extent_int(0);
  assert(dl.extent_int(0) == nrow);
  assert(du.extent_int(0) == nrow);
  const auto f = [&] () {
    for (int i = 1; i < nrow; ++i) {
      dl(i) /= d(i-1);
      d (i) -= dl(i) * du(i-1);
    }
  };
  Kokkos::single(Kokkos::PerTeam(team), f);
}

template <typename TeamMember, typename TridiagDiag, typename DataArray>
KOKKOS_INLINE_FUNCTION
void thomas_solve (const TeamMember& team,
                   TridiagDiag dl, TridiagDiag d, TridiagDiag du,
                   DataArray X) {
  const int nrow = d.extent_int(0);
  const int nrhs = X.extent_int(1);
  assert(X.extent_int(0) == nrow);
  assert(dl.extent_int(0) == nrow);
  assert(du.extent_int(0) == nrow);
  const auto f = [&] (const int& j) {
    const auto g = [&] () {
      for (int i = 1; i < nrow; ++i)
        X(i,j) -= dl(i) * X(i-1,j);
      X(nrow-1,j) /= d(nrow-1);
      for (int i = nrow-1; i > 0; --i)
        X(i-1,j) = (X(i-1,j) - du(i-1) * X(i,j)) / d(i-1);
    };
    Kokkos::single(Kokkos::PerThread(team), g);
  };
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nrhs), f);
}

template <typename DT, typename XT>
KOKKOS_INLINE_FUNCTION
void thomas_a1x1 (DT* const dl, DT* d, DT* const du, XT* X, const int nrow) {
  for (int i = 1; i < nrow; ++i) {
    const auto dli = dl[i] / d[i-1];
    d[i] -= dli * du[i-1];
    X[i] -= dli * X[i-1];
  }
  X[nrow-1] /= d[nrow-1];
  for (int i = nrow-1; i > 0; --i)
    X[i-1] = (X[i-1] - du[i-1] * X[i]) / d[i-1];
}

template <typename DT, typename XT>
KOKKOS_INLINE_FUNCTION
void thomas_a1xm (DT* const dl, DT* d, DT* const du, XT* X,
                  const int nrow, const int nrhs) {
  for (int i = 1; i < nrow; ++i) {
    const auto dli = dl[i] / d[i-1];
    d[i] -= dli * du[i-1];
    auto* const xim1 = X + (i-1)*nrhs;
    auto* const xi = X + i*nrhs;
    for (int j = 0; j < nrhs; ++j)
      xi[j] -= dli * xim1[j];
  }
  {
    auto* const xi = X + (nrow-1)*nrhs;
    for (int j = 0; j < nrhs; ++j)
      xi[j] /= d[nrow-1];
  }
  for (int i = nrow-1; i > 0; --i) {
    auto* const xim1 = X + (i-1)*nrhs;
    auto* const xi = X + i*nrhs;
    for (int j = 0; j < nrhs; ++j)
      xim1[j] = (xim1[j] - du[i-1] * xi[j]) / d[i-1];
  }
}

template <typename DT, typename XT>
KOKKOS_INLINE_FUNCTION
void thomas_amxm (DT* const dl, DT* d, DT* const du, XT* X,
                  const int nrow, const int nrhs) {
  for (int i = 1; i < nrow; ++i) {
    const int ios = i*nrhs;
    const int im1os = (i-1)*nrhs;
    auto* const dli = dl + ios;
    auto* const di = d + ios;
    auto* const dim1 = d + im1os;
    auto* const duim1 = du + im1os;
    auto* const xim1 = X + im1os;
    auto* const xi = X + ios;
    for (int j = 0; j < nrhs; ++j) {
      const auto dlij = dli[j] / dim1[j];
      di[j] -= dlij * duim1[j];
      xi[j] -= dlij * xim1[j];
    }
  }
  {
    const int ios = (nrow-1)*nrhs;
    auto* const di = d + ios;
    auto* const xi = X + ios;
    for (int j = 0; j < nrhs; ++j)
      xi[j] /= di[j];
  }
  for (int i = nrow-1; i > 0; --i) {
    const int ios = i*nrhs;
    const int im1os = (i-1)*nrhs;
    auto* const dim1 = d + im1os;
    auto* const duim1 = du + im1os;
    auto* const xim1 = X + im1os;
    auto* const xi = X + ios;
    for (int j = 0; j < nrhs; ++j)
      xim1[j] = (xim1[j] - duim1[j] * xi[j]) / dim1[j];
  }
}

} // namespace impl

template <typename TeamMember, typename TridiagDiag, typename DataArray>
KOKKOS_INLINE_FUNCTION
void thomas (const TeamMember& team,
             TridiagDiag dl, TridiagDiag d, TridiagDiag du, DataArray X,
             typename std::enable_if<TridiagDiag::rank == 1>::type* = 0) {
  const int nrow = d.extent_int(0);
  assert(X.extent_int(0) == nrow);
  assert(dl.extent_int(0) == nrow);
  assert(du.extent_int(0) == nrow);
  impl::thomas_factorize(team, dl, d, du);
  team.team_barrier();
  impl::thomas_solve(team, dl, d, du, X);
}

template <typename TridiagDiag, typename DataArray>
KOKKOS_INLINE_FUNCTION
void thomas (TridiagDiag dl, TridiagDiag d, TridiagDiag du, DataArray X,
             typename std::enable_if<TridiagDiag::rank == 1>::type* = 0,
             typename std::enable_if<DataArray::rank == 1>::type* = 0,
             impl::EnableIfCanUsePointer<TridiagDiag>* = 0,
             impl::EnableIfCanUsePointer<DataArray>* = 0) {
  const int nrow = d.extent_int(0);
  assert( X.extent_int(0) == nrow);
  assert(dl.extent_int(0) == nrow);
  assert(du.extent_int(0) == nrow);
  impl::thomas_a1x1(dl.data(), d.data(), du.data(), X.data(), nrow);
}

template <typename TridiagDiag, typename DataArray>
KOKKOS_INLINE_FUNCTION
void thomas (TridiagDiag dl, TridiagDiag d, TridiagDiag du, DataArray X,
             typename std::enable_if<TridiagDiag::rank == 1>::type* = 0,
             typename std::enable_if<DataArray::rank == 2>::type* = 0,
             impl::EnableIfCanUsePointer<TridiagDiag>* = 0,
             impl::EnableIfCanUsePointer<DataArray>* = 0) {
  const int nrow = d.extent_int(0);
  const int nrhs = X.extent_int(1);
  assert( X.extent_int(0) == nrow);
  assert(dl.extent_int(0) == nrow);
  assert(du.extent_int(0) == nrow);
  impl::thomas_a1xm(dl.data(), d.data(), du.data(), X.data(), nrow, nrhs);
}

template <typename TridiagDiag, typename DataArray>
KOKKOS_INLINE_FUNCTION
void thomas (TridiagDiag dl, TridiagDiag d, TridiagDiag du, DataArray X,
             typename std::enable_if<TridiagDiag::rank == 2>::type* = 0,
             typename std::enable_if<DataArray::rank == 2>::type* = 0,
             impl::EnableIfCanUsePointer<TridiagDiag>* = 0,
             impl::EnableIfCanUsePointer<DataArray>* = 0) {
  const int nrow = d.extent_int(0);
  const int nrhs = X.extent_int(1);
  assert(X .extent_int(0) == nrow);
  assert(dl.extent_int(0) == nrow);
  assert(du.extent_int(0) == nrow);
  assert(dl.extent_int(1) == nrhs);
  assert(d .extent_int(1) == nrhs);
  assert(du.extent_int(1) == nrhs);
  impl::thomas_amxm(dl.data(), d.data(), du.data(), X.data(), nrow, nrhs);
}

// Cyclic reduction at the Kokkos team level. Any (thread, vector)
// parameterization is intended to work.
template <typename TeamMember, typename TridiagDiag, typename DataArray>
KOKKOS_INLINE_FUNCTION
void cr (const TeamMember& team,
         TridiagDiag dl, TridiagDiag d, TridiagDiag du, DataArray X,
         typename std::enable_if<TridiagDiag::rank == 1>::type* = 0,
         typename std::enable_if<DataArray::rank == 1>::type* = 0) {
  using Scalar = typename TridiagDiag::non_const_value_type;
  const int nrow = d.extent_int(0);
  assert(dl.extent_int(0) == nrow);
  assert(d. extent_int(0) == nrow);
  assert(du.extent_int(0) == nrow);
  assert(X. extent_int(0) == nrow);
  const int team_id = impl::get_thread_id_within_team(team);
  const int nteam = impl::get_team_nthr(team);
  int os = 1, stride;
  // Go down reduction.
  while ((stride = (os << 1)) < nrow) {
    const int inc = stride*nteam;
    for (int i = stride*team_id; i < nrow; i += inc) {
      int im = i - os;
      int ip = i + os;
      // GPU does well with ternary ?: op. Use it throughout this
      // impl. It requires the trick noted in a few lines.
      const auto f1 = im >= 0   ? -dl(i)/d(im) : 0;
      const auto f2 = ip < nrow ? -du(i)/d(ip) : 0;
      // Trick to keep im, ip in bounds; the index is modified only
      // when the corresponding f is 0, so the resulting invalid
      // value is multipled by 0.
      im = im >= 0   ? im : i;
      ip = ip < nrow ? ip : i;
      dl(i)  = f1*dl(im);
      du(i)  = f2*du(ip);
      d (i) += f1*du(im) + f2*dl(ip);
      X (i) += f1*X (im) + f2*X (ip);
    }
    os <<= 1;
    // Go down in cyclic reduction level only when this level is complete.
    team.team_barrier();
  }
  // Bottom 1 or 2 levels of the reduction. This could be folded into
  // the previous loop, but it's a slight opt to handle these cases
  // separately.
  if (team_id == 0) {
    if (os >= nrow) {
      X(0) /= d(0);
    } else {
      const auto
        det = d(0)*d(os) - du(0)*dl(os),
        x0 = X(0), x1 = X(os);
      X( 0) = (d(os)*x0 - du( 0)*x1)/det;
      X(os) = (d( 0)*x1 - dl(os)*x0)/det;
    }
  }
  team.team_barrier();
  os >>= 1;
  assert(os < nrow);
  // Go up reduction.
  while (os) {
    stride = os << 1;
    const int inc = stride*nteam;
    for (int i = stride*team_id + os; i < nrow; i += inc) {
      const int im = i - os;
      const int ip = i + os;
      assert(im >= 0 || ip < nrow);
      Scalar f = 0;
      f += im >=   0 ? dl(i)*X(im) : 0;
      f += ip < nrow ? du(i)*X(ip) : 0;
      X(i) = (X(i) - f)/d(i);
    }
    os >>= 1;
    // Go up in cyclic reduction level only when this level is complete.
    team.team_barrier();
  }
}

template <typename TeamMember, typename TridiagDiag, typename DataArray>
KOKKOS_INLINE_FUNCTION
void cr (const TeamMember& team,
         TridiagDiag dl, TridiagDiag d, TridiagDiag du, DataArray X,
         typename std::enable_if<TridiagDiag::rank == 1>::type* = 0,
         typename std::enable_if<DataArray::rank == 2>::type* = 0) {
  using Scalar = typename TridiagDiag::non_const_value_type;
  const int nrow = d.extent_int(0);
  assert(X.extent_int(0) == nrow);
  assert(dl.extent_int(0) == nrow);
  assert(du.extent_int(0) == nrow);
  const int nrhs = X.extent_int(1);
  const int tid = impl::get_thread_id_within_team(team);
  const int nthr = impl::get_team_nthr(team);
  const int team_size = util::min(nrhs, nthr);
  const int nteam = nthr / team_size;
  const int team_id = tid / team_size;
  const int team_tid = tid % team_size;
  const bool team_lead = tid % team_size == 0;
  int os = 1, stride;
  while ((stride = (os << 1)) < nrow) {
    const int inc = stride*nteam;
    const int nit = (nrow + inc - 1)/inc;
    for (int i = stride*team_id, it = 0; it < nit; i += inc, ++it) {
      int im = i - os;
      int ip = i + os;
      Scalar f1 = 0, f2 = 0;
      const bool run = team_id < nteam && i < nrow;
      assert(team_id != 0 || run);
      if (run) {
        f1 = im >= 0   ? -dl(i)/d(im) : 0;
        f2 = ip < nrow ? -du(i)/d(ip) : 0;
        im = im >= 0   ? im : i;
        ip = ip < nrow ? ip : i;
        for (int j = team_tid; j < nrhs; j += team_size)
          X(i,j) += f1*X(im,j) + f2*X(ip,j);
      }
      // Update A only after all threads are done using current values.
      team.team_barrier();
      if (team_lead && run) {
        dl(i)  = f1*dl(im);
        du(i)  = f2*du(ip);
        d (i) += f1*du(im) + f2*dl(ip);
      }
    }
    os <<= 1;
    team.team_barrier();
  }
  if (team_id == 0) {
    if (os >= nrow) {
      for (int j = team_tid; j < nrhs; j += team_size)
        X(0,j) /= d(0);
    } else {
      for (int j = team_tid; j < nrhs; j += team_size) {
        const auto
          det = d(0)*d(os) - du(0)*dl(os),
          x0 = X(0,j), x1 = X(os,j);
        X( 0,j) = (d(os)*x0 - du( 0)*x1)/det;
        X(os,j) = (d( 0)*x1 - dl(os)*x0)/det;
      }
    }
  }
  team.team_barrier();
  os >>= 1;
  assert(os < nrow);
  while (os) {
    if (team_id < nteam) {
      stride = os << 1;
      const int inc = stride*nteam;
      for (int i = stride*team_id + os; i < nrow; i += inc) {
        const int im = i - os;
        const int ip = i + os;
        assert(im >= 0 || ip < nrow);
        for (int j = team_tid; j < nrhs; j += team_size) {
          Scalar f = 0;
          f += im >=   0 ? dl(i)*X(im,j) : 0;
          f += ip < nrow ? du(i)*X(ip,j) : 0;
          X(i,j) = (X(i,j) - f)/d(i);
        }
      }
    }
    os >>= 1;
    team.team_barrier();
  }
}

template <typename TeamMember, typename TridiagDiag, typename DataArray>
KOKKOS_INLINE_FUNCTION
void cr (const TeamMember& team,
         TridiagDiag dl, TridiagDiag d, TridiagDiag du, DataArray X,
         typename std::enable_if<TridiagDiag::rank == 2>::type* = 0,
         typename std::enable_if<DataArray::rank == 2>::type* = 0) {
  using Scalar = typename TridiagDiag::non_const_value_type;
  const int nrow = d.extent_int(0);
  const int nrhs = X.extent_int(1);
  assert(dl.extent_int(1) == nrhs);
  assert(d. extent_int(1) == nrhs);
  assert(du.extent_int(1) == nrhs);
  assert(dl.extent_int(0) == nrow);
  assert(d. extent_int(0) == nrow);
  assert(du.extent_int(0) == nrow);
  assert(X. extent_int(0) == nrow);
  const int tid = impl::get_thread_id_within_team(team);
  const int nthr = impl::get_team_nthr(team);
  const int team_size = util::min(nrhs, nthr);
  const int nteam = nthr / team_size;
  const int team_id = tid / team_size;
  const int team_tid = tid % team_size;
  int os = 1, stride;
  while ((stride = (os << 1)) < nrow) {
    if (team_id < nteam) {
      const int inc = stride*nteam;
      for (int i = stride*team_id; i < nrow; i += inc) {
        int im = i - os;
        int ip = i + os;
        const bool im_ok = im >= 0;
        const bool ip_ok = ip < nrow;
        im = im_ok ? im : i;
        ip = ip_ok ? ip : i;
        for (int j = team_tid; j < nrhs; j += team_size) {
          const auto f1 = im_ok ? -dl(i,j)/d(im,j) : 0;
          const auto f2 = ip_ok ? -du(i,j)/d(ip,j) : 0;
          dl(i,j)  = f1*dl(im,j);
          du(i,j)  = f2*du(ip,j);
          d (i,j) += f1*du(im,j) + f2*dl(ip,j);
          X (i,j) += f1*X (im,j) + f2*X (ip,j);
        }
      }
    }
    os <<= 1;
    team.team_barrier();
  }
  if (team_id == 0) {
    if (os >= nrow) {
      for (int j = team_tid; j < nrhs; j += team_size)
        X(0,j) /= d(0,j);
    } else {
      for (int j = team_tid; j < nrhs; j += team_size) {
        const auto
          det = d(0,j)*d(os,j) - du(0,j)*dl(os,j),
          x0 = X(0,j), x1 = X(os,j);
        X( 0,j) = (d(os,j)*x0 - du( 0,j)*x1)/det;
        X(os,j) = (d( 0,j)*x1 - dl(os,j)*x0)/det;
      }
    }
  }
  team.team_barrier();
  os >>= 1;
  assert(os < nrow);
  while (os) {
    if (team_id < nteam) {
      stride = os << 1;
      const int inc = stride*nteam;
      for (int i = stride*team_id + os; i < nrow; i += inc) {
        const int im = i - os;
        const int ip = i + os;
        assert(im >= 0 || ip < nrow);
        const bool im_ok = im >= 0;
        const bool ip_ok = ip < nrow;
        for (int j = team_tid; j < nrhs; j += team_size) {
          Scalar f = 0;
          f += im_ok ? dl(i,j)*X(im,j) : 0;
          f += ip_ok ? du(i,j)*X(ip,j) : 0;
          X(i,j) = (X(i,j) - f)/d(i,j);
        }
      }
    }
    os >>= 1;
    team.team_barrier();
  }
}

} // namespace tridiag
} // namespace scream

#endif // INCLUDE_SCREAM_TRIDIAG
