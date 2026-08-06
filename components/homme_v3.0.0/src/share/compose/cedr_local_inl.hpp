// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

#ifndef INCLUDE_CEDR_LOCAL_INL_HPP
#define INCLUDE_CEDR_LOCAL_INL_HPP

#include "cedr_util.hpp"

namespace cedr {
namespace local {

namespace impl {
KOKKOS_INLINE_FUNCTION
Real calc_r_tol (const Real b, const Real* a, const Real* y, const Int n) {
  Real ab = std::abs(b);
  for (Int i = 0; i < n; ++i) ab = cedr::impl::max(ab, std::abs(a[i]*y[i]));
  return 1e1*cedr::impl::TypeTraits<Real>::epsilon*std::abs(ab);
}

// Eval r at end points to check for feasibility, and also possibly a quick exit
// on a common case. Return -1 if infeasible, 1 if a corner is a solution, 0 if
// feasible and a corner is not.
KOKKOS_INLINE_FUNCTION
Int check_lu (const Int n, const Real* a, const Real& b, const Real* xlo,
              const Real* xhi, const Real& r_tol, Real* x) {
  Real r = -b;
  for (Int i = 0; i < n; ++i) {
    x[i] = xlo[i];
    r += a[i]*x[i];
  }
  if (std::abs(r) <= r_tol) return 1;
  if (r > 0) return -1;
  r = -b;
  for (Int i = 0; i < n; ++i) {
    x[i] = xhi[i];
    r += a[i]*x[i];
  }
  if (std::abs(r) <= r_tol) return 1;
  if (r < 0) return -1;
  return 0;
}

KOKKOS_INLINE_FUNCTION
void calc_r (const Int n, const Real* w, const Real* a, const Real b,
             const Real* xlo, const Real* xhi,  const Real* y, const Real& lambda,
             Real* x, Real& r, Real& r_lambda) {
  r = 0;
  r_lambda = 0;
  for (Int i = 0; i < n; ++i) {
    const Real q = a[i]/w[i];
    const Real x_trial = y[i] + lambda*q;
    Real xtmp;
    if (x_trial < (xtmp = xlo[i]))
      x[i] = xtmp;
    else if (x_trial > (xtmp = xhi[i]))
      x[i] = xtmp;
    else {
      x[i] = x_trial;
      r_lambda += a[i]*q;
    }
    r += a[i]*x[i];
  }
  r -= b;
}
} // namespace impl

// 2D special case for efficiency.
KOKKOS_INLINE_FUNCTION
Int solve_1eq_bc_qp_2d (const Real* w, const Real* a, const Real b,
                        const Real* xlo, const Real* xhi, 
                        const Real* y, Real* x,
                        const bool clip, const bool early_exit_on_tol) {
  Int info;
  if (early_exit_on_tol) {
    const Real r_tol = impl::calc_r_tol(b, a, y, 2);
    Int info = impl::check_lu(2, a, b, xlo, xhi, r_tol, x);
    if (info == -1) return info;
  }

  { // Check if the optimal point ignoring bound constraints is in bounds.
    Real qmass = 0, dm = b;
    for (int i = 0; i < 2; ++i) {
      const Real qi = a[i]/w[i];
      qmass += a[i]*qi;
      dm -= a[i]*y[i];
    }
    const Real lambda = dm/qmass;
    bool ok = true;
    for (int i = 0; i < 2; ++i) {
      x[i] = y[i] + lambda*(a[i]/w[i]);
      if (x[i] < xlo[i] || x[i] > xhi[i]) {
        ok = false;
        break;
      }
    }
    if (ok) return 1;
  }

  // Solve for intersection of a'x = b, given by the parameterized line
  //     p(alpa) = x_base + alpha x_dir,
  // with a bounding line.

  // Get parameterized line.
  Real x_base[2];
  for (int i = 0; i < 2; ++i)
    x_base[i] = 0.5*b/a[i];
  Real x_dir[] = {-a[1], a[0]};

  // Get the 4 alpha values.
  Real alphas[4];
  alphas[0] = (xlo[1] - x_base[1])/x_dir[1]; // bottom
  alphas[1] = (xhi[0] - x_base[0])/x_dir[0]; // right
  alphas[2] = (xhi[1] - x_base[1])/x_dir[1]; // top
  alphas[3] = (xlo[0] - x_base[0])/x_dir[0]; // left

  // Find the middle two in the sorted alphas.
  Real min = alphas[0], max = min;
  Int imin = 0, imax = 0;
  for (Int i = 1; i < 4; ++i) {
    const Real alpha = alphas[i];
    if (alpha < min) { min = alpha; imin = i; }
    if (alpha > max) { max = alpha; imax = i; }
  }
  Int ais[2];
  Int cnt = 0;
  for (Int i = 0; i < 4; ++i)
    if (i != imin && i != imax) {
      ais[cnt++] = i;
      if (cnt == 2) break;
    }

  Real objs[2];
  Real alpha_mid = 0;
  for (Int j = 0; j < 2; ++j) {
    const Real alpha = alphas[ais[j]];
    alpha_mid += alpha;
    Real obj = 0;
    for (Int i = 0; i < 2; ++i) {
      x[i] = x_base[i] + alpha*x_dir[i];
      obj += w[i]*cedr::util::square(y[i] - x[i]);
    }
    objs[j] = obj;
  }

  const Int ai = ais[objs[0] <= objs[1] ? 0 : 1];

  info = 1;
  Int i0 = 0;
  switch (ai) {
  case 0: case 2:
    x[1] = ai == 0 ? xlo[1] : xhi[1];
    i0 = 1;
    break;
  case 1: case 3:
    x[0] = ai == 1 ? xhi[0] : xlo[0];
    i0 = 0;
    break;
  default: cedr_kernel_assert(0); info = -2;
  }
  const Int i1 = (i0 + 1) % 2;
  x[i1] = (b - a[i0]*x[i0])/a[i1];
  if (clip)
    x[i1] = cedr::impl::min(xhi[i1], cedr::impl::max(xlo[i1], x[i1]));
  return info;
}

KOKKOS_INLINE_FUNCTION
Int solve_1eq_bc_qp (const Int n, const Real* w, const Real* a, const Real b,
                     const Real* xlo, const Real* xhi, const Real* y, Real* x,
                     const Int max_its) {
  const Real r_tol = impl::calc_r_tol(b, a, y, n);
  Int info = impl::check_lu(n, a, b, xlo, xhi, r_tol, x);
  if (info != 0) return info;

  for (int i = 0; i < n; ++i)
    if (x[i] != y[i]) {
      info = 1;
      x[i] = y[i];
    }

  // In our use case, the caller has already checked (more cheaply) for a quick
  // exit.
#if 0
  { // Check for a quick exit.
    bool all_in = true;
    Real r = 0;
    for (Int i = 0; i < n; ++i) {
      if (x[i] < xlo[i] || x[i] > xhi[i]) {
        all_in = false;
        break;
      }
      r += a[i]*x[i];
    }
    if (all_in) {
      r -= b;
      if (std::abs(r) <= r_tol)
        return info;
    }
  }
#endif

  const Real wall_dist = 1e-3;

  // Get lambda endpoints.
  Real lamlo = 0, lamhi = 0;
  for (Int i = 0; i < n; ++i) {
    const Real rq = w[i]/a[i];
    const Real lamlo_i = rq*(xlo[i] - y[i]);
    const Real lamhi_i = rq*(xhi[i] - y[i]);
    if (i == 0) {
      lamlo = lamlo_i;
      lamhi = lamhi_i;
    } else {
      lamlo = cedr::impl::min(lamlo, lamlo_i);
      lamhi = cedr::impl::max(lamhi, lamhi_i);
    }
  }
  const Real lamlo_feas = lamlo, lamhi_feas = lamhi;
  Real lambda = lamlo <= 0 && lamhi >= 0 ? 0 : lamlo;

  // Bisection-safeguarded Newton iteration for r(lambda) = 0.
  bool prev_step_bisect = false;
  Int nbisect = 0;
  info = -2;
  for (Int iteration = 0; iteration < max_its; ++iteration) {
    // Compute x, r, r_lambda.
    Real r, r_lambda;
    impl::calc_r(n, w, a, b, xlo, xhi, y, lambda, x, r, r_lambda);
    // Is r(lambda) - b sufficiently == 0?
    if (std::abs(r) <= r_tol) {
      info = 1;
      break;
    }
    // Check if the lambda bounds are too close.
    if (nbisect > 64) {
      if (lamhi == lamhi_feas || lamlo == lamlo_feas) {
        // r isn't small enough and one lambda bound is on the feasibility
        // limit. The QP must not be feasible.
        info = -1;
        break;
      }
      info = 1;
      break;
    }
    // Adjust lambda bounds.
    if (r > 0)
      lamhi = lambda;
    else
      lamlo = lambda;
    if (r_lambda != 0) {
      // Newton step.
      lambda -= r/r_lambda;
    } else {
      // Force bisection.
      lambda = lamlo;
    }
    // Safeguard. The wall distance check assures progress, but use it only
    // every other potential bisection.
    const Real D = prev_step_bisect ? 0 : wall_dist*(lamhi - lamlo);
    if (lambda - lamlo < D || lamhi - lambda < D) {
      lambda = 0.5*(lamlo + lamhi);
      ++nbisect;
      prev_step_bisect = true;
    } else {
      prev_step_bisect = false;
    }
  }

  return info;
}

KOKKOS_INLINE_FUNCTION
void caas (const Int n, const Real* a, const Real b,
           const Real* xlo, const Real* xhi,
           const Real* y, Real* x,
           const bool clip) {
  Real dm = b;
  for (Int i = 0; i < n; ++i) {
    x[i] = cedr::impl::max(xlo[i], cedr::impl::min(xhi[i], y[i]));
    dm -= a[i]*x[i];
  }
  if (dm == 0) return;
  if (dm > 0) {
    Real fac = 0;
    for (Int i = 0; i < n; ++i)
      fac += a[i]*(xhi[i] - x[i]);
    if (fac > 0) {
      fac = dm/fac;
      for (Int i = 0; i < n; ++i)
        x[i] += fac*(xhi[i] - x[i]);
    }
  } else if (dm < 0) {
    Real fac = 0;
    for (Int i = 0; i < n; ++i)
      fac += a[i]*(x[i] - xlo[i]);
    if (fac > 0) {
      fac = dm/fac;
      for (Int i = 0; i < n; ++i)
        x[i] += fac*(x[i] - xlo[i]);
    }
  }
  // Clip again for numerics.
  if (clip)
    for (Int i = 0; i < n; ++i)
      x[i] = cedr::impl::max(xlo[i], cedr::impl::min(xhi[i], x[i]));
}

KOKKOS_INLINE_FUNCTION
Int solve_1eq_nonneg (const Int n, const Real* a, const Real b, const Real* y, Real* x,
                      const Real* w,  const Method::Enum method) {
  cedr_kernel_assert(n <= 16);
  if (b < 0) return -1;

  const Real zero[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  // Set the upper bound to the value that implies that just one slot gets all
  // of the mass.
  Real xhi[16];
  for (int i = 0; i < n; ++i)
    xhi[i] = b/a[i];

  if (method == Method::caas) {
    caas(n, a, b, zero, xhi, y, x);
    return 1;
  } else {
    if (n == 2)
      return solve_1eq_bc_qp_2d(w, a, b, zero, xhi, y, x);
    else
      return solve_1eq_bc_qp(n, w, a, b, zero, xhi, y, x);
  }
}

} // namespace local
} // namespace cedr

#endif
