// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

#include "cedr_local.hpp"
#include "cedr_local_inl.hpp"

namespace cedr {
namespace local {
namespace test {
// Check the first-order optimality conditions. Return true if OK, false
// otherwise. If quiet, don't print anything.
bool check_1eq_bc_qp_foc (
  const char* label, const Int n, const Real* w, const Real* a, const Real b,
  const Real* xlo, const Real* xhi, const Real* y, const Real* x, const bool verbose)
{
  auto& os = std::cout;
  bool ok = true;
  Real xtmp;
  // Check the bound constraints.
  for (Int i = 0; i < n; ++i)
    if (x[i] < (xtmp = xlo[i])) {
      if (verbose)
        os << "x[" << i << "] = " << x[i]
           << " but x[i] - xlo[i] = " << (x[i] - xtmp) << "\n";
      ok = false;
    }
  for (Int i = 0; i < n; ++i)
    if (x[i] > (xtmp = xhi[i])) {
      if (verbose)
        os << "x[" << i << "] = " << x[i]
           << " but xhi[i] - x[i] = " << (xtmp - x[i]) << "\n";
      ok = false;
    }
  // Check the equality constraint.
  Real r = 0;
  for (Int i = 0; i < n; ++i)
    r += a[i]*x[i];
  r -= b;
  if (std::abs(r) > impl::calc_r_tol(b, a, y, n)) {
    if (verbose)
      os << "r = " << r << "\n";
    ok = false;
  }
  // Check the gradient is 0 when projected into the constraints. Compute
  //     g = W (x - y)
  //     g_reduced = g - C ((C'C) \ (C'g))
  // where
  //     IA = I(:,A)
  //     C = [IA a],
  // and A is the active set.
  const Real padtol = 1e5*std::numeric_limits<Real>::epsilon();
  Real lambda = 0, den = 0;
  for (Int i = 0; i < n; ++i) {
    const Real pad = padtol*(xhi[i] - xlo[i]);
    if (xlo[i] + pad <= x[i] && x[i] <= xhi[i] - pad) {
      const Real gi = w[i]*(x[i] - y[i]);
      lambda += a[i]*gi;
      den += a[i]*a[i];
    }
  }
  lambda /= den;
  Real normg = 0, normy = 0;
  for (Int i = 0; i < n; ++i) {
    normy += cedr::util::square(y[i]);
    const Real pad = padtol*(xhi[i] - xlo[i]);
    if (xlo[i] + pad <= x[i] && x[i] <= xhi[i] - pad)
      normg += cedr::util::square(w[i]*(x[i] - y[i]) - a[i]*lambda);
  }
  normy = std::sqrt(normy);
  normg = std::sqrt(normg);
  const Real gtol = 1e4*std::numeric_limits<Real>::epsilon()*normy;
  if (normg > gtol) {
    if (verbose)
      os << "norm(g) = " << normg << " gtol = " << gtol << "\n";
    ok = false;
  }
  // Check the gradient at the active boundaries.
  for (Int i = 0; i < n; ++i) {
    const bool onlo = x[i] == xlo[i];
    const bool onhi = onlo ? false : x[i] == xhi[i];
    if (onlo || onhi) {
      const Real rg = w[i]*(x[i] - y[i]) - a[i]*lambda;
      if (onlo && rg < -gtol) {
        if (verbose)
          os << "onlo but rg = " << rg << "\n";
        ok = false;
      } else if (onhi && rg > gtol) {
        if (verbose)
          os << "onhi but rg = " << rg << "\n";
        ok = false;
      }
    }
  }
  if ( ! ok && verbose)
    os << "label: " << label << "\n";
  return ok;
}

Int test_1eq_bc_qp () {
  bool verbose = true;
  Int nerr = 0;

  Int n;
  static const Int N = 16;
  Real w[N], a[N], b, xlo[N], xhi[N], y[N], x[N], al, au;

  auto run = [&] () {
    const Int info = solve_1eq_bc_qp(n, w, a, b, xlo, xhi, y, x);
    const bool ok = test::check_1eq_bc_qp_foc(
      "unittest", n, w, a, b, xlo, xhi, y, x, verbose);
    if ( ! ok) ++nerr;

    if (n == 2) {
      // This version never returns 0.
      Real x2[2];
      const Int info2 = solve_1eq_bc_qp_2d(w, a, b, xlo, xhi, y, x2);
      if (info2 != 1 && (info == 0 || info == 1)) {
        if (verbose) pr(puf(info) pu(info2));
        ++nerr;
      }
      const Real rd = cedr::util::reldif(x, x2, 2);
      if (rd > 1e4*std::numeric_limits<Real>::epsilon()) {
        if (verbose)
          printf("%1.1e | y %1.15e %1.15e | x %1.15e %1.15e | "
                 "x2 %1.15e %1.15e | l %1.15e %1.15e | u %1.15e %1.15e\n",
                 rd, y[0], y[1], x[0], x[1], x2[0], x2[1],
                 xlo[0], xlo[1], xhi[0], xhi[1]);
        ++nerr;
      }
    }

    caas(n, a, b, xlo, xhi, y, x);
    Real m = 0, den = 0;
    for (Int i = 0; i < n; ++i) {
      m += a[i]*x[i];
      den += std::abs(a[i]*x[i]);
      if (x[i] < xlo[i]) ++nerr;
      else if (x[i] > xhi[i]) ++nerr;
    }
    const Real rd = std::abs(b - m)/den;
    if (rd > 1e3*std::numeric_limits<Real>::epsilon()) {
      if (verbose) pr(puf(rd) pu(n) pu(b) pu(m));
      ++nerr;
    }
  };

  auto gena = [&] () {
    for (Int i = 0; i < n; ++i)
      a[i] = 0.1 + cedr::util::urand();
  };
  auto genw = [&] () {
    for (Int i = 0; i < n; ++i)
      w[i] = 0.1 + cedr::util::urand();
  };
  auto genbnds = [&] () {
    al = au = 0;
    for (Int i = 0; i < n; ++i) {
      xlo[i] = cedr::util::urand() - 0.5;
      al += a[i]*xlo[i];
      xhi[i] = xlo[i] + cedr::util::urand();
      au += a[i]*xhi[i];
    }
  };
  auto genb = [&] (const bool in) {
    if (in) {
      const Real alpha = cedr::util::urand();
      b = alpha*al + (1 - alpha)*au;
    } else {
      if (cedr::util::urand() > 0.5)
        b = au + 0.01 + cedr::util::urand();
      else
        b = al - 0.01 - cedr::util::urand();
    }
  };
  auto geny = [&] (const bool in) {
    if (in) {
      for (Int i = 0; i < n; ++i) {
        const Real alpha = cedr::util::urand();
        y[i] = alpha*xlo[i] + (1 - alpha)*xhi[i];
      }
    } else if (cedr::util::urand() > 0.2) {
      for (Int i = 1; i < n; i += 2) {
        const Real alpha = cedr::util::urand();
        y[i] = alpha*xlo[i] + (1 - alpha)*xhi[i];
        cedr_assert(y[i] >= xlo[i] && y[i] <= xhi[i]);
      }      
      for (Int i = 0; i < n; i += 4)
        y[i] = xlo[i] - cedr::util::urand();
      for (Int i = 2; i < n; i += 4)
        y[i] = xhi[i] + cedr::util::urand();
    } else {
      for (Int i = 0; i < n; i += 2)
        y[i] = xlo[i] - cedr::util::urand();
      for (Int i = 1; i < n; i += 2)
        y[i] = xhi[i] + cedr::util::urand();
    }
  };
  auto b4y = [&] () {
    b = 0;
    for (Int i = 0; i < n; ++i)
      b += a[i]*y[i];
  };

  for (n = 2; n <= N; ++n) {
    const Int count = n == 2 ? 100 : 10;
    for (Int i = 0; i < count; ++i) {
      gena();
      genw();
      genbnds();
      genb(true);
      geny(true);
      run();
      b4y();
      run();
      genb(true);
      geny(false);
      run();
    }
  }

  return  nerr;
}

Int test_1eq_nonneg () {
  using cedr::util::urand;
  using cedr::util::reldif;

  bool verbose = true;
  Int nerr = 0;

  Int n;
  static const Int N = 16;
  Real w[N], a[N], b, xlo[N], xhi[N], y[N], x_ls[N], x_caas[N], x1_ls[N], x1_caas[N];

  for (n = 2; n <= 2; ++n) {
    const Int count = 20;
    for (Int trial = 0; trial < count; ++trial) {
      b = 0.5*n*urand();
      for (Int i = 0; i < n; ++i) {
        w[i] = 0.1 + urand();
        a[i] = 0.1 + urand();
        xlo[i] = 0;
        xhi[i] = b/a[i];
        y[i] = urand();
        if (urand() > 0.8) y[i] *= -1;
        x1_caas[i] = urand() > 0.5 ? y[i] : -1;
        x1_ls[i] = urand() > 0.5 ? y[i] : -1;
      }
      solve_1eq_nonneg(n, a, b, y, x1_caas, w, Method::caas);
      caas(n, a, b, xlo, xhi, y, x_caas);
      solve_1eq_nonneg(n, a, b, y, x1_ls, w, Method::least_squares);
      solve_1eq_bc_qp(n, w, a, b, xlo, xhi, y, x_ls);
      const Real rd_caas = reldif(x_caas, x1_caas, 2);
      const Real rd_ls = reldif(x_ls, x1_ls, 2);
      if (rd_ls > 5e1*std::numeric_limits<Real>::epsilon() ||
          rd_caas > 1e1*std::numeric_limits<Real>::epsilon()) {
        pr(puf(rd_ls) pu(rd_caas));
        if (verbose) {
          using cedr::util::prarr;
          prarr("w", w, n);
          prarr("a", a, n);
          prarr("xhi", xhi, n);
          prarr("y", y, n);
          prarr("x_ls", x_ls, n);
          prarr("x1_ls", x1_ls, n);
          prarr("x_caas", x_caas, n);
          prarr("x1_caas", x1_caas, n);
          prc(b);
          Real mass = 0;
          for (Int i = 0; i < n; ++i) mass += a[i]*x_ls[i];
          prc(mass);
          mass = 0;
          for (Int i = 0; i < n; ++i) mass += a[i]*x_ls[i];
          prc(mass);
          mass = 0;
          for (Int i = 0; i < n; ++i) mass += a[i]*x_caas[i];
          prc(mass);
        }
        ++nerr;
      }
    }
  }

  return nerr;
}

} // namespace test

Int unittest () {
  Int nerr = 0;
  nerr += test::test_1eq_bc_qp();
  nerr += test::test_1eq_nonneg();
  return nerr;
}

} // namespace local
} // namespace cedr
