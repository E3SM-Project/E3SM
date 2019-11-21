#include "scream_tridiag.hpp"
#include "Types.hpp"

using Real = Homme::Real;

extern "C" {
void tridiag_diagdom_bfb_a1x1 (int n, Real* dl, Real* d, Real* du, Real* x) {
  Kokkos::View<Real*, Kokkos::DefaultHostExecutionSpace>
    dlv("dlv", n), dv(d, n), duv(du, n), xv(x, n);
  dlv(0) = 0;
  for (int i = 1; i < n; ++i) dlv(i) = dl[i-1];
  scream::tridiag::impl::bfb_thomas_factorize(dlv, dv, duv);
  scream::tridiag::impl::bfb_thomas_solve(dlv, dv, duv, xv);
}
}
