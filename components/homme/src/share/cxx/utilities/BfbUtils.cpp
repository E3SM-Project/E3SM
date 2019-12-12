#ifdef HAVE_CONFIG_H
#include "config.h.c"
#endif

#include "BfbUtils.hpp"
#include "scream_tridiag.hpp"

template<typename RealType>
void tridiag_diagdom_bfb_a1x1_impl (int n, RealType* dl, RealType* d, RealType* du, RealType* x) {
  Kokkos::View<RealType*, Kokkos::DefaultHostExecutionSpace>
    dlv("dlv", n), dv(d, n), duv(du, n), xv(x, n);
  dlv(0) = 0;
  for (int i = 1; i < n; ++i) dlv(i) = dl[i-1];
  scream::tridiag::impl::bfb_thomas_factorize(dlv, dv, duv);
  scream::tridiag::impl::bfb_thomas_solve(dlv, dv, duv, xv);
}

// For F90.
extern "C" {
void czeroulpn(int a_len, double* a, int nbit, double* replace) {
  for (int i = 0; i < a_len; ++i)
    a[i] = Homme::zeroulpn(a[i], nbit, replace[i]);
}

void cbfb_pow(int a_len, double* a, const double e) {
  for (int i = 0; i < a_len; ++i) {
    a[i] = std::pow(a[i], e);
#ifdef CUDA_BUILD
    const auto replace = a[i];
    a[i] = Homme::zeroulpn(a[i], 25, replace);
#endif
  }
}

void tridiag_diagdom_bfb_a1x1 (int n, void* dl, void* d, void* du, void* x, int real_size) {
  assert (real_size==4 || real_size==8);
  if (real_size==4) {
    tridiag_diagdom_bfb_a1x1_impl<float>(n,
      reinterpret_cast<float*>(dl),
      reinterpret_cast<float*>(d),
      reinterpret_cast<float*>(du),
      reinterpret_cast<float*>(x));
  } else if (real_size==8) {
    tridiag_diagdom_bfb_a1x1_impl<double>(n,
      reinterpret_cast<double*>(dl),
      reinterpret_cast<double*>(d),
      reinterpret_cast<double*>(du),
      reinterpret_cast<double*>(x));
  }
}

} // extern C
