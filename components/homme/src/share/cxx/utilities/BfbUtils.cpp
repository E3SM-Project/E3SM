#ifdef HAVE_CONFIG_H
#include "config.h.c"
#endif

#include "BfbUtils.hpp"
#include "ExecSpaceDefs.hpp"
#include "scream_tridiag.hpp"

extern "C" {
// These can be used by fortran targets to initialize/finalize kokkos
// when running bfb end-to-end runs.
void initialize_kokkos_f90 () {
  // Use this rather than Kokkos::initialize. See comment in ExecSpaceDefs.hpp
  Homme::initialize_kokkos();
}

void finalize_kokkos_f90 () {
  Kokkos::finalize();
}

} // extern "C"

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

// Allows F90 to call c++ log() on device for BFB testing
extern "C" {
double cxx_log(double input)
{
#ifdef HOMMEXX_ENABLE_GPU
  double result;
  Kokkos::RangePolicy<Homme::ExecSpace> policy(0,1);
  Kokkos::parallel_reduce(policy, KOKKOS_LAMBDA(const int&, double& value) {
    value = log(input);
  }, result);
  return result;
#else
  return log(input);
#endif
}
} // extern C
