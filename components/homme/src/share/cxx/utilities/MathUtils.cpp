#include "MathUtils.hpp"

// For F90.
extern "C" {
void czeroulpn(int a_len, double* a, int nbit, double* replace) {
  for (int i = 0; i < a_len; ++i)
    a[i] = Homme::zeroulpn(a[i], nbit, replace[i]);
}

void cbfb_pow(int a_len, double* a, const double e) {
  for (int i = 0; i < a_len; ++i) {
    const auto replace = a[i];
    a[i] = std::pow(a[i], e);
    if (Homme::OnGpu<Homme::ExecSpace>::value)
      a[i] = Homme::zeroulpn(a[i], 25, replace);
  }
}
} // extern C
