#include "MathUtils.hpp"

// For F90.
extern "C" {
void czeroulpn(int a_len, double* a, int nbit) {
  for (int i = 0; i < a_len; ++i)
    a[i] = Homme::zeroulpn(a[i], nbit);
}

void cbfb_pow(int a_len, double* a, const double e) {
  for (int i = 0; i < a_len; ++i)
    a[i] = std::pow(a[i], e);
  if (Homme::OnGpu<Homme::ExecSpace>::value)
    czeroulpn(a_len, a, 20);
}
} // extern C
