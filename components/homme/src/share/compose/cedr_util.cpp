// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

#include "cedr_util.hpp"

#include <random>

namespace cedr {
namespace util {

bool eq (const std::string& a, const char* const b1, const char* const b2) {
  return (a == std::string(b1) || (b2 && a == std::string(b2)) ||
          a == std::string("-") + std::string(b1));
}

namespace {
std::mt19937_64& urand_engine () {
  static std::mt19937_64 engine(1);
  return engine;
}
}

Real urand () {
  std::uniform_real_distribution<Real> pdf(0, 1);
  return pdf(urand_engine());
}

void set_urand_seed (const unsigned int seed) { urand_engine().seed(seed); }

Real reldif (const Real* a, const Real* b, const Int n) {
  Real num = 0, den = 0;
  for (Int i = 0; i < n; ++i) {
    num += std::abs(a[i] - b[i]);
    den += std::abs(a[i]);
  }
  return num/den;
}

}
}
