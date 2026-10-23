#ifndef SCREAM_COSP_RNG_IMPL_HPP
#define SCREAM_COSP_RNG_IMPL_HPP

#include "physics/cosp/cosp_cxx_functions.hpp"

namespace scream {
namespace cosp {

/*
 * Marsaglia CONG generator, mirroring get_rng_1 in mo_rng.F90. To be BFB:
 *  - The Fortran update 69069*seed+1234567 overflows a 32-bit integer and wraps. Signed
 *    overflow is UB in C++, so do it in uint32_t and convert back (modular since C++20).
 *  - Fortran mod() and C++ % both take the sign of the dividend.
 *  - The scaling to [0,1) is done in single precision, with a single precision literal.
 *  - Default integers are 32 bit, so the r = r+1; r = r-int(r) fold is always applied.
 */
template <typename S, typename D>
KOKKOS_INLINE_FUNCTION
S Functions<S,D>::get_rng(RngState& seed)
{
  const std::uint32_t u = 69069u*static_cast<std::uint32_t>(seed) + 1234567u;
  seed = static_cast<std::int32_t>(u) % (std::int32_t(1) << 30);

  float r = static_cast<float>(seed)*0.931322574615479E-09f;
  r = r + 1.0f;
  r = r - static_cast<float>(static_cast<int>(r));
  return static_cast<S>(r);
}

} // namespace cosp
} // namespace scream

#endif // SCREAM_COSP_RNG_IMPL_HPP
