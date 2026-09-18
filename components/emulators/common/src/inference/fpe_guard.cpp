/**
 * @file fpe_guard.cpp
 * @brief FpeGuard implementation.
 */

#include "fpe_guard.hpp"

#include <cfenv>

namespace emulator {
namespace inference {

FpeGuard::FpeGuard() {
#ifdef __GLIBC__
  m_saved_excepts = fegetexcept();
  if (m_saved_excepts > 0) {
    fedisableexcept(m_saved_excepts);
  }
#endif
}

FpeGuard::~FpeGuard() {
#ifdef __GLIBC__
  if (m_saved_excepts > 0) {
    // Clear flags raised while guarded, or re-enabling would trap on them.
    feclearexcept(FE_ALL_EXCEPT);
    feenableexcept(m_saved_excepts);
  }
#endif
}

} // namespace inference
} // namespace emulator
