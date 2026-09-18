/**
 * @file fpe_guard.hpp
 * @brief RAII suspension of floating-point exception traps.
 */

#ifndef E3SM_EMULATOR_FPE_GUARD_HPP
#define E3SM_EMULATOR_FPE_GUARD_HPP

namespace emulator {
namespace inference {

/**
 * @brief Disable FPE traps for the lifetime of the object, then restore
 *        the set that was enabled before.
 *
 * E3SM debug builds trap on invalid, divide-by-zero and overflow, which ML
 * library kernels raise benignly. EAMxx's PySession does the same.
 *
 * A no-op where feenableexcept does not exist (it is a glibc extension).
 */
class FpeGuard {
public:
  FpeGuard();
  ~FpeGuard();
  FpeGuard(const FpeGuard &) = delete;
  FpeGuard &operator=(const FpeGuard &) = delete;

private:
  int m_saved_excepts = 0;
};

} // namespace inference
} // namespace emulator

#endif // E3SM_EMULATOR_FPE_GUARD_HPP
