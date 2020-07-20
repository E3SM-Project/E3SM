#ifndef EKAT_ASSERT_HPP
#define EKAT_ASSERT_HPP

#include <sstream>
#include <exception>
#include <assert.h>
#include <stdexcept>  // For std::logic_error
#include <cfenv>

/*
 * Asserts and error checking for Scream.
 *
 * scream_k* are for error checking within kokkos kernels.
 *
 * Any check with "assert" in the name is disabled for release builds
 *
 * For _msg checks, the msg argument can contain '<<' if not a kernel check.
 */

// Internal do not call directly
#define impl_throw(condition, msg, exception_type)                      \
  do {                                                                  \
    if ( ! (condition) ) {                                              \
      std::stringstream _ss_;                                           \
      _ss_ << __FILE__ << ":" << __LINE__ << ": FAIL:\n" << #condition; \
      _ss_ << "\n" << msg;                                              \
      throw exception_type(_ss_.str());                                 \
    }                                                                   \
  } while(0)

#define impl_kthrow(condition, msg)             \
  do {                                          \
    if ( ! (condition) )                        \
      Kokkos::abort(#condition "\n" msg);       \
  } while (0)

#ifndef NDEBUG
#define scream_assert(condition)           impl_throw(condition, "", std::logic_error)
#define scream_kassert(condition)          impl_kthrow(condition, "")
#define scream_assert_msg(condition, msg)  impl_throw(condition, msg, std::logic_error)
#define scream_kassert_msg(condition, msg) impl_kthrow(condition, msg)
#else
#define scream_assert(condition)  ((void) (0))
#define scream_kassert(condition) ((void) (0))
#define scream_assert_msg(condition, msg)  ((void) (0))
#define scream_kassert_msg(condition, msg) ((void) (0))
#endif

#define scream_require(condition)           impl_throw(condition, "", std::logic_error)
#define scream_krequire(condition)          impl_kthrow(condition, "")
#define scream_require_msg(condition, msg)  impl_throw(condition, msg, std::logic_error)
#define scream_krequire_msg(condition, msg) impl_kthrow(condition, msg)

#define scream_error_msg(msg)               scream_require_msg(false, msg)
#define scream_kerror_msg(msg)              scream_krequire_msg(false, msg)

// Macros to do asserts inside constexpr functions
// This is not necessary with C++14, where constexpr functions are "regular"
// functions (they just can be evaluated at compile time). But in C++11 you can only have
// 'return blah;' in a constexpr function, and user-defined constexpr constructors must
// have an empty body (only initializer list syntax allowed).
// NOTE: doing `return assert(check), blah;` seems the right solution, but it seems
//       older versions of gcc may not like this option even if check evaluates to true.
// As soon as we can use C++14, this macro can be removed, and you can simply use
// `assert(check);` in the body of the constexpr function.
#if defined(NDEBUG) || !defined(EKAT_CONSTEXPR_ASSERT)
#define CONSTEXPR_ASSERT(CHECK) void(CHECK)
#else
namespace impl {
struct assert_failure {
  explicit assert_failure(const char *sz)
  {
    std::fprintf(stderr, "Assertion failure: %s\n", sz);
  }
};
}
#define CONSTEXPR_ASSERT(CHECK) \
    ( (CHECK) || (throw impl::assert_failure(#CHECK), false) )
#endif

namespace scream {
namespace error {

void runtime_check(bool cond, const std::string& message, int code = -1);
void runtime_abort(const std::string& message, int code = -1);

} // namespace error

/*
 * Routines to activate/deactivate floating point exceptions.
 * These routines are only meaningful if EKAT_FPE is defined.
 * The last two functions activate/deactivate a predefined set
 * of FPEs: FE_DIVBYZERO, FE_INVALID, and FE_OVERFLOW.
 * If you need to temporarily enable/disable a specific exception,
 * you should use the first set of functions, which accept
 * a specific fpe mask.
 *
 * WARNING: be *very* careful when using these functions, as they
 *          effectively change the FPE environment for the rest
 *          of the translation unit. If you only need to turn off
 *          FPEs for a few lines (perhaps you do some dirty trick
 *          that is correct, but would throw an fpe), don't forget
 *          to re-enable them after you're done.
 */

#ifdef EKAT_FPE
// #pragma STDC FENV_ACCESS ON
#endif

static unsigned int constexpr scream_default_fpes =
#ifdef EKAT_FPE
  FE_DIVBYZERO |
  FE_INVALID   |
  FE_OVERFLOW;
#else
  0;
#endif

void enable_fpes (const int mask);
void disable_fpes (const int mask);

inline void enable_default_fpes () {
  enable_fpes(scream_default_fpes);
}
inline void disable_default_fpes () {
  disable_fpes(scream_default_fpes);
}

int get_enabled_fpes ();
void disable_all_fpes ();

} // namespace scream

#endif // EKAT_ASSERT_HPP
