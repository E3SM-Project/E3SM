#ifndef SCREAM_ASSERT_HPP
#define SCREAM_ASSERT_HPP

#include <sstream>
#include <exception>
#include <stdexcept>  // For std::logic_error

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

namespace scream {
namespace error {

void runtime_check(bool cond, const std::string& message, int code = -1);
void runtime_abort(const std::string& message, int code = -1);

} // namespace error
} // namespace scream

#endif // SCREAM_ASSERT_HPP
