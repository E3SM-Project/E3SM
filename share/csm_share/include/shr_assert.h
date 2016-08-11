#ifdef NDEBUG
#define SHR_ASSERT(assert, msg)
#define SHR_ASSERT_ALL(assert, msg)
#define SHR_ASSERT_ANY(assert, msg)
#else
#define SHR_ASSERT(assert, msg) call shr_assert(assert, msg)
#define SHR_ASSERT_ALL(assert, msg) call shr_assert_all(assert, msg)
#define SHR_ASSERT_ANY(assert, msg) call shr_assert_any(assert, msg)
#endif
use shr_assert_mod
