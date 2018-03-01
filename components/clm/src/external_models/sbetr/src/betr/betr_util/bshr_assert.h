#ifdef NDEBUG
#define SHR_ASSERT(assert, msg, bstatus)
#define SHR_ASSERT_ALL(assert, msg, bstatus)
#define SHR_ASSERT_ANY(assert, msg, bstatus)
#else
#define SHR_ASSERT(assert, msg, bstatus) call shr_assert(assert, msg, bstatus)
#define SHR_ASSERT_ALL(assert, msg, bstatus) call shr_assert_all(assert, msg, bstatus)
#define SHR_ASSERT_ANY(assert, msg, bstatus) call shr_assert_any(assert, msg, bstatus)
#endif
use bshr_assert_mod, only : shr_assert
use bshr_assert_mod, only : shr_assert_all
use bshr_assert_mod, only : shr_assert_any
