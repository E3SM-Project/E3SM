#ifdef NDEBUG
#define SHR_ASSERT(assert, msg)
#define SHR_ASSERT_FL(assert, file, line)
#define SHR_ASSERT_MFL(assert, msg, file, line)
#define SHR_ASSERT_ALL(assert, msg)
#define SHR_ASSERT_ALL_FL(assert, file, line)
#define SHR_ASSERT_ALL_MFL(assert, msg, file, line)
#define SHR_ASSERT_ANY(assert, msg)
#define SHR_ASSERT_ANY_FL(assert, file, line)
#define SHR_ASSERT_ANY_MFL(assert, msg, file, line)
#else
#define SHR_ASSERT(assert, my_msg) call shr_assert(assert, msg=my_msg)
#define SHR_ASSERT_FL(assert, my_file, my_line) call shr_assert(assert, file=my_file, line=my_line)
#define SHR_ASSERT_MFL(assert, my_msg, my_file, my_line) call shr_assert(assert, msg=my_msg, file=my_file, line=my_line)
#define SHR_ASSERT_ALL(assert, my_msg) call shr_assert_all(assert, msg=my_msg)
#define SHR_ASSERT_ALL_FL(assert, my_file, my_line) call shr_assert_all(assert, file=my_file, line=my_line)
#define SHR_ASSERT_ALL_MFL(assert, my_msg, my_file, my_line) call shr_assert_all(assert, msg=my_msg, file=my_file, line=my_line)
#define SHR_ASSERT_ANY(assert, my_msg) call shr_assert_any(assert, msg=my_msg)
#define SHR_ASSERT_ANY_FL(assert, my_file, my_line) call shr_assert_any(assert, file=my_file, line=my_line)
#define SHR_ASSERT_ANY_MFL(assert, my_msg, my_file, my_line) call shr_assert_any(assert, msg=my_msg, file=my_file, line=my_line)
#endif
use shr_assert_mod
