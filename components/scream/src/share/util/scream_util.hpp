#ifndef INCLUDE_SCREAM_UTIL
#define INCLUDE_SCREAM_UTIL

#include <cstdio>
#include <sstream>
#include <memory>
#include <exception>

#ifndef NDEBUG
#define scream_assert(condition) do {                                   \
    if ( ! (condition)) {                                               \
      std::stringstream _ss_;                                           \
      _ss_ << __FILE__ << ":" << __LINE__ << ": FAIL:\n" << #condition  \
        << "\n";                                                        \
        throw std::logic_error(_ss_.str());                             \
    }                                                                   \
  } while (0)
#define scream_kernel_assert(condition) do {    \
    if ( ! (condition))                         \
      Kokkos::abort(#condition);                \
  } while (0)
#else
#define scream_assert(condition)
#define scream_kernel_assert(condition)
#endif
#define scream_throw_if(condition, message) do {                        \
    if (condition) {                                                    \
      std::stringstream _ss_;                                           \
      _ss_ << __FILE__ << ":" << __LINE__ << ": The condition:\n"       \
           << #condition "\nled to the exception\n" << message << "\n"; \
        throw std::logic_error(_ss_.str());                             \
    }                                                                   \
  } while (0)
#define scream_kernel_throw_if(condition, message) do {             \
    if (condition)                                                  \
      Kokkos::abort(#condition " led to the exception\n" message);  \
  } while (0)

#if defined __INTEL_COMPILER
# define vector_ivdep _Pragma("ivdep")
# ifdef _OPENMP
#  define vector_simd _Pragma("omp simd")
# else
#  define vector_simd _Pragma("simd")
# endif
#elif defined __GNUG__
# define vector_ivdep _Pragma("GCC ivdep")
# define vector_simd _Pragma("GCC ivdep")
# define restrict __restrict__
#else
# define vector_ivdep
# define vector_simd
# define restrict
#endif

namespace scream {
namespace util {

struct FILECloser { void operator() (FILE* fh) { fclose(fh); } };
using FILEPtr = std::unique_ptr<FILE, FILECloser>;

template<typename T>
void write (const T* v, size_t sz, const FILEPtr& fid) {
  size_t nwrite = fwrite(v, sizeof(T), sz, fid.get());
  scream_throw_if(nwrite != sz, "write: nwrite = " << nwrite << " sz = " << sz);
}

template<typename T>
void read (T* v, size_t sz, const FILEPtr& fid) {
  size_t nread = fread(v, sizeof(T), sz, fid.get());
  scream_throw_if(nread != sz, "read: nread = " << nread << " sz = " << sz);
}

template <typename Real> struct is_single_precision {};
template <> struct is_single_precision<float> { enum : bool { value = true }; };
template <> struct is_single_precision<double> { enum : bool { value = false }; };

bool eq(const std::string& a, const char* const b1, const char* const b2 = 0);

void activate_floating_point_exceptions_if_enabled();

std::string active_avx_string();
std::string config_string();

template <typename Real>
Real reldif (const Real& a, const Real& b) {
  return std::abs(b - a)/std::abs(a);
}

} // namespace util
} // namespace scream

#endif
