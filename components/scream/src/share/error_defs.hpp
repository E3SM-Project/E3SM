#ifndef SCREAM_ERROR_DEFS_HPP
#define SCREAM_ERROR_DEFS_HPP

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
#endif // NDEBUG

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

#endif // SCREAM_ERROR_DEFS_HPP
