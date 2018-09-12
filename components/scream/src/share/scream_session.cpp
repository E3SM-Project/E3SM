#include "scream_session.hpp"
#include "share/util/scream_utils.hpp"

#ifdef SCREAM_FPE
# include <xmmintrin.h>
#endif

#include <Kokkos_Core.hpp>

namespace scream {

void activate_floating_point_exceptions_if_enabled () {
#ifdef SCREAM_FPE
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() &
                         ~( _MM_MASK_INVALID |
                            _MM_MASK_DIV_ZERO |
                            _MM_MASK_OVERFLOW |
                            _MM_MASK_UNDERFLOW ));
#endif
}

void initialize (int argc, char **argv) {
  util::activate_floating_point_exceptions_if_enabled();
  Kokkos::initialize(argc, argv);
  std::cout << util::config_string() << "\n";
}

void finalize () {
  Kokkos::finalize();
}

} // namespace scream
