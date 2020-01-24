#include "share/scream_session.hpp"
#include "share/util/scream_arch.hpp"

#ifdef SCREAM_FPE
# include <xmmintrin.h>
#endif

#include <Kokkos_Core.hpp>

namespace scream {

void initialize_scream_session () {
  util::activate_floating_point_exceptions_if_enabled();
  Kokkos::initialize();
  std::cout << util::config_string() << "\n";
}

void initialize_scream_session (int argc, char **argv) {
  util::activate_floating_point_exceptions_if_enabled();
  Kokkos::initialize(argc, argv);
  std::cout << util::config_string() << "\n";
}

void finalize_scream_session () {
  Kokkos::finalize();
}

} // namespace scream
