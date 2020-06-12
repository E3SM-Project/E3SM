#include "ekat/scream_session.hpp"
#include "ekat/scream_kokkos.hpp"
#include "ekat/scream_assert.hpp"
#include "ekat/util/scream_arch.hpp"

namespace scream {

void initialize_scream_session () {
  enable_default_fpes ();
  Kokkos::initialize();
  std::cout << util::config_string() << "\n";
}

void initialize_scream_session (int argc, char **argv) {
  enable_default_fpes ();
  Kokkos::initialize(argc, argv);
  std::cout << util::config_string() << "\n";
}

void finalize_scream_session () {
  Kokkos::finalize();
}

} // namespace scream
