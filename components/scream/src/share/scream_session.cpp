#include "scream_session.hpp"
#include "scream_config.hpp"

#include "ekat/ekat_session.hpp"

#include <iostream>

namespace scream {

void initialize_scream_session (bool print_config) {
  ekat::initialize_ekat_session(false);

  if (print_config) 
    std::cout << scream_config_string() << "\n";
}

void initialize_scream_session (int argc, char **argv, bool print_config) {
  ekat::initialize_ekat_session(argc,argv,false);
  if (print_config) 
    std::cout << scream_config_string() << "\n";
}

extern "C" {
void finalize_scream_session () {
  ekat::finalize_ekat_session();
}
} // extern "C"

} // namespace scream
