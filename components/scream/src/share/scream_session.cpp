#include "scream_session.hpp"
#include "scream_config.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/ekat_session.hpp"

#include <iostream>
#include <cfenv>

namespace scream {

int get_default_fpes () {
#ifdef SCREAM_FPE
  return (FE_DIVBYZERO |
          FE_INVALID   |
          FE_OVERFLOW);
#else
  return 0;
#endif
}

void initialize_scream_session (bool print_config) {
  ekat::initialize_ekat_session(print_config);

  // Make sure scream only has its FPEs
  ekat::disable_all_fpes();
  ekat::enable_fpes(get_default_fpes());

  if (print_config) 
    std::cout << scream_config_string() << "\n";
}

void initialize_scream_session (int argc, char **argv, bool print_config) {
  ekat::initialize_ekat_session(argc,argv,print_config);

  // Make sure scream only has its FPEs
  ekat::disable_all_fpes();
  ekat::enable_fpes(get_default_fpes());

  if (print_config) 
    std::cout << scream_config_string() << "\n";
}

extern "C" {
void finalize_scream_session () {
  ekat::finalize_ekat_session();
}
} // extern "C"

} // namespace scream
