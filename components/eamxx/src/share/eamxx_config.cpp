#include "eamxx_config.hpp"
#include "eamxx_session.hpp"
#include "eamxx_types.hpp"

#include "ekat/util/ekat_arch.hpp"
#include "ekat/ekat_assert.hpp"

namespace scream {

std::string eamxx_config_string() {
  std::string config = "\n-------- EKAT CONFIGS --------\n\n";
  config += ekat::ekat_config_string();
  config += "\n-------- SCREAM CONFIGS --------\n\n";
  config += " sizeof(Real) = " + std::to_string(sizeof(Real)) + "\n";
  config += " default pack size = " + std::to_string(SCREAM_PACK_SIZE) + "\n";
  config += " default FPE mask: " +
      ( get_default_fpes() == 0 ? "0 (NONE) \n" :
        std::to_string(get_default_fpes()) + " (FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW) \n");
  config += "-------------------------------\n";

  return config;
}

bool& use_leap_year_impl () {
#ifdef SCREAM_HAS_LEAP_YEAR
  static bool use_leap = true;
#else
  static bool use_leap = false;
#endif
  return use_leap;
}

bool use_leap_year () {
  return use_leap_year_impl ();
}

void set_use_leap_year (const bool use_leap) {
  use_leap_year_impl () = use_leap;
}

} // namespace scream
