#include "eamxx_config.hpp"
#include "eamxx_session.hpp"
#include "eamxx_types.hpp"
#include "eamxx_version.h"

#include <ekat_kokkos_session.hpp>
#include <ekat_arch.hpp>
#include <ekat_fpe.hpp>
#include <ekat_assert.hpp>

namespace scream {

std::string eamxx_version () {
  return EAMXX_VERSION;
}

std::string eamxx_git_version () {
  return EAMXX_GIT_VERSION;
}

std::string eamxx_config_string() {
  std::string config = "\n-------- EKAT CONFIGS --------\n\n";
  config += "Active AVX settings: " + ekat::active_avx_string () + "\n";
  config += "Compiler Id: " + ekat::compiler_id_string () + "\n";
  config += ekat::kokkos_config_string();
  config += "\n-------- EAMXX CONFIGS --------\n\n";
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
