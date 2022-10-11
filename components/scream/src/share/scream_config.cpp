#include "scream_config.hpp"
#include "scream_session.hpp"
#include "scream_types.hpp"

#include "ekat/util/ekat_arch.hpp"
#include "ekat/ekat_assert.hpp"

namespace scream {

std::string scream_config_string() {
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

} // namespace scream
