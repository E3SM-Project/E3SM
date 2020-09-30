#include "scream_config.hpp"
#include "ekat/util/ekat_arch.hpp"
#include "scream_types.hpp"

namespace scream {

std::string scream_config_string() {
  auto config = ekat::ekat_config_string();

  config += "sizeof(Real) = " + std::to_string(sizeof(Real)) + "\n";
  config += "default pack size = " + std::to_string(SCREAM_PACK_SIZE) + "\n";

  return config;
}

} // namespace scream
