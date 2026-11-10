/**
 * @file inference_config.cpp
 * @brief InferenceConfig implementation.
 */

#include "inference_config.hpp"

#include "inference_error.hpp"

#include <cstddef>
#include <stdexcept>

namespace emulator {
namespace inference {

std::string InferenceConfig::get(const std::string &key,
                                 const std::string &fallback) const {
  const auto it = options.find(key);
  return it == options.end() ? fallback : it->second;
}

int InferenceConfig::get_int(const std::string &key, int fallback) const {
  const auto it = options.find(key);
  if (it == options.end()) {
    return fallback;
  }
  const std::string &value = it->second;
  std::size_t parsed = 0;
  int result = 0;
  try {
    result = std::stoi(value, &parsed);
  } catch (const std::logic_error &) {
    parsed = 0;
  }
  EMULATOR_INFER_REQUIRE(!value.empty() && parsed == value.size(),
                         "Option '" << key << "' must be an integer, got '"
                                    << value << "'.");
  return result;
}

bool InferenceConfig::get_bool(const std::string &key, bool fallback) const {
  const auto it = options.find(key);
  if (it == options.end()) {
    return fallback;
  }
  const std::string &value = it->second;
  if (value == "true" || value == "1") {
    return true;
  }
  EMULATOR_INFER_REQUIRE(value == "false" || value == "0",
                         "Option '" << key << "' must be true or false, got '"
                                    << value << "'.");
  return false;
}

} // namespace inference
} // namespace emulator
