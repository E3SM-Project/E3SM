/**
 * @file inference_config.hpp
 * @brief Configuration for inference backends.
 */

#ifndef E3SM_EMULATOR_INFERENCE_CONFIG_HPP
#define E3SM_EMULATOR_INFERENCE_CONFIG_HPP

#include <map>
#include <string>

namespace emulator {
namespace inference {

/**
 * @brief Configuration for inference backends.
 *
 * Settings common to all backends are named fields; backend-specific ones
 * go in `options` (each backend documents the keys it reads).
 */
struct InferenceConfig {
  std::string model_path;  ///< Path to the model file
  int input_channels = 0;  ///< Number of input features per grid point
  int output_channels = 0; ///< Number of output features per grid point
  bool verbose = false;    ///< Enable verbose output (for debugging)

  std::map<std::string, std::string> options; ///< Backend-specific settings

  void set(const std::string &key, const std::string &value) {
    options[key] = value;
  }

  /// @brief Option lookup; returns the fallback if the key is absent.
  std::string get(const std::string &key,
                  const std::string &fallback = "") const;

  /// @throws InferenceError if the value is not an integer
  int get_int(const std::string &key, int fallback) const;

  /// @throws InferenceError if the value is not true/false/1/0
  bool get_bool(const std::string &key, bool fallback) const;
};

} // namespace inference
} // namespace emulator

#endif // E3SM_EMULATOR_INFERENCE_CONFIG_HPP
