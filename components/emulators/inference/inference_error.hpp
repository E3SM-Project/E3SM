/**
 * @file inference_error.hpp
 * @brief Exception type used throughout the inference layer.
 */

#ifndef E3SM_EMULATOR_INFERENCE_ERROR_HPP
#define E3SM_EMULATOR_INFERENCE_ERROR_HPP

#include <sstream>
#include <stdexcept>
#include <string>

namespace e3sm::emulator::inference {

/**
 * @brief Any failure raised by the inference layer.
 */
class InferenceError : public std::runtime_error {
public:
  explicit InferenceError(const std::string &what)
      : std::runtime_error("[emulator::inference] " + what) {}
};

} // namespace e3sm::emulator::inference

/**
 * @brief Throw an InferenceError unless a condition holds.
 *
 * The message is streamed, so it can carry names and shapes.
 */
#define EMULATOR_INFER_REQUIRE(cond, msg)                                      \
  do {                                                                         \
    if (!(cond)) {                                                             \
      std::ostringstream _emulator_infer_oss;                                  \
      _emulator_infer_oss << msg;                                              \
      throw ::e3sm::emulator::inference::InferenceError(_emulator_infer_oss.str());  \
    }                                                                          \
  } while (false)

#endif // E3SM_EMULATOR_INFERENCE_ERROR_HPP
