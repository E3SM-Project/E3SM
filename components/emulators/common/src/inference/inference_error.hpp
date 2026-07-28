/**
 * @file inference_error.hpp
 * @brief The exception type used throughout the inference layer.
 */

#ifndef E3SM_EMULATOR_INFERENCE_ERROR_HPP
#define E3SM_EMULATOR_INFERENCE_ERROR_HPP

#include <sstream>
#include <stdexcept>
#include <string>

namespace emulator {
namespace inference {

/// @brief Any failure raised by the inference layer.
class InferenceError : public std::runtime_error {
public:
  explicit InferenceError(const std::string &what)
      : std::runtime_error("[emulator::inference] " + what) {}
};

} // namespace inference
} // namespace emulator

/**
 * @brief Throw an InferenceError unless a condition holds.
 *
 * The message is streamed, so it can carry shapes and names:
 *   EMULATOR_INFER_REQUIRE(t.rank() == 2, "Expected a 2-D tensor, got "
 *                                             << t.to_string() << ".");
 */
#define EMULATOR_INFER_REQUIRE(cond, msg)                                      \
  do {                                                                         \
    if (!(cond)) {                                                             \
      std::ostringstream _emulator_infer_oss;                                  \
      _emulator_infer_oss << msg;                                              \
      throw ::emulator::inference::InferenceError(_emulator_infer_oss.str());  \
    }                                                                          \
  } while (false)

#endif // E3SM_EMULATOR_INFERENCE_ERROR_HPP
