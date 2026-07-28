/**
 * @file inference_config.hpp
 * @brief Settings for an inference backend.
 */

#ifndef E3SM_EMULATOR_INFERENCE_CONFIG_HPP
#define E3SM_EMULATOR_INFERENCE_CONFIG_HPP

#include <map>
#include <string>
#include <vector>

namespace emulator {
namespace inference {

/**
 * @brief Configuration for an inference backend.
 *
 * A handful of named fields that the C++ layer itself acts on, plus a
 * free-form option map that is passed through to the backend (and, for the
 * python backend, straight on to the model).  Anything the C++ layer does
 * not recognise lands in `options` rather than being an error, so a Python
 * emulator can grow settings without touching this file.
 */
struct InferenceConfig {
  /// @brief Which backend to build: "stub" or "python".
  std::string backend = "stub";

  /// @brief Path to the model / checkpoint, handed to the backend verbatim.
  std::string model_path;

  /// @brief Ordered names of the tensors the model consumes.
  std::vector<std::string> inputs;

  /// @brief Ordered names of the tensors the model produces.
  std::vector<std::string> outputs;

  /// @brief Features per column, used only by the flat-array infer() overload.
  int input_channels = 0;

  /// @brief Features per column, used only by the flat-array infer() overload.
  int output_channels = 0;

  /// @brief Chatty initialization and per-step reporting.
  bool verbose = false;

  /// @brief Everything else, passed through to the backend.
  std::map<std::string, std::string> options;

  /**
   * @brief Apply one `key: value` setting.
   *
   * Recognised keys set the named fields above; `input`/`output` append to
   * the ordered lists.  Any other key is stored in `options`, with a leading
   * `option.` stripped if present.
   */
  void apply(const std::string &key, const std::string &value);

  /// @brief Set an option directly.
  void set(const std::string &key, const std::string &value) {
    options[key] = value;
  }

  /// @brief Option lookup with a fallback.
  std::string get(const std::string &key,
                  const std::string &fallback = "") const;
  int get_int(const std::string &key, int fallback) const;
  bool get_bool(const std::string &key, bool fallback) const;

  /**
   * @brief Read settings from a line-oriented `key: value` file.
   *
   * Deliberately the same shape as the parsing the emulator components
   * already do on their `*_in` namelists, so the inference settings can live
   * in `atm_in` next to everything else.  `#` and `!` start a comment.
   *
   * @param path   File to read; a missing file yields default settings.
   * @param prefix Only keys starting with this are used, and the prefix is
   *               stripped before `apply`.  Pass "inference." to read out of
   *               a component namelist; pass "" to read a dedicated file.
   */
  static InferenceConfig from_file(const std::string &path,
                                   const std::string &prefix = "");
};

} // namespace inference
} // namespace emulator

#endif // E3SM_EMULATOR_INFERENCE_CONFIG_HPP
