#ifndef E3SM_EMULATOR_CONFIG_HPP
#define E3SM_EMULATOR_CONFIG_HPP

#include <cstddef>
#include <field_registry.hpp>
#include <filesystem>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace emulator {

// [n_samples, n_features] or [n_features, n_samples]
enum class TensorLayout { SampleMajor, FeatureMajor };

/**
 * @brief Enumeration of emulator types in E3SM.
 */
enum class EmulatorType {
  ATM = 0,     ///< Atmosphere component emulator
  OCN = 1,     ///< Ocean component emulator
  ICE = 2,     ///< Sea ice component emulator
  LND = 3,     ///< Land component emulator
  Process = 4, ///< Process emulator
  Component = 5
};

enum class ModelFieldRole {
  Input = 0,
  Output = 1,
};

/**
 * @brief: Struct to hold information relevant to NN to be used
 * @fields:
 *  - path
 *  - input_layout
 *  - output_layout
 */
struct ModelConfig {
  std::filesystem::path path;
  TensorLayout input_layout = TensorLayout::SampleMajor;
  TensorLayout output_layout = TensorLayout::SampleMajor;
};

/**
 * @brief: Struct to hold information needed to allocate an emulator's
 * input/output field 
 *
 * @fields:
   - name
   - e3sm::coupler::RegisteredFieldAttributes attributes;
   - dimensions
 */
struct FieldConfig {
  // NOTE: maybe we will want some scaling or packing information
  std::string name;
  e3sm::coupler::RegisteredFieldAttributes attributes;
  std::vector<std::string> dimensions;
};

/**
 * @brief Struct to hold emulator configuration parsed from yaml file
 * @fields:
 *  - name
 *  - type
 *  - model
 *  - dimensions
 *  - inputs
 *  - outputs
 **/
struct EmulatorConfig {
  std::string name;
  std::string type;

  ModelConfig model;
  std::vector<std::string> dimensions;

  std::vector<FieldConfig> inputs;
  std::vector<FieldConfig> outputs;
};

EmulatorConfig read_emulator_config(const std::filesystem::path& cfg);

std::string to_string(const FieldConfig& field_cfg);
std::string to_string(const ModelConfig& model_cfg);
std::string to_string(const EmulatorConfig& cfg);

inline std::ostream& operator<<(std::ostream& os, const EmulatorConfig& cfg) {
  return os << to_string(cfg) << '\n';
}

} // namespace emulator

#endif
