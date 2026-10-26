#ifndef E3SM_EMULATOR_CONFIG_HPP
#define E3SM_EMULATOR_CONFIG_HPP

#include "inference_backend.hpp"
#include <field_registry.hpp>

#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>

namespace e3sm::emulator {

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
 *  - backend_type
 *  - input_layout
 *  - output_layout
 */
struct ModelConfig {
  std::filesystem::path path;
  inference::BackendType backend_type = inference::BackendType::STUB;
  TensorLayout input_layout = TensorLayout::SampleMajor;
  TensorLayout output_layout = TensorLayout::SampleMajor;
  std::map<std::string, std::string> options;
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

  std::vector<std::string> input_order;
  std::vector<std::string> output_order;

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
