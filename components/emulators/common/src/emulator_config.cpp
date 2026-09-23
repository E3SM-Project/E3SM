#include <ekat_yaml.hpp>
#include <emulator_config.hpp>
#include <field_registry.hpp>
#include <filesystem>
#include <sstream>
#include <stdexcept>
#include <string>

namespace emulator {
namespace {

constexpr std::string TABS = "  ";
EmulatorType get_type_from_string(std::string type, std::string name) {
  if (name == "process") {
    return EmulatorType::Process;
  } else {
    return EmulatorType::Component;
  }
}

TensorLayout get_layout_from_string(std::string layout) {
  if (layout == "SampleMajor") {
    return TensorLayout::SampleMajor;
  } else if (layout == "FeatureMajor") {
    return TensorLayout::FeatureMajor;
  }
  throw std::runtime_error("Unknown Layout requested: " + layout);
}

std::string to_string(TensorLayout layout) {
  if (layout == TensorLayout::FeatureMajor) {
    return "FeatureMajor";
  } else if (layout == TensorLayout::SampleMajor) {
    return "SampleMajor";
  } else {
    throw std::runtime_error("Unknown Layout");
  }
}

std::vector<FieldConfig> read_fields(const ekat::ParameterList& params,
                                     std::string label) {
  const auto& input_list = params.sublist(label);
  const auto& field_names = input_list.sublist_names();
  std::vector<FieldConfig> fields;
  fields.reserve(field_names.size());

  for (const auto& field_name : field_names) {
    const auto& field_list = input_list.sublist(field_name);
    const auto& attr_list = field_list.sublist("attributes");
    const auto attributes = e3sm::coupler::read_attributes(attr_list);
    auto dimensions = field_list.get<std::vector<std::string>>("dimensions");

    fields.emplace_back(FieldConfig{.name = field_name,
                                    .attributes = std::move(attributes),
                                    .dimensions = std::move(dimensions)});
  }

  return fields;
}
} // namespace

/**
 * @brief Parses a emulator_config.yaml file into
 */
EmulatorConfig read_emulator_config(const std::filesystem::path& filename) {
  if (!std::filesystem::exists(filename)) {
    throw std::runtime_error("Coupling configuration file not found: " +
                             filename.string());
  }

  const ekat::ParameterList params = ekat::parse_yaml_file(filename.string());
  EmulatorConfig cfg;

  const auto& attrs = params.sublist("emulator");
  cfg.name = attrs.get<std::string>("name");
  cfg.type = attrs.get<std::string>("type");

  // Get Model network configuration
  const auto& model_attrs = attrs.sublist("model");
  cfg.model.path = model_attrs.get<std::string>("path");
  cfg.model.output_layout =
      get_layout_from_string(model_attrs.get<std::string>("output_layout"));
  cfg.model.input_layout =
      get_layout_from_string(model_attrs.get<std::string>("input_layout"));

  // Get Sample Space Information
  const auto& sample_info = params.sublist("sample_space");
  cfg.dimensions = sample_info.get<std::vector<std::string>>("dimensions");

  // Get inputs
  cfg.inputs = read_fields(params, "inputs");
  // Get outputs
  cfg.outputs = read_fields(params, "outputs");
  return cfg;
}

std::string to_string(const ModelConfig& cfg) {
  std::ostringstream out_str;
  out_str << "| Model Config: \n";
  out_str << '|' << TABS << "path: " << cfg.path << '\n';
  out_str << '|' << TABS << "Input Layout: " << to_string(cfg.input_layout)
          << ", Output Layout: " << to_string(cfg.output_layout) << '\n';

  return out_str.str();
}

std::string to_string(const FieldConfig& cfg) {
  std::ostringstream out_str;
  out_str << '|' << TABS << cfg.name << "[ ";
  for (auto& dim : cfg.dimensions) {
    out_str << dim << ", ";
  }
  out_str << "]\n";
  out_str << '|' + TABS + "attributes:\n"
          << to_string(cfg.attributes, '|' + TABS + TABS);
  return out_str.str();
}

std::string to_string(const EmulatorConfig& cfg) {
  std::ostringstream out_string;
  out_string << "| Emulator Configuration for: " << cfg.name << '\n';
  out_string << "| type: " << cfg.type << '\n';
  out_string << to_string(cfg.model);
  out_string << "| Sample Space dimensions: [ ";
  for (const auto& dim : cfg.dimensions) {
    out_string << dim << ", ";
  }
  out_string << "]\n";

  out_string << "| Inputs: \n";
  for (const auto& field : cfg.inputs) {
    out_string << to_string(field);
  }
  out_string << "| Outputs: \n";
  for (const auto& field : cfg.outputs) {
    out_string << to_string(field);
  }

  return out_string.str();
}
} // namespace emulator
