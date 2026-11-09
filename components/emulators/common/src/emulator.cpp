/**
 * @file emulator.cpp
 * @brief Implementation of the Emulator base class.
 */

#include "create_inference_backend.hpp"
#include "inference_backend.hpp"
#include "inference_config.hpp"
#include "moab_adapter.hpp"
#include <coupler_driver.hpp>
#include <emulator.hpp>
#include <emulator_config.hpp>
#include <field_registry.hpp>
#include <memory>
#include <process.hpp>

#include <ekat_yaml.hpp>
#include <moab/Types.hpp>
#include <moab/iMOAB.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string_view>

namespace e3sm::emulator {

RegisteredField
Emulator::make_registered_field(const EmulatorFieldDesc& field) {
  const auto role = field.role == ModelFieldRole::Input ? FieldRole::Import
                                                        : FieldRole::Export;
  return RegisteredField{
      .role = role,
      .component = name_,
      .name = field.config.name,
      .size = field.size,
  };
}

std::string to_string(const DimensionDesc& desc) {
  std::ostringstream out;
  out << "Dimension: " << desc.name
      << " Local/Global sizes: " << desc.local_extent << ", "
      << desc.global_extent << "\n";
  return out.str();
}

const DimensionDesc& Emulator::dimension(std::string_view name) const {
  const auto it =
      std::ranges::find(grid_.dimensions, name, &DimensionDesc::name);
  if (it == grid_.dimensions.end()) {
    std::string known;
    for (const auto& dim : grid_.dimensions) {
      known += (known.empty() ? "" : ", ") + dim.name;
    }
    throw std::runtime_error("Emulator '" + name_ + "': dimension '" +
                             std::string(name) +
                             "' is not defined by its GridDesc. Defined "
                             "dimensions: " +
                             (known.empty() ? "none" : known));
  }
  return *it;
}

const EmulatorFieldDesc&
Emulator::get_field_descriptor(std::string_view name) const {

  const auto has_name = [name](const EmulatorFieldDesc& field) {
    return field.config.name == name;
  };

  if (std::ranges::count_if(field_descriptors_, has_name) > 1) {
    throw std::runtime_error("Emulator '" + name_ + "': field '" +
                             std::string(name) +
                             "' is both an input and an output; specify a "
                             "ModelFieldRole");
  }
  const auto it = std::ranges::find_if(field_descriptors_, has_name);
  if (it == field_descriptors_.end()) {
    throw std::runtime_error("Emulator '" + name_ + "': no field named '" +
                             std::string(name) + "'");
  }
  return *it;
}

const EmulatorFieldDesc&
Emulator::get_field_descriptor(std::string_view name,
                               ModelFieldRole role) const {

  const auto it = std::ranges::find_if(
      field_descriptors_, [name, role](const EmulatorFieldDesc& field) {
        return field.role == role && field.config.name == name;
      });
  if (it == field_descriptors_.end()) {
    throw std::runtime_error(
        "Emulator '" + name_ + "': no " +
        (role == ModelFieldRole::Input ? "input" : "output") +
        " field named '" + std::string(name) + "'");
  }
  return *it;
}

std::span<const EmulatorFieldDesc> Emulator::field_descriptors() const {
  return field_descriptors_;
}

std::span<const double> Emulator::get_field(std::string_view name) {
  const auto& desc = get_field_descriptor(name);
  return field_view(desc);
}

std::span<const double> Emulator::get_field(std::string_view name,
                                            ModelFieldRole role) {
  const auto& desc = get_field_descriptor(name, role);
  return field_view(desc);
}

std::size_t Emulator::local_field_size(const FieldConfig& field) const {
  std::size_t size = 1;

  if (field.dimensions.empty()) {
    throw std::runtime_error("Emulator '" + name_ + "': field '" + field.name +
                             "' does not declare any dimensions");
  }

  for (const auto& dim_name : field.dimensions) {
    const auto& desc = dimension(dim_name);
    if (desc.local_extent == 0) {
      throw std::runtime_error(
          "Emulator '" + name_ + "': dimension '" + dim_name +
          "' of field '" + field.name + "' has zero local extent");
    }
    if (size > std::numeric_limits<std::size_t>::max() / desc.local_extent) {
      throw std::runtime_error("Inappropriate dimensions for field:\n" +
                               to_string(field) + to_string(desc));
    }
    size *= desc.local_extent;
  }
  return size;
}

std::size_t Emulator::append_fields(std::span<const std::string> order,
                                    std::span<const FieldConfig> fields,
                                    ModelFieldRole role) {

  std::size_t next_offset = 0;

  for (const std::string& field_name : order) {
    auto field = std::ranges::find(fields, field_name, &FieldConfig::name);
    if (field == fields.end()) {
      const auto* order_key = role == ModelFieldRole::Input ? "input" : "output";
      throw std::runtime_error("Emulator '" + name_ + "': field '" + field_name +
                               "' is listed in " + order_key +
                               "_order but not defined in " + order_key + "s");
    }
    const auto field_size = local_field_size(*field);

    if (next_offset > std::numeric_limits<std::size_t>::max() - field_size) {
      throw std::overflow_error("Emulator field-buffer size overflow");
    }

    field_descriptors_.push_back(EmulatorFieldDesc{
        .role = role,
        .config = *field,
        .offset = next_offset,
        .size = field_size,
    });

    next_offset += field_size;
  }

  // this should be total size needed for the flat buffer
  return next_offset;
}

std::span<double> Emulator::field_view(const EmulatorFieldDesc& field) {
  auto& buffer =
      field.role == ModelFieldRole::Input ? input_fields_ : output_fields_;

  return std::span<double>{buffer.data() + field.offset, field.size};
}

std::span<double> Emulator::field_view(std::string_view name) {
  const auto& desc = get_field_descriptor(name);
  return field_view(desc);
}

Emulator::Emulator(const GridDesc grid, const EmulatorConfig emu_cfg)
    : config_path_("temp"), grid_(grid), emulator_config_(emu_cfg),
      component_id_(grid.id), moab_app_id_(-1), m_initialized(false),
      m_step_count(0) {

  if (grid.use_moab) {
    moab_.emplace(emu_cfg.name, component_id_,
                   MoabPartition{
                       .comm = grid_.comm,
                       .global_ids = std::span<int>{grid_.col_gids},
                       .latitudes = std::span<double>{grid_.lat},
                       .longitudes = std::span<double>{grid_.lon},
                   });
  }
  initialize();
}

void Emulator::apply_config(const EmulatorConfig& config) {
  name_ = config.name;

  const auto input_size =
      append_fields(config.input_order, config.inputs, ModelFieldRole::Input);
  const auto output_size = append_fields(config.output_order, config.outputs,
                                         ModelFieldRole::Output);

  input_fields_.resize(input_size);
  output_fields_.resize(output_size);

  inference::InferenceConfig infer_config;
  infer_config.model_path = config.model.path;
  infer_config.options = config.model.options;
  infer_config.input_channels = config.inputs.size();
  infer_config.output_channels = config.outputs.size();

  try {
    backend_ = inference::create_backend(emulator_config_.model.backend_type,
                                         infer_config);
  } catch (const std::exception& e) {
    throw std::runtime_error("Emulator '" + name_ +
                             "': failed to create the inference backend: " +
                             e.what());
  }
}

void Emulator::initialize() {
  if (m_finalized) {
    throw std::runtime_error("Emulator '" + name_ +
                             "': cannot initialize after finalize()");
  }
  if (m_initialized) {
    throw std::runtime_error("Emulator '" + name_ + "': already initialized");
  }
  apply_config(emulator_config_);
  coupling_.set_local_gids(grid_.col_gids);
  if (moab_) {
    register_with_moab();
  }
  m_initialized = true;
}

void Emulator::populate_registry(FieldRegistry& registry) {

  if (!m_initialized) {
    throw std::logic_error("Emulator '" + name_ +
                           "': cannot populate the registry after finalize()");
  }

  for (const auto& field : field_descriptors_) {
    coupling_.register_coupled_field(registry, make_registered_field(field),
                                      field_view(field));
  }
}

void Emulator::run(const ProcessRunOpts opts) {
  if (!m_initialized) {
    throw std::runtime_error("Emulator '" + name_ +
                             "': run() called after finalize()");
  }

  // wrap input/output fields in tensor maps
  inference::TensorMap inputs = {}, outputs = {};
  for (const EmulatorFieldDesc& field : field_descriptors()) {
    auto values = field_view(field);
    std::vector<std::int64_t> dims;
    for (const auto& dim_name : field.config.dimensions) {
      const DimensionDesc& dim_desc = dimension(dim_name);
      dims.push_back(static_cast<std::int64_t>(dim_desc.local_extent));
    }
    if (field.role == ModelFieldRole::Input) {
      const double* data = values.data();
      inputs.wrap(field.config.name, data, dims);
    } else { // field_desc.role == ModelFieldRole::Output
      outputs.wrap(field.config.name, values.data(), dims);
    }
  }
  try {
    if (!backend_->infer(inputs, outputs)) {
      throw std::runtime_error("the backend reported a failure");
    }
  } catch (const std::exception& e) {
    throw std::runtime_error("Emulator '" + name_ + "': inference failed at step " +
                             std::to_string(m_step_count) + ": " + e.what());
  }

  m_step_count++;
}

void Emulator::finalize() {
  if (!m_initialized) {
    return; // Already finalized
  }
  m_initialized = false;
  m_finalized = true;
  backend_->finalize();

  if (moab_) {
    moab_.reset();
  }
}

MoabTagDescriptor Emulator::make_moab_tag(const EmulatorFieldDesc& desc) const {
  if (desc.config.dimensions.empty()) {
    throw std::runtime_error("Emulator '" + name_ + "': field '" +
                             desc.config.name + "' has no dimensions");
  }
  if (grid_.num_local_columns == 0) {
    throw std::runtime_error("Emulator '" + name_ +
                             "': GridDesc has no local columns");
  }
  return MoabTagDescriptor{
      .name = desc.config.name,
      .role = (desc.role == ModelFieldRole::Input) ? FieldRole::Import
                                                   : FieldRole::Export,
      .values_per_vertex = desc.size / grid_.num_local_columns,
  };
}

void Emulator::print_info(std::ostream& os) const {
  os << "Emulator '" << name_ << "'\n";
  os << "  id            : " << component_id_ << "\n";
  os << "  initialized   : " << std::boolalpha << m_initialized << "\n";
  os << "  step_count    : " << m_step_count << "\n";

  // Grid summary
  int nx = get_nx();
  int ny = get_ny();
  int nloc = get_num_local_cols();
  int nglob = get_num_global_cols();

  os << "  grid          : nx=" << nx << " ny=" << ny
     << " num_local_cols=" << nloc << " num_global_cols=" << nglob << "\n";
}
void Emulator::register_with_moab() {
  for (const auto& field : field_descriptors()) {
    moab_->define_tag(make_moab_tag(field));
  }
}

} // namespace e3sm::emulator
