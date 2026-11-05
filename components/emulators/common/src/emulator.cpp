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
      .component = m_name,
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
      std::ranges::find(m_grid.dimensions, name, &DimensionDesc::name);
  if (it == m_grid.dimensions.end()) {
    throw std::runtime_error("Requested invalid dimension name: " +
                             std::string(name));
  }
  return *it;
}

const EmulatorFieldDesc&
Emulator::get_field_descriptor(std::string_view name) const {
  const auto it =
      std::ranges::find(field_descriptors_, name,
                        [](const EmulatorFieldDesc& field) -> std::string_view {
                          return field.config.name;
                        });
  if (it == field_descriptors_.end()) {
    throw std::runtime_error("Requested field not found: " + std::string(name));
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

std::size_t Emulator::local_field_size(const FieldConfig& field) const {
  std::size_t size = 1;

  if (field.dimensions.empty()) {
    throw std::runtime_error("Field '" + field.name +
                             "' does not declare any dimensions");
  }

  for (const auto& dim_name : field.dimensions) {
    const auto& desc = dimension(dim_name);
    if (desc.local_extent == 0) {
      throw std::runtime_error(
          "Dimension '" + dim_name +
          "' has zero local extent for field: " + field.name);
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
      throw std::runtime_error("Field " + field_name + " Not in fields");
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
    : m_config_path("temp"), m_grid(grid), m_emulator_config(emu_cfg),
      m_component_id(grid.id), m_moab_app_id(-1), m_initialized(false),
      m_step_count(0) {

  if (grid.use_moab) {
    m_moab.emplace(emu_cfg.name, m_component_id,
                   MoabPartition{
                       .comm = m_grid.comm,
                       .global_ids = std::span<int>{m_grid.col_gids},
                       .latitudes = std::span<double>{m_grid.lat},
                       .longitudes = std::span<double>{m_grid.lon},
                   });
  }
  initialize();
}

void Emulator::apply_config(const EmulatorConfig& config) {
  m_name = config.name;
  if (config.type == "process") {
    m_type = EmulatorType::Process;
  } else {
    m_type = EmulatorType::Component;
  }

  const auto input_size =
      append_fields(config.input_order, config.inputs, ModelFieldRole::Input);
  const auto output_size = append_fields(config.output_order, config.outputs,
                                         ModelFieldRole::Output);

  // for (const auto& field : field_descriptors()) {
  //   std::cout << to_string(field) << std::endl;
  // }
  input_fields_.resize(input_size);
  output_fields_.resize(output_size);

  inference::InferenceConfig infer_config;
  infer_config.model_path = config.model.path;
  infer_config.options = config.model.options;
  infer_config.input_channels = config.inputs.size();
  infer_config.output_channels = config.outputs.size();

  backend_ = inference::create_backend(m_emulator_config.model.backend_type,
                                       infer_config);
}

void Emulator::initialize() {
  if (m_initialized) {
    throw std::runtime_error("Emulator already initialized");
  }
  apply_config(m_emulator_config);
  m_coupling.set_local_gids(m_grid.col_gids);
  if (m_moab) {
    register_with_moab();
  }
  m_initialized = true;
}

void Emulator::populate_registry(FieldRegistry& registry) {

  if (!m_initialized) {
    throw std::logic_error(
        "Cannot populate registry before emulator initialization");
  }

  for (const auto& field : field_descriptors_) {
    m_coupling.register_coupled_field(registry, make_registered_field(field),
                                      field_view(field));
  }
}

void Emulator::run(const ProcessRunOpts opts) {
  if (!m_initialized) {
    throw std::runtime_error("Emulator::run() called before initialize()");
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
  backend_->infer(inputs, outputs);

  m_step_count++;
}

void Emulator::finalize() {
  if (!m_initialized) {
    return; // Already finalized or never initialized
  }
  backend_->finalize();

  if (m_moab) {
    m_moab.reset();
  }
  m_initialized = false;
}

MoabTagDescriptor Emulator::make_moab_tag(const EmulatorFieldDesc& desc) const {
  if (desc.config.dimensions.empty()) {
    throw std::runtime_error("No dimensions for field " + desc.config.name);
  }
  if (m_grid.num_local_columns == 0) {
    throw std::runtime_error("No Local columns for this emulator!!");
  }
  return MoabTagDescriptor{
      .name = desc.config.name,
      .role = (desc.role == ModelFieldRole::Input) ? FieldRole::Import
                                                   : FieldRole::Export,
      .values_per_vertex = desc.size / m_grid.num_local_columns,
  };
}

static std::string to_string(EmulatorType t) {
  switch (t) {
  case EmulatorType::ATM:
    return "ATM";
  case EmulatorType::OCN:
    return "OCN";
  case EmulatorType::ICE:
    return "ICE";
  case EmulatorType::LND:
    return "LND";
  default:
    return "UNKNOWN";
  }
}

void Emulator::print_info(std::ostream& os) const {
  os << "Emulator '" << m_name << "'\n";
  os << "  type          : " << to_string(m_type) << "\n";
  os << "  id            : " << m_component_id << "\n";
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
    m_moab->define_tag(make_moab_tag(field));
  }
}

} // namespace e3sm::emulator
