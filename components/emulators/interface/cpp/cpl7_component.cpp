#include "field_registry.hpp"
#include <cpl7_component.hpp>
#include <stdexcept>

namespace e3sm::emulator {

void CPL7Component::attach_fields(
    std::span<const ExternalField> requested_fields,
    std::span<double> import_storage, std::span<double> export_storage) {
  if (initialized_) {
    throw std::runtime_error("Re-attaching fields not allowed");
  }
  std::vector<BoundField> bindings;
  bindings.reserve(requested_fields.size());

  for (const auto& ext_field : requested_fields) {
    std::span<double> buffer =
        (ext_field.role == FieldRole::Import) ? import_storage : export_storage;

    if (ext_field.extent > buffer.size() - ext_field.offset) {
      throw std::runtime_error("Invalid slice requested for " +
                               ext_field.label);
    }

    bindings.push_back(BoundField{
        .field =
            RegisteredField{
                .role = ext_field.role,
                .component = std::string(name()),
                .name = ext_field.label,
                .size = ext_field.extent,
            },
        .storage = buffer.subspan(ext_field.offset, ext_field.extent)});
  }
  fields_ = std::move(bindings);
  initialized_ = true;
}

void CPL7Component::populate_registry(FieldRegistry& registry) {
  if (!initialized_) {
    throw std::runtime_error(std::string(name()) +
                             "Must initialize before populating registry");
  }
  for (auto& field : fields_) {
    coupling_.register_coupled_field(registry, field.field, field.storage);
  }
}
} // namespace e3sm::emulator
