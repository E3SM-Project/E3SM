#include "field_registry.hpp"
#include <cpl7_component.hpp>
#include <stdexcept>

namespace e3sm::emulator {
void CPL7Component::attach_fields(ExternalBindings bindings) {
  if (initialized_) {
    throw std::runtime_error("Re-attaching fields not allowed");
  }
  bindings_ = std::move(bindings);
  initialized_ = true;
}

void CPL7Component::populate_registry(FieldRegistry& registry) {
  if (!initialized_) {
    throw std::runtime_error(std::string(name()) +
                             "Must initialize before populating registry");
  }
  for (auto& field : bindings_.to_external) {
    coupling_.register_coupled_field(
        registry,
        RegisteredField{.role = FieldRole::Import,
                        .component = std::string(name()),
                        .name = field.name,
                        .size = field.data.size()},
        field.data);
  }
  for (auto& field : bindings_.from_external) {
    coupling_.register_coupled_field(
        registry,
        RegisteredField{.role = FieldRole::Export,
                        .component = std::string(name()),
                        .name = field.name,
      .size=field.data.size()},
        field.data);
  }
}
} // namespace e3sm::emulator
