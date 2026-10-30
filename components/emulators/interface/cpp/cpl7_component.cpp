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
                        .size = field.data.size()},
        field.data);
  }
}

void CPL7Component::print_binding(std::ostream& os, ExternalField field) {

  os << "      " << "::" << field.name << " ["
     << " size=" << field.data.size()
     << ", addr=" << static_cast<const void*>(field.data.data()) << "]\n";

  os << "       [ ";
  for (std::size_t i = 0; i < field.data.size(); i++) {
    os << field.data[i] << ", ";
  }
  os << "]\n";
}

void CPL7Component::print_bindings() {
  std::cout << " ---- to_external ----- " << std::endl;
  for (auto& binding : bindings_.to_external) {
    print_binding(std::cout, binding);
  }

  std::cout << " ----- from_external ------ " << std::endl;
  for (auto& binding : bindings_.from_external) {
    print_binding(std::cout, binding);
  }
  std::cout << std::endl;
}

} // namespace e3sm::emulator
