#include "field_registry.hpp"
#include <cpl7_component.hpp>
#include <stdexcept>

namespace e3sm::emulator {
void CPL7Component::attach_fields(ExternalBindings bindings) {
  if (initialized_) {
    throw std::runtime_error("CPL7Component: fields were already attached");
  }
  bindings_ = std::move(bindings);
  initialized_ = true;
}

void CPL7Component::populate_registry(FieldRegistry& registry) {
  if (!initialized_) {
    throw std::runtime_error(
        "CPL7Component: attach fields before populating the registry");
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

void CPL7Component::print_binding(std::ostream& os,
                                  const ExternalField& field) const {

  os << "      " << "::" << field.name << " ["
     << " size=" << field.data.size()
     << ", addr=" << static_cast<const void*>(field.data.data()) << "]\n";

  os << "       [ ";
  for (std::size_t i = 0; i < field.data.size(); i++) {
    os << field.data[i] << ", ";
  }
  os << "]\n";
}

void CPL7Component::print_bindings(std::ostream& os) const {
  os << " ---- to_external ----- \n";
  for (const auto& binding : bindings_.to_external) {
    print_binding(os, binding);
  }

  os << " ----- from_external ------ \n";
  for (const auto& binding : bindings_.from_external) {
    print_binding(os, binding);
  }
  os << '\n';
}

} // namespace e3sm::emulator
