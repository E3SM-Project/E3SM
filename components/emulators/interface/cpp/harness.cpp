#include "harness.hpp"
#include "cpl7_component.hpp"
#include "emulator.hpp"
#include "field_registry.hpp"
#include <cstddef>
#include <limits>
#include <memory>
#include <ranges>
#include <stdexcept>

namespace e3sm::emulator {

void Harness::transfer_bindings(ExternalBindings bindings) {

  auto external = std::make_unique<CPL7Component>();
  CPL7Component* view = external.get();

  external->attach_fields(bindings);
  driver_.add_component(std::move(external));

  external_component_ = view;
}

void Harness::initialize(const std::filesystem::path& plan) {
  if (initialized_) {
    throw std::runtime_error("Harness already initialized");
  }
  driver_.initialize(plan);
}

void Harness::add_process(ProcessDescriptor desc) {

  driver_.add_component(
      std::make_unique<Emulator>(desc.grid, desc.emulator_config));
  processes_[desc.emulator_config.name] = std::move(desc.run_opts);
}

void Harness::run_process(std::string_view name) {
  ProcessMap::iterator it = processes_.find(std::string(name));

  if (it == processes_.end()) {
    throw std::runtime_error("Invalid name for BoundProcess " +
                             std::string(name));
  }
  driver_.run_component(name, it->second);
}

} // namespace e3sm::emulator
