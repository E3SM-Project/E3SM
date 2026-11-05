#include "harness.hpp"
#include "coupler.hpp"
#include "cpl7_component.hpp"
#include "emulator.hpp"
#include "process.hpp"
#include <cstddef>
#include <memory>
#include <stdexcept>

namespace e3sm::emulator {
static_assert(CoupledProcess<CPL7Component>);

void Harness::transfer_bindings(ExternalBindings bindings) {
  if (external_component_ != nullptr) {
    throw std::runtime_error(
        "Harness::transfer_bindings: bindings were already set");
  }

  auto external = std::make_unique<CPL7Component>();
  CPL7Component* view = external.get();

  external->attach_fields(std::move(bindings));
  driver_.add_component(std::move(external));

  external_component_ = view;
}

void Harness::print_bindings(std::ostream& os) const {
  if (external_component_ == nullptr) {
    throw std::runtime_error(
        "Harness::print_bindings: call transfer_bindings() first");
  }
  external_component_->print_bindings(os);
}

void Harness::initialize(const std::filesystem::path& plan,RouteBufferMap external) {
  if (initialized_) {
    throw std::runtime_error("Harness::initialize: already initialized");
  }
  driver_.initialize(plan,external);
  initialized_ = true;
}

void Harness::add_process(ProcessDescriptor desc) {

  driver_.add_component(
      std::make_unique<Emulator>(desc.grid, desc.emulator_config));
  processes_[desc.emulator_config.name] = std::move(desc.run_opts);
}

void Harness::run_process(std::string_view name) {
  if (!initialized_) {
    throw std::runtime_error("Harness: cannot run process '" +
                             std::string(name) + "' before initialize()");
  }

  ProcessMap::iterator it = processes_.find(std::string(name));

  if (it == processes_.end()) {
    std::string known;
    for (const auto& [process, opts] : processes_) {
      known += (known.empty() ? "" : ", ") + process;
    }
    throw std::runtime_error("Harness: no process named '" + std::string(name) +
                             "'. Processes added: " +
                             (known.empty() ? "none" : known));
  }
  driver_.run_component(name, it->second);
}

void Harness::import_to_external(std::string_view name) {
  driver_.import_component(name);
}
void Harness::export_from_external(std::string_view name) {
  driver_.export_component(name);
}

} // namespace e3sm::emulator
