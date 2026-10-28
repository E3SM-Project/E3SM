#include "harness.hpp"
#include "emulator.hpp"
#include <memory>
#include <ranges>
#include <stdexcept>

namespace e3sm::emulator {

void Harness::add_process(ProcessDescriptor desc) {

  driver_.add_component(
      std::make_unique<Emulator>(desc.grid, desc.emulator_config));

  processes_[desc.emulator_config.name] =
      BoundProcess{.descriptor = std::move(desc)};
}

void Harness::run_process(std::string_view name) {
  ProcessMap::iterator it = processes_.find(std::string(name));

  if (it == processes_.end()) {
    throw std::runtime_error("Invalid name for BoundProcess " +
                             std::string(name));
  }
  driver_.run_component(name, it->second.descriptor.run_opts);
}

} // namespace e3sm::emulator
