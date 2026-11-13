#include "harness.hpp"
#include "cpl7_component.hpp"

#include <coupler.hpp>
#include <emulator.hpp>
#include <process.hpp>

#include <cmath>
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

void Harness::initialize(const std::filesystem::path& plan,
                         RouteBufferMap external) {
  if (initialized_) {
    throw std::runtime_error("Harness::initialize: already initialized");
  }

  // Every process's time step interval must be a whole number of the same base step
  ComponentIntervals intervals;
  const std::string* base_process = nullptr;

  double base_step = 0.0;
  for (const auto& [name, state] : processes_) {
    const auto interval = state.run_opts.coupling_interval;
    if (interval < 1) {
      throw std::runtime_error(
          "Harness: process '" + name + "' has a coupling_interval of " +
          std::to_string(interval) + "; it must be at least 1");
    }
    intervals.emplace(name, interval);

    const double step = static_cast<double>(state.run_opts.dt) / interval;
    if (base_process == nullptr) {
      base_process = &name;
      base_step = step;
    } else if (std::abs(step - base_step) > 1e-6 * std::abs(base_step)) {
      throw std::runtime_error(
          "Harness: dt / coupling_interval must be the same for every process, "
          "but it's " +
          std::to_string(base_step) + " for '" + *base_process + "' and " +
          std::to_string(step) + " for '" + name + "'");
    }
  }

  driver_.initialize(plan, external, intervals);
  initialized_ = true;
}

void Harness::add_process(ProcessDescriptor desc) {

  driver_.add_component(
      std::make_unique<Emulator>(desc.grid, desc.emulator_config));

  processes_[desc.emulator_config.name] =
      ProcessState{.run_opts = desc.run_opts};
}

bool Harness::run_process(std::string_view name) {
  if (!initialized_) {
    throw std::runtime_error("Harness: cannot run process '" +
                             std::string(name) + "' before initialize()");
  }

  ProcessMap::iterator it = processes_.find(std::string(name));

  if (it == processes_.end()) {
    std::string known;
    for (const auto& [process, state] : processes_) {
      known += (known.empty() ? "" : ", ") + process;
    }
    throw std::runtime_error(
        "Harness: no process named '" + std::string(name) +
        "'. Processes added: " + (known.empty() ? "none" : known));
  }

  auto& state = it->second;
  if (++state.calls_since_run < state.run_opts.coupling_interval) {
    return false;
  }
  state.calls_since_run = 0;
  driver_.run_component(name, state.run_opts);
  return true;
}

void Harness::import_to_external(std::string_view name) {
  driver_.import_component(name);
}
void Harness::export_from_external(std::string_view name) {
  driver_.export_component(name);
}

} // namespace e3sm::emulator
