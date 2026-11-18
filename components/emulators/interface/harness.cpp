#include "harness.hpp"
#include "external_process.hpp"

#include <coupler.hpp>
#include <emulator.hpp>
#include <process.hpp>

#include <cmath>
#include <cstddef>
#include <memory>
#include <stdexcept>

namespace e3sm::emulator {
static_assert(CoupledProcess<ExternalProcess>);

namespace {

coupler::ComponentIntervals calc_coupling_intervals() {}

} // namespace

void Harness::add_bindings(ExternalBindings bindings) {
  if (initialized_) {
    throw std::runtime_error(
        "Harness::add_bindings: cannot attach fields after initialize()");
  }
  // Create an external process if none exists
  if (external_component_ == nullptr) {
    auto external = std::make_unique<ExternalProcess>();
    ExternalProcess* view = external.get();
    external_component_ = view;
    driver_.add_component(std::move(external));
  }
  external_component_->attach_fields(std::move(bindings));
}

void Harness::print_bindings(std::ostream& os) const {
  if (external_component_ == nullptr) {
    throw std::runtime_error(
        "Harness::print_bindings: call add_bindings() first");
  }
  external_component_->print_bindings(os);
}

void Harness::initialize(const std::vector<CoupledFieldEntry>& coupled_fields,
                         RouteBufferMap route_map) {
  if (initialized_) {
    throw std::runtime_error("Harness::initialize: already initialized");
  }

  // Every process's time step interval must be a whole number of the same base
  // step
  coupler::ComponentIntervals intervals;
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

  driver_.initialize(coupled_fields, route_map, intervals);
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

void Harness::import_to_external() {
  driver_.import_component(external_component_->name());
}
void Harness::export_from_external() {
  driver_.export_component(external_component_->name());
}

} // namespace e3sm::emulator
