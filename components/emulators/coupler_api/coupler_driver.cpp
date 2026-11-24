#include "coupler.hpp"
#include "coupler_driver.hpp"
#include "process.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>

namespace e3sm::emulator::coupler {
void CouplerDriver::initialize(const std::string& filename,
                               RouteBufferMap external_buffers,
                               const ComponentIntervals& intervals) {
  if (initialized_) {
    throw std::runtime_error("CouplerDriver::initialize: already initialized");
  }
  for (auto& component : components_) {
    component.populate_registry(coupler_.registry());
  }

  coupler_.build_routes(filename);
  coupler_.build_buffers(external_buffers, intervals);

  for (auto& component : components_) {
    const auto name = component.name();
    const auto& plan = coupler_.coupling_plan(std::string(name));
    component.configure_coupling(plan);
  }

  initialized_ = true;
}

AnyComponent& CouplerDriver::get_component(std::string_view name) {
  if (auto result = std::ranges::find_if(
          components_,
          [name](AnyComponent& comp) { return name == comp.name(); });
      result != components_.end())
    return *result;
  else {
    std::string known;
    for (const auto& component : components_) {
      known += (known.empty() ? "" : ", ") + std::string(component.name());
    }
    throw std::runtime_error("CouplerDriver: no component named '" +
                             std::string(name) + "'. Components added: " +
                             (known.empty() ? "none" : known));
  }
}

void CouplerDriver::run_component(std::string_view name,
                                  const ProcessRunOpts opts) {
  auto& comp = get_component(name);
  import_component(name);
  comp.run(opts);
  export_component(name);
}

void CouplerDriver::export_component(std::string_view name) {
  require_initialized("export", name);
  auto& comp = get_component(name);
  comp.export_fields(coupler_.export_buffers(name));
  coupler_.accumulate_exports(name);
}

void CouplerDriver::import_component(std::string_view name) {
  require_initialized("import", name);
  auto& comp = get_component(name);
  coupler_.prepare_imports(name);
  comp.import_fields(coupler_.import_buffers(name));
}

void CouplerDriver::require_initialized(std::string_view action,
                                        std::string_view name) const {
  if (!initialized_) {
    throw std::runtime_error("CouplerDriver: cannot " + std::string(action) +
                             " '" + std::string(name) +
                             "' before initialize()");
  }
}

} // namespace e3sm::emulator::coupler
