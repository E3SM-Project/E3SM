#include <algorithm>
#include <coupler_driver.hpp>
#include <iostream>
#include <stdexcept>

namespace e3sm::coupler {
void CouplerDriver::initialize(const std::string& filename) {
  for (auto& component : components_) {
    component.populate_registry(coupler_.registry());
  }
  coupler_.build_routes(filename);

  for (auto& component : components_) {
    const auto name = component.name();
    const auto& plan = coupler_.coupling_plan(std::string(name));
    component.configure_coupling(plan);
  }

  // coupler_.print_coupling_plan(std::cout);
}

AnyComponent& CouplerDriver::get_component(std::string_view name) {
  if (auto result = std::ranges::find_if(
          components_,
          [name](AnyComponent& comp) { return name == comp.name(); });
      result != components_.end())
    return *result;
  else
    throw std::runtime_error("Invalid component requested" + std::string(name));
}

void CouplerDriver::run_component(std::string_view name) {
  auto& comp = get_component(name);
  import_component(name);
  comp.run();
  export_component(name);
}

void CouplerDriver::export_component(std::string_view name) {
  auto& comp = get_component(name);
  const auto& buffers = coupler_.export_buffers(name);
  comp.export_fields(buffers);

  std::cout << "Export buffers for " << name << " after export\n";
  for (const auto& buf : buffers) {
    coupler_.print_field_buffer(std::cout, buf.id, buf.data, 10);
  }
  std::cout << std::endl;
}
void CouplerDriver::import_component(std::string_view name) {
  auto& comp = get_component(name);
  const auto& buffers = coupler_.import_buffers(name);

  std::cout << "Import buffers for " << name << " before import\n";
  std::cout << std::endl;
  for (const auto& buf : buffers) {
    coupler_.print_field_buffer(std::cout, buf.id, buf.data,10);
  }
  comp.import_fields(buffers);
}

} // namespace e3sm::coupler
