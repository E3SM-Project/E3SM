#include <algorithm>
#include <coupler_driver.hpp>
#include <stdexcept>

namespace e3sm::coupler {
void CouplerDriver::initialize(const std::string& filename) {
  for (auto& component : components_) {
    component.populate_registry(coupler_.registry());
  }
  coupler_.build_routes(filename);

}

AnyComponent& CouplerDriver::get_component(std::string_view name) {
  if (auto result = std::ranges::find_if(
          components_,
          [name](AnyComponent& comp) { return name == comp.name(); });
      result != components_.end())
    return *result;
  else
    std::runtime_error("Invalid component requested" + std::string(name));
}



} // namespace e3sm::coupler
