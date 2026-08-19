#ifndef E3SM_COUPLER_DRIVER_HPP
#define E3SM_COUPLER_DRIVER_HPP

#include <components.hpp>
#include <coupler.hpp>

namespace e3sm::coupler {


/**
 */
class CouplerDriver {
public:
  template <E3SMComponent C> void add_component(C& component) {
    components_.emplace_back(component);
  }

  void initialize(const std::string& filename);

  void run() {
    for (auto& component : components_) {
      component.run();
    }
  }
  AnyComponent& get_component(std::string_view name);

  Coupler& coupler() noexcept { return coupler_; }

private:
  Coupler coupler_;
  std::vector<AnyComponent> components_;
};

}

#endif
