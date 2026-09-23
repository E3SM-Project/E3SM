#ifndef E3SM_COUPLER_DRIVER_HPP
#define E3SM_COUPLER_DRIVER_HPP

#include <coupler.hpp>
#include <iostream>
#include <process.hpp>

namespace e3sm::coupler {

/**
 */
class CouplerDriver {
public:
  template <CoupledProcess C> void add_component(C& component) {
    components_.emplace_back(component);
  }

  void initialize(const std::string& filename);

  void run(const ProcessRunOpts run_opts) {
    for (auto& component : components_) {
      component.run(run_opts);
    }
  }

  void run_component(std::string_view name, const ProcessRunOpts opts);

  void export_component(std::string_view name);
  void import_component(std::string_view name);

  AnyComponent& get_component(std::string_view name);

  Coupler& coupler() noexcept { return coupler_; }

private:
  Coupler coupler_;
  std::vector<AnyComponent> components_;
};

} // namespace e3sm::coupler

#endif
