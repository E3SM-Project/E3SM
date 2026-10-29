#ifndef E3SM_COUPLER_DRIVER_HPP
#define E3SM_COUPLER_DRIVER_HPP

#include <coupler.hpp>
#include <memory>
#include <process.hpp>

namespace e3sm::coupler {

/**
 */
class CouplerDriver {
public:
  CouplerDriver() = default;
  // neither copyable nor movable
  CouplerDriver(CouplerDriver&) = delete;
  CouplerDriver(const CouplerDriver&) = delete;
  CouplerDriver& operator=(const CouplerDriver&) = delete;
  CouplerDriver& operator=(const CouplerDriver&&) = delete;

  template <CoupledProcess C> void add_component(std::unique_ptr<C> component) {
    if (initialized_) {
      throw std::runtime_error(
          "CouplerDriver: Can't add component after initialization!");
    }
    if (!component) {
      throw std::invalid_argument("CouplerDriver: Null component");
    }
    components_.emplace_back(std::move(component));
  }

  void initialize(const std::string& filename);
  bool is_initialized() const { return initialized_; }

  void run(const ProcessRunOpts run_opts) {
    for (auto& component : components_) {
      component.run(run_opts);
    }
  }

  void run_component(std::string_view name, const ProcessRunOpts opts);

  void export_component(std::string_view name);
  void import_component(std::string_view name);

  AnyComponent& get_component(std::string_view name);

  template <CoupledProcess C> C& get_component_as(std::string_view name) {
    return get_component(name).as<C>();
  }

  Coupler& coupler() noexcept { return coupler_; }

private:
  std::vector<AnyComponent> components_;
  Coupler coupler_;
  bool initialized_ = false;
};

} // namespace e3sm::coupler

#endif
