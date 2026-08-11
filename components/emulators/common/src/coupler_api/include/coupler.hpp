#ifndef E3SM_COUPLER_HPP
#define E3SM_COUPLER_HPP

#include <coupler_types.hpp>
#include <field_registry.hpp>

namespace e3sm::coupler {

class Coupler {
public:
  FieldRegistry& registry() noexcept { return registry_; }
  const FieldRegistry& registry() const noexcept { return registry_; }

  void build_routes();
  // void exchange();

  const std::vector<CoupledFieldEntry>
  read_coupling_fields_from_yaml(const std::string& filename);

private:
  FieldRegistry registry_;
  std::vector<CouplingRoute> routes_;
};
} // namespace e3sm::coupler
#endif
