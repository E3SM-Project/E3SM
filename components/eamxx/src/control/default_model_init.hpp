#ifndef SCREAM_DEFAULT_MODEL_INIT_HPP
#define SCREAM_DEFAULT_MODEL_INIT_HPP

#include "share/data_managers/model_init.hpp"
#include "share/field/field_identifier.hpp"

namespace scream {

/*
 * Default implementation of ModelInit.
 *
 * Handles the standard restart and initial-conditions workflows,
 * including the fvphyshack pg2/gll logic, IOP data, and topography.
 * Registers under the "default" (and "homme") keys in ModelInitFactory.
 */
class DefaultModelInit : public ModelInit {
public:
  explicit DefaultModelInit (const ekat::ParameterList& params);

  void run () override;

protected:
  void restart_model () override;
  void set_initial_conditions () override;

private:
  using strvec_t = std::vector<std::string>;

#ifndef KOKKOS_ENABLE_CUDA
  // Cuda requires methods enclosing __device__ lambda's to be public
public:
#endif
  void initialize_constant_field (const FieldIdentifier& fid,
                                   const ekat::ParameterList& ic_pl);
};

std::shared_ptr<ModelInit>
create_default_model_init (const ekat::ParameterList& params);

void register_default_model_init ();

} // namespace scream

#endif // SCREAM_DEFAULT_MODEL_INIT_HPP
