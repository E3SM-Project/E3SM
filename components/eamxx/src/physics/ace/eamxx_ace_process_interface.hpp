#ifndef EAMXX_ACE_PROCESS_INTERFACE_HPP
#define EAMXX_ACE_PROCESS_INTERFACE_HPP

#include <torch/script.h>

#include "share/atm_process/atmosphere_process.hpp"

namespace scream {

class ACE : public AtmosphereProcess {
 public:
  ACE(const ekat::Comm &comm, const ekat::ParameterList &params);

  AtmosphereProcessType type() const override {
    return AtmosphereProcessType::Physics;
  }

  std::string name() const override { return "ace"; }

  void set_grids(
      const std::shared_ptr<const GridsManager> grids_manager) override;

#ifndef KOKKOS_ENABLE_CUDA
 protected:
#endif
  void initialize_impl(const RunType run_type) override;
  void run_impl(const double dt) override;
  void finalize_impl() override;

  std::shared_ptr<const AbstractGrid> m_grid;
  int m_forward_steps = 1;
  int m_forward_steps_in_memory = 1;
  int m_batch = 1;
  int m_in_ch = 39;
  int m_out_ch = 44;
  std::string m_checkpoint_path;
  // torch module containing the model
  torch::jit::script::Module m_module;
  // a guard for inference-only mode
  c10::InferenceMode m_guard;
};

}  // namespace scream

#endif  // EAMXX_ACE_PROCESS_INTERFACE_HPP
