#ifndef EAMXX_ACE_PROCESS_INTERFACE_HPP
#define EAMXX_ACE_PROCESS_INTERFACE_HPP

#include <torch/script.h>

#include "share/atm_process/atmosphere_process.hpp"
#include "share/grid/remap/abstract_remapper.hpp"

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

  // Need three grids (latlon, and pointgrid)
  std::shared_ptr<const AbstractGrid> m_ll_grid, m_pt_grid, m_sc_grid;

  // Need tensors for input and output
  torch::Tensor m_input_tensor, m_output_tensor;
  // Tensor geometries
  int m_batch = 1;
  int m_in_ch = 39;
  int m_out_ch = 44;
  int m_height = 180;
  int m_width = 360;
  // Levels
  int m_levels = 8;

  // Checkpoint path from the namelist 
  std::string m_checkpoint_path;
  // torch module containing the model
  torch::jit::script::Module m_module;
  // a guard for inference-only mode
  c10::InferenceMode m_guard;

  // Need three remappers (two horiz, one vert)
  std::shared_ptr<scream::AbstractRemapper> m_sc2ace_remapper, m_ace2sc_remapper, m_vert_remapper;

  // Helper functions
  void preprocess();
  void postprocess();

  // Store ak and bk values in a map {level, {ak, bk}}
  //  p_i = ak_i + bk_i * PRESsfc
  std::map<int, std::pair<Real, Real>> m_ak_bk_values = {
      {0, {64.247, 0.0}},
      {1, {5167.14603, 0.0}},
      {2, {12905.42546, 0.01755}},
      {3, {13982.4677, 0.11746}},
      {4, {12165.28766, 0.2896}},
      {5, {8910.07678, 0.49806}},
      {6, {4955.72632, 0.72625}},
      {7, {2155.78385, 0.88192}}
  };
};;

}  // namespace scream

#endif  // EAMXX_ACE_PROCESS_INTERFACE_HPP
