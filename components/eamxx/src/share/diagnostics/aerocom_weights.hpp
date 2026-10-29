#ifndef EAMXX_AEROCOM_WEIGHTS_DIAG_HPP
#define EAMXX_AEROCOM_WEIGHTS_DIAG_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"

namespace scream {

/*
 * The AeroCom cloud-top (or cloud-bottom) weights of each layer: the probability
 * that the layer is the first cloudy one seen from the model top (or bottom),
 * under the maximum-random overlap assumption, times the fraction of the cloud
 * condensate in the requested phase (if any). The weights are 0 at layers that
 * do not pass the cloud thresholds below.
 *
 * Params:
 *  - position: 'top' or 'bot'
 *  - phase: 'liq', 'ice', or 'none' (default)
 *
 * AeroComCld contracts a field against these weights. They are a diag of their
 * own so that all fields sampled at the same position and phase share them.
 * See docs/technical/aerocom_cldtop.md for the algorithm.
 */

class AeroComWeights : public AbstractDiagnostic {
 public:
  // A layer is cloudy if qc+qi and cldfrac_tot exceed these
  static constexpr Real q_threshold = 0;
  static constexpr Real cldfrac_tot_threshold = 0.001;

  AeroComWeights(const ekat::Comm &comm, const ekat::ParameterList &params,
                 const std::shared_ptr<const AbstractGrid>& grid);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name() const override { return "AeroComWeights"; }

  // The name of the output field: the expression that builds this diag, and
  // the name AeroComCld asks for. Throws on an invalid position or phase.
  static std::string expression (const std::string& position, const std::string& phase);

 protected:
#ifdef KOKKOS_ENABLE_CUDA
 public:
#endif
  void compute_impl() override;

  void initialize_impl() override;

  std::string m_position;
  std::string m_phase;
  std::string m_diag_name;
};

}  // namespace scream

#endif  // EAMXX_AEROCOM_WEIGHTS_DIAG_HPP
