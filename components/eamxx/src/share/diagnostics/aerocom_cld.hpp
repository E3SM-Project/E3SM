#ifndef EAMXX_AEROCOMCLD_DIAG
#define EAMXX_AEROCOMCLD_DIAG

#include "share/diagnostics/abstract_diagnostic.hpp"

namespace scream {

/*
 * A field X sampled at cloud top (or bottom), following the AeroCom
 * recommendation: per column, x = sum_k X(k)*W(k), where W are the
 * AeroComWeights for the same position and phase. With field_name='mask',
 * X is 1, and x is the cloud fraction seen from the top (or bottom).
 *
 * Params:
 *  - field_name: a Real (COL,LEV) field, or 'mask'
 *  - position: 'top' or 'bot'
 *  - phase: 'liq', 'ice', or 'none' (default)
 *
 * Levels with zero weight are skipped, so X need not be finite there.
 * If X is masked, its invalid levels are skipped too.
 */

class AeroComCld : public AbstractDiagnostic {
 public:
  AeroComCld(const ekat::Comm &comm, const ekat::ParameterList &params,
             const std::shared_ptr<const AbstractGrid>& grid);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name() const override { return "AeroComCld"; }

 protected:
#ifdef KOKKOS_ENABLE_CUDA
 public:
#endif
  void compute_impl() override;

  void initialize_impl() override;

  std::string m_field_name;
  std::string m_weights_name;
  std::string m_diag_name;
};

}  // namespace scream

#endif  // EAMXX_AEROCOMCLD_DIAG
