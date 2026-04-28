#ifndef EAMXX_AEROCOMCLD_DIAG
#define EAMXX_AEROCOMCLD_DIAG

#include "share/diagnostics/abstract_diagnostic.hpp"

namespace scream {

/*
 * This diagnostic will compute the AeroCom diagnostics.
 */

class AeroComCld : public AbstractDiagnostic {
 public:
  // Constructors
  AeroComCld(const ekat::Comm &comm,
             const ekat::ParameterList &params,
             const std::shared_ptr<const AbstractGrid>& grid);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name() const override { return "AeroComCld"; }

 protected:
#ifdef KOKKOS_ENABLE_CUDA
 public:
#endif
  void compute_diagnostic_impl() override;

  // How many diags we have
  int m_ndiag;

  // Bot or Top
  std::string m_topbot;

  // Attribute maps for self-documentation
  std::map<std::string, int> m_index_map;
  std::map<std::string, std::string> m_units_map;

  // Set a field for dz
  Field m_dz;
};

}  // namespace scream

#endif  // EAMXX_AEROCOMCLD_DIAG
