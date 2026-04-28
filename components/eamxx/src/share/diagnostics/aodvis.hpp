#ifndef EAMXX_AODVIS_DIAG
#define EAMXX_AODVIS_DIAG

#include "share/diagnostics/abstract_diagnostic.hpp"

namespace scream {

/*
 * This diagnostic will compute the visible aerosol optical depth.
 */

class AODVis : public AbstractDiagnostic {
 public:
  // Constructors
  AODVis(const ekat::Comm &comm, const ekat::ParameterList &params,
         const std::shared_ptr<const AbstractGrid>& grid);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name() const override { return "AerosolOpticalDepth550nm"; }

 protected:
#ifdef KOKKOS_ENABLE_CUDA
 public:
#endif
  void initialize_impl() override;
  void compute_diagnostic_impl() override;

  int m_swbands;
  int m_vis_bnd;
};

}  // namespace scream

#endif  // EAMXX_AODVIS_DIAG
