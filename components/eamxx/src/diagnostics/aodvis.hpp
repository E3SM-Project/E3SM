#ifndef EAMXX_AODVIS_DIAG
#define EAMXX_AODVIS_DIAG

#include "share/atm_process/atmosphere_diagnostic.hpp"
#include "share/util/eamxx_utils.hpp"

namespace scream {

/*
 * This diagnostic will compute the visible aerosol optical depth.
 */

class AODVis : public AtmosphereDiagnostic {
 public:
  // Constructors
  AODVis(const ekat::Comm &comm, const ekat::ParameterList &params);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name() const override { return "AerosolOpticalDepth550nm"; }

  // Set the grid
  void set_grids(
      const std::shared_ptr<const GridsManager> grids_manager) override;

 protected:
#ifdef KOKKOS_ENABLE_CUDA
 public:
#endif
  void compute_diagnostic_impl() override;

  int m_ncols;
  int m_nlevs;

  int m_swbands = eamxx_swbands();
  int m_vis_bnd = eamxx_vis_swband_idx();
};

}  // namespace scream

#endif  // EAMXX_AODVIS_DIAG
