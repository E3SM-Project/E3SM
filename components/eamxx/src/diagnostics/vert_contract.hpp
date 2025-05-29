
#ifndef EAMXX_VERT_CONTRACT_HPP
#define EAMXX_VERT_CONTRACT_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"

namespace scream {

/*
 * This diagnostic will calculate the dp- or dz-weighted average of a field
 * across the LEV tag dimension, producing an N-1 dimensional field
 * that is a weighted average of the input field.
 */

class VertContractDiag : public AtmosphereDiagnostic {
 public:
  // Constructors
  VertContractDiag(const ekat::Comm &comm, const ekat::ParameterList &params);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name() const { return "VertContractDiag"; }

  // Set the grid
  void set_grids(const std::shared_ptr<const GridsManager> grids_manager);

 protected:
#ifdef KOKKOS_ENABLE_CUDA
 public:
#endif
  void compute_diagnostic_impl();
  void initialize_impl(const RunType /*run_type*/);
  // Additional function to scale the weights
  void scale_wts(Field &wts, const Field &wts_sum);


  // Name of each field (because the diagnostic impl is generic)
  std::string m_diag_name;
  // Name of contraction method (avg, sum)
  std::string m_contract_method;
  // Name of weighting method (dp, dz, none)
  std::string m_weighting_method;

  // Need some weighting, if unweighted, we will make it 1
  Field m_weighting;
  // Need a weighting field set to all ones
  Field m_weighting_one;
  // Need weighting summed vertically
  Field m_weighting_sum;
};

}  // namespace scream

#endif  // EAMXX_VERT_CONTRACT_HPP
