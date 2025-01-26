#ifndef EAMXX_NUMBER_PATH_DIAGNOSTIC_HPP
#define EAMXX_NUMBER_PATH_DIAGNOSTIC_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"

namespace scream {

/*
 * This diagnostic will produce the "number" path.
 */

class NumberPathDiagnostic : public AtmosphereDiagnostic {
 public:
  // Constructors
  NumberPathDiagnostic(const ekat::Comm &comm,
                       const ekat::ParameterList &params);

  // The name of the diagnostic
  std::string name() const override { return "NumberPath"; }

  // Set the grid
  void set_grids(const std::shared_ptr<const GridsManager> grids_manager);

 protected:
#ifdef KOKKOS_ENABLE_CUDA
 public:
#endif
  void compute_diagnostic_impl();

 protected:
  // Keep track of field dimensions
  int m_num_cols;
  int m_num_levs;

  // Keep track of which quantities to compute
  std::map<std::string, std::string> m_qnames;
  std::map<std::string, std::string> m_nnames;
  std::vector<std::string> m_kinds;
};  // class NumberPathDiagnostic

}  // namespace scream

#endif  // EAMXX_NUMBER_PATH_DIAGNOSTIC_HPP
