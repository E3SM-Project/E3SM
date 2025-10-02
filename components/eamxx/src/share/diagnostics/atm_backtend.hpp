#ifndef EAMXX_ATM_BACKTEND_DIAG_HPP
#define EAMXX_ATM_BACKTEND_DIAG_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"
#include "share/util/eamxx_time_stamp.hpp"

namespace scream {

/*
 * This diagnostic will back out the atmosphere tendency of a given field.
 */

class AtmBackTendDiag : public AtmosphereDiagnostic {
 public:
  // Constructors
  AtmBackTendDiag(const ekat::Comm &comm, const ekat::ParameterList &params);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name() const { return "AtmBackTendDiag"; }

  // Set the grid
  void set_grids(const std::shared_ptr<const GridsManager> grids_manager);

 protected:
#ifdef KOKKOS_ENABLE_CUDA
 public:
#endif
  void compute_diagnostic_impl();

  // Let's override the init time step method
  void init_timestep(const util::TimeStamp &start_of_step) override;

  // Let's override the initialize method to set the fields below
  void initialize_impl(const RunType /*run_type*/) override;

  // Keep track of field dimensions
  int m_num_cols;
  int m_num_levs;

  // The tendency of what?
  std::string m_name;

  // Store the previous field
  Field m_f_prev;

};  // class AtmBackTendDiag

}  // namespace scream

#endif  // EAMXX_ATM_BACKTEND_DIAG_HPP
