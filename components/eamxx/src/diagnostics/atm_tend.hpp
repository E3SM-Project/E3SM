#ifndef EAMXX_ATM_TEND_DIAG_HPP
#define EAMXX_ATM_TEND_DIAG_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"
#include "share/util/scream_time_stamp.hpp"

namespace scream {

/*
 * This diagnostic will produce the atmospheric tendency.
 */

class AtmTendDiag : public AtmosphereDiagnostic {
 public:
  // Constructors
  AtmTendDiag(const ekat::Comm &comm, const ekat::ParameterList &params);

  // The name of the diagnostic
  std::string name() const;

  // Set the grid
  void set_grids(const std::shared_ptr<const GridsManager> grids_manager);

 protected:
#ifdef KOKKOS_ENABLE_CUDA
 public:
#endif
  void compute_diagnostic_impl();

 protected:
  void initialize_impl(const RunType /*run_type*/);

  // Keep track of field dimensions
  int m_num_cols;
  int m_num_levs;

  // The tendency of what?
  std::string m_name;

  // Store the previous field
  Field m_field_prev;

  // Store a time stamp
  util::TimeStamp m_ts;

};  // class AtmTendDiag

}  // namespace scream

#endif  // EAMXX_ATM_TEND_DIAG_HPP
