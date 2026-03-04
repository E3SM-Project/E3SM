#ifndef EAMXX_FIELD_PREV_DIAG_HPP
#define EAMXX_FIELD_PREV_DIAG_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"
#include "share/util/eamxx_time_stamp.hpp"

namespace scream {

/*
 * This diagnostic stores the value of a given field at the beginning
 * of the previous timestep. Users can request X_prev for any field X
 * to capture the start-of-timestep value before model updates occur.
 */

class FieldPrevDiag : public AtmosphereDiagnostic {
 public:
  // Constructors
  FieldPrevDiag(const ekat::Comm &comm, const ekat::ParameterList &params);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name() const { return "FieldPrevDiag"; }

  // Set the grid
  void create_requests();

 protected:
#ifdef KOKKOS_ENABLE_CUDA
 public:
#endif
  void compute_diagnostic_impl();

  // Store the field at the start of each timestep
  void init_timestep(const util::TimeStamp &start_of_step) override;

  // Set up the output field
  void initialize_impl(const RunType /*run_type*/) override;

  // The name of the field to track
  std::string m_name;

  // Store the field value from the beginning of the timestep
  Field m_f_prev;

};  // class FieldPrevDiag

}  // namespace scream

#endif  // EAMXX_FIELD_PREV_DIAG_HPP
