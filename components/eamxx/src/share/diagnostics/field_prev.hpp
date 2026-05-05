#ifndef EAMXX_FIELD_PREV_HPP
#define EAMXX_FIELD_PREV_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"

namespace scream {

/*
 * This diagnostic stores the value of a given field at the beginning
 * of the current timestep (i.e., at the end of the previous timestep).
 * Users can request X_prev for any field X to capture the start-of-timestep
 * value before model updates occur.
 */

class FieldPrev : public AbstractDiagnostic {
 public:
  // Constructors
  FieldPrev(const ekat::Comm &comm, const ekat::ParameterList &params,
            const std::shared_ptr<const AbstractGrid>& grid);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name() const override { return "FieldPrev"; }

 protected:
#ifdef KOKKOS_ENABLE_CUDA
 public:
#endif
  void compute_impl() override;

  // Store the field at the start of each timestep
  void init_timestep(const util::TimeStamp &start_of_step) override;

  // Set up the output field
  void initialize_impl() override;

  // The name of the field to track
  std::string m_name;

};

}  // namespace scream

#endif  // EAMXX_FIELD_PREV_HPP
