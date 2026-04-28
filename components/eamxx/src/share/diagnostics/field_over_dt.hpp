#ifndef EAMXX_FIELD_OVER_DT_DIAG_HPP
#define EAMXX_FIELD_OVER_DT_DIAG_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"
#include "share/util/eamxx_time_stamp.hpp"

namespace scream {

/*
 * This diagnostic divides any input field by the model timestep dt,
 * where dt = current_timestamp - start_of_step_timestamp.
 * Users can request X_over_dt for any field X.
 */

class FieldOverDt : public AbstractDiagnostic {
 public:
  FieldOverDt(const ekat::Comm &comm, const ekat::ParameterList &params,
              const std::shared_ptr<const AbstractGrid>& grid);

  std::string name() const override { return "FieldOverDt"; }

 protected:
#ifdef KOKKOS_ENABLE_CUDA
 public:
#endif
  void compute_diagnostic_impl() override;

  void init_timestep(const util::TimeStamp &start_of_step) override;

  void initialize_impl() override;

  std::string      m_name;
  util::TimeStamp  m_start_ts;  // saved at beginning of each timestep
};

}  // namespace scream

#endif  // EAMXX_FIELD_OVER_DT_DIAG_HPP
