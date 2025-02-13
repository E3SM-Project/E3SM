#ifndef EAMXX_POMPEI_PROCESS_INTERFACE_HPP
#define EAMXX_POMPEI_PROCESS_INTERFACE_HPP

#include "share/atm_process/atmosphere_process.hpp"

namespace scream {

// Practicing On Manually Programming EAMxx Interfaces (POMPEI)
// Inject a tracer (e.g. ash) in the atmosphere

class PompeiEruption : public AtmosphereProcess {
 public:
  PompeiEruption(const ekat::Comm &comm, const ekat::ParameterList &params);

  AtmosphereProcessType type() const override {
    return AtmosphereProcessType::Physics;
  }

  std::string name() const override { return "POMPEI"; }

  void set_grids(
      const std::shared_ptr<const GridsManager> grids_manager) override;

 protected:
  void initialize_impl(const RunType run_type) override;
  void run_impl(const double dt) override;
  void finalize_impl() override;

  Field m_emission_mask;

  util::TimeStamp m_eruption_start;
};

}  // namespace scream

#endif  // EAMXX_POMPEI_PROCESS_INTERFACE_HPP
