// ============================================================
//  atm_backtend2.hpp   --  second finite difference in time
//  δ²φ(t) = BT(t) - BT(t-1)
//           = φ(t) - 2·φ(t-1) + φ(t-2)
// ============================================================
#ifndef SCREAM_ATM_BACKTEND2_HPP
#define SCREAM_ATM_BACKTEND2_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"
#include "share/util/eamxx_time_stamp.hpp"

namespace scream {

class AtmBacktend2 : public AtmosphereDiagnostic {
public:
  AtmBacktend2(const ekat::Comm& comm, const ekat::ParameterList& params);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name() const { return "AtmBacktend2"; }

  void create_requests ();

protected:
  void initialize_impl(const RunType run_type) override;
  void compute_diagnostic_impl() override;

private:
  // The tendency of what?
  std::string m_name;
  std::string m_diag_name;

  // Carry-forward state (never reset between output intervals)
  Field m_f_pre;    // φ(t-1)
  Field m_bt_pre;   // BT(t-1) = φ(t-1) - φ(t-2)
};

} // namespace scream
#endif // SCREAM_ATM_BACKTEND2_HPP
