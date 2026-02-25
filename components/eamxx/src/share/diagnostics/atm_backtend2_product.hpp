// ============================================================
//  atm_backtend2_product.hpp
//
//  Signed BT product:  P(t) = BT(t) · BT(t-1)
//
//  For a 2Δt oscillation, consecutive backtends are large and
//  opposite in sign so P(t) ≪ 0.  For smooth evolution P ≈ 0.
// ============================================================
#ifndef SCREAM_ATM_BACKTEND2_PRODUCT_HPP
#define SCREAM_ATM_BACKTEND2_PRODUCT_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"
#include "share/util/eamxx_time_stamp.hpp"

namespace scream {

class AtmBacktend2Product : public AtmosphereDiagnostic {
public:
  AtmBacktend2Product(const ekat::Comm& comm,
                     const ekat::ParameterList& params);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name() const { return "AtmBacktend2Product"; }

  void create_requests ();

protected:
  void initialize_impl(const RunType run_type) override;
  void compute_diagnostic_impl() override;

private:
  std::string m_name;
  std::string m_diag_name;

  // Carry-forward state
  Field m_f_pre;    // φ(t-1)
  Field m_bt_pre;   // BT(t-1) = φ(t-1) - φ(t-2)
};

} // namespace scream
#endif // SCREAM_ATM_BACKTEND2_PRODUCT_HPP
