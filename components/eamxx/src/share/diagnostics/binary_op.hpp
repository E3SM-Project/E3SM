#ifndef EAMXX_BINARY_OP_DIAG_HPP
#define EAMXX_BINARY_OP_DIAG_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"

#include <limits>

namespace scream {

/*
 * This diagnostic will perform binary ops
 * like +, -, *, ÷ on two inputs, which can be fields or known physical constants.
 *
 * It is requested as an expression, e.g.:
 *
 *  - T_mid/p_mid
 *  - qc+qr
 *  - Rgas*T_mid
 *  - cp-Rgas
 */

class BinaryOp : public AbstractDiagnostic {
 public:
  BinaryOp (const ekat::Comm &comm, const ekat::ParameterList &params,
            const std::shared_ptr<const AbstractGrid>& grid);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name() const  override{ return "BinaryOp"; }

 protected:
  void compute_impl() override;

  void initialize_impl() override;

  std::string m_arg1_name;
  std::string m_arg2_name;
  std::string m_binary_op_str;
  int         m_binary_op_code;

  bool m_arg1_is_field;
  bool m_arg2_is_field;

  // Value/units of a non-field operand, resolved once in the constructor.
  // Meaningless (and NaN) when the corresponding operand IS a field.
  Real m_arg1_value = std::numeric_limits<Real>::quiet_NaN();
  Real m_arg2_value = std::numeric_limits<Real>::quiet_NaN();
  ekat::units::Units m_arg1_units = ekat::units::Units::invalid();
  ekat::units::Units m_arg2_units = ekat::units::Units::invalid();

  bool m_arg1_has_mask = false;
  bool m_arg2_has_mask = false;
};

}  // namespace scream

#endif  // EAMXX_BINARY_OP_DIAG_HPP
