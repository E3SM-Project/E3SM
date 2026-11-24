#ifndef EAMXX_BINARY_OP_DIAG_HPP
#define EAMXX_BINARY_OP_DIAG_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"

#include "share/util/eamxx_quantity.hpp"

#include <optional>

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

  // Value and units (invalid for a literal) of an operand that is not a field. Empty for a field.
  std::optional<Quantity<Real>> m_arg1_scalar;
  std::optional<Quantity<Real>> m_arg2_scalar;

  bool m_arg1_has_mask = false;
  bool m_arg2_has_mask = false;
};

}  // namespace scream

#endif  // EAMXX_BINARY_OP_DIAG_HPP
