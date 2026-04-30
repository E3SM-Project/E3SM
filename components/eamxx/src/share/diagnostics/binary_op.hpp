#ifndef EAMXX_BINARY_OP_DIAG_HPP
#define EAMXX_BINARY_OP_DIAG_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"

namespace scream {

/*
 * This diagnostic will perform binary ops
 * like +, -, *, ÷ on two input fields or on an input field and a known physical constant.
 */

class BinaryOp : public AbstractDiagnostic {
 public:
  BinaryOp(const ekat::Comm &comm, const ekat::ParameterList &params,
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

  bool m_arg1_has_mask = false;
  bool m_arg2_has_mask = false;
};

}  // namespace scream

#endif  // EAMXX_BINARY_OP_DIAG_HPP
