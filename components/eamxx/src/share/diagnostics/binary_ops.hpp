#ifndef EAMXX_BINARY_OPS_DIAG_HPP
#define EAMXX_BINARY_OPS_DIAG_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"

namespace scream {

/*
 * This diagnostic will perform binary ops
 * like +, -, *, ÷ on two input fields or on an input field and a known physical constant.
 */

class BinaryOpsDiag : public AtmosphereDiagnostic {
 public:
  // Constructors
  BinaryOpsDiag(const ekat::Comm &comm, const ekat::ParameterList &params);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name() const  override{ return "BinaryOpsDiag"; }

  void create_requests () override;

 protected:
  void compute_diagnostic_impl() override;

  void initialize_impl(const RunType /*run_type*/) override;

  std::string m_arg1_name;
  std::string m_arg2_name;
  std::string m_binary_op_str;
  int         m_binary_op_code;

  bool m_arg1_is_field;
  bool m_arg2_is_field;
};  // class BinaryOpsDiag

}  // namespace scream

#endif  // EAMXX_BINARY_OPS_DIAG_HPP
