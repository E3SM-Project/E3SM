#ifndef EAMXX_UNARY_MATH_DIAG_HPP
#define EAMXX_UNARY_MATH_DIAG_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"

namespace scream {

/*
 * Applies an elementwise math function to a field.
 *
 * Parameters:
 *  - field_name: the input field X (Real data type, any layout)
 *  - op:         one of
 *                  log   natural log      invalid where X<=0
 *                  exp   exponential
 *                  sqrt  square root      invalid where X<0
 *                  abs   absolute value
 *                  pow   X^exponent       invalid where X<0 and the exponent is
 *                                         not an integer, or X==0 and the
 *                                         exponent is negative
 *  - exponent:   (pow only) a number, as a string (e.g., "2", "0.5", "-1")
 *  - invalid:    (log/sqrt/pow only) what an invalid entry becomes:
 *                  nan   (default) whatever the math library returns (NaN or
 *                        -inf). With floating point exceptions enabled, this
 *                        traps: use it when X is not expected to be invalid.
 *                  mask  the entry is masked out (and set to the fill value)
 *  - units:      (log/exp/pow only)
 *                  check  (default) log and exp require a dimensionless X, and
 *                         pow requires a rational exponent (n/d, with d<=12)
 *                         unless X is dimensionless
 *                  ignore X's units are not checked, and the output is unitless
 *
 * Output units: log/exp: none; sqrt: sqrt(u); abs: u; pow: u^exponent.
 *
 * If X has a valid mask, the output has one too, and masked entries of X are
 * not computed (they may hold the fill value, whose power may overflow).
 */

class UnaryMath : public AbstractDiagnostic {
public:
  enum class Op { Log, Exp, Sqrt, Abs, Pow };

  UnaryMath (const ekat::Comm &comm, const ekat::ParameterList &params,
             const std::shared_ptr<const AbstractGrid>& grid);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name() const override { return "UnaryMath"; }

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_impl() override;

  void initialize_impl() override;

  std::string m_field_name;
  std::string m_op_str;
  Op          m_op;
  Real        m_exponent = 0;
  bool        m_mask_invalid = false;
  bool        m_ignore_units = false;
  bool        m_out_masked   = false;
};

} // namespace scream

#endif // EAMXX_UNARY_MATH_DIAG_HPP
