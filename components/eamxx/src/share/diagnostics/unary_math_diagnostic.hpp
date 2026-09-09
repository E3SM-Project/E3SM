#ifndef EAMXX_UNARY_MATH_DIAGNOSTIC_HPP
#define EAMXX_UNARY_MATH_DIAGNOSTIC_HPP

#include "abstract_diagnostic.hpp"

namespace scream {

template <typename Op>
class UnaryMathDiagnostic : public AbstractDiagnostic {
public:
  UnaryMathDiagnostic(const ekat::Comm& comm,
                      const ekat::ParameterList& params,
                      const std::shared_ptr<const AbstractGrid>& grid);

  std::string name() const override;

protected:
  void initialize_impl() override;

  void compute_impl() override;
};

#include "unary_mth_ops.hpp"

// Decalre the ones that will be available
#define EAMXX_ETI_UNARY_DIAG(OpName) \
  template class UnaryMathDiagnostic<math_ops::OpName>; \
  using OpName##Diagnostic = UnaryMathDiagnostic<math_ops::OpName>;

EAMXX_ETI_UNARY_DIAG(Sqrt)
EAMXX_ETI_UNARY_DIAG(Log)
EAMXX_ETI_UNARY_DIAG(Exp)
EAMXX_ETI_UNARY_DIAG(Sin)
EAMXX_ETI_UNARY_DIAG(Cos)
EAMXX_ETI_UNARY_DIAG(Tan)
EAMXX_ETI_UNARY_DIAG(Atan)

} // namespace scream

#endif // EAMXX_UNARY_MATH_DIAGNOSTIC_HPP
