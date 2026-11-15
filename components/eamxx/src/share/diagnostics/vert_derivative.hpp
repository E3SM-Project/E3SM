#ifndef EAMXX_VERT_DERIVATIVE_HPP
#define EAMXX_VERT_DERIVATIVE_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"

namespace scream {

/*
 * This diagnostic will calculate the derivative of a field X in
 * the vertical direction such that dX/dp or dX/dz
 */

class VertDerivative : public AbstractDiagnostic {
public:
  // Constructors
  VertDerivative(const ekat::Comm &comm, const ekat::ParameterList &params,
                 const std::shared_ptr<const AbstractGrid> &grid);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name() const { return "VertDerivative"; }

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_impl();
  void initialize_impl();

  std::string m_field_name; // Input field name

  // Name of derivative method (differential dp or dz)
  std::string m_derivative_method;
};

} // namespace scream

#endif // EAMXX_VERT_DERIVATIVE_HPP
