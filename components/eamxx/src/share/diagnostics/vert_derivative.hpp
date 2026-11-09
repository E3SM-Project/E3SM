#ifndef EAMXX_VERT_DERIVATIVE_HPP
#define EAMXX_VERT_DERIVATIVE_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"

namespace scream {

/*
 * The derivative of a field X along the vertical, with respect to a coordinate
 * field on the same layout (e.g., p_mid or z_mid), computed in each column the
 * way xarray's differentiate (numpy's gradient, edge_order=1) computes it:
 * second-order central differences on the coordinate's (uneven) spacing at
 * interior levels, and one-sided first-order differences at the end levels.
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
  std::string m_coord_name; // Coordinate field name (e.g., p_mid or z_mid)
};

} // namespace scream

#endif // EAMXX_VERT_DERIVATIVE_HPP
