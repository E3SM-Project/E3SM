#ifndef EAMXX_PRESSURE_LEVEL_INDEX_HPP
#define EAMXX_PRESSURE_LEVEL_INDEX_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"

namespace scream
{

/*
 * This diagnostic computes, for a given target pressure and a given
 * vertical layer (mid or int), the pair of level indices that bracket
 * the target pressure in each column, ready to be used for a linear
 * interpolation. The output field has layout (ncol,2) and IntType data.
 *
 * For columns where the target pressure falls outside the column's
 * pressure range, both indices are set to -1, signaling that the
 * bracket is not meaningful (rather than attaching a separate valid
 * mask field, which would have to share this field's own (ncol,2)
 * layout, and hence could not represent a single flag per column).
 *
 * This diagnostic exists so that multiple diagnostics/fields requesting
 * the same "value at pressure level" (e.g. via FieldAtPressureLevel) can
 * share a single bracket-index computation, rather than each one
 * performing its own binary search.
 */

class PressureLevelIndex : public AbstractDiagnostic
{
public:

  // Constructors
  PressureLevelIndex (const ekat::Comm& comm, const ekat::ParameterList& params,
                      const std::shared_ptr<const AbstractGrid>& grid);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name () const { return "PressureLevelIndex"; }

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_impl ();
protected:
  void initialize_impl ();

  std::string         m_pressure_name;
  std::string         m_diag_name;

  Real                m_pressure_level;
};

} //namespace scream

#endif // EAMXX_PRESSURE_LEVEL_INDEX_HPP
