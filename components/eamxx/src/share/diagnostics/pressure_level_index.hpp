#ifndef EAMXX_PRESSURE_LEVEL_INDEX_HPP
#define EAMXX_PRESSURE_LEVEL_INDEX_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"

namespace scream
{

/*
 * This diagnostic computes, for a given target pressure and a given
 * vertical layer (mid or int), the index of the first level whose pressure
 * is strictly greater than the target, in each column. This is numpy's
 * searchsorted (with side='right') on the pressure of the column, as std::upper_bound.
 * The output field has layout (ncol) and IntType data, with values in [0,nlevs].
 *
 * The output is masked: its valid mask is 1
 * where the target pressure falls within the column's pressure range (ends included),
 * and 0 otherwise. The index alone cannot tell, since nlevs is returned both when
 * the target is above the last level and when it is exactly equal to it.
 * The mask is the same for any field, since it only depends on the
 * pressure profile, so it has no dimension besides the columns, and it can be
 * broadcast to the layout of whichever field is interpolated.
 *
 * To linearly interpolate in a column, use the levels (k-1,k), with k clamped
 * to [1,nlevs-1]. This is always possible when the mask is 1, and the clamp
 * also covers a target exactly equal to the last entry. Where the mask is 0, the
 * clamped indices are still valid levels, so one can access data without
 * out-of-bounds accesses, and then ignore it.
 *
 * This diagnostic exists so that multiple diagnostics/fields requesting
 * the same "value at pressure level" (e.g. via FieldAtPressureLevel) can
 * share a single search, rather than each one performing its own
 * binary search, and can share the same mask.
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
  std::string         m_mask_name;
  std::string         m_diag_name;

  Real                m_pressure_level;
};

} //namespace scream

#endif // EAMXX_PRESSURE_LEVEL_INDEX_HPP
