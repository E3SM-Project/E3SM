#ifndef EAMXX_WATER_TRACERS_FUNCTIONS_HPP
#define EAMXX_WATER_TRACERS_FUNCTIONS_HPP

#include "share/core/eamxx_types.hpp"  // for scream::Real and scream::sp()

#include <ekat_pack.hpp>
#include <ekat_pack_utils.hpp>  // ekat::impl::merge for conditional selection

namespace scream {
namespace water_tracers {

/*
 * Utility functions for water tracer computations.
 *
 * These are pure, device-callable functions for computing tracer ratios
 * and related quantities. Ported from the Fortran module water_tracers.F90
 * (author David Noone). Templated on ScalarT so they work for both a plain
 * Real and an ekat::Pack<Real,N>.
 */

// Water tracer species enumeration
// This is a simplified version - the actual species would be defined
// based on which tracers are active in the simulation
enum WaterTracerSpecies {
  H2O_BASE = 0  // Base water (H2O-16)
};

struct WaterTracerFunctions
{
  using Real = scream::Real;

  // Minimum water amount for stable ratio computation.
  // From water_tracers.F90 (David Noone): smaller makes scheme more accurate.
  static constexpr Real wtrc_qmin = sp(1.e-22);

  // -----------------------------------------------------------------------
  // Standard tracer ratio by species.
  //
  // Returns the model standard tracer ratio R_std for a given species.
  // For simple water tracers (non-fractionating), this is unity.
  // -----------------------------------------------------------------------
  KOKKOS_INLINE_FUNCTION
  static Real wtrc_get_rstd(const WaterTracerSpecies species)
  {
    // For non-fractionating water tracers, standard ratio is 1.0
    (void)species;  // Unused parameter in this implementation
    return sp(1.0);
  }

  // -----------------------------------------------------------------------
  // Tracer ratio from masses with numerical checks.
  //
  // Computes the tracer ratio R = qtrc / qtot with safeguards against
  // division by tiny or zero denominators. When abs(qtot) < wtrc_qmin,
  // returns the standard ratio R_std instead of computing the division.
  //
  // Port of water_tracers.F90::wtrc_ratio (David Noone).
  //
  // Parameters:
  //   species - tracer species identifier
  //   qtrc    - tracer water mass [kg/kg]
  //   qtot    - "base" water mass [kg/kg]
  //
  // Returns:
  //   Tracer ratio R (dimensionless)
  //
  // Templated on ScalarT to work with both plain Real and ekat::Pack<Real,N>.
  // -----------------------------------------------------------------------
  template <typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static ScalarT wtrc_ratio(const WaterTracerSpecies species,
                            const ScalarT& qtrc,
                            const ScalarT& qtot)
  {
    using ekat::abs;  // ADL-compatible abs for Pack<Real,N>

    const ScalarT abs_qtot = abs(qtot);
    const ScalarT qmin = sp(wtrc_qmin);

    // Mask: use standard ratio when total water is below threshold.
    const auto use_standard = (abs_qtot < qmin);

    const ScalarT ratio_computed = qtrc / qtot;
    const ScalarT ratio_standard = sp(wtrc_get_rstd(species));

    // Blend based on mask (works for both Real and Pack).
    return ekat::impl::merge(use_standard, ratio_standard, ratio_computed);
  }
};

} // namespace water_tracers
} // namespace scream

#endif // EAMXX_WATER_TRACERS_FUNCTIONS_HPP
