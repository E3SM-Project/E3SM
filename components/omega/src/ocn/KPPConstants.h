#ifndef OMEGA_KPP_CONSTANTS_H
#define OMEGA_KPP_CONSTANTS_H
//===-- ocn/KPPConstants.h - KPP Constants and Profiles --------*- C++ -*-===//
//
/// \file
/// \brief KPP-specific constants, parameters, and profile functions
///
/// This header defines constants, parameters, and inline functions for the
/// K-Profile Parameterization boundary layer mixing scheme. Includes
/// non-dimensional profile functions used in KPP coefficient calculations.
//
//===----------------------------------------------------------------------===//

#include "GlobalConstants.h"
#include "OmegaKokkos.h"

namespace OMEGA::KPP {

// ==========================================================================
// Monin-Obukhov stability function parameters (Large et al. 1994, App. B)
//
// The stability coordinate is zeta = d/L, where d is depth below the surface
// and L is the Monin-Obukhov length. zeta > 0 is stable (surface warming or
// salinification), zeta < 0 is unstable (convective). The values below set
// where the weakly-unstable branch hands off to the strongly-unstable branch.
// ==========================================================================

/// Transition zeta for momentum: below this value, strongly-unstable formula
/// is used. Default: -0.2 (CVMix default).
constexpr Real ZetaM = -0.2_Real;

/// Transition zeta for scalars: below this value, strongly-unstable formula
/// is used. Default: -1.0 (CVMix default).
constexpr Real ZetaS = -1.0_Real;

// The four constants below are fixed by requiring the strongly-unstable
// branch to match the weakly-unstable branch in value at the transition zeta.
// These constants should be treated as fixed, they are fits to previous
// field data for the surface layer under varying forcing conditions
//
/// a_m = (1-16*ZetaM)^{-0.25} * (1 - 4*ZetaM)
constexpr Real AMoM = 1.2573615702_Real;
/// c_m = (1-16*ZetaM)^{-0.25} * 12
constexpr Real CMoM = 8.3824104679_Real;

/// a_s = sqrt(1-16*ZetaS) * (1 + 8*ZetaS)  (can be negative)
constexpr Real AMoS = -28.8617393793_Real;
/// c_s = 24 * sqrt(1-16*ZetaS)
constexpr Real CMoS = 98.9545350148_Real;

/// Surface value of the momentum shape function (dimensionless)
constexpr Real HuOn = 0.03_Real;

/// Floor on friction velocity (m/s). Keeps the turbulent velocity scales and
/// the Monin-Obukhov length finite in near-calm conditions.
constexpr Real MinUStar = 1.0e-4_Real;

// ===========================================================================
// Numerical Safeguards
// ===========================================================================

/// Generic tolerance for zero checks and degenerate geometric intervals
constexpr Real NumericalTolerance = 1.0e-12_Real;

/// Effectively nonzero denominator used where no physical floor is intended
constexpr Real Tiny = 1.0e-20_Real;

/// Minimum unresolved turbulent shear squared (m^2/s^2) in the bulk Ri
/// denominator
constexpr Real MinUnresolvedShearSq = 1.0e-10_Real;

/// Minimum Stokes drift (m/s) used to keep the Langmuir number finite
constexpr Real MinStokesDrift = 1.0e-8_Real;

// ==========================================================================
// Ocean Surface Boundary Layer (OSBL) Parameters
// ==========================================================================

/// The OSBL is the near-surface ocean layer where KPP represents turbulent
/// mixing driven by surface momentum and buoyancy forcing.
///
/// Critical bulk Richardson number defining the OSBL base
/// REFERENCES: Large et al. (1994) Eq. (21)
constexpr Real CriticalRi = 0.25_Real;

/// Minimum OSBL depth under sea ice (m), applied above IceSuppressThresh
constexpr Real MinOSBLUnderIce = 5.0_Real;

/// Ice fraction above which Langmuir enhancement is disabled
constexpr Real IceFracThresh = 0.05_Real;

/// Ice fraction above which the minimum OSBL depth is enforced
constexpr Real IceSuppressThresh = 0.15_Real;

// ==========================================================================
// Surface Layer Parameters
// ==========================================================================

/// Surface layer extent as a fraction of the trial OSBL depth (epsilon in
/// Large et al. 1994). Reference values entering the bulk Richardson number
/// are averaged over the top SurfaceLayerExtent * d of the column.
constexpr Real SurfaceLayerExtent = 0.1_Real;

/// Empirical convective velocity coefficient in the turbulent velocity scale
/// w_s = (u*^3 + ConvectiveVelFac * max(-B_0,0) * h)^(1/3)
constexpr Real ConvectiveVelFac = 0.35_Real;

// ==========================================================================
// KPP Shape and Stability Functions
//
// Vertical position inside the OSBL is expressed two ways:
//   Sigma    in [-1,0], the Omega convention: 0 at the surface, -1 at the
//            OSBL base (it follows the sign of the z coordinate).
//   SigmaMu  in [ 0,1], the CVMix/Large et al. convention: SigmaMu = -Sigma,
//            so 0 at the surface and 1 at the OSBL base.
// The shape functions below take Sigma and convert internally, so callers
// never need to flip signs.
// ==========================================================================

/// @brief Momentum gradient shape function, multiplied by h and the turbulent
/// velocity scale to give the KPP viscosity: Kx = h * w_m * G(sigma).
/// REFERENCES: Large et al. (1994) Eq. (11)
///
/// @param Sigma Normalized vertical position (-z/h)
/// @return Dimensionless shape value
KOKKOS_INLINE_FUNCTION
Real kppShapeMomentum(Real Sigma) {
   Sigma = Kokkos::fmax(-1.0_Real, Kokkos::fmin(0.0_Real, Sigma));

   const Real SigmaMu = -Sigma;
   return SigmaMu * (1.0_Real - SigmaMu) * (1.0_Real - SigmaMu);
}

/// @brief Matched KPP gradient shape function.
///
/// Uses the SimpleShapes gradient shape plus a smooth correction that is
/// zero at the surface and equals ShapeAtBase at the OSBL base. This lets
/// MatchBoth profiles meet pre-existing interior mixing at the OSBL base while
/// preserving SimpleShapes behavior when ShapeAtBase is zero.
KOKKOS_INLINE_FUNCTION
Real kppShapeMatched(Real Sigma, Real ShapeAtBase) {
   Sigma = Kokkos::fmax(-1.0_Real, Kokkos::fmin(0.0_Real, Sigma));

   const Real SigmaMu = -Sigma;
   const Real Simple  = SigmaMu * (1.0_Real - SigmaMu) * (1.0_Real - SigmaMu);
   const Real Smooth  = SigmaMu * SigmaMu * (3.0_Real - 2.0_Real * SigmaMu);
   return Simple + ShapeAtBase * Smooth;
}

/// @brief phi_m^{-1}(zeta) - Inverse momentum Monin-Obukhov stability function
/// Multiplied by von Karman constant and friction velocity to give turbulent
/// momentum velocity scale: w_m = kappa * u* * phi_m^{-1}(zeta)
/// Three-regime formulation per Large et al. (1994) Appendix B and CVMix.
///
/// @param Zeta Monin-Obukhov stability coordinate d/L (dimensionless)
/// @return phi_m^{-1} (dimensionless, > 0)
KOKKOS_INLINE_FUNCTION
Real kppPhiInvMomentum(Real Zeta) {
   if (Zeta >= 0.0_Real) {
      // Stable regime
      return 1.0_Real / (1.0_Real + 5.0_Real * Zeta);
   } else if (Zeta >= ZetaM) {
      // Weakly unstable: (1 - 16*zeta)^{1/4}
      return Kokkos::pow(1.0_Real - 16.0_Real * Zeta, 0.25_Real);
   } else {
      // Strongly unstable (convective): (a_m - c_m*zeta)^{1/3}
      return Kokkos::pow(AMoM - CMoM * Zeta, 1.0_Real / 3.0_Real);
   }
}

/// @brief Scalar gradient shape function, multiplied by h and the turbulent
/// velocity scale to give the KPP diffusivity: Kx = h * w_s * G(sigma).
/// The non-local flux always uses this shape, scaled by C_s instead of h*w_s,
/// even under MatchBoth, so that gamma vanishes at the OSBL base.
/// REFERENCES: Large et al. (1994) Eq. (11), Eq. (12)-(13), Large et al. (1997)
///
/// @param Sigma Normalized vertical position (-z/h)
/// @return Dimensionless shape value
KOKKOS_INLINE_FUNCTION
Real kppShapeScalar(Real Sigma) {
   Sigma = Kokkos::fmax(-1.0_Real, Kokkos::fmin(0.0_Real, Sigma));

   const Real SigmaMu = -Sigma;
   return SigmaMu * (1.0_Real - SigmaMu) * (1.0_Real - SigmaMu);
}

/// @brief phi_s^{-1}(zeta) - Inverse scalar Monin-Obukhov stability function
/// Multiplied by von Karman constant and friction velocity to give turbulent
/// scalar velocity scale: w_s = kappa * u* * phi_s^{-1}(zeta)
/// Three-regime formulation per Large et al. (1994) Appendix B and CVMix.
/// Scalars mix more efficiently than momentum in unstable conditions, which
/// is why the weakly-unstable exponent is 1/2 here and 1/4 for momentum.
///
/// @param Zeta Monin-Obukhov stability coordinate d/L (dimensionless)
/// @return phi_s^{-1} (dimensionless, > 0)
KOKKOS_INLINE_FUNCTION
Real kppPhiInvScalar(Real Zeta) {
   if (Zeta >= 0.0_Real) {
      // Stable regime
      return 1.0_Real / (1.0_Real + 5.0_Real * Zeta);
   } else if (Zeta >= ZetaS) {
      // Weakly unstable: (1 - 16*zeta)^{1/2}
      return Kokkos::sqrt(1.0_Real - 16.0_Real * Zeta);
   } else {
      // Strongly unstable (convective): (a_s - c_s*zeta)^{1/3}
      return Kokkos::pow(AMoS - CMoS * Zeta, 1.0_Real / 3.0_Real);
   }
}

/// @brief Hu(sigma) - Momentum surface value scaling
/// Sets surface boundary condition for momentum mixing
/// REFERENCES: Large et al. (1994)
///
/// @param Sigma Normalized vertical position
/// @return Normalized shape value
KOKKOS_INLINE_FUNCTION
Real kppSurfaceMomentumScale(Real Sigma) {
   // Linear decay from HuOn at the surface to zero at the OSBL base
   return HuOn * (1.0_Real + Sigma);
}

// ==========================================================================
// Langmuir Enhancement Factor (Theory-based Wave Model)
//
// Langmuir circulations are wind-and-wave driven counter-rotating vortices
// that deepen and strengthen boundary layer mixing beyond what wind stress
// alone produces. Their strength is measured by the turbulent Langmuir
// number La = sqrt(u* / u_stokes): small La means wave forcing dominates.
// With no wave model coupled, the Stokes drift is estimated from the 10 m
// wind following Li et al. (2016).
// ==========================================================================

/// @brief Estimate the surface-layer Stokes drift velocity scale from wind
/// REFERENCES: Li et al. 2016, cvmix_kpp_ustokes_SL_model
///
/// TODO: this is a placeholder empirical fit, not the full surface-layer
/// averaged Stokes drift of the reference model, which needs the wave
/// spectrum (or a wave component) to evaluate
///   u_s,BL = (U10/362) * sqrt(2*alpha*Cd) * lambda/HBL
/// with alpha = 0.84, Cd ~ 1.2e-3 and lambda ~ 2*pi*g/omega^2. Until those
/// inputs are available, HBL is accepted but unused.
///
/// @param Wind10m Wind speed at 10 m height (m/s)
/// @param HBL Boundary layer depth (m)
/// @return Stokes drift velocity scale (m/s)
KOKKOS_INLINE_FUNCTION
Real estimateStokesDriftSL(Real Wind10m, Real HBL) {
   Wind10m = Kokkos::fmax(0.0_Real, Wind10m);
   HBL     = Kokkos::fmax(1.0_Real, HBL);

   const Real UStokes = 0.016_Real * Wind10m;

   return Kokkos::fmax(0.0_Real, UStokes);
}

/// @brief Turbulent Langmuir number from friction velocity and Stokes drift
/// REFERENCES: Large et al. 2015 Eq. 6
///
/// @param UStar Friction velocity (m/s)
/// @param UStokes Stokes drift at surface (m/s)
/// @return Langmuir number (dimensionless)
KOKKOS_INLINE_FUNCTION
Real computeLangmuirNumber(Real UStar, Real UStokes) {
   UStar   = Kokkos::fmax(MinUStar, UStar);
   UStokes = Kokkos::fmax(MinStokesDrift, UStokes);

   return Kokkos::sqrt(UStar / UStokes);
}

/// @brief Langmuir enhancement factor applied to the KPP velocity scales
/// REFERENCES: Li et al. (2016) Eq. 1-3
///
/// @param Wind10m Wind speed at 10 m (m/s)
/// @param UStar Friction velocity (m/s)
/// @param HBL Boundary layer depth (m)
/// @return Enhancement factor R_L (dimensionless, > 1.0 enhances mixing)
KOKKOS_INLINE_FUNCTION
Real computeLangmuirEnhancement(Real Wind10m, Real UStar, Real HBL) {
   UStar   = Kokkos::fmax(MinUStar, UStar);
   Wind10m = Kokkos::fmax(0.0_Real, Wind10m);
   HBL     = Kokkos::fmax(1.0_Real, HBL);

   const Real UStokes = estimateStokesDriftSL(Wind10m, HBL);
   const Real La      = computeLangmuirNumber(UStar, UStokes);

   // R_L = sqrt(1 + 0.5/La^2); La is floored at 0.5 so that the weak-wave
   // limit returns the unenhanced scales rather than diverging.
   const Real LaInv = 1.0_Real / Kokkos::fmax(0.5_Real, La);
   const Real RL    = Kokkos::sqrt(1.0_Real + 0.5_Real * LaInv * LaInv);

   return Kokkos::fmin(2.0_Real, Kokkos::fmax(1.0_Real, RL));
}

// ==========================================================================
// Utility Functions for OSBL Depth Computation
// ==========================================================================

/// @brief Check if a point should be suppressed (e.g., under ice)
/// Sets OSBL to minimum if ice coverage or land ice present
///
/// @param IceFrac Sea ice coverage (0-1)
/// @param LandIceMask Land ice mask (0=ocean, non-zero=ice)
/// @return True if suppression applies
KOKKOS_INLINE_FUNCTION
bool shouldSuppressOSBL(Real IceFrac, I4 LandIceMask) {
   return (LandIceMask != 0) || (IceFrac > IceSuppressThresh);
}

/// @brief Apply OSBL depth constraints based on column properties
///
/// @param HOSBL Current OSBL depth (m)
/// @param LayerThickness Surface layer thickness (m)
/// @param WaterDepth Total water depth (m)
/// @param IceFrac Sea ice coverage (0-1)
/// @return Constrained OSBL depth (m)
KOKKOS_INLINE_FUNCTION
Real constrainOSBLDepth(Real HOSBL, Real LayerThickness, Real WaterDepth,
                        Real IceFrac) {
   // Lower bound: at least the full surface layer thickness
   HOSBL = Kokkos::fmax(HOSBL, LayerThickness);

   // Enforce minimum under ice
   if (IceFrac > IceSuppressThresh) {
      HOSBL = Kokkos::fmax(HOSBL, MinOSBLUnderIce);
   }

   // Upper bound: cannot exceed water depth
   HOSBL = Kokkos::fmin(HOSBL, WaterDepth * 0.95_Real);

   return HOSBL;
}

// ==========================================================================
// Turbulent Velocity Scale Computation
// ==========================================================================

/// @brief Compute the depth-independent turbulent velocity scale
/// Blends the shear-driven scale u* with the convective scale so that the
/// result stays finite in both the wind-driven and free-convection limits.
/// REFERENCES: Large et al. (1994) Eq. (9)-(10)
///
/// @param UStar Friction velocity (m/s)
/// @param BuoyFlux Surface buoyancy flux (m^2/s^3), negative when convective
/// @param HOBL Boundary layer depth (m)
/// @return Turbulent velocity scale w_s (m/s)
KOKKOS_INLINE_FUNCTION
Real computeTurbVelocityScale(Real UStar, Real BuoyFlux, Real HOBL) {
   UStar = Kokkos::fmax(0.0_Real, UStar);
   HOBL  = Kokkos::fmax(0.0_Real, HOBL);

   // Momentum contribution
   const Real WMom = UStar * UStar * UStar;

   // Buoyancy contribution for unstable (cooling/densifying) forcing.
   // In this sign convention, free convection corresponds to BuoyFlux < 0.
   const Real WBuoy =
       ConvectiveVelFac * Kokkos::fmax(0.0_Real, -BuoyFlux) * HOBL;

   return Kokkos::pow(Kokkos::fmax(0.0_Real, WMom + WBuoy),
                      1.0_Real / 3.0_Real);
}

/// @brief Momentum and scalar turbulent velocity scales at a point in the OSBL
///
/// CVMix-style scales: w = kappa * u* * phi^{-1}(zeta) in general, with
/// explicit free-convection limits when u* vanishes. Both scales are returned
/// together because they share the stability coordinate zeta.
/// REFERENCES: Large et al. (1994) Eq. (13), Appendix B
///
/// @param UStar Friction velocity (m/s)
/// @param BuoyFlux Surface buoyancy flux (m^2/s^3), negative when convective
/// @param HOBL Boundary layer depth (m)
/// @param SigmaLoc Normalized depth [0,1], capped at SurfaceLayerExtent
/// @param Kappa von Karman constant
/// @param WMTurb [out] Momentum turbulent velocity scale w_m (m/s)
/// @param WSTurb [out] Scalar turbulent velocity scale w_s (m/s)
KOKKOS_INLINE_FUNCTION
void kppTurbScales(Real UStar, Real BuoyFlux, Real HOBL, Real SigmaLoc,
                   Real Kappa, Real &WMTurb, Real &WSTurb) {
   WMTurb = 0.0_Real;
   WSTurb = 0.0_Real;

   if (UStar > 0.0_Real) {
      const Real U3 = UStar * UStar * UStar;
      const Real Zeta =
          SigmaLoc * HOBL * BuoyFlux * Kappa / Kokkos::max(U3, Real(Tiny));

      // These return phi^{-1}, which is positive over their whole domain;
      // do not invert or floor again.
      WMTurb = Kappa * UStar * kppPhiInvMomentum(Zeta);
      WSTurb = Kappa * UStar * kppPhiInvScalar(Zeta);
   } else if (BuoyFlux < 0.0_Real) {
      // Free-convection edge case (u*=0, unstable forcing).
      const Real WM3 = -CMoM * SigmaLoc * HOBL * Kappa * BuoyFlux;
      const Real WS3 = -CMoS * SigmaLoc * HOBL * Kappa * BuoyFlux;
      WMTurb =
          Kappa * Kokkos::pow(Kokkos::max(0.0_Real, WM3), 1.0_Real / 3.0_Real);
      WSTurb =
          Kappa * Kokkos::pow(Kokkos::max(0.0_Real, WS3), 1.0_Real / 3.0_Real);
   }
}

/// @brief Shape value the KPP profile must reach at the OSBL base so that it
/// joins the pre-existing interior coefficient there (MatchBoth only).
///
/// Callers must skip this when no interior mixing is supplied, since
/// InteriorCoeff is then read from an unallocated array.
///
/// @param InteriorCoeff Interior diffusivity or viscosity at the OSBL base
/// @param HOBL Boundary layer depth (m)
/// @param W Turbulent velocity scale matching InteriorCoeff (m/s)
/// @return Shape value at the OSBL base (dimensionless)
KOKKOS_INLINE_FUNCTION
Real kppMatchShape(Real InteriorCoeff, Real HOBL, Real W) {
   if (HOBL <= 0.0_Real || W <= 0.0_Real) {
      return 0.0_Real;
   }

   return InteriorCoeff / (HOBL * W);
}

/// @brief Non-local flux normalization constant
/// C_s = C* * kappa * (c_s * kappa * epsilon)^(1/3), with C* = 10.
/// Evaluates to roughly 6.33 with the default constants.
/// REFERENCES: Large et al. (1994) Eq. (20)
///
/// @param Kappa von Karman constant
/// @param SurfLayerExtent Surface layer extent epsilon (dimensionless)
/// @return C_s (dimensionless)
KOKKOS_INLINE_FUNCTION
Real kppNonLocalCs(Real Kappa, Real SurfLayerExtent) {
   return 10.0_Real * Kappa *
          Kokkos::pow(CMoS * Kappa * SurfLayerExtent, 1.0_Real / 3.0_Real);
}

/// @brief Clamp a trial OSBL depth to the range supported by the column
///
/// @param OSBLDepth Trial OSBL depth (m)
/// @param MinOSBLDepth Lower bound, typically the top layer thickness (m)
/// @param MaxOSBLDepth Upper bound, typically the deepest cell center (m)
/// @param ApplyIceMinimum Whether the sea-ice minimum depth applies
/// @param MinOSBLUnderIce Minimum OSBL depth under sea ice (m)
/// @return Clamped OSBL depth (m)
KOKKOS_INLINE_FUNCTION
Real kppClampOSBLDepth(Real OSBLDepth, Real MinOSBLDepth, Real MaxOSBLDepth,
                       bool ApplyIceMinimum, Real MinOSBLUnderIce) {
   OSBLDepth = Kokkos::fmax(OSBLDepth, MinOSBLDepth);

   if (ApplyIceMinimum) {
      OSBLDepth = Kokkos::fmax(OSBLDepth, MinOSBLUnderIce);
   }

   return Kokkos::fmin(OSBLDepth, MaxOSBLDepth);
}

/// @brief Index of the cell layer containing a given OSBL depth
///
/// @param ZInterface Geometric height of layer interfaces (m)
/// @param ICell Cell index
/// @param KMin Index of the top active layer
/// @param KMax Index of the bottom active layer
/// @param Ssh Sea surface height (m), since depths are measured below it
/// @param OSBLDepth OSBL depth (m)
/// @return Layer index bracketing OSBLDepth, or KMax if none does
KOKKOS_INLINE_FUNCTION
I4 kppOSBLIndex(const Array2DReal &ZInterface, I4 ICell, I4 KMin, I4 KMax,
                Real Ssh, Real OSBLDepth) {
   for (I4 K = KMin; K < KMax; ++K) {
      const Real ZAbove = Ssh - ZInterface(ICell, K);
      const Real ZBelow = Ssh - ZInterface(ICell, K + 1);
      if (OSBLDepth >= ZAbove && OSBLDepth <= ZBelow) {
         return K;
      }
   }

   return KMax;
}

} // namespace OMEGA::KPP

#endif // OMEGA_KPP_CONSTANTS_H
