#ifndef OMEGA_PGRAD_RECON_H
#define OMEGA_PGRAD_RECON_H
//===-- ocn/PGradRecon.h - Pressure Gradient Reconstruction ---*- C++ -*-===//
///
/// Mean-preserving vertical reconstruction of conservative temperature and
/// absolute salinity in pressure, used by the FiniteVolume pressure gradient.
///
/// Within layer K of a column the prognostic layer mean is supplemented by a
/// deviation that integrates to zero over the layer,
///
///    Theta(p) = Theta_k + Slope_k * (p - PressureMid_k)
///
/// The deviations are linear, which makes the pressure gradient exact for
/// profiles that vary linearly with pressure. There is no limiter: the
/// reconstruction feeds an integral rather than an advective flux, and a
/// limiter active on smooth data would break the cancellation the scheme rests
/// on precisely where the profile is well resolved.
///
/// Two properties carry the weight.
///
/// Mean-preserving. PressureMid is the exact arithmetic midpoint of the two
/// interface pressures, so the integral of (p - PressureMid) over the layer is
/// zero and the deviation cannot move the layer mean. This is also what makes
/// VertCoord's midpoint rule the exact layer integral of a linear
/// reconstruction, and hence what spares VertCoord any change.
///
/// Exact on a non-uniform grid. For a profile linear in pressure the layer
/// means lie exactly on the line as a function of mid-layer pressure, because
/// the mean of a linear function over a layer is its value at the layer's
/// midpoint. A centered difference of the layer means with respect to
/// PressureMid therefore recovers the slope exactly, for any distribution of
/// layer thicknesses. This is why the estimator differences against the actual
/// mid-layer pressures and not against layer index: a formula that assumed
/// uniform thickness would pass a uniform-thickness test and fail on Omega's
/// grid.
///
/// The functions here take scalars so that they can be unit tested without a
/// mesh.
//
//===----------------------------------------------------------------------===//

#include "DataTypes.h"
#include "OmegaKokkos.h"

namespace OMEGA {

/// The two layers the linear slope estimator differences, for layer K of a
/// column whose valid layers run from KMin to KMax. The difference is centered
/// in the interior and one-sided in the shallowest and deepest valid layer.
/// The one-sided branches are not an edge case: the deepest valid layer is a
/// partial cell, which carries the whole signal where bathymetry steps. Where
/// a column has a single valid layer both indices are K, which the slope
/// estimator turns into a zero slope -- a constant being the only
/// mean-preserving reconstruction available there.
KOKKOS_INLINE_FUNCTION void
linearReconStencil(const I4 K,    ///< [in] layer to reconstruct
                   const I4 KMin, ///< [in] shallowest valid layer
                   const I4 KMax, ///< [in] deepest valid layer
                   I4 &KLo,       ///< [out] shallower layer to difference
                   I4 &KHi        ///< [out] deeper layer to difference
) {
   KLo = (K > KMin) ? K - 1 : K;
   KHi = (K < KMax) ? K + 1 : K;
}

/// The slope, per Pa, of the mean-preserving linear reconstruction, from the
/// layer means and mid-layer pressures of the two layers linearReconStencil
/// selects. Exact for a profile linear in pressure on an arbitrary grid and
/// second-order accurate otherwise. Returns zero where the two mid-layer
/// pressures coincide, which happens only when the column has a single valid
/// layer.
KOKKOS_INLINE_FUNCTION Real linearReconSlope(
    const Real ValueLo,    ///< [in] layer mean at KLo
    const Real ValueHi,    ///< [in] layer mean at KHi
    const Real PressMidLo, ///< [in] mid-layer pressure at KLo
    const Real PressMidHi  ///< [in] mid-layer pressure at KHi
) {
   const Real DeltaPress = PressMidHi - PressMidLo;
   return (DeltaPress != 0.0_Real) ? (ValueHi - ValueLo) / DeltaPress
                                   : 0.0_Real;
}

/// The deviation of the reconstruction from the layer mean at a pressure.
/// Integrates to zero over the layer because PressMid is the exact arithmetic
/// midpoint of the two interface pressures. The pressure need not lie within
/// the layer: near the sea floor a column's deepest reconstruction is
/// evaluated slightly below its own floor, where this extrapolates.
KOKKOS_INLINE_FUNCTION Real linearReconDeviation(
    const Real Slope,    ///< [in] slope from
                         ///< linearReconSlope
    const Real PressMid, ///< [in] mid-layer pressure
    const Real Press     ///< [in] pressure to evaluate at
) {
   return Slope * (Press - PressMid);
}

/// The reconstruction evaluated at a pressure. Equal to the layer mean at
/// PressMid, by construction.
KOKKOS_INLINE_FUNCTION Real linearReconEval(
    const Real Value,    ///< [in] prognostic layer mean
    const Real Slope,    ///< [in] slope from linearReconSlope
    const Real PressMid, ///< [in] mid-layer pressure
    const Real Press     ///< [in] pressure to evaluate at
) {
   return Value + linearReconDeviation(Slope, PressMid, Press);
}

/// The layer of column ICell whose interfaces bracket the pressure Press,
/// searched from the hint KHint and clamped to the column's valid range
/// [KMin, KMax].
///
/// This is the lookup that makes the scheme a fixed-pressure comparison rather
/// than a fixed-layer-index one, and it is the single place where that
/// distinction lives. At a pressure in edge layer K, each column supplies its
/// own temperature and salinity from whichever of *its own* layers contains
/// that pressure, which under tilt is generally not layer K -- at a coordinate
/// tilt of 50 m/km with 64 m layers the two columns' layer K are offset by
/// nearly three layer thicknesses and do not overlap in pressure at all.
/// Looking a column's state up by layer index instead would silently
/// reintroduce the very error the scheme exists to remove, and would still
/// pass every exactness and accuracy check available, which is why the
/// property tests on this function are not optional.
///
/// Interface pressures increase downward, so the search walks up while Press
/// lies above the layer's top interface and down while it lies below the
/// bottom one. Within the column scan the answer advances monotonically with
/// K, so passing the previous answer as KHint makes this a pair of incremented
/// cursors rather than a search; the result does not depend on the hint.
///
/// Where Press lies outside the column altogether the outermost valid layer is
/// returned and its reconstruction is extrapolated. This is the rule at the
/// top and bottom of the column: the edge control volume's pressure range is
/// in general neither column's own, so near the sea floor a column's deepest
/// reconstruction must be evaluated slightly below its own floor. Exactness
/// survives extrapolation, because on the exact set an extrapolated
/// reconstruction still reproduces the true profile.
///
/// Templated on the array type so that the same code is exercised by the
/// device kernels and by the host-side property tests that pin it.
template <class ArrayType>
KOKKOS_INLINE_FUNCTION I4
findLayerForPress(const ArrayType &PressInterface, ///< [in] interface pressures
                  const I4 ICell,                  ///< [in] cell to search
                  const I4 KMin,    ///< [in] shallowest valid layer
                  const I4 KMax,    ///< [in] deepest valid layer
                  const Real Press, ///< [in] pressure to locate
                  const I4 KHint    ///< [in] starting guess
) {
   I4 K = KHint;
   if (K < KMin)
      K = KMin;
   if (K > KMax)
      K = KMax;

   while (K > KMin && Press < PressInterface(ICell, K))
      --K;
   while (K < KMax && Press > PressInterface(ICell, K + 1))
      ++K;

   return K;
}

} // namespace OMEGA

//===----------------------------------------------------------------------===//
#endif
