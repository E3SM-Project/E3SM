#ifndef OMEGA_PGRAD_FINITE_VOLUME_H
#define OMEGA_PGRAD_FINITE_VOLUME_H
//===-- ocn/PGradFiniteVolume.h - Finite-Volume PGrad ---------*- C++ -*-===//
///
/// The equation-of-state expansion shared across an edge, used by the
/// FiniteVolume pressure gradient.
///
/// Specific volume is never integrated directly. For each cell and layer it is
/// expanded to first order about a reference state,
///
///    alpha(Ct, Sa, p) = alpha0 + alphaCt (Ct - Ct0) + alphaSa (Sa - Sa0)
///                              + alphaP  (p  - p0)
///
/// with the four coefficients coming from a single equation-of-state
/// evaluation per cell per layer -- the SpecVol, SpecVolDCt, SpecVolDSa and
/// SpecVolDP fields Eos already computes. Because temperature and salinity are
/// reconstructed as low-order polynomials in pressure and the expansion is
/// linear in them, the expanded specific volume is a low-order polynomial in
/// pressure and every integral the scheme needs is available in closed form.
/// No equation-of-state evaluation occurs inside any integral, which is what
/// bounds the cost.
///
/// The expansion point and the coefficients are shared across each edge:
/// averaged from the edge's two cells and used for *both* columns. This is
/// what the robustness property rests on. Give each column its own expansion
/// point and the two columns describe the same water with two slightly
/// different approximations to the same equation of state, their alpha0 and
/// alphaP terms no longer cancel in the matched-pressure difference, and the
/// scheme generates spurious flow out of nothing but its own equation-of-state
/// approximation.
///
/// Note what is and is not load-bearing here. That *one* set multiplies both
/// columns is essential. *Which* set it is -- selected by edge layer, or by
/// some other rule -- is an ordinary accuracy question and cannot break
/// exactness, because the coefficients end up multiplying a quantity that is
/// identically zero on the exact set.
//
//===----------------------------------------------------------------------===//

#include "DataTypes.h"
#include "OmegaKokkos.h"

namespace OMEGA {

/// Largest number of quadrature points supported within an edge layer
inline constexpr I4 MaxPGradQuadPoints = 4;

/// The equation-of-state expansion shared by both columns of an edge layer.
/// One instance multiplies both columns' contributions.
struct PGradEdgeEos {
   Real SpecVol0;    ///< specific volume at the shared expansion state
   Real SpecVolDCt;  ///< d(alpha)/d(ConservTemp) at the shared state
   Real SpecVolDSa;  ///< d(alpha)/d(AbsSalinity) at the shared state
   Real SpecVolDP;   ///< d(alpha)/d(Pressure) at the shared state
   Real ConservTemp; ///< shared expansion state, conservative temperature
   Real AbsSalinity; ///< shared expansion state, absolute salinity
   Real Press;       ///< shared expansion state, pressure
};

/// Build a shared expansion by averaging the coefficients and reference states
/// of layer K0 of the first cell and layer K1 of the second. Within an edge
/// layer the two are the same index, through the overload below. The sea-floor
/// anchor pairs the two layers that meet where the columns' floor heights are
/// compared, which need not share an index once the floors are more than a
/// layer apart.
///
/// Templated on the array type so that the same code is exercised by the
/// device kernels and by the host-side reference that tests them.
template <class ArrayType>
KOKKOS_INLINE_FUNCTION PGradEdgeEos
buildEdgeEos(const ArrayType &SpecVol,     ///< [in] specific volume
             const ArrayType &SpecVolDCt,  ///< [in] d(alpha)/d(ConservTemp)
             const ArrayType &SpecVolDSa,  ///< [in] d(alpha)/d(AbsSalinity)
             const ArrayType &SpecVolDP,   ///< [in] d(alpha)/d(Pressure)
             const ArrayType &ConservTemp, ///< [in] layer-mean temperature
             const ArrayType &AbsSalinity, ///< [in] layer-mean salinity
             const ArrayType &PressureMid, ///< [in] mid-layer pressure
             const I4 ICell0,              ///< [in] first cell on the edge
             const I4 K0,                  ///< [in] layer of the first cell
             const I4 ICell1,              ///< [in] second cell on the edge
             const I4 K1                   ///< [in] layer of the second cell
) {
   PGradEdgeEos Eos;
   Eos.SpecVol0 = 0.5_Real * (SpecVol(ICell0, K0) + SpecVol(ICell1, K1));
   Eos.SpecVolDCt =
       0.5_Real * (SpecVolDCt(ICell0, K0) + SpecVolDCt(ICell1, K1));
   Eos.SpecVolDSa =
       0.5_Real * (SpecVolDSa(ICell0, K0) + SpecVolDSa(ICell1, K1));
   Eos.SpecVolDP = 0.5_Real * (SpecVolDP(ICell0, K0) + SpecVolDP(ICell1, K1));
   Eos.ConservTemp =
       0.5_Real * (ConservTemp(ICell0, K0) + ConservTemp(ICell1, K1));
   Eos.AbsSalinity =
       0.5_Real * (AbsSalinity(ICell0, K0) + AbsSalinity(ICell1, K1));
   Eos.Press = 0.5_Real * (PressureMid(ICell0, K0) + PressureMid(ICell1, K1));
   return Eos;
}

/// Build the shared expansion for one edge layer by averaging the two
/// adjacent cells' coefficients and reference states. The pressure of the
/// shared state is the edge average of the two columns' mid-layer pressures,
/// which is exactly the mid-layer pressure of the edge control volume when
/// that volume is taken as the average of the two columns' interface
/// pressures.
template <class ArrayType>
KOKKOS_INLINE_FUNCTION PGradEdgeEos
buildEdgeEos(const ArrayType &SpecVol,     ///< [in] specific volume
             const ArrayType &SpecVolDCt,  ///< [in] d(alpha)/d(ConservTemp)
             const ArrayType &SpecVolDSa,  ///< [in] d(alpha)/d(AbsSalinity)
             const ArrayType &SpecVolDP,   ///< [in] d(alpha)/d(Pressure)
             const ArrayType &ConservTemp, ///< [in] layer-mean temperature
             const ArrayType &AbsSalinity, ///< [in] layer-mean salinity
             const ArrayType &PressureMid, ///< [in] mid-layer pressure
             const I4 ICell0,              ///< [in] first cell on the edge
             const I4 ICell1,              ///< [in] second cell on the edge
             const I4 K                    ///< [in] edge layer
) {
   return buildEdgeEos(SpecVol, SpecVolDCt, SpecVolDSa, SpecVolDP, ConservTemp,
                       AbsSalinity, PressureMid, ICell0, K, ICell1, K);
}

/// One column's expanded specific volume at a pressure, using the edge-shared
/// expansion. Needed only where a single column's specific volume is
/// integrated on its own, which is only the anchor of the column
/// scan; the interior of the scheme differences two columns at matched
/// pressure, where the SpecVol0 and SpecVolDP terms cancel and never have to
/// be formed.
KOKKOS_INLINE_FUNCTION Real edgeSpecVol(
    const PGradEdgeEos &Eos, ///< [in] edge-shared expansion
    const Real ConservTemp,  ///< [in] this column's temperature at p
    const Real AbsSalinity,  ///< [in] this column's salinity at p
    const Real Press         ///< [in] pressure
) {
   return Eos.SpecVol0 + Eos.SpecVolDCt * (ConservTemp - Eos.ConservTemp) +
          Eos.SpecVolDSa * (AbsSalinity - Eos.AbsSalinity) +
          Eos.SpecVolDP * (Press - Eos.Press);
}

/// The matched-pressure difference in specific volume between the two columns
/// of an edge, evaluated at one pressure with one shared expansion.
///
/// This is the central quantity of the scheme. Because both columns use the
/// same expansion, the SpecVol0 and SpecVolDP terms are identical in the two
/// and cancel in the difference, leaving a coefficient times the horizontal
/// contrast in reconstructed temperature and salinity *at matched pressure*.
/// That contrast is identically zero, pointwise, whenever the two columns'
/// reconstructions describe the same water -- which is why exactness depends
/// neither on the values of the coefficients, nor on the quadrature, nor on
/// the two columns' interfaces lining up.
///
/// Compressibility drops out entirely: SpecVolDP does not appear. That is
/// correct physics rather than an approximation. If temperature and salinity
/// are horizontally uniform then so is specific volume as a function of
/// pressure, and a horizontally uniform compressibility exerts no horizontal
/// pressure gradient.
KOKKOS_INLINE_FUNCTION Real matchedPressSpecVolDiff(
    const PGradEdgeEos &Eos, ///< [in] edge-shared expansion
    const Real ConservTemp0, ///< [in] first column's temperature at p
    const Real AbsSalinity0, ///< [in] first column's salinity at p
    const Real ConservTemp1, ///< [in] second column's temperature at p
    const Real AbsSalinity1  ///< [in] second column's salinity at p
) {
   return Eos.SpecVolDCt * (ConservTemp1 - ConservTemp0) +
          Eos.SpecVolDSa * (AbsSalinity1 - AbsSalinity0);
}

/// Gauss-Legendre nodes and weights on [-1, 1] for NPoints points, NPoints
/// running from 1 to MaxPGradQuadPoints.
///
/// The number of points is an accuracy setting only. The integrand is zero at
/// every point for any profile the reconstruction resolves exactly, so no rule
/// can break the robustness property; the choice trades cost against accuracy
/// off the exact set with nothing else at stake. Two points is exact for the
/// integrand within a sub-interval and is the default. More is worth
/// considering only where the two columns' interfaces are strongly offset,
/// since the integrand is piecewise linear with breakpoints at the union of
/// the two columns' interfaces and a fixed rule does not resolve those breaks.
KOKKOS_INLINE_FUNCTION void
gaussLegendreRule(const I4 NPoints,                 ///< [in] number of points
                  Real Nodes[MaxPGradQuadPoints],   ///< [out] nodes on [-1,1]
                  Real Weights[MaxPGradQuadPoints]) ///< [out] weights
{
   switch (NPoints) {
   case 1:
      Nodes[0]   = 0.0_Real;
      Weights[0] = 2.0_Real;
      break;
   case 3:
      Nodes[0]   = -0.77459666924148337704_Real;
      Nodes[1]   = 0.0_Real;
      Nodes[2]   = 0.77459666924148337704_Real;
      Weights[0] = 0.55555555555555555556_Real;
      Weights[1] = 0.88888888888888888889_Real;
      Weights[2] = 0.55555555555555555556_Real;
      break;
   case 4:
      Nodes[0]   = -0.86113631159405257522_Real;
      Nodes[1]   = -0.33998104358485626480_Real;
      Nodes[2]   = 0.33998104358485626480_Real;
      Nodes[3]   = 0.86113631159405257522_Real;
      Weights[0] = 0.34785484513745385737_Real;
      Weights[1] = 0.65214515486254614263_Real;
      Weights[2] = 0.65214515486254614263_Real;
      Weights[3] = 0.34785484513745385737_Real;
      break;
   default: // two points, the default
      Nodes[0]   = -0.57735026918962576451_Real;
      Nodes[1]   = 0.57735026918962576451_Real;
      Weights[0] = 1.0_Real;
      Weights[1] = 1.0_Real;
      break;
   }
}

} // namespace OMEGA

//===----------------------------------------------------------------------===//
#endif
