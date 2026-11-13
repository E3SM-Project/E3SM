//===-- Test driver for OMEGA Pressure Gradient (PGrad) --------------*-
// C++-*-===/
//
/// \file
/// \brief Test driver for PressureGrad module
//
//===----------------------------------------------------------------------===/

#include "PGrad.h"

#include "DataTypes.h"
#include "Decomp.h"
#include "Dimension.h"
#include "Eos.h"
#include "Error.h"
#include "Field.h"
#include "GlobalConstants.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "IOStream.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanState.h"
#include "OmegaKokkos.h"
#include "PGrad.h"
#include "PGradFiniteVolume.h"
#include "PGradRecon.h"
#include "Pacer.h"
#include "TimeStepper.h"
#include "Tracers.h"
#include "VertCoord.h"
#include "mpi.h"
#include <limits>
#include <vector>

using namespace OMEGA;

void initPGradTest() {

   Error Err;
   int Err1;

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   // Initialize the Logging system
   initLogging(DefEnv);

   // Read default config if present
   Config("Omega");
   Config::readAll("omega.yml");

   // The default instance is the centered reference that the finite-volume
   // scheme is checked against, so select it here rather than relying on the
   // type in the default config. Note the explicit std::string: a string
   // literal would select the boolean overload of Config::set.
   Config *OmegaConfig = Config::getOmegaConfig();
   Config PGradConfig("PressureGrad");
   OmegaConfig->get(PGradConfig);
   PGradConfig.set("PressureGradType", std::string("Centered"));

   // First step of time stepper initialization needed for IOstream
   TimeStepper::init1();
   TimeStepper *DefStepper = TimeStepper::getDefault();
   Clock *ModelClock       = DefStepper->getClock();

   // Initialize the IO system
   IO::init(DefComm);

   // Create the default decomposition (initializes the decomposition)
   Decomp::init();

   // Initialize Field infrastructure
   Field::init(ModelClock);

   // Initialize IOStreams - this does not yet validate the contents
   // of each file, only creates streams from Config
   IOStream::init(ModelClock);

   // Initialize the default halo
   Err1 = Halo::init();
   if (Err1 != 0) {
      LOG_ERROR("PGrad: error initializing default halo");
      Err += Error(ErrorCode::Fail, "PGrad: error initializing default halo");
   }

   // Initialize the default mesh
   HorzMesh::init(ModelClock);

   // Initialize the default vertical coordinate
   VertCoord::init();

   // Initialize the equation of state
   Eos::init();

   // Initialize ocean state
   OceanState::init();

   // Initialize tracers
   Tracers::init();

   CHECK_ERROR_ABORT(Err, "PGrad: error during initialization");
}

// Check the mean-preserving linear reconstruction estimator on a single
// column, without a mesh.
//
// Layer means are sampled from a profile that is exactly linear in pressure,
// so the recovered slope must match the profile's slope to round-off in every
// layer, including the two that use a one-sided difference. The layer
// thicknesses are deliberately non-uniform, which is the case this test exists
// for: a formula that assumed uniform thickness passes the uniform control
// below and fails here. Both are run, so a failure says which of the two it
// is.
//
// Also asserted are the two properties the scheme leans on: that the
// reconstruction reproduces the layer mean at mid-layer pressure, and that its
// deviation integrates to zero over the layer.
int testReconstruction() {

   int Err = 0;

   const Real Eps = std::numeric_limits<Real>::epsilon();
   const Real Tol = 16.0_Real;

   // a profile linear in pressure: Theta(p) = Value0 + Slope0 * p
   const Real Value0 = 12.5_Real;
   const Real Slope0 = -1.5e-6_Real;

   // deliberately non-uniform thicknesses in Pa, then a uniform control
   const std::vector<std::vector<Real>> ThicknessSets = {
       {1.0e5, 3.0e5, 2.0e5, 7.0e5, 1.5e5, 9.0e5, 4.0e5, 2.5e5},
       {3.0e5, 3.0e5, 3.0e5, 3.0e5, 3.0e5, 3.0e5, 3.0e5, 3.0e5}};
   const std::vector<std::string> SetNames = {"non-uniform", "uniform"};

   for (int ISet = 0; ISet < ThicknessSets.size(); ++ISet) {

      const std::vector<Real> &DeltaPress = ThicknessSets[ISet];
      const I4 NLayers                    = DeltaPress.size();
      const I4 KMin                       = 0;
      const I4 KMax                       = NLayers - 1;

      // interface pressures accumulated from the surface, mid-layer pressures
      // as their exact arithmetic midpoint, and layer means sampled from the
      // profile at mid-layer pressure -- which is the exact layer average of a
      // linear profile
      std::vector<Real> PressInterface(NLayers + 1);
      std::vector<Real> PressMid(NLayers);
      std::vector<Real> Value(NLayers);
      PressInterface[0] = 0.0_Real;
      for (int K = 0; K < NLayers; ++K) {
         PressInterface[K + 1] = PressInterface[K] + DeltaPress[K];
         PressMid[K] = 0.5_Real * (PressInterface[K] + PressInterface[K + 1]);
         Value[K]    = Value0 + Slope0 * PressMid[K];
      }

      Real MaxSlopeErr = 0.0_Real;
      Real MaxMeanErr  = 0.0_Real;
      Real MaxIntegral = 0.0_Real;

      for (int K = 0; K < NLayers; ++K) {

         I4 KLo, KHi;
         linearReconStencil(K, KMin, KMax, KLo, KHi);
         const Real Slope = linearReconSlope(Value[KLo], Value[KHi],
                                             PressMid[KLo], PressMid[KHi]);

         // round-off in the difference is set by the size of the layer means,
         // divided by the pressure interval differenced over
         const Real ValueScale =
             std::max(std::abs(Value[KLo]), std::abs(Value[KHi]));
         const Real SlopeBound =
             Eps * ValueScale / std::abs(PressMid[KHi] - PressMid[KLo]);
         MaxSlopeErr =
             std::max(MaxSlopeErr, std::abs(Slope - Slope0) / SlopeBound);

         // the reconstruction reproduces the layer mean at mid-layer pressure
         const Real AtMid =
             linearReconEval(Value[K], Slope, PressMid[K], PressMid[K]);
         MaxMeanErr = std::max(MaxMeanErr,
                               std::abs(AtMid - Value[K]) / (Eps * ValueScale));

         // the deviation integrates to zero over the layer, by two-point
         // Gauss quadrature, which is exact for a linear integrand
         const Real HalfWidth = 0.5_Real * DeltaPress[K];
         const Real Offset    = HalfWidth / std::sqrt(3.0_Real);
         const Real DevLo =
             linearReconDeviation(Slope, PressMid[K], PressMid[K] - Offset);
         const Real DevHi =
             linearReconDeviation(Slope, PressMid[K], PressMid[K] + Offset);
         const Real Integral = HalfWidth * (DevLo + DevHi);
         const Real DevScale = HalfWidth * std::abs(Slope) * HalfWidth;
         if (DevScale > 0.0_Real)
            MaxIntegral =
                std::max(MaxIntegral, std::abs(Integral) / (Eps * DevScale));
      }

      LOG_INFO("PGradTest: reconstruction on {} layers: slope error {} eps, "
               "layer-mean error {} eps, layer integral {} eps",
               SetNames[ISet], MaxSlopeErr, MaxMeanErr, MaxIntegral);

      if (MaxSlopeErr > Tol) {
         LOG_ERROR("PGradTest: reconstruction slope FAIL on {} layers: {} eps",
                   SetNames[ISet], MaxSlopeErr);
         ++Err;
      }
      if (MaxMeanErr > Tol) {
         LOG_ERROR("PGradTest: reconstruction layer mean FAIL on {} layers: "
                   "{} eps",
                   SetNames[ISet], MaxMeanErr);
         ++Err;
      }
      if (MaxIntegral > Tol) {
         LOG_ERROR("PGradTest: reconstruction layer integral FAIL on {} "
                   "layers: {} eps",
                   SetNames[ISet], MaxIntegral);
         ++Err;
      }
   }

   // a column with a single valid layer has nothing to difference against, so
   // a constant is the only mean-preserving reconstruction available
   I4 KLo, KHi;
   linearReconStencil(3, 3, 3, KLo, KHi);
   const Real SingleSlope =
       linearReconSlope(Value0, Value0, 1.0e5_Real, 1.0e5_Real);
   if (KLo == 3 && KHi == 3 && SingleSlope == 0.0_Real) {
      LOG_INFO("PGradTest: reconstruction single-layer column PASS");
   } else {
      LOG_ERROR("PGradTest: reconstruction single-layer column FAIL");
      ++Err;
   }

   if (Err == 0)
      LOG_INFO("PGradTest: reconstruction estimator PASS");

   return Err;

} // end testReconstruction

// Pin the per-column pressure lookup with direct property tests, on a
// fabricated pair of columns and without a mesh.
//
// These tests are not optional, and they are the only thing standing between a
// correct implementation and a silently wrong one. No answer-level check
// anywhere in the test plan can distinguish looking a column's state up by
// pressure from looking it up by layer index: on a profile linear in pressure
// every layer's mean-preserving reconstruction is that same line, so looking
// up the wrong layer costs nothing, and on a curved profile the two rules
// differ by less than a factor of two. An implementation that indexes by layer
// therefore passes every exactness and accuracy gate specified.
//
// Two columns are built sharing a surface and a bottom pressure, with the
// second column's interfaces bulging away from the first in the middle of the
// column so that the two are strongly offset there -- the situation coordinate
// tilt produces. Four properties are asserted:
//
//   - the returned layer's interfaces actually bracket the pressure;
//   - the returned layer differs from the edge layer index by at least two
//     somewhere, so the lookup is demonstrably not an index lookup;
//   - the answer does not depend on the starting hint, so the incremented
//     cursors used in the column scan cannot disagree with a search;
//   - a pressure outside the column clamps to the outermost valid layer, which
//     is the rule where the edge control volume extends past a column's own
//     floor.
int testPressureLookup() {

   int Err = 0;

   const I4 NLayers = 32;
   const I4 KMin    = 0;
   const I4 KMax    = NLayers - 1;

   // Column 0 is uniform in pressure; column 1 has the same surface and
   // bottom pressure but bulges away from it in between, by up to three layer
   // thicknesses. The bulge is small enough that column 1's interfaces stay
   // monotonic.
   const Real DeltaP = 1.25e6_Real;
   const Real Bulge  = 6.0_Real * DeltaP;

   const Real Pi = 3.14159265358979323846_Real;

   HostArray2DReal PressInterface("PressInterface", 2, NLayers + 1);
   for (int K = 0; K <= NLayers; ++K) {
      const Real Uniform   = K * DeltaP;
      PressInterface(0, K) = Uniform;
      PressInterface(1, K) = Uniform + Bulge * std::sin(Pi * K / NLayers);
   }

   // check that column 1 stayed monotonic, or the test is not testing what it
   // means to
   for (int K = 0; K < NLayers; ++K) {
      if (PressInterface(1, K + 1) <= PressInterface(1, K)) {
         LOG_ERROR("PGradTest: pressure lookup setup FAIL: column 1 is not "
                   "monotonic at layer {}",
                   K);
         ++Err;
      }
   }

   I4 MaxOffset     = 0;
   int NBracketFail = 0;
   int NHintFail    = 0;

   for (int K = 0; K < NLayers; ++K) {

      // the edge layer's mid pressure, as the edge average of the two
      // columns' interface pressures
      const Real PressTop =
          0.5_Real * (PressInterface(0, K) + PressInterface(1, K));
      const Real PressBot =
          0.5_Real * (PressInterface(0, K + 1) + PressInterface(1, K + 1));
      const Real Press = 0.5_Real * (PressTop + PressBot);

      for (int ICell = 0; ICell < 2; ++ICell) {

         const I4 KFound =
             findLayerForPress(PressInterface, ICell, KMin, KMax, Press, K);

         // the interfaces of the returned layer must bracket the pressure,
         // unless the search clamped at an end of the column
         const bool Above =
             (KFound == KMin) && (Press < PressInterface(ICell, KMin));
         const bool Below =
             (KFound == KMax) && (Press > PressInterface(ICell, KMax + 1));
         const bool Brackets = Press >= PressInterface(ICell, KFound) &&
                               Press <= PressInterface(ICell, KFound + 1);
         if (!Brackets && !Above && !Below)
            ++NBracketFail;

         const I4 Offset = std::abs(KFound - K);
         if (Offset > MaxOffset)
            MaxOffset = Offset;

         // the answer must not depend on where the search started
         for (I4 KHint = KMin; KHint <= KMax; ++KHint) {
            if (findLayerForPress(PressInterface, ICell, KMin, KMax, Press,
                                  KHint) != KFound)
               ++NHintFail;
         }
      }
   }

   if (NBracketFail == 0) {
      LOG_INFO("PGradTest: pressure lookup bracketing PASS");
   } else {
      LOG_ERROR("PGradTest: pressure lookup bracketing FAIL: {} pressures not "
                "bracketed by the returned layer",
                NBracketFail);
      ++Err;
   }

   // Under tilt the layer containing a pressure is generally not the edge
   // layer index. If this never differed by much, the test configuration
   // would not be exercising the distinction the lookup exists for.
   if (MaxOffset >= 2) {
      LOG_INFO("PGradTest: pressure lookup is not an index lookup PASS: "
               "returned layer differs from the edge layer by up to {}",
               MaxOffset);
   } else {
      LOG_ERROR("PGradTest: pressure lookup is not an index lookup FAIL: "
                "returned layer never differs from the edge layer by more "
                "than {}",
                MaxOffset);
      ++Err;
   }

   if (NHintFail == 0) {
      LOG_INFO("PGradTest: pressure lookup hint independence PASS");
   } else {
      LOG_ERROR("PGradTest: pressure lookup hint independence FAIL: {} "
                "disagreements between starting hints",
                NHintFail);
      ++Err;
   }

   // Pressures outside the column clamp to the outermost valid layer, whose
   // reconstruction is then extrapolated
   const Real AboveTop = PressInterface(0, KMin) - DeltaP;
   const Real BelowBot = PressInterface(0, KMax + 1) + DeltaP;
   const I4 KAbove =
       findLayerForPress(PressInterface, 0, KMin, KMax, AboveTop, KMax);
   const I4 KBelow =
       findLayerForPress(PressInterface, 0, KMin, KMax, BelowBot, KMin);
   if (KAbove == KMin && KBelow == KMax) {
      LOG_INFO("PGradTest: pressure lookup out-of-column clamping PASS");
   } else {
      LOG_ERROR("PGradTest: pressure lookup out-of-column clamping FAIL: got "
                "{} above the column and {} below it",
                KAbove, KBelow);
      ++Err;
   }

   return Err;

} // end testPressureLookup

// The exact layer average over [PressLo, PressHi] of a profile that is
// quadratic in pressure, used to initialize layer means the way the exactness
// gate requires: as the exact layer averages of a prescribed continuous
// profile, which under tilt come out different in the two columns.
KOKKOS_INLINE_FUNCTION Real exactLayerAvg(
    const Real Coeff0,  ///< [in] constant term
    const Real Coeff1,  ///< [in] coefficient of p
    const Real Coeff2,  ///< [in] coefficient of p squared
    const Real Coeff3,  ///< [in] coefficient of p cubed
    const Real PressLo, ///< [in] top interface pressure
    const Real PressHi  ///< [in] bottom interface pressure
) {
   const Real PressMid = 0.5_Real * (PressLo + PressHi);
   const Real HalfSq   = 0.25_Real * (PressHi - PressLo) * (PressHi - PressLo);
   return Coeff0 + Coeff1 * PressMid +
          Coeff2 * (PressMid * PressMid + HalfSq / 3.0_Real) +
          Coeff3 * (PressMid * PressMid * PressMid + PressMid * HalfSq);
}

// A prescribed continuous profile of temperature and salinity in pressure, up
// to cubic. Linear is the exact set of the linear reconstruction; quadratic
// and cubic lie outside it, where the residual should shrink like the square
// of the layer thickness.
struct PGradProfile {
   Real Temp0, Temp1, Temp2, Temp3;
   Real Salt0, Salt1, Salt2, Salt3;
};

// Check the matched-pressure integrand on a fabricated column pair, without a
// mesh. This exercises the reconstruction, the pressure lookup and the
// edge-shared expansion together, which is the combination the whole scheme
// rests on.
//
// On a profile linear in pressure the integrand must be zero at *every*
// quadrature point, not merely in the integral. That is the sharpest available
// form of the check: it localizes a failure to the layer that caused it, and
// it is what a per-layer-index sharing rule fails. Both columns' layer means
// are the exact layer averages of one prescribed continuous profile, so under
// the offset between the two columns those means genuinely differ -- a
// configuration built by copying identical means into both columns would not
// exercise the property at all.
//
// A quadratic profile runs as a control. The reconstruction does not resolve
// it, so the integrand is nonzero there, and the gap between the two cases is
// what makes the linear result meaningful rather than a statement that
// everything in the test is zero.
int testMatchedPressIntegrand() {

   int Err = 0;

   const Real Eps = std::numeric_limits<Real>::epsilon();
   const Real Tol = 16.0_Real;

   const I4 NLayers = 32;
   const I4 KMin    = 0;
   const I4 KMax    = NLayers - 1;

   const Real DeltaP = 1.25e6_Real;
   const Real Bulge  = 6.0_Real * DeltaP;
   const Real Pi     = 3.14159265358979323846_Real;

   // representative shared edge coefficients; their values cannot affect the
   // exact-set result, since they multiply a quantity that is identically zero
   PGradEdgeEos EdgeEos;
   EdgeEos.SpecVol0    = 9.7e-4_Real;
   EdgeEos.SpecVolDCt  = -2.5e-7_Real;
   EdgeEos.SpecVolDSa  = -7.5e-7_Real;
   EdgeEos.SpecVolDP   = -4.0e-13_Real;
   EdgeEos.ConservTemp = 8.0_Real;
   EdgeEos.AbsSalinity = 34.9_Real;
   EdgeEos.Press       = 2.0e7_Real;

   // the prescribed continuous profile, linear then quadratic in pressure
   const Real Temp0 = 12.0_Real, Temp1 = -2.0e-7_Real;
   const Real Salt0 = 34.5_Real, Salt1 = 2.5e-8_Real;

   Real MaxLinear       = 0.0_Real;
   Real MaxQuadratic    = 0.0_Real;
   Real MaxAbsLinear    = 0.0_Real;
   Real MaxAbsQuadratic = 0.0_Real;

   for (int ICase = 0; ICase < 2; ++ICase) {

      const bool Curved = (ICase == 1);
      const Real Temp2  = Curved ? 3.0e-15_Real : 0.0_Real;
      const Real Salt2  = Curved ? -4.0e-16_Real : 0.0_Real;

      // Column 0 is uniform in pressure; column 1 shares its surface and
      // bottom pressure but is offset from it by up to three layer
      // thicknesses in between, as coordinate tilt would make it.
      HostArray2DReal PressInterface("PressInterface", 2, NLayers + 1);
      HostArray2DReal PressMid("PressMid", 2, NLayers);
      HostArray2DReal ConservTemp("ConservTemp", 2, NLayers);
      HostArray2DReal AbsSalinity("AbsSalinity", 2, NLayers);
      HostArray2DReal SlopeCt("SlopeCt", 2, NLayers);
      HostArray2DReal SlopeSa("SlopeSa", 2, NLayers);

      for (int ICell = 0; ICell < 2; ++ICell) {
         for (int K = 0; K <= NLayers; ++K) {
            const Real Uniform = K * DeltaP;
            PressInterface(ICell, K) =
                (ICell == 1) ? Uniform + Bulge * std::sin(Pi * K / NLayers)
                             : Uniform;
         }
         for (int K = 0; K < NLayers; ++K) {
            PressMid(ICell, K)    = 0.5_Real * (PressInterface(ICell, K) +
                                             PressInterface(ICell, K + 1));
            ConservTemp(ICell, K) = exactLayerAvg(Temp0, Temp1, Temp2, 0.0_Real,
                                                  PressInterface(ICell, K),
                                                  PressInterface(ICell, K + 1));
            AbsSalinity(ICell, K) = exactLayerAvg(Salt0, Salt1, Salt2, 0.0_Real,
                                                  PressInterface(ICell, K),
                                                  PressInterface(ICell, K + 1));
         }
         for (int K = 0; K < NLayers; ++K) {
            I4 KLo, KHi;
            linearReconStencil(K, KMin, KMax, KLo, KHi);
            SlopeCt(ICell, K) = linearReconSlope(
                ConservTemp(ICell, KLo), ConservTemp(ICell, KHi),
                PressMid(ICell, KLo), PressMid(ICell, KHi));
            SlopeSa(ICell, K) = linearReconSlope(
                AbsSalinity(ICell, KLo), AbsSalinity(ICell, KHi),
                PressMid(ICell, KLo), PressMid(ICell, KHi));
         }
      }

      Real Nodes[MaxPGradQuadPoints];
      Real Weights[MaxPGradQuadPoints];
      const I4 NQuad = 2;
      gaussLegendreRule(NQuad, Nodes, Weights);

      Real MaxScaled = 0.0_Real;
      Real MaxAbs    = 0.0_Real;

      for (int K = 0; K < NLayers; ++K) {

         // the edge layer spans the average of the two columns' interfaces
         const Real EdgeTop =
             0.5_Real * (PressInterface(0, K) + PressInterface(1, K));
         const Real EdgeBot =
             0.5_Real * (PressInterface(0, K + 1) + PressInterface(1, K + 1));
         const Real EdgeMid  = 0.5_Real * (EdgeTop + EdgeBot);
         const Real EdgeHalf = 0.5_Real * (EdgeBot - EdgeTop);

         for (int IQuad = 0; IQuad < NQuad; ++IQuad) {

            const Real Press = EdgeMid + EdgeHalf * Nodes[IQuad];

            Real TempAt[2];
            Real SaltAt[2];
            for (int ICell = 0; ICell < 2; ++ICell) {
               // each column's own layer containing this pressure, which
               // under the offset is generally not layer K
               const I4 KFound = findLayerForPress(PressInterface, ICell, KMin,
                                                   KMax, Press, K);
               TempAt[ICell]   = linearReconEval(ConservTemp(ICell, KFound),
                                                 SlopeCt(ICell, KFound),
                                                 PressMid(ICell, KFound), Press);
               SaltAt[ICell]   = linearReconEval(AbsSalinity(ICell, KFound),
                                                 SlopeSa(ICell, KFound),
                                                 PressMid(ICell, KFound), Press);
            }

            const Real SpecVolDiff = matchedPressSpecVolDiff(
                EdgeEos, TempAt[0], SaltAt[0], TempAt[1], SaltAt[1]);

            // scale by the size of the terms that had to cancel, not by an
            // absolute constant
            const Real Scale = std::abs(EdgeEos.SpecVolDCt * TempAt[0]) +
                               std::abs(EdgeEos.SpecVolDSa * SaltAt[0]);

            MaxAbs = std::max(MaxAbs, std::abs(SpecVolDiff));
            MaxScaled =
                std::max(MaxScaled, std::abs(SpecVolDiff) / (Eps * Scale));
         }
      }

      LOG_INFO("PGradTest: matched-pressure integrand, {} profile: max "
               "|dSpecVol| = {} m3/kg = {} eps of the cancelling terms",
               Curved ? "quadratic" : "linear", MaxAbs, MaxScaled);

      if (Curved) {
         MaxQuadratic    = MaxScaled;
         MaxAbsQuadratic = MaxAbs;
      } else {
         MaxLinear    = MaxScaled;
         MaxAbsLinear = MaxAbs;
      }
   }

   if (MaxLinear <= Tol) {
      LOG_INFO("PGradTest: matched-pressure integrand PASS: {} eps on a "
               "profile linear in pressure",
               MaxLinear);
   } else {
      LOG_ERROR("PGradTest: matched-pressure integrand FAIL: {} eps on a "
                "profile linear in pressure exceeds {} eps",
                MaxLinear, Tol);
      ++Err;
   }

   // The linear result means nothing unless the same machinery gives a large
   // answer on a profile the reconstruction does not resolve. The comparison
   // is a ratio of the two absolute values rather than a count of epsilons,
   // so that it means the same thing in a single-precision build, where the
   // linear result rises with the round-off floor while the quadratic one, a
   // truncation error, does not.
   const Real ControlRatio =
       (MaxAbsLinear > 0.0_Real) ? MaxAbsQuadratic / MaxAbsLinear : 0.0_Real;
   if (ControlRatio > 100.0_Real) {
      LOG_INFO("PGradTest: matched-pressure integrand control PASS: the "
               "quadratic profile is {} times the linear one ({} eps)",
               ControlRatio, MaxQuadratic);
   } else {
      LOG_ERROR("PGradTest: matched-pressure integrand control FAIL: the "
                "quadratic profile is only {} times the linear one, so the "
                "linear result may be trivial",
                ControlRatio);
      ++Err;
   }

   return Err;

} // end testMatchedPressIntegrand

// Build the two-column test state. Every edge joins cells 0 and 1 and is
// DC apart. Layer interfaces are staggered between the two columns by
// TiltFactor, which tilts them relative to surfaces of constant pressure;
// TiltFactor = 0.5 gives two identical columns and no tilt. Temperature and
// salinity vary linearly with geometric height. Pseudo-thickness, specific
// volume and pressure are iterated to consistency, after which the VertCoord
// pressure and geometric height fields are filled.
void setupTwoColumnState(HorzMesh *Mesh,    ///< [in] Horizontal mesh
                         VertCoord *VCoord, ///< [inout] Vertical coordinate
                         OceanState *State, ///< [inout] Ocean state
                         Eos *EqState,      ///< [inout] Equation of state
                         I4 NVertLayers,    ///< [in] Number of layers
                         Real DC,           ///< [in] Distance between cells
                         Real TiltFactor,   ///< [in] Interface tilt, 0.5 = none
                         Real ZBottom       ///< [in] Depth of the sea floor
) {

   const I4 NCellsAll = Mesh->NCellsAll;
   const I4 NEdgesAll = Mesh->NEdgesAll;

   VCoord->NVertLayers   = NVertLayers;
   VCoord->NVertLayersP1 = NVertLayers + 1;

   Array2DReal SpecVolOld("SpecVolOld", Mesh->NCellsSize, NVertLayers);
   Array1DReal SurfacePressure("SurfacePressure", Mesh->NCellsSize);

   auto &MinLayerCell = VCoord->MinLayerCell;
   auto &MaxLayerCell = VCoord->MaxLayerCell;
   parallelFor(
       {NCellsAll}, KOKKOS_LAMBDA(int i) {
          MinLayerCell(i) = 0;
          MaxLayerCell(i) = NVertLayers - 1;
       });

   auto &MinLayerEdgeBot = VCoord->MinLayerEdgeBot;
   auto &MaxLayerEdgeTop = VCoord->MaxLayerEdgeTop;
   parallelFor(
       {NEdgesAll}, KOKKOS_LAMBDA(int i) {
          MinLayerEdgeBot(i) = 0;
          MaxLayerEdgeTop(i) = NVertLayers - 1;
       });

   auto &CellsOnEdge = Mesh->CellsOnEdge;
   auto &DcEdge      = Mesh->DcEdge;
   parallelFor(
       {NEdgesAll}, KOKKOS_LAMBDA(int i) {
          CellsOnEdge(i, 0) = 0;
          CellsOnEdge(i, 1) = 1;
          DcEdge(i)         = DC;
       });

   // Fetch reference density from Config
   const Real Density0 = RhoSw;

   const I4 TimeLevel = 0;

   // get state and tracer arrays
   Array2DReal PseudoThick = State->getPseudoThickness(TimeLevel);
   Array2DReal Temp        = Tracers::getByName(TimeLevel, "Temperature");
   Array2DReal Salinity    = Tracers::getByName(TimeLevel, "Salinity");

   // set Z interface and mid-point locations
   const Real DZ         = 2.0_Real * (-ZBottom / NVertLayers);
   auto &BottomGeomDepth = VCoord->BottomGeomDepth;
   auto &GeomZInterface  = VCoord->GeomZInterface;
   auto &GeomZMid        = VCoord->GeomZMid;
   parallelFor(
       {NCellsAll}, KOKKOS_LAMBDA(int i) {
          GeomZInterface(i, NVertLayers) = ZBottom;
          SurfacePressure(i)             = 0.0_Real;
          BottomGeomDepth(i)             = 0.0_Real;
          for (int k = NVertLayers - 1; k >= 0; --k) {
             Real X = (k + i) % 2;
             Real Dz =
                 (2.0_Real * TiltFactor - 1.0_Real) * X * DZ +
                 (1.0_Real - TiltFactor) * DZ; // staggered pseudo-thickness
             GeomZInterface(i, k) = GeomZInterface(i, k + 1) + Dz;
             PseudoThick(i, k) =
                 GeomZInterface(i, k) - GeomZInterface(i, k + 1);
             GeomZMid(i, k) =
                 0.5_Real * (GeomZInterface(i, k) + GeomZInterface(i, k + 1));
             BottomGeomDepth(i) += Dz;
          }
       });

   // set simple temperature and salinity profiles
   auto &SpecVol = EqState->SpecVol;
   parallelFor(
       {NCellsAll, NVertLayers}, KOKKOS_LAMBDA(int i, int k) {
          Real T0 = 30.0;
          Real TB = 5.0;
          Real S0 = 30.0;
          Real SB = 40.0;

          Real Phi0 = (GeomZMid(i, k) - ZBottom) / (-ZBottom);
          Real PhiB = 1.0_Real - Phi0;

          Temp(i, k)       = T0 * Phi0 + TB * PhiB;
          Salinity(i, k)   = S0 * Phi0 + SB * PhiB;
          SpecVol(i, k)    = 1.0_Real / Density0;
          SpecVolOld(i, k) = SpecVol(i, k);
       });

   // Iterate to converge PseudoThick, SpecVol, PressureMid
   auto &PressureMid = VCoord->PressureMid;
   VCoord->computePressure(PseudoThick, SurfacePressure);
   for (int Iteration = 0; Iteration < 15; ++Iteration) {

      // compute specific volume from EOS. This fills SpecVol only, not the
      // derivatives the FiniteVolume scheme needs; a test running that scheme
      // on this state must fill them itself.
      VCoord->computePressure(PseudoThick, SurfacePressure);
      EqState->computeSpecVol(Temp, Salinity, PressureMid);

      // compute psuedo thickness from specific volume
      parallelFor(
          {NCellsAll, NVertLayers}, KOKKOS_LAMBDA(int i, int k) {
             PseudoThick(i, k) =
                 1.0_Real / (SpecVol(i, k) * Density0) *
                 (GeomZInterface(i, k) - GeomZInterface(i, k + 1));
          });

      // compute difference from previous iteration
      Real MaxValue = 0.0_Real;
      parallelReduce(
          {NCellsAll, NVertLayers},
          KOKKOS_LAMBDA(int i, int k, Real &max) {
             Real Diff = Kokkos::abs(SpecVol(i, k) - SpecVolOld(i, k));
             if (Diff > max)
                max = Diff;
          },
          Kokkos::Max<Real>(MaxValue));

      // check convergence
      if (MaxValue < 1e-12_Real) {
         LOG_INFO("converged: max diff = {}", MaxValue);
         break;
      } else {
         parallelFor(
             {NCellsAll, NVertLayers},
             KOKKOS_LAMBDA(int i, int k) { SpecVolOld(i, k) = SpecVol(i, k); });
      }
   }

   // compute pressure once more with converged PseudoThick
   VCoord->computePressure(PseudoThick, SurfacePressure);

   // compute z levels
   VCoord->computeGeomZHeight(PseudoThick, SpecVol);

} // end setupTwoColumnState

// Build a two-column state whose vertical profile is a prescribed continuous
// function of pressure, with each column's layer means set to the *exact layer
// averages* of that one profile over that column's own layers.
//
// Under the offset between the two columns those averages come out different
// in the two, and that is the point. A configuration built by copying identical
// layer means into both columns would not exercise the property at all, because
// the condition that matters is a property of the reconstructed profiles rather
// than of the layer means.
//
// The two columns sit on a flat floor and differ only in how pressure is
// distributed between their interfaces, which is what an ALE coordinate does.
// SurfPressDiff offsets the second column's surface pressure, and hence its
// bottom pressure; the exactness gate must be run with it zero, since where
// surface pressure varies the state carries a real fixed-pressure height
// difference and zero is the wrong expectation. It is nonzero only for the
// anchor guard, which needs the two columns' end pressures to differ.
//
// Returns the pressure at the bottom of the columns, so tests can scale their
// tolerances by the hydrostatic terms.
Real setupProfileState(HorzMesh *Mesh,    ///< [in] Horizontal mesh
                       VertCoord *VCoord, ///< [inout] Vertical coordinate
                       OceanState *State, ///< [inout] Ocean state
                       Eos *EqState,      ///< [inout] Equation of state
                       I4 NVertLayers,    ///< [in] Number of layers
                       Real DC,           ///< [in] Distance between cells
                       Real BulgePress,   ///< [in] interface offset, in Pa
                       PGradProfile Prof, ///< [in] prescribed profile
                       Real SurfPressDiff ///< [in] surface pressure contrast
) {

   const I4 NCellsAll = Mesh->NCellsAll;
   const I4 NEdgesAll = Mesh->NEdgesAll;

   VCoord->NVertLayers   = NVertLayers;
   VCoord->NVertLayersP1 = NVertLayers + 1;

   const Real PressBot = 4.0e7_Real;
   const Real Pi       = 3.14159265358979323846_Real;
   const Real DeltaP   = PressBot / NVertLayers;

   auto &MinLayerCell = VCoord->MinLayerCell;
   auto &MaxLayerCell = VCoord->MaxLayerCell;
   parallelFor(
       {NCellsAll}, KOKKOS_LAMBDA(int i) {
          MinLayerCell(i) = 0;
          MaxLayerCell(i) = NVertLayers - 1;
       });

   auto &MinLayerEdgeBot = VCoord->MinLayerEdgeBot;
   auto &MaxLayerEdgeTop = VCoord->MaxLayerEdgeTop;
   parallelFor(
       {NEdgesAll}, KOKKOS_LAMBDA(int i) {
          MinLayerEdgeBot(i) = 0;
          MaxLayerEdgeTop(i) = NVertLayers - 1;
       });

   auto &CellsOnEdge = Mesh->CellsOnEdge;
   auto &DcEdge      = Mesh->DcEdge;
   parallelFor(
       {NEdgesAll}, KOKKOS_LAMBDA(int i) {
          CellsOnEdge(i, 0) = 0;
          CellsOnEdge(i, 1) = 1;
          DcEdge(i)         = DC;
       });

   Array1DReal SurfacePressure("SurfacePressure", Mesh->NCellsSize);
   Array2DReal PseudoThick = State->getPseudoThickness(0);
   Array2DReal Temp        = Tracers::getByName(0, "Temperature");
   Array2DReal Salinity    = Tracers::getByName(0, "Salinity");

   // Column 0 is uniform in pressure; column 1 bulges away from it in the
   // interior while spanning the same pressure range. Pseudo-thickness and
   // pressure thickness are proportional, so setting one sets the other. The
   // bulge stays below PressBot/pi, which keeps column 1's interfaces
   // monotonic.
   auto &BottomGeomDepth = VCoord->BottomGeomDepth;
   parallelFor(
       {NCellsAll}, KOKKOS_LAMBDA(int i) {
          SurfacePressure(i) = (i == 1) ? SurfPressDiff : 0.0_Real;
          // a flat floor: at the sea-floor anchor the height difference is
          // then exact input and vanishes identically
          BottomGeomDepth(i) = 4000.0_Real;
          for (int k = 0; k < NVertLayers; ++k) {
             const Real Shape = (i == 1) ? 1.0_Real : 0.0_Real;
             const Real PTop =
                 k * DeltaP +
                 Shape * BulgePress * Kokkos::sin(Pi * k / NVertLayers);
             const Real PBot =
                 (k + 1) * DeltaP +
                 Shape * BulgePress * Kokkos::sin(Pi * (k + 1) / NVertLayers);
             PseudoThick(i, k) = (PBot - PTop) / (Gravity * RhoSw);
          }
       });

   // pressure from the pseudo-thickness the model will actually use
   VCoord->computePressure(PseudoThick, SurfacePressure);

   // layer means as the exact layer averages of the prescribed profile, taken
   // over the model's own interface pressures so the two are consistent
   const auto &PressureInterface = VCoord->PressureInterface;
   parallelFor(
       {NCellsAll, NVertLayers}, KOKKOS_LAMBDA(int i, int k) {
          Temp(i, k)     = exactLayerAvg(Prof.Temp0, Prof.Temp1, Prof.Temp2,
                                         Prof.Temp3, PressureInterface(i, k),
                                         PressureInterface(i, k + 1));
          Salinity(i, k) = exactLayerAvg(Prof.Salt0, Prof.Salt1, Prof.Salt2,
                                         Prof.Salt3, PressureInterface(i, k),
                                         PressureInterface(i, k + 1));
       });

   // specific volume and its derivatives, from one equation-of-state pass
   EqState->computeSpecVolAndDerivs(Temp, Salinity, VCoord->PressureMid);

   VCoord->computeGeomZHeight(PseudoThick, EqState->SpecVol);

   return PressBot;

} // end setupProfileState

// Build a two-column resting state on terrain-following layers over a sloping
// floor. Both columns hold the same prescribed profile in the same number of
// layers, each spanning its whole column, and the second column's floor is
// FloorPressDiff deeper in pressure. Each column's bathymetry is the height its
// own layers stack to under VertCoord's rule, so the sea surface VertCoord
// derives is flat and the state is at rest.
//
// Where FloorPressDiff exceeds a layer thickness the two floors are more than a
// layer apart in pressure, which is what terrain-following layers over a steep
// slope produce and what the anchor's comparison against the deeper column's
// own interface is for.
void setupSlopingFloorState(HorzMesh *Mesh,    ///< [in] Horizontal mesh
                            VertCoord *VCoord, ///< [inout] Vertical coord
                            OceanState *State, ///< [inout] Ocean state
                            Eos *EqState,      ///< [inout] Equation of state
                            I4 NVertLayers,    ///< [in] Number of layers
                            Real DC,           ///< [in] Distance between cells
                            Real FloorPressDiff, ///< [in] floor offset, in Pa
                            PGradProfile Prof    ///< [in] prescribed profile
) {

   const I4 NCellsAll = Mesh->NCellsAll;
   const I4 NEdgesAll = Mesh->NEdgesAll;

   VCoord->NVertLayers   = NVertLayers;
   VCoord->NVertLayersP1 = NVertLayers + 1;

   const Real PressBot = 4.0e7_Real;

   auto &MinLayerCell = VCoord->MinLayerCell;
   auto &MaxLayerCell = VCoord->MaxLayerCell;
   parallelFor(
       {NCellsAll}, KOKKOS_LAMBDA(int i) {
          MinLayerCell(i) = 0;
          MaxLayerCell(i) = NVertLayers - 1;
       });

   auto &MinLayerEdgeBot = VCoord->MinLayerEdgeBot;
   auto &MaxLayerEdgeTop = VCoord->MaxLayerEdgeTop;
   parallelFor(
       {NEdgesAll}, KOKKOS_LAMBDA(int i) {
          MinLayerEdgeBot(i) = 0;
          MaxLayerEdgeTop(i) = NVertLayers - 1;
       });

   auto &CellsOnEdge = Mesh->CellsOnEdge;
   auto &DcEdge      = Mesh->DcEdge;
   parallelFor(
       {NEdgesAll}, KOKKOS_LAMBDA(int i) {
          CellsOnEdge(i, 0) = 0;
          CellsOnEdge(i, 1) = 1;
          DcEdge(i)         = DC;
       });

   Array1DReal SurfacePressure("SurfacePressure", Mesh->NCellsSize);
   Array2DReal PseudoThick = State->getPseudoThickness(0);
   Array2DReal Temp        = Tracers::getByName(0, "Temperature");
   Array2DReal Salinity    = Tracers::getByName(0, "Salinity");

   // equal layers spanning each column, the second column's deeper
   parallelFor(
       {NCellsAll}, KOKKOS_LAMBDA(int i) {
          SurfacePressure(i) = 0.0_Real;
          const Real ColPress =
              PressBot + ((i == 1) ? FloorPressDiff : 0.0_Real);
          for (int k = 0; k < NVertLayers; ++k)
             PseudoThick(i, k) = ColPress / NVertLayers / (Gravity * RhoSw);
       });

   VCoord->computePressure(PseudoThick, SurfacePressure);

   const auto &PressureInterface = VCoord->PressureInterface;
   parallelFor(
       {NCellsAll, NVertLayers}, KOKKOS_LAMBDA(int i, int k) {
          Temp(i, k)     = exactLayerAvg(Prof.Temp0, Prof.Temp1, Prof.Temp2,
                                         Prof.Temp3, PressureInterface(i, k),
                                         PressureInterface(i, k + 1));
          Salinity(i, k) = exactLayerAvg(Prof.Salt0, Prof.Salt1, Prof.Salt2,
                                         Prof.Salt3, PressureInterface(i, k),
                                         PressureInterface(i, k + 1));
       });

   EqState->computeSpecVolAndDerivs(Temp, Salinity, VCoord->PressureMid);

   // the bathymetry each column's own layers stack to, so that the sea surface
   // VertCoord derives from it is flat
   auto &BottomGeomDepth = VCoord->BottomGeomDepth;
   const auto &SpecVol   = EqState->SpecVol;
   parallelFor(
       {NCellsAll}, KOKKOS_LAMBDA(int i) {
          Real Depth = 0.0_Real;
          for (int k = NVertLayers - 1; k >= 0; --k)
             Depth += RhoSw * SpecVol(i, k) * PseudoThick(i, k);
          BottomGeomDepth(i) = Depth;
       });

   VCoord->computeGeomZHeight(PseudoThick, EqState->SpecVol);

} // end setupSlopingFloorState

// Which of the scheme's rules the reference assembly below applies. Turning
// one off is how the guards of design section 5.2, and the sloping-floor
// guard, are made to fire.
struct PGradGuardRules {
   bool SharedExpansion = true; ///< one EOS expansion per edge layer, shared
                                ///< by both columns, rather than one per cell
   bool PressureLookup = true;  ///< each column's state from the layer
                                ///< containing the pressure, not from layer K
   bool ShiftedAnchor = true;   ///< the anchor shifted to a common pressure,
                                ///< rather than the raw height difference
   bool AnchorWithinOneLayer = true; ///< the shallower floor compared with the
                                     ///< deeper column's own interface just
                                     ///< below it, rather than both floors
                                     ///< shifted to their mean pressure
};

// Assemble the column scan and the tendency on the host for the single edge
// joining cells 0 and 1, out of the same helper functions the kernel uses, with
// each of the scheme's rules switchable.
//
// With all three rules on this must reproduce the kernel to round-off, and that
// is checked before any guard is trusted. A guard that cannot fire is worse
// than no guard, because it looks like protection.
void referenceScan(const HostArray2DReal &PressInterface, ///< [in]
                   const HostArray2DReal &PressMid,       ///< [in]
                   const HostArray2DReal &Temp,           ///< [in]
                   const HostArray2DReal &Salt,           ///< [in]
                   const HostArray2DReal &SpecVol,        ///< [in]
                   const HostArray2DReal &SpecVolDCt,     ///< [in]
                   const HostArray2DReal &SpecVolDSa,     ///< [in]
                   const HostArray2DReal &SpecVolDP,      ///< [in]
                   const HostArray2DReal &GeomZ,          ///< [in]
                   const I4 NLayers,                      ///< [in]
                   const Real Dc,                         ///< [in]
                   const I4 NQuad,                        ///< [in]
                   const PGradGuardRules &Rules,          ///< [in]
                   std::vector<Real> &DeltaZ,             ///< [out] NLayers+1
                   std::vector<Real> &Tend                ///< [out] NLayers
) {

   const I4 KMin = 0;
   const I4 KMax = NLayers - 1;

   DeltaZ.assign(NLayers + 1, 0.0_Real);
   Tend.assign(NLayers, 0.0_Real);

   // the per-cell reconstruction slopes
   std::vector<std::vector<Real>> SlopeCt(2, std::vector<Real>(NLayers));
   std::vector<std::vector<Real>> SlopeSa(2, std::vector<Real>(NLayers));
   for (int ICell = 0; ICell < 2; ++ICell) {
      for (int K = 0; K < NLayers; ++K) {
         I4 KLo, KHi;
         linearReconStencil(K, KMin, KMax, KLo, KHi);
         SlopeCt[ICell][K] =
             linearReconSlope(Temp(ICell, KLo), Temp(ICell, KHi),
                              PressMid(ICell, KLo), PressMid(ICell, KHi));
         SlopeSa[ICell][K] =
             linearReconSlope(Salt(ICell, KLo), Salt(ICell, KHi),
                              PressMid(ICell, KLo), PressMid(ICell, KHi));
      }
   }

   Real Nodes[MaxPGradQuadPoints];
   Real Weights[MaxPGradQuadPoints];
   gaussLegendreRule(NQuad, Nodes, Weights);

   const Real InvGravity = 1.0_Real / Gravity;

   // one expansion per cell, used when the shared-expansion rule is off
   auto cellEos = [&](int ICell, int K) {
      PGradEdgeEos Eos;
      Eos.SpecVol0    = SpecVol(ICell, K);
      Eos.SpecVolDCt  = SpecVolDCt(ICell, K);
      Eos.SpecVolDSa  = SpecVolDSa(ICell, K);
      Eos.SpecVolDP   = SpecVolDP(ICell, K);
      Eos.ConservTemp = Temp(ICell, K);
      Eos.AbsSalinity = Salt(ICell, K);
      Eos.Press       = PressMid(ICell, K);
      return Eos;
   };

   auto lookup = [&](int ICell, Real Press, I4 K) {
      return Rules.PressureLookup ? findLayerForPress(PressInterface, ICell,
                                                      KMin, KMax, Press, K)
                                  : K;
   };

   std::vector<Real> Incr(NLayers, 0.0_Real);
   std::vector<Real> Moment(NLayers, 0.0_Real);

   for (int K = 0; K < NLayers; ++K) {

      const PGradEdgeEos EdgeEos =
          buildEdgeEos(SpecVol, SpecVolDCt, SpecVolDSa, SpecVolDP, Temp, Salt,
                       PressMid, 0, 1, K);

      const Real EdgeTop =
          0.5_Real * (PressInterface(0, K) + PressInterface(1, K));
      const Real EdgeBot =
          0.5_Real * (PressInterface(0, K + 1) + PressInterface(1, K + 1));
      const Real EdgeMid  = 0.5_Real * (EdgeTop + EdgeBot);
      const Real EdgeHalf = 0.5_Real * (EdgeBot - EdgeTop);

      for (int IQuad = 0; IQuad < NQuad; ++IQuad) {
         const Real Press  = EdgeMid + EdgeHalf * Nodes[IQuad];
         const Real Weight = EdgeHalf * Weights[IQuad];

         const I4 KF0 = lookup(0, Press, K);
         const I4 KF1 = lookup(1, Press, K);

         const Real T0 = linearReconEval(Temp(0, KF0), SlopeCt[0][KF0],
                                         PressMid(0, KF0), Press);
         const Real S0 = linearReconEval(Salt(0, KF0), SlopeSa[0][KF0],
                                         PressMid(0, KF0), Press);
         const Real T1 = linearReconEval(Temp(1, KF1), SlopeCt[1][KF1],
                                         PressMid(1, KF1), Press);
         const Real S1 = linearReconEval(Salt(1, KF1), SlopeSa[1][KF1],
                                         PressMid(1, KF1), Press);

         Real Diff;
         if (Rules.SharedExpansion) {
            Diff = matchedPressSpecVolDiff(EdgeEos, T0, S0, T1, S1);
         } else {
            // two expansion points mean the SpecVol0 and SpecVolDP terms no
            // longer cancel, and the full specific volume must be formed for
            // each column separately
            Diff = edgeSpecVol(cellEos(1, KF1), T1, S1, Press) -
                   edgeSpecVol(cellEos(0, KF0), T0, S0, Press);
         }

         Incr[K] += Weight * Diff;
         Moment[K] += Weight * (Press - EdgeTop) * Diff;
      }
      Incr[K] *= InvGravity;
      Moment[K] *= InvGravity;
   }

   // the anchor, at the sea floor: the shallower floor against the deeper
   // column's own interface at or just below that floor's pressure
   const Real FloorPress0 = PressInterface(0, KMax + 1);
   const Real FloorPress1 = PressInterface(1, KMax + 1);
   const bool Shallow0    = FloorPress0 <= FloorPress1;
   I4 KAnchor[2]          = {KMax, KMax};
   if (Rules.AnchorWithinOneLayer)
      KAnchor[Shallow0 ? 1 : 0] =
          findLayerForPress(PressInterface, Shallow0 ? 1 : 0, KMin, KMax,
                            Shallow0 ? FloorPress0 : FloorPress1, KMax);

   Real Anchor = GeomZ(1, KAnchor[1] + 1) - GeomZ(0, KAnchor[0] + 1);
   if (Rules.ShiftedAnchor) {
      const Real AnchorPress = 0.5_Real * (PressInterface(0, KAnchor[0] + 1) +
                                           PressInterface(1, KAnchor[1] + 1));
      const PGradEdgeEos AnchorEos =
          buildEdgeEos(SpecVol, SpecVolDCt, SpecVolDSa, SpecVolDP, Temp, Salt,
                       PressMid, 0, KAnchor[0], 1, KAnchor[1]);
      for (int ISide = 0; ISide < 2; ++ISide) {
         const int ICell     = ISide;
         const Real Sign     = (ISide == 0) ? -1.0_Real : 1.0_Real;
         const Real ColPress = PressInterface(ICell, KAnchor[ICell] + 1);
         const Real Mid      = 0.5_Real * (AnchorPress + ColPress);
         const Real Half     = 0.5_Real * (ColPress - AnchorPress);
         Real Integral       = 0.0_Real;
         for (int IQuad = 0; IQuad < NQuad; ++IQuad) {
            const Real Press  = Mid + Half * Nodes[IQuad];
            const Real Weight = Half * Weights[IQuad];
            const I4 KF       = lookup(ICell, Press, KAnchor[ICell]);
            const Real T = linearReconEval(Temp(ICell, KF), SlopeCt[ICell][KF],
                                           PressMid(ICell, KF), Press);
            const Real S = linearReconEval(Salt(ICell, KF), SlopeSa[ICell][KF],
                                           PressMid(ICell, KF), Press);
            Integral += Weight * edgeSpecVol(AnchorEos, T, S, Press);
         }
         Anchor += Sign * InvGravity * Integral;
      }

      // carried down to the edge's deepest interface where it sits above it
      if (KAnchor[0] != KMax || KAnchor[1] != KMax) {
         const PGradEdgeEos FloorEos =
             buildEdgeEos(SpecVol, SpecVolDCt, SpecVolDSa, SpecVolDP, Temp,
                          Salt, PressMid, 0, 1, KMax);
         const Real FloorPress = 0.5_Real * (FloorPress0 + FloorPress1);
         const Real Mid        = 0.5_Real * (AnchorPress + FloorPress);
         const Real Half       = 0.5_Real * (FloorPress - AnchorPress);
         Real Integral         = 0.0_Real;
         for (int IQuad = 0; IQuad < NQuad; ++IQuad) {
            const Real Press  = Mid + Half * Nodes[IQuad];
            const Real Weight = Half * Weights[IQuad];
            const I4 KF0      = lookup(0, Press, KMax);
            const I4 KF1      = lookup(1, Press, KMax);
            const Real T0     = linearReconEval(Temp(0, KF0), SlopeCt[0][KF0],
                                                PressMid(0, KF0), Press);
            const Real S0     = linearReconEval(Salt(0, KF0), SlopeSa[0][KF0],
                                                PressMid(0, KF0), Press);
            const Real T1     = linearReconEval(Temp(1, KF1), SlopeCt[1][KF1],
                                                PressMid(1, KF1), Press);
            const Real S1     = linearReconEval(Salt(1, KF1), SlopeSa[1][KF1],
                                                PressMid(1, KF1), Press);
            Integral +=
                Weight * matchedPressSpecVolDiff(FloorEos, T0, S0, T1, S1);
         }
         Anchor -= InvGravity * Integral;
      }
   }

   DeltaZ[KMax + 1] = Anchor;
   for (int K = KMax; K >= KMin; --K)
      DeltaZ[K] = DeltaZ[K + 1] + Incr[K];

   for (int K = 0; K < NLayers; ++K) {
      const Real DeltaPress =
          0.5_Real * ((PressInterface(0, K + 1) - PressInterface(0, K)) +
                      (PressInterface(1, K + 1) - PressInterface(1, K)));
      const Real LayerMean = DeltaZ[K + 1] + Moment[K] / DeltaPress;
      Tend[K]              = -Gravity / Dc * LayerMean;
   }

} // end referenceScan

// Copy the state the reference assembly needs to the host.
struct PGradHostState {
   HostArray2DReal PressInterface, PressMid, Temp, Salt;
   HostArray2DReal SpecVol, SpecVolDCt, SpecVolDSa, SpecVolDP, GeomZ;
};

PGradHostState copyStateToHost(VertCoord *VCoord, Eos *EqState) {
   PGradHostState H;
   H.PressInterface = createHostMirrorCopy(VCoord->PressureInterface);
   H.PressMid       = createHostMirrorCopy(VCoord->PressureMid);
   H.Temp       = createHostMirrorCopy(Tracers::getByName(0, "Temperature"));
   H.Salt       = createHostMirrorCopy(Tracers::getByName(0, "Salinity"));
   H.SpecVol    = createHostMirrorCopy(EqState->SpecVol);
   H.SpecVolDCt = createHostMirrorCopy(EqState->SpecVolDCt);
   H.SpecVolDSa = createHostMirrorCopy(EqState->SpecVolDSa);
   H.SpecVolDP  = createHostMirrorCopy(EqState->SpecVolDP);
   H.GeomZ      = createHostMirrorCopy(VCoord->GeomZInterface);
   return H;
}

// Run the FiniteVolume scheme on the current state and return the largest
// tendency magnitude, along with its root-mean-square over the active edges
// and layers and the largest fixed-pressure height difference the column scan
// produced.
//
// The exactness gate uses the maximum, since it asserts a zero at every edge
// and layer. The convergence rates use the RMS: a maximum is set by whichever
// layer happens to be worst, and that layer moves under refinement, which
// makes the measured rate noisy and non-monotone. RMS is also the norm the
// design's measured rates were taken in, so the two are comparable.
Real runFiniteVolume(HorzMesh *Mesh, VertCoord *VCoord, OceanState *State,
                     Eos *EqState, I4 NVertLayers, Real &MaxDeltaZ,
                     Real &MaxDeltaZUpper, Real &MaxDeltaZLower,
                     Real &RmsTend) {

   Config *Options = Config::getOmegaConfig();
   Config PGradConfig("PressureGrad");
   Options->get(PGradConfig);
   PGradConfig.set("PressureGradType", std::string("FiniteVolume"));
   PressureGrad *FVPGrad =
       PressureGrad::create("TestFV", Mesh, VCoord, Options);
   PGradConfig.set("PressureGradType", std::string("Centered"));
   // create returns null if an instance of this name already exists, and
   // dereferencing that gives a segmentation fault rather than a message
   CHECK_ERROR_ABORT(
       Error(FVPGrad ? ErrorCode::Success : ErrorCode::Fail, ""),
       "PGradTest: could not create the TestFV pressure gradient");

   Array2DReal Tend("TendFV", Mesh->NEdgesSize, NVertLayers);
   deepCopy(Tend, 0.0_Real);

   Array2DReal PseudoThick = State->getPseudoThickness(0);
   Array2DReal Temp        = Tracers::getByName(0, "Temperature");
   Array2DReal Salinity    = Tracers::getByName(0, "Salinity");

   FVPGrad->computePressureGrad(
       Tend, VCoord->PressureMid, VCoord->PressureInterface, EqState->SpecVol,
       VCoord->GeomZInterface, PseudoThick, Temp, Salinity, EqState);

   const auto &EdgeMask = VCoord->EdgeMask;
   Real MaxTend         = 0.0_Real;
   Real SumSq           = 0.0_Real;
   I4 NActive           = 0;
   parallelReduce(
       {Mesh->NEdgesAll, NVertLayers},
       KOKKOS_LAMBDA(int IEdge, int K, Real &MaxV, Real &LSum, I4 &LCount) {
          if (EdgeMask(IEdge, K) <= 0.0_Real)
             return;
          const Real V = Kokkos::abs(Tend(IEdge, K));
          if (V > MaxV)
             MaxV = V;
          LSum += V * V;
          ++LCount;
       },
       Kokkos::Max<Real>(MaxTend), Kokkos::Sum<Real>(SumSq),
       Kokkos::Sum<I4>(NActive));

   RmsTend = (NActive > 0) ? std::sqrt(SumSq / NActive) : 0.0_Real;

   const auto &DeltaZFixedP = FVPGrad->getDeltaZFixedP();
   MaxDeltaZ                = 0.0_Real;
   MaxDeltaZUpper           = 0.0_Real;
   MaxDeltaZLower           = 0.0_Real;
   parallelReduce(
       {Mesh->NEdgesAll, NVertLayers + 1},
       KOKKOS_LAMBDA(int IEdge, int K, Real &MaxA, Real &MaxU, Real &MaxL) {
          const Real V = Kokkos::abs(DeltaZFixedP(IEdge, K));
          if (V > MaxA)
             MaxA = V;
          if (K <= NVertLayers / 2 && V > MaxU)
             MaxU = V;
          if (K > NVertLayers / 2 && V > MaxL)
             MaxL = V;
       },
       Kokkos::Max<Real>(MaxDeltaZ), Kokkos::Max<Real>(MaxDeltaZUpper),
       Kokkos::Max<Real>(MaxDeltaZLower));

   PressureGrad::erase("TestFV");

   return MaxTend;

} // end runFiniteVolume

// Run PressureGradCentered on the current state and return the largest
// tendency magnitude.
Real runCentered(HorzMesh *Mesh, VertCoord *VCoord, OceanState *State,
                 Eos *EqState, I4 NVertLayers) {

   Array2DReal Tend("TendCentered", Mesh->NEdgesSize, NVertLayers);
   deepCopy(Tend, 0.0_Real);

   Array2DReal PseudoThick = State->getPseudoThickness(0);
   Array2DReal Temp        = Tracers::getByName(0, "Temperature");
   Array2DReal Salinity    = Tracers::getByName(0, "Salinity");

   PressureGrad::getDefault()->computePressureGrad(
       Tend, VCoord->PressureMid, VCoord->PressureInterface, EqState->SpecVol,
       VCoord->GeomZInterface, PseudoThick, Temp, Salinity, EqState);

   const auto &EdgeMask = VCoord->EdgeMask;
   Real MaxTend         = 0.0_Real;
   parallelReduce(
       {Mesh->NEdgesAll, NVertLayers},
       KOKKOS_LAMBDA(int IEdge, int K, Real &MaxV) {
          if (EdgeMask(IEdge, K) <= 0.0_Real)
             return;
          const Real V = Kokkos::abs(Tend(IEdge, K));
          if (V > MaxV)
             MaxV = V;
       },
       Kokkos::Max<Real>(MaxTend));

   return MaxTend;

} // end runCentered

// The gating test of design section 5.2: exactness on the exact set, the
// convergence of the residual off it, and the guards.
int testExactnessAndGuards(HorzMesh *Mesh,    ///< [in] Horizontal mesh
                           VertCoord *VCoord, ///< [inout] Vertical coordinate
                           OceanState *State, ///< [inout] Ocean state
                           Eos *EqState       ///< [inout] Equation of state
) {

   int Err = 0;

   const Real Eps   = std::numeric_limits<Real>::epsilon();
   const Real DC    = 4000.0_Real;
   const Real Bulge = 7.5e6_Real;
   const I4 NLayers = 32;
   const I4 NQuad   = 2;

   // linear in pressure: the exact set of the reconstruction
   const PGradProfile Linear = {12.0_Real, -2.0e-7_Real, 0.0_Real, 0.0_Real,
                                34.0_Real, 2.5e-8_Real,  0.0_Real, 0.0_Real};
   // quadratic and cubic: outside it
   const PGradProfile Quadratic = {12.0_Real,     -2.0e-7_Real, 3.0e-15_Real,
                                   0.0_Real,      34.0_Real,    2.5e-8_Real,
                                   -4.0e-16_Real, 0.0_Real};
   const PGradProfile Cubic     = {12.0_Real,     -2.0e-7_Real, 3.0e-15_Real,
                                   -5.0e-23_Real, 34.0_Real,    2.5e-8_Real,
                                   -4.0e-16_Real, 8.0e-24_Real};

   //
   // Group one: the exact set. The tendency must be zero to machine precision
   // at every edge and layer, for any tilt or thickness, over the flat floor
   // setupProfileState builds.
   //
   const Real PressBot = setupProfileState(
       Mesh, VCoord, State, EqState, NLayers, DC, Bulge, Linear, 0.0_Real);

   // the size of the terms that had to cancel: the hydrostatic contribution of
   // the temperature and salinity structure over the column
   const auto &SpecVolDCt = EqState->SpecVolDCt;
   const auto &SpecVolDSa = EqState->SpecVolDSa;
   Array2DReal Temp       = Tracers::getByName(0, "Temperature");
   Array2DReal Salinity   = Tracers::getByName(0, "Salinity");
   Real CancelScale       = 0.0_Real;
   parallelReduce(
       {2, NLayers},
       KOKKOS_LAMBDA(int i, int k, Real &MaxV) {
          const Real V = (Kokkos::abs(SpecVolDCt(i, k) * Temp(i, k)) +
                          Kokkos::abs(SpecVolDSa(i, k) * Salinity(i, k))) *
                         PressBot / Gravity;
          if (V > MaxV)
             MaxV = V;
       },
       Kokkos::Max<Real>(CancelScale));

   Real MaxDeltaZ, MaxUpper, MaxLower, RmsFV;
   const Real MaxFV = runFiniteVolume(Mesh, VCoord, State, EqState, NLayers,
                                      MaxDeltaZ, MaxUpper, MaxLower, RmsFV);
   const Real MaxCentered = runCentered(Mesh, VCoord, State, EqState, NLayers);

   // The tendency scale is the height scale divided by the cell spacing and
   // multiplied by gravity, so that the gate tracks Real's epsilon and the
   // size of the hydrostatic terms rather than a physical tolerance.
   const Real HeightScale = CancelScale;
   const Real TendScale   = Gravity * CancelScale / DC;

   const Real ScaledDeltaZ = MaxDeltaZ / (Eps * HeightScale);
   const Real ScaledTend   = MaxFV / (Eps * TendScale);

   LOG_INFO("PGradTest: exact set ({} precision): FiniteVolume max |Tend| = "
            "{} m/s2 ({} eps), Centered max |Tend| = {} m/s2, ratio {}",
            (Eps < 1.0e-10_Real) ? "double" : "single", MaxFV, ScaledTend,
            MaxCentered, MaxFV / MaxCentered);
   LOG_INFO("PGradTest: exact set: max |DeltaZFixedP| = {} m ({} eps); upper "
            "column {} m, lower column {} m -- the residual grows away from "
            "the sea-floor anchor",
            MaxDeltaZ, ScaledDeltaZ, MaxUpper, MaxLower);

   // Gates set from the measured values: 0.027 eps for both in double
   // precision, so four eps leaves well over a hundredfold margin while still
   // catching drift. A single-precision build has fewer guard digits in the
   // accumulation down the column, so the same measurement in epsilons is
   // expected to be somewhat larger; the gate is loosened to match rather
   // than being turned into a different kind of claim.
   const bool DoublePrec = Eps < 1.0e-10_Real;
   const Real ExactTol   = DoublePrec ? 4.0_Real : 16.0_Real;
   if (ScaledTend <= ExactTol && ScaledDeltaZ <= ExactTol) {
      LOG_INFO("PGradTest: exactness gate PASS");
   } else {
      LOG_ERROR("PGradTest: exactness gate FAIL: tendency {} eps, height "
                "difference {} eps, against {} eps",
                ScaledTend, ScaledDeltaZ, ExactTol);
      ++Err;
   }

   // Guard (a), tilt sensitivity. This is the only guard that can fire on a
   // configuration where every other check passes: a bug that zeroed the tilt
   // response would satisfy all of them perfectly.
   Real PrevCentered  = 0.0_Real;
   bool CenteredGrows = true;
   for (int ITilt = 1; ITilt <= 3; ++ITilt) {
      setupProfileState(Mesh, VCoord, State, EqState, NLayers, DC,
                        Bulge * ITilt / 3.0_Real, Linear, 0.0_Real);
      const Real C = runCentered(Mesh, VCoord, State, EqState, NLayers);
      if (C <= PrevCentered)
         CenteredGrows = false;
      PrevCentered = C;
   }
   if (CenteredGrows && MaxCentered > 0.0_Real && MaxFV < MaxCentered) {
      LOG_INFO("PGradTest: guard (a) tilt sensitivity PASS: Centered grows "
               "with tilt and the two schemes differ");
   } else {
      LOG_ERROR("PGradTest: guard (a) tilt sensitivity FAIL");
      ++Err;
   }

   // Everything from here on measures a truncation error rather than a
   // cancellation, and a single-precision build cannot resolve it: the
   // round-off floor of the tendency there is about 3e-8 m/s2, which swamps
   // the quadratic profile's residual at 60 layers and leaves only a factor
   // of ten between guard (b)'s signal and the noise. The exactness gate
   // above is the check design section 5.2 asks to be run in both
   // precisions, since it is the one that tests whether the scheme forms
   // large quantities; the rest is a statement about the discretization and
   // is precision-independent, so measuring it once in double precision is
   // enough.
   if (!DoublePrec) {
      LOG_INFO("PGradTest: single precision: skipping the convergence rates "
               "and guards (b) to (d), which measure truncation errors below "
               "this build's round-off floor");
      return Err;
   }

   //
   // The guards below need an assembly whose rules can be switched, so they
   // run against a host reference built from the same helper functions. It is
   // only trustworthy if it reproduces the kernel first.
   //
   setupProfileState(Mesh, VCoord, State, EqState, NLayers, DC, Bulge, Linear,
                     0.0_Real);
   Real KernelDeltaZ, KernelUpper, KernelLower, KernelRms;
   const Real KernelTend =
       runFiniteVolume(Mesh, VCoord, State, EqState, NLayers, KernelDeltaZ,
                       KernelUpper, KernelLower, KernelRms);

   PGradHostState H = copyStateToHost(VCoord, EqState);
   std::vector<Real> RefDeltaZ, RefTend;
   PGradGuardRules Rules;
   referenceScan(H.PressInterface, H.PressMid, H.Temp, H.Salt, H.SpecVol,
                 H.SpecVolDCt, H.SpecVolDSa, H.SpecVolDP, H.GeomZ, NLayers, DC,
                 NQuad, Rules, RefDeltaZ, RefTend);

   Real RefMax = 0.0_Real;
   for (int K = 0; K < NLayers; ++K)
      RefMax = std::max(RefMax, std::abs(RefTend[K]));

   if (std::abs(RefMax - KernelTend) <= 8.0_Real * Eps * TendScale) {
      LOG_INFO("PGradTest: guard harness fidelity on the exact set PASS: "
               "reference {} m/s2 against kernel {} m/s2",
               RefMax, KernelTend);
   } else {
      LOG_ERROR("PGradTest: guard harness fidelity on the exact set FAIL: "
                "reference {} m/s2 against kernel {} m/s2; the guards below "
                "say nothing",
                RefMax, KernelTend);
      ++Err;
   }

   // Agreeing at 1e-18 on the exact set says little, since a harness that
   // computed nothing at all would also agree. The check that has content is
   // on a curved profile, where both sides are large.
   setupProfileState(Mesh, VCoord, State, EqState, NLayers, DC, Bulge,
                     Quadratic, 0.0_Real);
   Real CurvedDeltaZ, CurvedUpper, CurvedLower, CurvedRms;
   const Real CurvedKernel =
       runFiniteVolume(Mesh, VCoord, State, EqState, NLayers, CurvedDeltaZ,
                       CurvedUpper, CurvedLower, CurvedRms);

   PGradHostState HC = copyStateToHost(VCoord, EqState);
   std::vector<Real> CurvedRefDeltaZ, CurvedRefTend;
   Rules = PGradGuardRules();
   referenceScan(HC.PressInterface, HC.PressMid, HC.Temp, HC.Salt, HC.SpecVol,
                 HC.SpecVolDCt, HC.SpecVolDSa, HC.SpecVolDP, HC.GeomZ, NLayers,
                 DC, NQuad, Rules, CurvedRefDeltaZ, CurvedRefTend);
   Real CurvedRefMax = 0.0_Real;
   for (int K = 0; K < NLayers; ++K)
      CurvedRefMax = std::max(CurvedRefMax, std::abs(CurvedRefTend[K]));

   const Real CurvedRelDiff =
       (CurvedKernel > 0.0_Real)
           ? std::abs(CurvedRefMax - CurvedKernel) / CurvedKernel
           : 1.0_Real;
   LOG_INFO("PGradTest: guard harness fidelity on a curved profile: reference "
            "{} m/s2 against kernel {} m/s2, relative difference {}",
            CurvedRefMax, CurvedKernel, CurvedRelDiff);
   if (CurvedRelDiff <= 1.0e-10_Real) {
      LOG_INFO("PGradTest: guard harness fidelity on a curved profile PASS");
   } else {
      LOG_ERROR("PGradTest: guard harness fidelity on a curved profile FAIL: "
                "the reference assembly does not reproduce the kernel, so the "
                "guards below say nothing");
      ++Err;
   }

   // back to the exact-set state for the guards
   setupProfileState(Mesh, VCoord, State, EqState, NLayers, DC, Bulge, Linear,
                     0.0_Real);
   H = copyStateToHost(VCoord, EqState);

   // Guard (b), a cell-local expansion point in place of the edge-shared one.
   // Two expansion points mean the SpecVol0 and SpecVolDP terms no longer
   // cancel, so this must fire.
   Rules                 = PGradGuardRules();
   Rules.SharedExpansion = false;
   referenceScan(H.PressInterface, H.PressMid, H.Temp, H.Salt, H.SpecVol,
                 H.SpecVolDCt, H.SpecVolDSa, H.SpecVolDP, H.GeomZ, NLayers, DC,
                 NQuad, Rules, RefDeltaZ, RefTend);
   Real GuardB = 0.0_Real;
   for (int K = 0; K < NLayers; ++K)
      GuardB = std::max(GuardB, std::abs(RefTend[K]));

   LOG_INFO("PGradTest: guard (b) cell-local expansion point gives {} m/s2, "
            "against {} m/s2 with the shared point",
            GuardB, RefMax);
   if (GuardB > 1.0e4_Real * Eps * TendScale) {
      LOG_INFO("PGradTest: guard (b) PASS: the edge-shared expansion point is "
               "load-bearing");
   } else {
      LOG_ERROR("PGradTest: guard (b) FAIL: replacing the shared expansion "
                "point with a cell-local one changed nothing, so the sharing "
                "is not being exercised");
      ++Err;
   }

   // Guard (c), each column's state taken from its own layer K rather than
   // from the layer containing the pressure. This *cannot* fire on the exact
   // set: the profile is a single line in pressure, so every layer's
   // reconstruction is that same line and looking up the wrong layer costs
   // nothing. The check here records that, because the alternative -- writing
   // it as a guard that must fire -- would be asserting something false, and
   // because it is the reason the property tests on the lookup exist at all.
   Rules                = PGradGuardRules();
   Rules.PressureLookup = false;
   referenceScan(H.PressInterface, H.PressMid, H.Temp, H.Salt, H.SpecVol,
                 H.SpecVolDCt, H.SpecVolDSa, H.SpecVolDP, H.GeomZ, NLayers, DC,
                 NQuad, Rules, RefDeltaZ, RefTend);
   Real GuardC = 0.0_Real;
   for (int K = 0; K < NLayers; ++K)
      GuardC = std::max(GuardC, std::abs(RefTend[K]));

   LOG_INFO("PGradTest: guard (c) layer-index lookup gives {} m/s2, against {} "
            "m/s2 with the pressure lookup -- this guard cannot fire on the "
            "exact set, which is why the lookup is pinned by property tests",
            GuardC, RefMax);
   if (GuardC <= ExactTol * Eps * TendScale) {
      LOG_INFO("PGradTest: guard (c) behaves as documented: no answer-level "
               "check can distinguish the two lookups here");
   } else {
      LOG_INFO("PGradTest: guard (c) unexpectedly fired at {} m/s2; an "
               "answer-level check on the lookup may now be available",
               GuardC);
   }

   // Guard (d), the anchor taken as the raw height difference, dropping the
   // short integrals that shift it to a common pressure. It fires only where
   // the two columns' end pressures differ, and it is flat with depth, which
   // is what distinguishes it from guard (b).
   setupProfileState(Mesh, VCoord, State, EqState, NLayers, DC, Bulge, Linear,
                     2.0e4_Real);
   PGradHostState HD = copyStateToHost(VCoord, EqState);

   std::vector<Real> DZCorrect, TendCorrect, DZBroken, TendBroken;
   Rules = PGradGuardRules();
   referenceScan(HD.PressInterface, HD.PressMid, HD.Temp, HD.Salt, HD.SpecVol,
                 HD.SpecVolDCt, HD.SpecVolDSa, HD.SpecVolDP, HD.GeomZ, NLayers,
                 DC, NQuad, Rules, DZCorrect, TendCorrect);
   Rules               = PGradGuardRules();
   Rules.ShiftedAnchor = false;
   referenceScan(HD.PressInterface, HD.PressMid, HD.Temp, HD.Salt, HD.SpecVol,
                 HD.SpecVolDCt, HD.SpecVolDSa, HD.SpecVolDP, HD.GeomZ, NLayers,
                 DC, NQuad, Rules, DZBroken, TendBroken);

   Real GuardD    = 0.0_Real;
   Real OffsetMin = std::numeric_limits<Real>::max();
   Real OffsetMax = 0.0_Real;
   for (int K = 0; K <= NLayers; ++K) {
      const Real Offset = std::abs(DZBroken[K] - DZCorrect[K]);
      OffsetMin         = std::min(OffsetMin, Offset);
      OffsetMax         = std::max(OffsetMax, Offset);
   }
   for (int K = 0; K < NLayers; ++K)
      GuardD = std::max(GuardD, std::abs(TendBroken[K] - TendCorrect[K]));

   LOG_INFO("PGradTest: guard (d) unshifted anchor changes the tendency by {} "
            "m/s2; the offset in the height difference runs from {} m to {} m "
            "over the column",
            GuardD, OffsetMin, OffsetMax);
   const bool Fires = GuardD > 1.0e4_Real * Eps * TendScale;
   const bool Flat  = OffsetMax > 0.0_Real &&
                     (OffsetMax - OffsetMin) <= 1.0e-6_Real * OffsetMax;
   if (Fires && Flat) {
      LOG_INFO("PGradTest: guard (d) PASS: the anchor shift is load-bearing "
               "and its omission is flat with depth");
   } else {
      LOG_ERROR("PGradTest: guard (d) FAIL: fires {}, flat with depth {}",
                Fires, Flat);
      ++Err;
   }

   //
   // Group two: profiles the reconstruction does not resolve. The residual
   // must shrink like the square of the layer thickness under vertical
   // refinement at fixed tilt, matching the design's table. A residual that
   // does not shrink at that rate means one of the scheme's conditions has
   // been broken; it is a bug to find, not a tolerance to widen.
   //
   const PGradProfile Profiles[2] = {Quadratic, Cubic};
   const char *ProfileNames[2]    = {"quadratic", "cubic"};

   for (int IProf = 0; IProf < 2; ++IProf) {
      Real PrevErr   = 0.0_Real;
      Real WorstRate = 1.0e30_Real;
      // 15, 30, 60 layers: a refinement factor of two, staying within the
      // vertical dimension the mesh file provides
      const I4 RefLayers[3] = {15, 30, 60};
      for (int IRef = 0; IRef < 3; ++IRef) {
         const I4 NRef = RefLayers[IRef];
         setupProfileState(Mesh, VCoord, State, EqState, NRef, DC, Bulge,
                           Profiles[IProf], 0.0_Real);
         Real DZ, Up, Low, Rms;
         const Real MaxE = runFiniteVolume(Mesh, VCoord, State, EqState, NRef,
                                           DZ, Up, Low, Rms);
         if (IRef > 0 && Rms > 0.0_Real) {
            const Real Rate = std::log2(PrevErr / Rms);
            WorstRate       = std::min(WorstRate, Rate);
            LOG_INFO("PGradTest: {} profile, {} layers: RMS |Tend| = {} m/s2 "
                     "(max {}), rate {}",
                     ProfileNames[IProf], NRef, Rms, MaxE, Rate);
         } else {
            LOG_INFO("PGradTest: {} profile, {} layers: RMS |Tend| = {} m/s2 "
                     "(max {})",
                     ProfileNames[IProf], NRef, Rms, MaxE);
         }
         PrevErr = Rms;
      }
      if (WorstRate >= 1.7_Real) {
         LOG_INFO("PGradTest: {} profile convergence PASS: worst rate {}",
                  ProfileNames[IProf], WorstRate);
      } else {
         LOG_ERROR("PGradTest: {} profile convergence FAIL: worst rate {} is "
                   "below the second order the design predicts",
                   ProfileNames[IProf], WorstRate);
         ++Err;
      }
   }

   return Err;

} // end testExactnessAndGuards

// The anchor where the two columns' floors are more than a layer apart in
// pressure.
//
// Terrain-following layers over a steep slope put one column's floor several
// layers deeper than its neighbour's. Shifting both floors to their mean
// pressure with one first-order expansion leaves an error that grows with the
// cube of the offset and is the same at every depth; comparing the shallower
// floor with the deeper column's own interface just below it keeps every
// shift within one layer. The state is at rest and linear in pressure, so
// what remains is the difference between the two columns' VertCoord heights,
// which converges with layer thickness rather than growing with the slope.
int testSlopingFloorAnchor(HorzMesh *Mesh,    ///< [in] Horizontal mesh
                           VertCoord *VCoord, ///< [inout] Vertical coordinate
                           OceanState *State, ///< [inout] Ocean state
                           Eos *EqState       ///< [inout] Equation of state
) {

   int Err = 0;

   const Real Eps        = std::numeric_limits<Real>::epsilon();
   const bool DoublePrec = Eps < 1.0e-10_Real;
   const Real DC         = 4000.0_Real;
   const I4 NLayers      = 32;
   const I4 NQuad        = 2;
   // 800 dbar, or 6.4 layers: a floor sloping at 200 m/km at this spacing
   const Real FloorPressDiff = 8.0e6_Real;

   const PGradProfile Linear = {12.0_Real, -2.0e-7_Real, 0.0_Real, 0.0_Real,
                                34.0_Real, 2.5e-8_Real,  0.0_Real, 0.0_Real};

   setupSlopingFloorState(Mesh, VCoord, State, EqState, NLayers, DC,
                          FloorPressDiff, Linear);

   Real MaxDeltaZ, MaxUpper, MaxLower, RmsFV;
   const Real MaxFV = runFiniteVolume(Mesh, VCoord, State, EqState, NLayers,
                                      MaxDeltaZ, MaxUpper, MaxLower, RmsFV);
   const Real MaxCentered = runCentered(Mesh, VCoord, State, EqState, NLayers);

   // the same assembly on the host, and again with both floors shifted to
   // their mean pressure
   PGradHostState H = copyStateToHost(VCoord, EqState);
   std::vector<Real> RefDeltaZ, RefTend, MeanDeltaZ, MeanTend;
   PGradGuardRules Rules;
   referenceScan(H.PressInterface, H.PressMid, H.Temp, H.Salt, H.SpecVol,
                 H.SpecVolDCt, H.SpecVolDSa, H.SpecVolDP, H.GeomZ, NLayers, DC,
                 NQuad, Rules, RefDeltaZ, RefTend);
   Rules.AnchorWithinOneLayer = false;
   referenceScan(H.PressInterface, H.PressMid, H.Temp, H.Salt, H.SpecVol,
                 H.SpecVolDCt, H.SpecVolDSa, H.SpecVolDP, H.GeomZ, NLayers, DC,
                 NQuad, Rules, MeanDeltaZ, MeanTend);

   Real RefMax  = 0.0_Real;
   Real MeanMax = 0.0_Real;
   for (int K = 0; K < NLayers; ++K) {
      RefMax  = std::max(RefMax, std::abs(RefTend[K]));
      MeanMax = std::max(MeanMax, std::abs(MeanTend[K]));
   }

   const Real Ratio     = MaxFV / MaxCentered;
   const Real MeanRatio = MeanMax / MaxCentered;
   LOG_INFO("PGradTest: sloping floor, {} Pa offset: FiniteVolume max |Tend| "
            "= {} m/s2, {} of Centered's {} m/s2; with both floors shifted to "
            "their mean pressure {} m/s2, {} of Centered's",
            FloorPressDiff, MaxFV, Ratio, MaxCentered, MeanMax, MeanRatio);

   // The anchor forms and cancels the height difference between the two
   // floors, so that is the scale of the round-off the kernel and the host
   // reference may differ by.
   const Real FloorStep   = std::abs(H.GeomZ(1, NLayers) - H.GeomZ(0, NLayers));
   const Real FidelityTol = 64.0_Real * Eps * Gravity * FloorStep / DC;
   if (std::abs(RefMax - MaxFV) <= FidelityTol) {
      LOG_INFO("PGradTest: sloping-floor harness fidelity PASS: reference {} "
               "m/s2 against kernel {} m/s2",
               RefMax, MaxFV);
   } else {
      LOG_ERROR("PGradTest: sloping-floor harness fidelity FAIL: reference {} "
                "m/s2 against kernel {} m/s2, tolerance {}",
                RefMax, MaxFV, FidelityTol);
      ++Err;
   }

   // Gate set from the measured 0.014 of Centered in double precision. In
   // single precision the floor heights near 4000 m carry round-off of order
   // 1e-4 m, a tendency of order 1e-7 m/s2 at this spacing, so the gate is
   // loosened to match while staying well below the mean-pressure anchor.
   const Real RatioTol = DoublePrec ? 0.05_Real : 0.15_Real;
   if (Ratio <= RatioTol) {
      LOG_INFO("PGradTest: sloping-floor anchor PASS");
   } else {
      LOG_ERROR("PGradTest: sloping-floor anchor FAIL: {} of Centered against "
                "{}",
                Ratio, RatioTol);
      ++Err;
   }

   // The gate is only worth having if the anchor it protects is what decides
   // it: shifting both floors to their mean pressure, measured at 0.45 of
   // Centered, must fail it by a wide margin.
   if (MeanRatio >= 0.25_Real) {
      LOG_INFO("PGradTest: sloping-floor anchor guard PASS: the mean-pressure "
               "anchor gives {} of Centered",
               MeanRatio);
   } else {
      LOG_ERROR("PGradTest: sloping-floor anchor guard FAIL: the mean-pressure "
                "anchor gives only {} of Centered, so the gate cannot tell the "
                "two anchors apart",
                MeanRatio);
      ++Err;
   }

   return Err;

} // end testSlopingFloorAnchor

// Assert the identity of design section 3.9 (test section 5.5): with the tidal
// and self-attraction-and-loading potentials zero, PressureGradCentered is
// exactly
//
//    -g / d_e * S_{e,k}
//    S_{e,k} = 1/2 (dZ_k + dZ_{k+1}) + alphaBar / (2 g) (dq_k + dq_{k+1})
//
// where dZ and dq are the cross-edge differences of GeomZInterface and
// PressureInterface and alphaBar is the edge average of SpecVol. S is the
// first-order conversion of a height difference taken at fixed layer index
// into one taken at fixed pressure, so this pins what error the centered
// scheme makes as well as what it computes.
//
// This is a permanent regression test rather than a transitional one. Because
// the two expressions read the mesh, VertCoord and Eos state through
// independently written code, their agreement tests the shared upstream state
// -- edge masks, interface indexing, VertCoord conventions -- and not just the
// pressure gradient arithmetic. Expect agreement to round-off rather than
// bit-for-bit: PressureMid may be formed as zBot + h/2 rather than as the mean
// of the two interface pressures, which is algebraically but not bitwise the
// same.
int testCenteredIdentity(HorzMesh *Mesh,    ///< [in] Horizontal mesh
                         VertCoord *VCoord, ///< [inout] Vertical coordinate
                         OceanState *State, ///< [inout] Ocean state
                         Eos *EqState       ///< [inout] Equation of state
) {

   int Err = 0;

   PressureGrad *DefPGrad = PressureGrad::getDefault();
   if (!DefPGrad || DefPGrad->getType() != PressureGradType::Centered) {
      LOG_ERROR("PGradTest: centered identity needs a Centered default");
      return 1;
   }

   const I4 NVertLayers = 60;
   const Real DC        = 30000.0_Real;
   const Real ZBottom   = -1000.0_Real;
   const Real Eps       = std::numeric_limits<Real>::epsilon();

   // Round-off in the centered functor is set by the size of the hydrostatic
   // terms it forms and cancels, so the tolerance is scaled by those rather
   // than by their differences. Measured at or below half an epsilon; the
   // factor here is headroom for compiler and precision differences.
   const Real Tol = 16.0_Real;

   // A sweep of tilts, with 0.5 the untilted case
   const std::vector<Real> Tilts = {0.5, 0.49, 0.45, 0.4, 0.3, 0.2, 0.1, 0.05};

   Real MaxResid = 0.0_Real;
   Real MinGuard = std::numeric_limits<Real>::max();

   for (Real TiltFactor : Tilts) {

      setupTwoColumnState(Mesh, VCoord, State, EqState, NVertLayers, DC,
                          TiltFactor, ZBottom);

      Array2DReal Tend("TendCentered", Mesh->NEdgesSize, NVertLayers);
      deepCopy(Tend, 0.0_Real);

      Array2DReal PseudoThick = State->getPseudoThickness(0);
      Array2DReal Temp        = Tracers::getByName(0, "Temperature");
      Array2DReal Salinity    = Tracers::getByName(0, "Salinity");

      const auto &PressureMid       = VCoord->PressureMid;
      const auto &PressureInterface = VCoord->PressureInterface;
      const auto &GeomZInterface    = VCoord->GeomZInterface;
      const auto &EdgeMask          = VCoord->EdgeMask;
      const auto &SpecVol           = EqState->SpecVol;
      const auto &CellsOnEdge       = Mesh->CellsOnEdge;
      const auto &DcEdge            = Mesh->DcEdge;

      DefPGrad->computePressureGrad(Tend, PressureMid, PressureInterface,
                                    SpecVol, GeomZInterface, PseudoThick, Temp,
                                    Salinity, EqState);

      // Residual of the identity, in units of eps times the hydrostatic
      // scale, and the same measure for a deliberately broken version of the
      // identity that drops the pressure half of S. The guard is what makes
      // the check non-vacuous: a mis-derivation must miss by far more than
      // the tolerance.
      Real LocResid = 0.0_Real;
      Real LocGuard = std::numeric_limits<Real>::max();
      parallelReduce(
          {Mesh->NEdgesAll, NVertLayers},
          KOKKOS_LAMBDA(int IEdge, int K, Real &MaxR, Real &MinG) {
             const I4 ICell0  = CellsOnEdge(IEdge, 0);
             const I4 ICell1  = CellsOnEdge(IEdge, 1);
             const Real InvDc = 1.0_Real / DcEdge(IEdge);

             const Real Z0K  = GeomZInterface(ICell0, K);
             const Real Z1K  = GeomZInterface(ICell1, K);
             const Real Z0K1 = GeomZInterface(ICell0, K + 1);
             const Real Z1K1 = GeomZInterface(ICell1, K + 1);

             const Real Q0K  = PressureInterface(ICell0, K);
             const Real Q1K  = PressureInterface(ICell1, K);
             const Real Q0K1 = PressureInterface(ICell0, K + 1);
             const Real Q1K1 = PressureInterface(ICell1, K + 1);

             const Real AlphaBar =
                 0.5_Real * (SpecVol(ICell0, K) + SpecVol(ICell1, K));

             const Real ShiftZ = 0.5_Real * ((Z1K - Z0K) + (Z1K1 - Z0K1));
             const Real ShiftQ = AlphaBar / (2.0_Real * Gravity) *
                                 ((Q1K - Q0K) + (Q1K1 - Q0K1));

             const Real Expected =
                 EdgeMask(IEdge, K) * (-Gravity * InvDc) * (ShiftZ + ShiftQ);
             const Real Broken =
                 EdgeMask(IEdge, K) * (-Gravity * InvDc) * ShiftZ;

             // the hydrostatic terms the centered functor forms and cancels
             const Real HydroScale =
                 InvDc * 0.25_Real *
                 (Gravity * (Kokkos::abs(Z0K) + Kokkos::abs(Z1K) +
                             Kokkos::abs(Z0K1) + Kokkos::abs(Z1K1)) +
                  AlphaBar * (Kokkos::abs(Q0K) + Kokkos::abs(Q1K) +
                              Kokkos::abs(Q0K1) + Kokkos::abs(Q1K1)));

             const Real Resid = Kokkos::abs(Tend(IEdge, K) - Expected);
             const Real Guard = Kokkos::abs(Tend(IEdge, K) - Broken);

             // masked edges carry no tendency, so neither expression has
             // anything to say there
             if (EdgeMask(IEdge, K) > 0.0_Real && HydroScale > 0.0_Real) {
                const Real Scale = Eps * HydroScale;
                if (Resid / Scale > MaxR)
                   MaxR = Resid / Scale;
                if (Guard / Scale < MinG)
                   MinG = Guard / Scale;
             }
          },
          Kokkos::Max<Real>(LocResid), Kokkos::Min<Real>(LocGuard));

      LOG_INFO("PGradTest: centered identity at tilt {}: residual {} eps, "
               "dropped-pressure-term guard {} eps",
               TiltFactor, LocResid, LocGuard);

      if (LocResid > MaxResid)
         MaxResid = LocResid;
      // the untilted case has nothing for the guard to detect
      if (TiltFactor != 0.5_Real && LocGuard < MinGuard)
         MinGuard = LocGuard;
   }

   if (MaxResid <= Tol) {
      LOG_INFO("PGradTest: centered identity PASS: max residual {} eps",
               MaxResid);
   } else {
      LOG_ERROR("PGradTest: centered identity FAIL: max residual {} eps "
                "exceeds {} eps",
                MaxResid, Tol);
      ++Err;
   }

   // A mis-derivation that drops the pressure half of S misses by many orders
   // of magnitude; if it does not, the identity is being satisfied trivially.
   //
   // The comparison is against the residual rather than against an absolute
   // count of epsilons. Both quantities are measured in epsilons of the same
   // scale, so their ratio is what the check is really about, and it means
   // the same thing in either precision: an absolute threshold of a thousand
   // epsilons happens to leave a comfortable margin in double precision and
   // only a factor of 1.4 in single, which is too near the edge for a check
   // whose whole purpose is to prove the identity is not vacuous.
   const Real GuardRatio =
       (MaxResid > 0.0_Real) ? MinGuard / MaxResid : 0.0_Real;
   if (GuardRatio > 100.0_Real) {
      LOG_INFO("PGradTest: centered identity guard PASS: dropping the "
               "pressure term misses by {} eps, {} times the residual",
               MinGuard, GuardRatio);
   } else {
      LOG_ERROR("PGradTest: centered identity guard FAIL: dropping the "
                "pressure term changes the answer by only {} eps, {} times "
                "the residual",
                MinGuard, GuardRatio);
      ++Err;
   }

   return Err;

} // end testCenteredIdentity

// Check that the PressureGrad configuration group is parsed as expected and
// that a FiniteVolume instance can be created and dispatched to. The
// FiniteVolume sub-options are read from the same group as PressureGradType
// and take their default values when a key is absent.
int testPGradConfig(const HorzMesh *Mesh,   ///< [in] Horizontal mesh
                    const VertCoord *VCoord ///< [in] Vertical coordinate
) {

   int Err = 0;

   // The default instance is built from the PressureGrad group in the
   // config, which initPGradTest set to the centered scheme
   PressureGrad *DefPGrad = PressureGrad::getDefault();

   if (DefPGrad->getType() == PressureGradType::Centered) {
      LOG_INFO("PGradTest: default PressureGradType parse PASS");
   } else {
      LOG_ERROR("PGradTest: default PressureGradType parse FAIL");
      ++Err;
   }

   // The only implemented sub-option values, whether read from the config
   // or taken from the defaults
   if (DefPGrad->getHorzOrder() == 2 &&
       DefPGrad->getVertRecon() == PressureGradVertRecon::Linear &&
       DefPGrad->getQuadraturePoints() >= 1) {
      LOG_INFO("PGradTest: FiniteVolume sub-option parse PASS");
   } else {
      LOG_ERROR("PGradTest: FiniteVolume sub-option parse FAIL: HorzOrder={} "
                "QuadraturePoints={}",
                DefPGrad->getHorzOrder(), DefPGrad->getQuadraturePoints());
      ++Err;
   }

   // Create a second instance that selects the FiniteVolume scheme. The
   // sub-config shares its node with the parent, so resetting the type here
   // is what the new instance sees. Note the explicit std::string: a string
   // literal would select the boolean overload of Config::set.
   Config *Options = Config::getOmegaConfig();
   Config PGradConfig("PressureGrad");
   Err += (Options->get(PGradConfig).isFail() ? 1 : 0);
   PGradConfig.set("PressureGradType", std::string("FiniteVolume"));

   PressureGrad *FVPGrad =
       PressureGrad::create("TestFiniteVolume", Mesh, VCoord, Options);
   CHECK_ERROR_ABORT(
       Error(FVPGrad ? ErrorCode::Success : ErrorCode::Fail, ""),
       "PGradTest: could not create the TestFiniteVolume pressure gradient");

   if (FVPGrad && FVPGrad->getType() == PressureGradType::FiniteVolume) {
      LOG_INFO("PGradTest: FiniteVolume PressureGradType parse PASS");
   } else {
      LOG_ERROR("PGradTest: FiniteVolume PressureGradType parse FAIL");
      ++Err;
   }

   // Dispatch check: computePressureGrad must take the FiniteVolume branch
   // and produce a tendency of the same order as the centered scheme's on the
   // same state.
   //
   // The specific volume derivatives are computed first. In a run they are
   // filled by AuxiliaryState, which calls computeSpecVolAndDerivs whenever
   // the FiniteVolume scheme is selected; here the state was built by a
   // helper that fills SpecVol alone, so without this the scheme would read
   // the fill value the fields were attached with. That is worth being
   // explicit about: checking only for finiteness let this go unnoticed in
   // double precision, where a fill value of 1e37 propagates through the
   // arithmetic without overflowing, and it surfaced only in a
   // single-precision build, where it does.
   OceanState *State       = OceanState::getDefault();
   Eos *EqState            = Eos::getInstance();
   Array2DReal PseudoThick = State->getPseudoThickness(0);
   Array2DReal TendFV("TendFiniteVolume", Mesh->NEdgesSize,
                      VCoord->NVertLayers);
   deepCopy(TendFV, 0.0_Real);

   Array2DReal ConservTemp = Tracers::getByName(0, "Temperature");
   Array2DReal AbsSalinity = Tracers::getByName(0, "Salinity");

   EqState->computeSpecVolAndDerivs(ConservTemp, AbsSalinity,
                                    VCoord->PressureMid);

   FVPGrad->computePressureGrad(
       TendFV, VCoord->PressureMid, VCoord->PressureInterface, EqState->SpecVol,
       VCoord->GeomZInterface, PseudoThick, ConservTemp, AbsSalinity, EqState);

   Array2DReal TendCtr("TendCentered", Mesh->NEdgesSize, VCoord->NVertLayers);
   deepCopy(TendCtr, 0.0_Real);
   PressureGrad::getDefault()->computePressureGrad(
       TendCtr, VCoord->PressureMid, VCoord->PressureInterface,
       EqState->SpecVol, VCoord->GeomZInterface, PseudoThick, ConservTemp,
       AbsSalinity, EqState);

   I4 NBad     = 0;
   Real MaxFV  = 0.0_Real;
   Real MaxCtr = 0.0_Real;
   parallelReduce(
       {Mesh->NEdgesAll, VCoord->NVertLayers},
       KOKKOS_LAMBDA(int IEdge, int K, I4 &LSum, Real &MaxF, Real &MaxC) {
          const Real F = Kokkos::abs(TendFV(IEdge, K));
          const Real C = Kokkos::abs(TendCtr(IEdge, K));
          if (!Kokkos::isfinite(TendFV(IEdge, K))) {
             ++LSum;
          } else if (F > MaxF) {
             MaxF = F;
          }
          if (C > MaxC)
             MaxC = C;
       },
       Kokkos::Sum<I4>(NBad), Kokkos::Max<Real>(MaxFV),
       Kokkos::Max<Real>(MaxCtr));

   // What can be asserted here is one-sided. The finite-volume tendency may
   // legitimately be very much *smaller* than the centered one -- that is the
   // whole point of the scheme, and this state has a profile close to linear
   // in pressure, so it is -- but it cannot be very much larger, since the two
   // agree to second order in the cross-edge pressure difference. An upper
   // bound plus a nonzero lower bound is therefore the honest check: it
   // catches a fill value, which is off by thirty orders, and it catches a
   // dispatch that quietly did nothing, without forbidding the scheme from
   // working.
   const Real Ratio = (MaxCtr > 0.0_Real) ? MaxFV / MaxCtr : 0.0_Real;

   LOG_INFO("PGradTest: FiniteVolume dispatch: max |Tend| = {} m/s2 against "
            "Centered {} m/s2, ratio {}",
            MaxFV, MaxCtr, Ratio);

   if (NBad == 0 && MaxFV > 0.0_Real && Ratio < 1.0e3_Real) {
      LOG_INFO("PGradTest: FiniteVolume dispatch PASS");
   } else {
      LOG_ERROR("PGradTest: FiniteVolume dispatch FAIL: {} non-finite values, "
                "max |Tend| {} m/s2, ratio to the centered scheme {}",
                NBad, MaxFV, Ratio);
      ++Err;
   }

   PGradConfig.set("PressureGradType", std::string("Centered"));
   PressureGrad::erase("TestFiniteVolume");

   return Err;

} // end testPGradConfig

int main(int argc, char *argv[]) {
   int RetVal = 0;

   MPI_Init(&argc, &argv);
   Kokkos::initialize();
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");

   {
      initPGradTest();

      // Flush every message rather than only those at warn and above. Omega's
      // default buffers info messages, so a run that aborts loses the whole
      // record of what it had measured, which is the case that most needs it.
      // This test is small enough that the cost does not matter.
      spdlog::flush_on(spdlog::level::info);

      // Initialize default PressureGrad
      PressureGrad::init();

      HorzMesh *DefMesh    = HorzMesh::getDefault();
      VertCoord *VCoord    = VertCoord::getDefault();
      OceanState *DefState = OceanState::getDefault();
      Eos *DefEos          = Eos::getInstance();

      // create arrays: Tend on edges, Pressure/Geopotential/SpecVol on cells
      Array2DReal Tend("Tend", DefMesh->NEdgesSize, VCoord->NVertLayers);

      I4 NEdgesAll = DefMesh->NEdgesAll;

      I4 NVertLayers  = 60;
      Real DC         = 30000.0_Real;
      Real TiltFactor = 0.495_Real;
      Real ZBottom    = -1000.0_Real;
      I4 NRefinements = 4;
      HostArray1DReal Rmse("Rmse", NRefinements);
      for (int Refinement = 0; Refinement < NRefinements; ++Refinement) {

         LOG_INFO("PGradTest: Starting refinement level {}", Refinement);

         setupTwoColumnState(DefMesh, VCoord, DefState, DefEos, NVertLayers, DC,
                             TiltFactor, ZBottom);

         Array2DReal PseudoThick = DefState->getPseudoThickness(0);
         Array2DReal Temp        = Tracers::getByName(0, "Temperature");
         Array2DReal Salinity    = Tracers::getByName(0, "Salinity");
         auto &SpecVol           = DefEos->SpecVol;
         auto &PressureMid       = VCoord->PressureMid;
         auto &GeomZInterface    = VCoord->GeomZInterface;

         LOG_INFO("NVertLayers = {}", NVertLayers);
         LOG_INFO("dC = {}", DC);
         DefState->copyToHost(0);
         HostArray2DReal PseudoThickH = DefState->getPseudoThicknessH(0);
         for (int i = 0; i < 2; ++i) {
            for (int k = 0; k < 2; ++k) {
               LOG_INFO("PseudoThick({}, {}) = {}", i, k, PseudoThickH(i, k));
            }
         }

         // get PressureGrad instance
         PressureGrad *DefPGrad = PressureGrad::getDefault();
         if (!DefPGrad) {
            LOG_INFO("PGrad: default instance not created by init");
         }

         // compute pressure gradient
         deepCopy(Tend, 0.0_Real);

         const auto &PressureInterface = VCoord->PressureInterface;
         DefPGrad->computePressureGrad(Tend, PressureMid, PressureInterface,
                                       SpecVol, GeomZInterface, PseudoThick,
                                       Temp, Salinity, DefEos);

         // compute errors
         Real MaxValue = 0.0_Real;
         parallelReduce(
             {NEdgesAll, NVertLayers - 2},
             KOKKOS_LAMBDA(int i, int k, Real &max) {
                Real Val = Kokkos::abs(Tend(i, k + 1));
                if (Val > max)
                   max = Val;
             },
             Kokkos::Max<Real>(MaxValue));
         Real SumValue = 0.0_Real;
         parallelReduce(
             {NEdgesAll, NVertLayers - 2},
             KOKKOS_LAMBDA(int i, int k, Real &LSum) {
                LSum += Tend(i, k + 1) * Tend(i, k + 1);
             },
             Kokkos::Sum<Real>(SumValue));
         Real RmseVal = std::sqrt(SumValue / (NEdgesAll * (NVertLayers - 2)));
         Rmse(Refinement) = RmseVal;

         LOG_INFO("refinement level {}: max |Tend| = {}, average Tend = {}",
                  Refinement, MaxValue, RmseVal);

         // coarsen for next iteration
         DC          = DC * 2.0_Real;
         NVertLayers = NVertLayers / 2;

      } // refinement loop

      // Test the reconstruction estimator on its own. It needs no mesh, and
      // it is the most likely place for the exactness gate to fail, so a
      // failure here localizes immediately.
      int ReconErr = testReconstruction();

      // Pin the per-column pressure lookup. No answer-level check anywhere
      // can distinguish this from a layer-index lookup, so these property
      // tests are the only protection against that failure mode.
      int LookupErr = testPressureLookup();

      // The matched-pressure integrand must be zero at every quadrature
      // point on a resolved profile, not merely in the integral.
      int IntegrandErr = testMatchedPressIntegrand();

      // The gating test: exactness on the exact set, the convergence of the
      // residual off it, and the guards.
      int ScanErr = testExactnessAndGuards(DefMesh, VCoord, DefState, DefEos);

      // The anchor where the two floors are more than a layer apart.
      int SlopeErr = testSlopingFloorAnchor(DefMesh, VCoord, DefState, DefEos);

      // The centered identity of design section 3.9 across a sweep of tilts.
      // This resets the two-column state, so it runs after the others.
      int IdentityErr = testCenteredIdentity(DefMesh, VCoord, DefState, DefEos);

      // Test parsing and dispatch of the PressureGrad configuration options
      int ConfigErr = testPGradConfig(DefMesh, VCoord);

      // Test for second order convergence of the centered scheme under
      // refinement; resolution (dC) increases in the refinement loop.
      //
      // Double precision only. This measures PressureGradCentered's
      // truncation error, and a single-precision build cannot resolve it at
      // the finest resolution: the centered scheme forms and cancels
      // Montgomery potentials of order 1e4 m2 s-2, so dividing by a 30 km cell
      // spacing leaves a round-off floor near 4e-8 m s-2. The finest
      // resolution sits below that, measuring round-off rather than
      // truncation -- 4.9e-8 in single precision against 7.1e-9 in double on
      // the same state -- so the convergence relation cannot hold there. That
      // is a property of the scheme being measured, not a defect: it is the
      // arithmetic the finite-volume scheme avoids by differencing the
      // integrand before integrating it.
      const bool DoublePrecision =
          std::numeric_limits<Real>::epsilon() < 1.0e-10_Real;

      RetVal = 0;
      if (DoublePrecision) {
         if (Rmse(0) >=
             Rmse(NRefinements - 1) / pow(4.0_Real, NRefinements - 1))
            RetVal = 1;
      } else {
         LOG_INFO("PGradTest: single precision: skipping the centered "
                  "refinement convergence check; its finest resolution gives "
                  "{}, which is the round-off floor of a scheme that cancels "
                  "large quantities rather than its truncation error",
                  Rmse(0));
      }

      RetVal += ReconErr + LookupErr + IntegrandErr + ScanErr + SlopeErr +
                IdentityErr + ConfigErr;

      // Flush before teardown. Omega logs at info level but only flushes at
      // warn and above, so a run that aborts during cleanup loses every
      // measurement it made -- which is exactly the case that most needs the
      // record.
      spdlog::default_logger()->flush();

      // cleanup
      PressureGrad::clear();
      IOStream::finalize();
      TimeStepper::clear();
      Tracers::clear();
      Eos::destroyInstance();
      OceanState::clear();
      VertCoord::clear();
      Field::clear();
      Dimension::clear();
      HorzMesh::clear();
      Halo::clear();
      Decomp::clear();
      MachEnv::removeAll();
   }
   Pacer::finalize();
   Kokkos::finalize();
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();

   if (RetVal >= 256)
      RetVal = 255;
   return RetVal;
}
